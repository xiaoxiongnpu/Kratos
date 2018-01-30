//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Vicente Mataix
//

#if !defined(KRATOS_RESIDUALBASED_RAMM_ARC_LENGTH_STRATEGY)
#define KRATOS_RESIDUALBASED_RAMM_ARC_LENGTH_STRATEGY

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes

namespace Kratos
{

template<class TSparseSpace,class TDenseSpace,class TLinearSolver>

class ResidualBasedRammArcLengthStrategy : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedRammArcLengthStrategy);
    
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> MotherType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef TSparseSpace SparseSpaceType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    using MotherType::mpConvergenceCriteria;
    using MotherType::mpScheme;
    using MotherType::mpBuilderAndSolver;
    using MotherType::mpA; //Tangent matrix
    using MotherType::mpb; //Residual vector of iteration i
    using MotherType::mpDx; //Delta x of iteration i
    using MotherType::mReformDofSetAtEachStep;
    using MotherType::mCalculateReactionsFlag;
    using MotherType::mSolutionStepIsInitialized;
    using MotherType::mMaxIterationNumber;
    using MotherType::mInitializeWasPerformed;
    using MotherType::mSubModelPartList;
    using MotherType::mVariableNames;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    ResidualBasedRammArcLengthStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters& rParameters,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {
            //only include validation with c++11 since raw_literals do not exist in c++03
            Parameters default_parameters( R"(
            {
                "desired_iterations": 4,
                "max_radius_factor": 20.0,
                "min_radius_factor": 0.5,
                "loads_sub_model_part_list": [],
                "loads_variable_list" : []
            }  )" );
            
            // Validate agains defaults -- this also ensures no type mismatch
            rParameters.ValidateAndAssignDefaults(default_parameters);
            
            mpParameters = &rParameters;
            
            // Set Load SubModelParts and Variable names
            if(rParameters["loads_sub_model_part_list"].size() > 0) {
                mSubModelPartList.resize(rParameters["loads_sub_model_part_list"].size());
                mVariableNames.resize(rParameters["loads_variable_list"].size());

                KRATOS_ERROR_IF(mSubModelPartList.size() != mVariableNames.size()) << "For each SubModelPart there must be a corresponding nodal Variable" << std::endl;

                for(unsigned int i = 0; i < mVariableNames.size(); i++) {
                    mSubModelPartList[i] = &( rModelPart.GetSubModelPart(rParameters["loads_sub_model_part_list"][i].GetString()) );
                    mVariableNames[i] = rParameters["loads_variable_list"][i].GetString();
                }
            }

            mDesiredIterations = rParameters["desired_iterations"].GetInt();
            mMaxRadiusFactor = rParameters["max_radius_factor"].GetDouble();
            mMinRadiusFactor = rParameters["min_radius_factor"].GetDouble();
        }

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~ResidualBasedRammArcLengthStrategy() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize() override
    {
        KRATOS_TRY

        if (mInitializeWasPerformed == false)
		{
            MotherType::Initialize();
            
            //set up the system
            if (mpBuilderAndSolver->GetDofSetIsInitializedFlag() == false)
            {
                //setting up the list of the DOFs to be solved
                mpBuilderAndSolver->SetUpDofSet(mpScheme, BaseType::GetModelPart());

                //shaping correctly the system
                mpBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());
            }
            
            // Compute initial radius (mRadius_0)
            mpBuilderAndSolver->ResizeAndInitializeVectors(mpScheme, mpA, mpDx, mpb, BaseType::GetModelPart().Elements(),
                                                            BaseType::GetModelPart().Conditions(), BaseType::GetModelPart().GetProcessInfo());
            TSystemMatrixType& mA = *mpA;
            TSystemVectorType& mDx = *mpDx;
            TSystemVectorType& mb = *mpb;
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);

            mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
            
            mRadius_0 = TSparseSpace::TwoNorm(mDx);
            mRadius = mRadius_0;

            // Compute vector of reference external force (mf)
            this->InitializeSystemVector(mpf);
            TSystemVectorType& mf = *mpf;
            TSparseSpace::SetToZero(mf);

            mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mf);
            
            //Initialize the loading factor Lambda
            mLambda = 0.0;
            mLambda_old = 1.0;
            
            // Initialize Norm of solution
            mNormxEquilibrium = 0.0;
            
            if (BaseType::mEchoLevel > 0)
              std::cout << "Ramm's Arc Length Strategy Initialized" << std::endl;

        }

        KRATOS_CATCH( "" )
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY

		if (mSolutionStepIsInitialized == false)
		{
            MotherType::InitializeSolutionStep();
            
            this->SaveInitializeSystemVector(mpf);
            this->InitializeSystemVector(mpDxf);
            this->InitializeSystemVector(mpDxb);
            this->InitializeSystemVector(mpDxPred);
            this->InitializeSystemVector(mpDxStep);
        }
        
        KRATOS_CATCH( "" )
    }

	bool SolveSolutionStep() override
	{
        // ********** Prediction phase **********
        if (BaseType::mEchoLevel > 0) 
            std::cout << "ARC-LENGTH RADIUS: " << mRadius/mRadius_0 << " X initial radius" << std::endl;
        
        // Initialize variables
		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;
        TSystemVectorType& mf = *mpf;
        TSystemVectorType& mDxb = *mpDxb;
        TSystemVectorType& mDxf = *mpDxf;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;
        
        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
                
        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mDxf);
        
        // Note: This is not so efficient, but I want to solve mA*mDxf=mf without losing mf
        this->BuildWithDirichlet(mA, mDxf, mb);
        noalias(mb) = mf;
        mpBuilderAndSolver->SystemSolve(mA, mDxf, mb);

        //update results
        double DLambda = mRadius/TSparseSpace::TwoNorm(mDxf);
        mDLambdaStep = DLambda;
        mLambda += DLambda;
        noalias(mDxPred) = DLambda*mDxf;
        noalias(mDxStep) = mDxPred;
        this->Update(rDofSet, mA, mDxPred, mb);

        //move the mesh if needed
        if(BaseType::MoveMeshFlag() == true)
            BaseType::MoveMesh();
        
        // ********** Correction phase (iteration cicle) **********

        //initializing the parameters of the iteration loop
        bool is_converged = false;
        mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), rDofSet, mA, mDxf, mb);
        if (mpConvergenceCriteria->GetActualizeRHSflag() == true) {
            TSparseSpace::SetToZero(mb);
            mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
        }
        is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDxf, mb);

        unsigned int iteration_number = 0;
        double NormDx;
        
        while (is_converged == false && iteration_number < mMaxIterationNumber) {
            // Setting the number of iteration
            iteration_number += 1;
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            
            mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDxf);
            
            // Note: This is not so efficient, but I want to solve mA*mDxf=mf without losing mf
            this->BuildWithDirichlet(mA, mDxf, mb);
            noalias(mb) = mf;
            mpBuilderAndSolver->SystemSolve(mA, mDxf, mb);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDxb);
            
            mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDxb, mb);
            
            DLambda = -TSparseSpace::Dot(mDxPred, mDxb)/TSparseSpace::Dot(mDxPred, mDxf);
            
            noalias(mDx) = mDxb + DLambda*mDxf;
            
            //Check solution before update
            const double tolerance = 1.0e-10;
            const double max_ratio = 1.0e3;
            if( mNormxEquilibrium > tolerance ) {
                NormDx = TSparseSpace::TwoNorm(mDx);
                
                if( (NormDx/mNormxEquilibrium) > max_ratio || (std::abs(DLambda)/std::abs(mLambda-mDLambdaStep)) > max_ratio )
                {
                    is_converged = false;
                    break;
                }
            }
            
            //update results
            mDLambdaStep += DLambda;
            mLambda += DLambda;
            noalias(mDxStep) += mDx;
            this->Update(rDofSet, mA, mDx, mb);

            //move the mesh if needed
            if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
            
            mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
            
            // *** Check Convergence ***
            
            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(mb);
                mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
            }
            is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        }//While
        
        // Check iteration_number 
        if (iteration_number >= mMaxIterationNumber) {
            is_converged = true;
            //plots a warning if the maximum number of iterations is exceeded
            if(BaseType::GetModelPart().GetCommunicator().MyPID() == 0) {
                this->MaxIterationsExceeded();
            }
        }
        
        //calculate reactions if required
        if (mCalculateReactionsFlag == true) {
            mpBuilderAndSolver->CalculateReactions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
        }
        
		return is_converged;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	void FinalizeSolutionStep() override
	{
		KRATOS_TRY
        
        unsigned int iteration_number = BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER];
        
        // Update the radius
        mRadius = mRadius*sqrt(double(mDesiredIterations)/double(iteration_number));
        
        DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        if (BaseType::GetModelPart().GetProcessInfo()[IS_CONVERGED] == true)
        {
            // Modify the radius to advance faster when convergence is achieved
            if (mRadius > mMaxRadiusFactor*mRadius_0)
                mRadius = mMaxRadiusFactor*mRadius_0;
            else if(mRadius < mMinRadiusFactor*mRadius_0)
                mRadius = mMinRadiusFactor*mRadius_0;
            
            // Update Norm of x
            mNormxEquilibrium = this->CalculateReferenceDofsNorm(rDofSet);
        }
        else
        {
            if (BaseType::mEchoLevel > 0)
               std::cout << "************ NO CONVERGENCE: restoring equilibrium path ************" << std::endl;
            
            TSystemVectorType& mDxStep = *mpDxStep;
            
            //update results
            mLambda -= mDLambdaStep;
            noalias(mDx) = -mDxStep;
            this->Update(rDofSet, mA, mDx, mb);

            //move the mesh if needed
            if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
        }

        BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_LAMBDA] = mLambda;
        BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_RADIUS_FACTOR] = mRadius/mRadius_0;

        mpScheme->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);
        mpBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(), mA, mDx, mb);

        //Cleaning memory after the solution
        mpScheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            this->ClearStep();
        }

		KRATOS_CATCH("")
	}
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Clear() override
    {
        KRATOS_TRY
        
        SparseSpaceType::Clear(mpf);
        SparseSpaceType::Clear(mpDxf);
        SparseSpaceType::Clear(mpDxb);
        SparseSpaceType::Clear(mpDxPred);
        SparseSpaceType::Clear(mpDxStep);
        
        TSystemVectorType& mf = *mpf;
        TSystemVectorType& mDxf = *mpDxf;
        TSystemVectorType& mDxb = *mpDxb;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;
        
        SparseSpaceType::Resize(mf, 0);
        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);
        
        MotherType::Clear();
        
        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void UpdateLoads()
    {
        KRATOS_TRY
                
        mLambda = BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_LAMBDA];
        mRadius = (BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_RADIUS_FACTOR])*mRadius_0;
        
        // Update External Loads        
        this->UpdateExternalLoads();
        
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    Parameters* mpParameters;
    std::vector<ModelPart*> mSubModelPartList; /// List of every SubModelPart associated to an external load
    std::vector<std::string> mVariableNames; /// Name of the nodal variable associated to every SubModelPart

    TSystemVectorPointerType mpf; /// Vector of reference external forces
    TSystemVectorPointerType mpDxf; /// Delta x of A*Dxf=f
    TSystemVectorPointerType mpDxb; /// Delta x of A*Dxb=b
    TSystemVectorPointerType mpDxPred; /// Delta x of prediction phase
    TSystemVectorPointerType mpDxStep; /// Delta x of the current step
    
    unsigned int mDesiredIterations; /// This is used to calculate the radius of the next step
    
    double mMaxRadiusFactor, mMinRadiusFactor; /// Used to limit the radius of the arc length strategy
    double mRadius, mRadius_0; /// Radius of the arc length strategy
    double mLambda, mLambda_old; /// Loading factor
    double mNormxEquilibrium; /// Norm of the solution vector in equilibrium
    double mDLambdaStep; /// Delta lambda of the current step

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check() override
    {
        KRATOS_TRY
        
        int ierr = MotherType::Check();
        if(ierr != 0) return ierr;
        
        KRATOS_ERROR_IF((ARC_LENGTH_LAMBDA.Key() == 0)) << "ARC_LENGTH_LAMBDA Key is 0. Check if all applications were correctly registered." << std::endl;
        KRATOS_ERROR_IF((ARC_LENGTH_RADIUS_FACTOR.Key() == 0)) <<"ARC_LENGTH_RADIUS_FACTOR Key is 0. Check if all applications were correctly registered." << std::endl;
        
        return ierr;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSystemVector(TSystemVectorPointerType& pv)
    {
        if (pv == NULL)
        {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            pv.swap(pNewv);
        }

        TSystemVectorType& v = *pv;

        if (v.size() != mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(mpBuilderAndSolver->GetEquationSystemSize(), false);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SaveInitializeSystemVector(TSystemVectorPointerType& pv)
    {
        if (pv == NULL)
        {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            pv.swap(pNewv);
        }

        TSystemVectorType& v = *pv;

        if (v.size() != mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(mpBuilderAndSolver->GetEquationSystemSize(), true);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void BuildWithDirichlet(TSystemMatrixType& mA, TSystemVectorType& mDx, TSystemVectorType& mb)
    {
        KRATOS_TRY

        mpBuilderAndSolver->Build(mpScheme, BaseType::GetModelPart(), mA, mb);
        mpBuilderAndSolver->ApplyDirichletConditions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    virtual void Update(DofsArrayType& rDofSet, TSystemMatrixType& mA, TSystemVectorType& mDx, TSystemVectorType& mb)
    {
        KRATOS_TRY
        
        // Update scheme
        mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        
        // Update External Loads
        this->UpdateExternalLoads();

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ClearStep()
    {
        KRATOS_TRY

        SparseSpaceType::Clear(mpDxf);
        SparseSpaceType::Clear(mpDxb);
        SparseSpaceType::Clear(mpDxPred);
        SparseSpaceType::Clear(mpDxStep);

        TSystemVectorType& mDxf = *mpDxf;
        TSystemVectorType& mDxb = *mpDxb;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;

        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);

        MotherType::Clear();

        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void UpdateExternalLoads()
    {
        // Update External Loads
        for(unsigned int i = 0; i < mVariableNames.size(); i++)
        {
            ModelPart& rSubModelPart = *(mSubModelPartList[i]);
            const std::string& VariableName = mVariableNames[i];
            
            if( KratosComponents< Variable<double> >::Has( VariableName ) )
            {
                Variable<double> var = KratosComponents< Variable<double> >::Get( VariableName );
                
                #pragma omp parallel
                {
                    ModelPart::NodeIterator NodesBegin;
                    ModelPart::NodeIterator NodesEnd;
                    OpenMPUtils::PartitionedIterators(rSubModelPart.Nodes(),NodesBegin,NodesEnd);
                    
                    for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
                    {
                        double& rvalue = itNode->FastGetSolutionStepValue(var);
                        rvalue *= (mLambda/mLambda_old);
                    }
                }
            }
            else if( KratosComponents< Variable<array_1d<double,3> > >::Has(VariableName) )
            {
                typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
                component_type varx = KratosComponents< component_type >::Get(VariableName+std::string("_X"));
                component_type vary = KratosComponents< component_type >::Get(VariableName+std::string("_Y"));
                component_type varz = KratosComponents< component_type >::Get(VariableName+std::string("_Z"));
                
                #pragma omp parallel
                {
                    ModelPart::NodeIterator NodesBegin;
                    ModelPart::NodeIterator NodesEnd;
                    OpenMPUtils::PartitionedIterators(rSubModelPart.Nodes(),NodesBegin,NodesEnd);
                    
                    for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
                    {
                        double& rvaluex = itNode->FastGetSolutionStepValue(varx);
                        rvaluex *= (mLambda/mLambda_old);
                        double& rvaluey = itNode->FastGetSolutionStepValue(vary);
                        rvaluey *= (mLambda/mLambda_old);
                        double& rvaluez = itNode->FastGetSolutionStepValue(varz);
                        rvaluez *= (mLambda/mLambda_old);
                    }
                }
            }
            else
            {
                KRATOS_THROW_ERROR( std::logic_error, "One variable of the applied loads has a non supported type. Variable: ", VariableName )
            }
        }
        
        // Save the applied Lambda factor
        mLambda_old = mLambda;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double CalculateReferenceDofsNorm(DofsArrayType& rDofSet)
    {
        double ReferenceDofsNorm = 0.0;

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector DofSetPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

        #pragma omp parallel reduction(+:ReferenceDofsNorm)
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofsBegin = rDofSet.begin() + DofSetPartition[k];
            typename DofsArrayType::iterator DofsEnd = rDofSet.begin() + DofSetPartition[k+1];
            
            for (typename DofsArrayType::iterator itDof = DofsBegin; itDof != DofsEnd; ++itDof)
            {                    
                if (itDof->IsFree())
                {
                    const double& temp = itDof->GetSolutionStepValue();
                    ReferenceDofsNorm += temp*temp;
                }
            }
        }
                
        return sqrt(ReferenceDofsNorm);
    }

}; // Class ResidualBasedRammArcLengthStrategy

} // namespace Kratos

#endif // KRATOS_RESIDUALBASED_RAMM_ARC_LENGTH_STRATEGY  defined
