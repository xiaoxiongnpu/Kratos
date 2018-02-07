//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Vicente Mataix
//

#if !defined(KRATOS_RESIDUALBASED_RAMM_ARC_LENGTH_STRATEGY)
#define KRATOS_RESIDUALBASED_RAMM_ARC_LENGTH_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
/** 
 * @class ResidualBasedRammArcLengthStrategy
 * @ingroup KratosCore
 * @brief This is a Arc-Length strategy based on Ramm algorithm
 * @details The Arc-Length method is a very efficient method in solving non-linear systems of equations when the problem under consideration exhibits one or more critical points. In terms of a simple mechanical loading-unloading problem, a critical point could be interpreted as the point at which the loaded body cannot support an increase of the external forces and an instability occurs.
 * @author Ignasi de Pouplana
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedRammArcLengthStrategy 
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

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
    
//     typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>            StrategyBaseType;
//     
//     typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
//     
//     typedef typename BaseType::TBuilderAndSolverType                        TBuilderAndSolverType;
// 
//     typedef typename BaseType::TDataType                                                TDataType;
// 
//     typedef TSparseSpace                                                          SparseSpaceType;
// 
//     typedef typename BaseType::TSchemeType                                            TSchemeType;
// 
//     typedef typename BaseType::DofsArrayType                                        DofsArrayType;
// 
//     typedef typename BaseType::TSystemMatrixType                                TSystemMatrixType;
// 
//     typedef typename BaseType::TSystemVectorType                                TSystemVectorType;
// 
//     typedef typename BaseType::LocalSystemVectorType                        LocalSystemVectorType;
// 
//     typedef typename BaseType::LocalSystemMatrixType                        LocalSystemMatrixType;
// 
//     typedef typename BaseType::TSystemMatrixPointerType                  TSystemMatrixPointerType;
//     
//     typedef typename BaseType::TSystemVectorPointerType                  TSystemVectorPointerType;
//     
//     typedef ModelPart::NodesContainerType                                          NodesArrayType;
//     
//     typedef ModelPart::ConditionsContainerType                                ConditionsArrayType;
    
    /// Pointer definition of TreeContactSearch
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedRammArcLengthStrategy);

    ///@}
    ///@name Life Cycle
    ///@{
    
    /**
     * Default constructor 
     * @param rModelPart The model part of the problem
     * @param pScheme The time integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param rParameters The configuration parameters
     * @param MaxIterationNumber The maximum number of iterations
     * @param CalculateReactions: The flag for the reaction calculation
     * @param ReformDofSetAtEachStep: The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag: The flag that allows to move the mesh
     */
    ResidualBasedRammArcLengthStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters& rParameters,
        int MaxIterationNumber = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterationNumber, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
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
    
    /**
     * Default constructor 
     * @param rModelPart The model part of the problem
     * @param pScheme The time integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param rParameters The configuration parameters
     * @param MaxIterationNumber The maximum number of iterations
     * @param CalculateReactions: The flag for the reaction calculation
     * @param ReformDofSetAtEachStep: The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag: The flag that allows to move the mesh
     */
    ResidualBasedRammArcLengthStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        Parameters& rParameters,
        int MaxIterationNumber = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, MaxIterationNumber, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {
//             //only include validation with c++11 since raw_literals do not exist in c++03
//             Parameters default_parameters( R"(
//             {
//                 "desired_iterations": 4,
//                 "max_radius_factor": 20.0,
//                 "min_radius_factor": 0.5,
//                 "loads_sub_model_part_list": [],
//                 "loads_variable_list" : []
//             }  )" );
//             
//             // Validate agains defaults -- this also ensures no type mismatch
//             rParameters.ValidateAndAssignDefaults(default_parameters);
//             
//             mpParameters = &rParameters;
//             
//             // Set Load SubModelParts and Variable names
//             if(rParameters["loads_sub_model_part_list"].size() > 0) {
//                 mSubModelPartList.resize(rParameters["loads_sub_model_part_list"].size());
//                 mVariableNames.resize(rParameters["loads_variable_list"].size());
// 
//                 KRATOS_ERROR_IF(mSubModelPartList.size() != mVariableNames.size()) << "For each SubModelPart there must be a corresponding nodal Variable" << std::endl;
// 
//                 for(unsigned int i = 0; i < mVariableNames.size(); i++) {
//                     mSubModelPartList[i] = &( rModelPart.GetSubModelPart(rParameters["loads_sub_model_part_list"][i].GetString()) );
//                     mVariableNames[i] = rParameters["loads_variable_list"][i].GetString();
//                 }
//             }
// 
//             mDesiredIterations = rParameters["desired_iterations"].GetInt();
//             mMaxRadiusFactor = rParameters["max_radius_factor"].GetDouble();
//             mMinRadiusFactor = rParameters["min_radius_factor"].GetDouble();
    }

    ///Destructor
    ~ResidualBasedRammArcLengthStrategy() {} override
    = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * @brief Initialization of member variables and prior operations
     */
        
    void Initialize() override
    {
        KRATOS_TRY

        if (mInitializeWasPerformed == false) {
            MotherType::Initialize();
            
            //set up the system
            if (mpBuilderAndSolver->GetDofSetIsInitializedFlag() == false) {
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

    /**
     * @brief Performs all the required operations that should be done (for each step) 
     * before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    
    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        if (mSolutionStepIsInitialized == false) {
            MotherType::InitializeSolutionStep();
            
            this->SaveInitializeSystemVector(mpf);
            this->InitializeSystemVector(mpDxf);
            this->InitializeSystemVector(mpDxb);
            this->InitializeSystemVector(mpDxPred);
            this->InitializeSystemVector(mpDxStep);
        }
        
        KRATOS_CATCH( "" )
    }

    /**
     * @brief Solves the current step. 
     * @details This function returns true if a solution has been found, false otherwise.
     */
        
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

    /**
     * @brief Performs all the required operations that should be done (for each step) 
     * after solving the solution step.
     */
    
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
    
    /**
     * @brief Clears the internal storage
     */
    
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

    /**
    * @brief This method updates the external loads and the lambda factors
    */
    
    virtual void UpdateLoads()
    {
        KRATOS_TRY
                
        mLambda = BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_LAMBDA];
        mRadius = (BaseType::GetModelPart().GetProcessInfo()[ARC_LENGTH_RADIUS_FACTOR])*mRadius_0;
        
        // Update External Loads        
        this->UpdateExternalLoads();
        
        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /************************************ GET INFO *************************************/
    /***********************************************************************************/
    
    virtual std::string Info() const
    {
        return "ResidualBasedRammArcLengthStrategy";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/
    
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
    
protected:
    
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

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

    ///@}
    ///@name Protected Operators
    ///@{
    
    ///@}
    ///@name Protected Operations
    ///@{
    
    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
     the
    int Check() override
    {
        KRATOS_TRY
        
        int ierr = MotherType::Check();
        if(ierr != 0) return ierr;
        
        // TODO: replace by the macros from checks.h
        KRATOS_ERROR_IF((ARC_LENGTH_LAMBDA.Key() == 0)) << "ARC_LENGTH_LAMBDA Key is 0. Check if all applications were correctly registered." << std::endl;
        KRATOS_ERROR_IF((ARC_LENGTH_RADIUS_FACTOR.Key() == 0)) <<"ARC_LENGTH_RADIUS_FACTOR Key is 0. Check if all applications were correctly registered." << std::endl;
        
        return ierr;

        KRATOS_CATCH( "" )
    }

    /**
     * @brief This method initializes the system 
     * @param rSystemVectorPointer The pointer that contains the system of equations
     */
    
    void InitializeSystemVector(TSystemVectorPointerType& rSystemVectorPointer)
    {
        if (rSystemVectorPointer == nullptr) {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            rSystemVectorPointer.swap(pNewv);
        }

        TSystemVectorType& v = *rSystemVectorPointer;

        if (v.size() != mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(mpBuilderAndSolver->GetEquationSystemSize(), false);
    }

    /**
     * @brief This method saves the initialized system 
     * @param rSystemVectorPointer The pointer that contains the system of equations
     */
    
    void SaveInitializeSystemVector(TSystemVectorPointerType& rSystemVectorPointer)
    {
        if (rSystemVectorPointer == nullptr) {
            TSystemVectorPointerType pNewv = TSystemVectorPointerType(new TSystemVectorType(0));
            rSystemVectorPointer.swap(pNewv);
        }

        TSystemVectorType& v = *rSystemVectorPointer;

        if (v.size() != mpBuilderAndSolver->GetEquationSystemSize())
            v.resize(mpBuilderAndSolver->GetEquationSystemSize(), true);
    }

    /**
    * @brief This method applies the boundary conditions of the system
    * @param mA The LHS of the system
    * @param mDx The increment of solution 
    * @param mb The RHS of the system
    */
    
    void BuildWithDirichlet( the
        TSystemMatrixType& mA, 
        TSystemVectorType& mDx, 
        TSystemVectorType& mb
        )
    {
        KRATOS_TRY

        mpBuilderAndSolver->Build(mpScheme, BaseType::GetModelPart(), mA, mb);
        mpBuilderAndSolver->ApplyDirichletConditions(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);

        KRATOS_CATCH( "" )
    }

    /**
     * @brief This method calll the update from the scheme and updates the external loads
     * @param rDofSet The set of degrees of freedom to consider
     * @param mA The LHS of the system
     * @param mDx The increment of solution 
     * @param mb The RHS of the system
     */
    
    virtual void Update(
        DofsArrayType& rDofSet, 
        TSystemMatrixType& mA, 
        TSystemVectorType& mDx, 
        TSystemVectorType& mb
        )
    {
        KRATOS_TRY
        
        // Update scheme
        mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        
        // Update External Loads
        this->UpdateExternalLoads();

        KRATOS_CATCH( "" )
    }

    /**
     * @brief This method clears partially the solution step
     */
    
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
 the
        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);

        MotherType::Clear();

        KRATOS_CATCH("");
    }

    /**
    * @brief This method updates the external loads multipliying by the lambda coefficient
    */
    
    void UpdateExternalLoads()
    {
        // Update External Loads
        for(unsigned int i = 0; i < mVariableNames.size(); i++) {
            ModelPart& rSubModelPart = *(mSubModelPartList[i]);
            const std::string& VariableName = mVariableNames[i];
            
            if( KratosComponents< Variable<double> >::Has( VariableName ) ) {
                Variable<double> var = KratosComponents< Variable<double> >::Get( VariableName );
                
                #pragma omp parallel
                {
                    ModelPart::NodeIterator NodesBegin;
                    ModelPart::NodeIterator NodesEnd;
                    OpenMPUtils::PartitionedIterators(rSubModelPart.Nodes(),NodesBegin,NodesEnd);
                    
                    for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode) {
                        double& rvalue = itNode->FastGetSolutionStepValue(var);
                        rvalue *= (mLambda/mLambda_old);
                    }
                }
            }
            else if( KratosComponents< Variable<array_1d<double,3> > >::Has(VariableName) ) {
                Variable<array_1d<double,3>> var = KratosComponents< Variable<array_1d<double,3>> >::Get( VariableName );                
                
                #pragma omp parallel
                {
                    ModelPart::NodeIterator NodesBegin;
                    ModelPart::NodeIterator NodesEnd;
                    OpenMPUtils::PartitionedIterators(rSubModelPart.Nodes(),NodesBegin,NodesEnd);
                    
                    for (ModelPart::NodeIterator it_node = NodesBegin; it_node != NodesEnd; ++it_node) {
                        array_1d<double, 3>& rvalue = it_node->FastGetSolutionStepValue(var);
                        rvalue *= (mLambda/mLambda_old);
                    }
                }
            }
            else {
                KARTOS_ERROR << "One variable of the applied loads has a non supported type. Variable: " << VariableName << std::endl;
            }
        }
        
        // Save the applied Lambda factor
        mLambda_old = mLambda;
    }

    /**
    * @brief This method computes the norm of reference
    * @param rDofSet The set of degrees of freedom to compute
    * @return The reference norm of the DoF studied 
    */
    
    double CalculateReferenceDofsNorm(DofsArrayType& rDofSet)
    {
        double reference_dofs_norm = 0.0;

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector DofSetPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

        #pragma omp parallel reduction(+:reference_dofs_norm)
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofsBegin = rDofSet.begin() + DofSetPartition[k];
            typename DofsArrayType::iterator DofsEnd = rDofSet.begin() + DofSetPartition[k+1];
            
            for (typename DofsArrayType::iterator it_dof = DofsBegin; it_dof != DofsEnd; ++it_dof) {                    
                if (it_dof->IsFree()) {
                    const double& temp = it_dof->GetSolutionStepValue();
                    reference_dofs_norm += temp*temp;
                }
            }
        }
                
        return std::sqrt(reference_dofs_norm);
    }
    
    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class ResidualBasedRammArcLengthStrategy
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

// /****************************** INPUT STREAM FUNCTION ******************************/
// /***********************************************************************************/
// 
// template<class TPointType, class TPointerType>
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ResidualBasedRammArcLengthStrategy& rThis);
// 
// /***************************** OUTPUT STREAM FUNCTION ******************************/
// /***********************************************************************************/
// 
// template<class TPointType, class TPointerType>
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ResidualBasedRammArcLengthStrategy& rThis)
// {
//     return rOStream;
// }

} // namespace Kratos

#endif // KRATOS_RESIDUALBASED_RAMM_ARC_LENGTH_STRATEGY  defined
