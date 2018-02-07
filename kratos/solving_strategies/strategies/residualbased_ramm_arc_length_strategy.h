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
#include "includes/checks.h"
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
    
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>            StrategyBaseType;
    
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    
    typedef typename BaseType::TBuilderAndSolverType                        TBuilderAndSolverType;

    typedef typename BaseType::TDataType                                                TDataType;

    typedef TSparseSpace                                                          SparseSpaceType;

    typedef typename BaseType::TSchemeType                                            TSchemeType;

    typedef typename BaseType::DofsArrayType                                        DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                                TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                                TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType                        LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType                        LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType                  TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType                  TSystemVectorPointerType;
    
    typedef ModelPart::NodesContainerType                                          NodesArrayType;
    
    typedef ModelPart::ElementsContainerType                                ElementsContainerType;
    
    typedef ModelPart::ConditionsContainerType                                ConditionsArrayType;
    
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
     * @param Parameters The configuration parameters
     * @param MaxIterationNumber The maximum number of iterations
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    ResidualBasedRammArcLengthStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters& Parameters,
        int MaxIterationNumber = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterationNumber, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag),
            mParameters(Parameters)
    {
        // Validate agains defaults -- this also ensures no type mismatch
        Parameters default_parameters = GetDefaultParameters();
        mParameters.ValidateAndAssignDefaults(default_parameters);
        
        // We create the list of submodelparts
        GenerateSubmodelpartList();
    }
    
    /**
     * Default constructor 
     * @param rModelPart The model part of the problem
     * @param pScheme The time integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param Parameters The configuration parameters
     * @param MaxIterationNumber The maximum number of iterations
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    ResidualBasedRammArcLengthStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        Parameters Parameters = Parameters(R"({})"),
        int MaxIterationNumber = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, MaxIterationNumber, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag),
            mParameters(Parameters)
    {
        // Validate agains defaults -- this also ensures no type mismatch
        Parameters default_parameters = GetDefaultParameters();
        mParameters.ValidateAndAssignDefaults(default_parameters);
        
        // We create the list of submodelparts
        GenerateSubmodelpartList();
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
            BaseType::Initialize();
            
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
            typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
            
            // Set up the system
            if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false) {
                // Setting up the list of the DOFs to be solved
                pBuilderAndSolver->SetUpDofSet(pScheme, StrategyBaseType::GetModelPart());

                // Shaping correctly the system
                pBuilderAndSolver->SetUpSystem(StrategyBaseType::GetModelPart());
            }
            
            // Compute initial radius (mRadius0)
            pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, BaseType::mpA, BaseType::mpDx, BaseType::mpb, StrategyBaseType::GetModelPart().Elements(),
                                                            StrategyBaseType::GetModelPart().Conditions(), StrategyBaseType::GetModelPart().GetProcessInfo());
            TSystemMatrixType& mA = *BaseType::mpA;
            TSystemVectorType& mDx = *BaseType::mpDx;
            TSystemVectorType& mb = *BaseType::mpb;
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildAndSolve(pScheme, StrategyBaseType::GetModelPart(), mA, mDx, mb);
            
            mRadius0 = TSparseSpace::TwoNorm(mDx);
            mRadius = mRadius0;

            // Compute vector of reference external force (mf)
            this->InitializeSystemVector(mpf);
            TSystemVectorType& mf = *mpf;
            TSparseSpace::SetToZero(mf);

            pBuilderAndSolver->BuildRHS(pScheme, StrategyBaseType::GetModelPart(), mf);
            
            //Initialize the loading factor Lambda
            mLambda = 0.0;
            mLambdaOld = 1.0;
            
            // Initialize Norm of solution
            mNormxEquilibrium = 0.0;
            
            if (StrategyBaseType::mEchoLevel > 0)
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
            BaseType::InitializeSolutionStep();
            
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
        if (StrategyBaseType::mEchoLevel > 0) 
            std::cout << "ARC-LENGTH RADIUS: " << mRadius/mRadius0 << " X initial radius" << std::endl;
        
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        
        // Initialize variables
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA = *BaseType::mpA;
        TSystemVectorType& mDx = *BaseType::mpDx;
        TSystemVectorType& mb = *BaseType::mpb;
        TSystemVectorType& mf = *mpf;
        TSystemVectorType& mDxb = *mpDxb;
        TSystemVectorType& mDxf = *mpDxf;
        TSystemVectorType& mDxPred = *mpDxPred;
        TSystemVectorType& mDxStep = *mpDxStep;
        
        pScheme->InitializeNonLinIteration(StrategyBaseType::GetModelPart(), mA, mDx, mb);
                
        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mDxf);
        
        // Note: This is not so efficient, but I want to solve mA*mDxf=mf without losing mf
        BuildWithDirichlet(mA, mDxf, mb);
        noalias(mb) = mf;
        pBuilderAndSolver->SystemSolve(mA, mDxf, mb);

        //update results
        double DLambda = mRadius/TSparseSpace::TwoNorm(mDxf);
        mDLambdaStep = DLambda;
        mLambda += DLambda;
        noalias(mDxPred) = DLambda*mDxf;
        noalias(mDxStep) = mDxPred;
        this->Update(rDofSet, mA, mDxPred, mb);

        //move the mesh if needed
        if(StrategyBaseType::MoveMeshFlag() == true)
            StrategyBaseType::MoveMesh();
        
        // ********** Correction phase (iteration cicle) **********

        //initializing the parameters of the iteration loop
        bool is_converged = false;
        mpConvergenceCriteria->InitializeSolutionStep(StrategyBaseType::GetModelPart(), rDofSet, mA, mDxf, mb);
        if (mpConvergenceCriteria->GetActualizeRHSflag() == true) {
            TSparseSpace::SetToZero(mb);
            pBuilderAndSolver->BuildRHS(pScheme, StrategyBaseType::GetModelPart(), mb);
        }
        is_converged = mpConvergenceCriteria->PostCriteria(StrategyBaseType::GetModelPart(), rDofSet, mA, mDxf, mb);

        unsigned int iteration_number = 0;
        double NormDx;
        
        while (is_converged == false && iteration_number < mMaxIterationNumber) {
            // Setting the number of iteration
            iteration_number += 1;
            StrategyBaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            
            pScheme->InitializeNonLinIteration(StrategyBaseType::GetModelPart(), mA, mDx, mb);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDxf);
            
            // Note: This is not so efficient, but I want to solve mA*mDxf=mf without losing mf
            BuildWithDirichlet(mA, mDxf, mb);
            noalias(mb) = mf;
            pBuilderAndSolver->SystemSolve(mA, mDxf, mb);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDxb);
            
            pBuilderAndSolver->BuildAndSolve(pScheme, StrategyBaseType::GetModelPart(), mA, mDxb, mb);
            
            DLambda = -TSparseSpace::Dot(mDxPred, mDxb)/TSparseSpace::Dot(mDxPred, mDxf);
            
            noalias(mDx) = mDxb + DLambda*mDxf;
            
            //Check solution before update
            const double tolerance = 1.0e-10;
            const double max_ratio = 1.0e3;
            if( mNormxEquilibrium > tolerance ) {
                NormDx = TSparseSpace::TwoNorm(mDx);
                
                if( (NormDx/mNormxEquilibrium) > max_ratio || (std::abs(DLambda)/std::abs(mLambda-mDLambdaStep)) > max_ratio ) {
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
            if(StrategyBaseType::MoveMeshFlag() == true) StrategyBaseType::MoveMesh();
            
            pScheme->FinalizeNonLinIteration(StrategyBaseType::GetModelPart(), mA, mDx, mb);
            
            // *** Check Convergence ***
            
            if (mpConvergenceCriteria->GetActualizeRHSflag() == true) {
                TSparseSpace::SetToZero(mb);
                pBuilderAndSolver->BuildRHS(pScheme, StrategyBaseType::GetModelPart(), mb);
            }
            is_converged = mpConvergenceCriteria->PostCriteria(StrategyBaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        }//While
        
        // Check iteration_number 
        if (iteration_number >= mMaxIterationNumber) {
            is_converged = true;
            //plots a warning if the maximum number of iterations is exceeded
            if(StrategyBaseType::GetModelPart().GetCommunicator().MyPID() == 0) {
                this->MaxIterationsExceeded();
            }
        }
        
        // calculate reactions if required
        if (mCalculateReactionsFlag == true)
            pBuilderAndSolver->CalculateReactions(pScheme, StrategyBaseType::GetModelPart(), mA, mDx, mb);
        
        return is_converged;
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) 
     * after solving the solution step.
     */
    
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY
        
        ModelPart& this_model_part = StrategyBaseType::GetModelPart();
        ProcessInfo& this_process_info = this_model_part.GetProcessInfo();
        
        const double iteration_number = static_cast<double>(this_process_info[NL_ITERATION_NUMBER]);
        
        // This is used to calculate the radius of the next step
        const double desired_iterations = static_cast<double>(mParameters["desired_iterations"].GetInt());
        // Update the radius
        mRadius *= std::sqrt(desired_iterations/iteration_number);
        
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA = *BaseType::mpA;
        TSystemVectorType& mDx = *BaseType::mpDx;
        TSystemVectorType& mb = *BaseType::mpb;

        if (this_process_info[IS_CONVERGED] == true) {
            // Used to limit the radius of the arc length strategy
            const double max_radius_factor = mParameters["max_radius_factor"].GetDouble();
            const double min_radius_factor = mParameters["min_radius_factor"].GetDouble();
            // Modify the radius to advance faster when convergence is achieved
            if (mRadius > max_radius_factor * mRadius0)
                mRadius = max_radius_factor * mRadius0;
            else if(mRadius < min_radius_factor * mRadius0)
                mRadius = min_radius_factor * mRadius0;
            
            // Update Norm of x
            mNormxEquilibrium = CalculateReferenceDofsNorm(rDofSet);
        } else {
            if (StrategyBaseType::mEchoLevel > 0)
               std::cout << "************ NO CONVERGENCE: restoring equilibrium path ************" << std::endl;
            
            TSystemVectorType& mDxStep = *mpDxStep;
            
            //update results
            mLambda -= mDLambdaStep;
            noalias(mDx) = -mDxStep;
            this->Update(rDofSet, mA, mDx, mb);

            //move the mesh if needed
            if(StrategyBaseType::MoveMeshFlag() == true) StrategyBaseType::MoveMesh();
        }

        this_process_info[ARC_LENGTH_LAMBDA] = mLambda;
        this_process_info[ARC_LENGTH_RADIUS_FACTOR] = mRadius/mRadius0;

        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        
        pScheme->FinalizeSolutionStep(this_model_part, mA, mDx, mb);
        pBuilderAndSolver->FinalizeSolutionStep(this_model_part, mA, mDx, mb);

        // Cleaning memory after the solution
        pScheme->Clean();

        // Reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
            this->ClearStep();

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
        
        BaseType::Clear();
        
        KRATOS_CATCH( "" )
    }

    /**
    * @brief This method updates the external loads and the lambda factors
    */
    
    virtual void UpdateLoads()
    {
        KRATOS_TRY
        
        const ProcessInfo& this_process_info = StrategyBaseType::GetModelPart().GetProcessInfo();
        mLambda = this_process_info[ARC_LENGTH_LAMBDA];
        mRadius = this_process_info[ARC_LENGTH_RADIUS_FACTOR] * mRadius0;
        
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

    Parameters& mParameters;
    std::vector<ModelPart*> mSubModelPartList; /// List of every SubModelPart associated to an external load
    std::vector<std::string> mVariableNames; /// Name of the nodal variable associated to every SubModelPart

    TSystemVectorPointerType mpf; /// Vector of reference external forces
    TSystemVectorPointerType mpDxf; /// Delta x of A*Dxf=f
    TSystemVectorPointerType mpDxb; /// Delta x of A*Dxb=b
    TSystemVectorPointerType mpDxPred; /// Delta x of prediction phase
    TSystemVectorPointerType mpDxStep; /// Delta x of the current step
    
    double mRadius, mRadius0; /// Radius of the arc length strategy
    double mLambda, mLambdaOld; /// Loading factor
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
    int Check() override
    {
        KRATOS_TRY
        
        int ierr = BaseType::Check();
        if(ierr != 0) return ierr;
        
        KRATOS_CHECK_VARIABLE_KEY(ARC_LENGTH_LAMBDA)
        KRATOS_CHECK_VARIABLE_KEY(ARC_LENGTH_RADIUS_FACTOR)
        
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

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        
        if (v.size() != pBuilderAndSolver->GetEquationSystemSize())
            v.resize(pBuilderAndSolver->GetEquationSystemSize(), false);
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

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        
        if (v.size() != pBuilderAndSolver->GetEquationSystemSize())
            v.resize(pBuilderAndSolver->GetEquationSystemSize(), true);
    }

    /**
    * @brief This method applies the boundary conditions of the system
    * @param mA The LHS of the system
    * @param mDx The increment of solution 
    * @param mb The RHS of the system
    */
    
    void BuildWithDirichlet(
        TSystemMatrixType& mA, 
        TSystemVectorType& mDx, 
        TSystemVectorType& mb
        )
    {
        KRATOS_TRY

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        
        pBuilderAndSolver->Build(pScheme, StrategyBaseType::GetModelPart(), mA, mb);
        pBuilderAndSolver->ApplyDirichletConditions(pScheme, StrategyBaseType::GetModelPart(), mA, mDx, mb);

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
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        pScheme->Update(StrategyBaseType::GetModelPart(), rDofSet, mA, mDx, mb);
        
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
        
        SparseSpaceType::Resize(mDxf, 0);
        SparseSpaceType::Resize(mDxb, 0);
        SparseSpaceType::Resize(mDxPred, 0);
        SparseSpaceType::Resize(mDxStep, 0);

        BaseType::Clear();

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
                        rvalue *= (mLambda/mLambdaOld);
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
                        rvalue *= (mLambda/mLambdaOld);
                    }
                }
            }
            else {
                KARTOS_ERROR << "One variable of the applied loads has a non supported type. Variable: " << VariableName << std::endl;
            }
        }
        
        // Save the applied Lambda factor
        mLambdaOld = mLambda;
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
    
    /**
     * @brief This method is defined in order to reduce the code duplication due to the default parameters
     */
    Parameters GetDefaultParameters()
    {
        // Only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
        {
            "desired_iterations": 4,
            "max_radius_factor": 20.0,
            "min_radius_factor": 0.5,
            "loads_sub_model_part_list": [],
            "loads_variable_list" : []
        }  )" );
        
        return default_parameters;
    }
    
    /**
     * @brief This method is created in order to reduce code duplication for the creation of the submodelpats list
     */
    void GenerateSubmodelpartList()
    {
        // Set Load SubModelParts and Variable names
        if(mParameters["loads_sub_model_part_list"].size() > 0) {
            mSubModelPartList.resize(mParameters["loads_sub_model_part_list"].size());
            mVariableNames.resize(mParameters["loads_variable_list"].size());

            KRATOS_ERROR_IF(mSubModelPartList.size() != mVariableNames.size()) << "For each SubModelPart there must be a corresponding nodal Variable" << std::endl;

            for(unsigned int i = 0; i < mVariableNames.size(); i++) {
                mSubModelPartList[i] = &( rModelPart.GetSubModelPart(mParameters["loads_sub_model_part_list"][i].GetString()) );
                mVariableNames[i] = mParameters["loads_variable_list"][i].GetString();
            }
        } else {
            KRATOS_ERROR << "Not submodelparts defined for apply the arc length" << std::endl;
        }
    }
    
    
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
