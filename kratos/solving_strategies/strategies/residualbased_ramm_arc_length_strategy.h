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
//  Refactor code:   Alejandro Cornejo Velazquez
//                   Vicente Mataix Ferrandiz
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

    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace>               TConvergenceCriteriaType;

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
     * @param ThisParameters The configuration parameters
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
        Parameters ThisParameters,
        int MaxIterationNumber = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterationNumber, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag),
            mParameters(ThisParameters)
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
     * @param ThisParameters The configuration parameters
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
        Parameters ThisParameters,
        int MaxIterationNumber = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, MaxIterationNumber, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag),
            mParameters(ThisParameters)
    {
        // Validate agains defaults -- this also ensures no type mismatch
        Parameters default_parameters = GetDefaultParameters();
        mParameters.ValidateAndAssignDefaults(default_parameters);

        // We create the list of submodelparts
        GenerateSubmodelpartList();
    }

    ///Destructor
    ~ResidualBasedRammArcLengthStrategy() override
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

        if (BaseType::mInitializeWasPerformed == false) {
            BaseType::Initialize();

            ModelPart& r_model_part = StrategyBaseType::GetModelPart();
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
            typename TSchemeType::Pointer pScheme = BaseType::GetScheme();

            // Set up the system
            if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false) {
                // Setting up the list of the DOFs to be solved
                pBuilderAndSolver->SetUpDofSet(pScheme, r_model_part);

                // Shaping correctly the system
                pBuilderAndSolver->SetUpSystem(r_model_part);
            }

            // Compute initial radius (mRadius0)
            pBuilderAndSolver->ResizeAndInitializeVectors(pScheme,
                                                          BaseType::mpA,
                                                          BaseType::mpDx,
                                                          BaseType::mpb,
                                                          r_model_part
                                                          );
            TSystemMatrixType& mA = *BaseType::mpA;
            TSystemVectorType& mDx = *BaseType::mpDx;
            TSystemVectorType& mb = *BaseType::mpb;
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);

            pBuilderAndSolver->BuildAndSolve(pScheme, r_model_part, mA, mDx, mb);

            mRadius0 = TSparseSpace::TwoNorm(mDx);
            mRadius = mRadius0;

            // Compute vector of reference external force (mf)
            InitializeSystemVector(mpf);
            TSystemVectorType& mf = *mpf;
            TSparseSpace::SetToZero(mf);

            ModelPart& r_external_forces_model_part = r_model_part.GetSubModelPart("EXTERNAL_FORCES_SUBMODELPART");
            
//             pBuilderAndSolver->BuildRHS(pScheme, StrategyBaseType::GetModelPart(), mf);
            pBuilderAndSolver->BuildRHS(pScheme, r_external_forces_model_part, mf);

            // Initialize the loading factor Lambda
            mLambda = 0.0;
            mLambdaOld = 1.0;

            // Initialize Norm of solution
            mNormxEquilibrium = 0.0;

            if (StrategyBaseType::mEchoLevel > 0)
              KRATOS_INFO("Ramm's Arc Length Strategy ") << "Initialized" << std::endl;
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

        if (BaseType::mSolutionStepIsInitialized == false) {
            BaseType::InitializeSolutionStep();

            SaveInitializeSystemVector(mpf);
            InitializeSystemVector(mpDxf);
            InitializeSystemVector(mpDxb);
            InitializeSystemVector(mpDxPred);
            InitializeSystemVector(mpDxStep);
        }

        KRATOS_CATCH( "" )
    }

    /**
     * @brief Solves the current step.
     * @details This function returns true if a solution has been found, false otherwise.
     */

    bool SolveSolutionStep() override
    {
        /* Prediction phase */
        if (StrategyBaseType::mEchoLevel > 0) // TODO: Replace for the new Kratos log when available
            KRATOS_INFO("ARC-LENGTH RADIUS: ") << mRadius/mRadius0 << " X initial radius" << std::endl;

        ModelPart& r_model_part = StrategyBaseType::GetModelPart();

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

        pScheme->InitializeNonLinIteration(r_model_part, mA, mDx, mb);

        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mDxf);

        // We solve mA*mDxf=mf
        BuildWithDirichlet(mA, mDxf, mb);
//         TSparseSpace::Assign(mb, 1.0, mf); // TODO: Remove this
        pBuilderAndSolver->SystemSolve(mA, mDxf, mf);

        // Update results
        double delta_lambda = mRadius/TSparseSpace::TwoNorm(mDxf);
        mDeltaLamdaStep = delta_lambda;
        mLambda += delta_lambda;
        TSparseSpace::Assign(mDxPred, delta_lambda, mDxf);
        TSparseSpace::Assign(mDxStep, 1.0, mDxPred);

        UpdateDatabase(mA, mDxPred, mb, StrategyBaseType::MoveMeshFlag());

        /* Correction phase (iteration cicle) */

        // Initializing the parameters of the iteration loop
        bool is_converged = false;
        BaseType::mpConvergenceCriteria->InitializeSolutionStep(r_model_part, rDofSet, mA, mDxf, mb);
        if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true) {
            TSparseSpace::SetToZero(mb);
            pBuilderAndSolver->BuildRHS(pScheme, r_model_part, mb);
        }
        is_converged = BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, rDofSet, mA, mDxf, mb);

        unsigned int iteration_number = 0;
        double norm_dx;

        while (is_converged == false && iteration_number < BaseType::mMaxIterationNumber) {
            // Setting the number of iteration
            iteration_number += 1;
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            pScheme->InitializeNonLinIteration(r_model_part, mA, mDx, mb);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDxf);

            // We solve mA*mDxf=mf
            BuildWithDirichlet(mA, mDxf, mb);
//             TSparseSpace::Assign(mb, 1.0, mf); // TODO: Remove this
            pBuilderAndSolver->SystemSolve(mA, mDxf, mf);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDxb);

            pBuilderAndSolver->BuildAndSolve(pScheme, r_model_part, mA, mDxb, mb);

            delta_lambda = -TSparseSpace::Dot(mDxPred, mDxb)/TSparseSpace::Dot(mDxPred, mDxf);

            if (std::abs(delta_lambda) < std::numeric_limits<double>::epsilon())
                KRATOS_WARNING("Zero Delta Lambda: ") << delta_lambda << std::endl;

            // Doing mDx = mDxb + delta_lambda * mDxf using spaces
            TSparseSpace::Assign(mDx, delta_lambda, mDxf);
            TSparseSpace::UnaliasedAdd(mDx, 1.0, mDxb);

            // Check solution before update
            const double tolerance = 1.0e-10;
            const double max_ratio = 1.0e3;
            if( mNormxEquilibrium > tolerance ) {
                norm_dx = TSparseSpace::TwoNorm(mDx);

                if( (norm_dx/mNormxEquilibrium) > max_ratio ||
                    (std::abs(delta_lambda)/std::abs(mLambda - mDeltaLamdaStep)) > max_ratio ) {
                    is_converged = false;
                    break;
                }
            }

            // Update results
            mDeltaLamdaStep += delta_lambda;
            mLambda += delta_lambda;
            TSparseSpace::UnaliasedAdd(mDxStep, 1.0, mDx);
            UpdateDatabase(mA, mDx, mb, StrategyBaseType::MoveMeshFlag());

            pScheme->FinalizeNonLinIteration(r_model_part, mA, mDx, mb);

            // Check convergence
            if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true) {
                TSparseSpace::SetToZero(mb);
                pBuilderAndSolver->BuildRHS(pScheme, r_model_part, mb);
            }
            is_converged = BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, rDofSet, mA, mDx, mb);
        }//While

        // Check iteration_number
        if (iteration_number >= BaseType::mMaxIterationNumber) {
            is_converged = true;
            // Plots a warning if the maximum number of iterations is exceeded
            if(r_model_part.GetCommunicator().MyPID() == 0)
                BaseType::MaxIterationsExceeded();
        }

        // Calculate reactions if required
        if (BaseType::mCalculateReactionsFlag == true)
            pBuilderAndSolver->CalculateReactions(pScheme, r_model_part, mA, mDx, mb);

        return is_converged;
    }

    /**
     * @brief Performs all the required operations that should be done (for each step)
     * after solving the solution step.
     */

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY

        ModelPart& r_model_part = StrategyBaseType::GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        const double iteration_number = static_cast<double>(r_process_info[NL_ITERATION_NUMBER]);

        // This is used to calculate the radius of the next step
        const double desired_iterations = static_cast<double>(mParameters["desired_iterations"].GetInt());
        // Update the radius
        mRadius *= std::sqrt(desired_iterations/iteration_number);

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA = *BaseType::mpA;
        TSystemVectorType& mDx = *BaseType::mpDx;
        TSystemVectorType& mb = *BaseType::mpb;

        if (r_process_info[IS_CONVERGED] == true) {
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
               KRATOS_WARNING("NO CONVERGENCE ACHIEVED: ") << "Restoring equilibrium path" << std::endl;

            TSystemVectorType& mDxStep = *mpDxStep;

            // Update results
            mLambda -= mDeltaLamdaStep;
            TSparseSpace::Assign(mDx, -1.0, mDxStep);
            UpdateDatabase(mA, mDx, mb, StrategyBaseType::MoveMeshFlag());
        }

        if (StrategyBaseType::mEchoLevel > 0)
            KRATOS_INFO("Arc length param:") << "\tARC_LENGTH_LAMBDA " << mLambda << "\tARC_LENGTH_RADIUS_FACTOR " << mRadius/mRadius0 << std::endl;

        r_process_info[ARC_LENGTH_LAMBDA] = mLambda;
        r_process_info[ARC_LENGTH_RADIUS_FACTOR] = mRadius/mRadius0;

        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();

        pScheme->FinalizeSolutionStep(r_model_part, mA, mDx, mb);
        pBuilderAndSolver->FinalizeSolutionStep(r_model_part, mA, mDx, mb);

        // Cleaning memory after the solution
        pScheme->Clean();

        // Reset flags for next step
        BaseType::mSolutionStepIsInitialized = false;

        if (BaseType::mReformDofSetAtEachStep == true) //deallocate the systemvectors
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

        const ProcessInfo& r_process_info = StrategyBaseType::GetModelPart().GetProcessInfo();
        mLambda = r_process_info[ARC_LENGTH_LAMBDA];
        mRadius = r_process_info[ARC_LENGTH_RADIUS_FACTOR] * mRadius0;

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

    std::string Info() const override
    {
        return "ResidualBasedRammArcLengthStrategy";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/

    void PrintInfo(std::ostream& rOStream) const override
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

    Parameters mParameters;                    /// The configuration parameters
    std::vector<ModelPart*> mSubModelPartList; /// List of every SubModelPart associated to an external load
    std::vector<std::string> mVariableNames;   /// Name of the nodal variable associated to every SubModelPart

    TSystemVectorPointerType mpf;      /// Vector of reference external forces
    TSystemVectorPointerType mpDxf;    /// Delta x of A*Dxf=f
    TSystemVectorPointerType mpDxb;    /// Delta x of A*Dxb=b
    TSystemVectorPointerType mpDxPred; /// Delta x of prediction phase
    TSystemVectorPointerType mpDxStep; /// Delta x of the current step

    double mRadius, mRadius0;   /// Radius of the arc length strategy
    double mLambda, mLambdaOld; /// Loading factor
    double mNormxEquilibrium;   /// Norm of the solution vector in equilibrium
    double mDeltaLamdaStep;     /// Delta lambda of the current step

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
            TSystemVectorPointerType pnew_v = TSystemVectorPointerType(new TSystemVectorType(0));
            rSystemVectorPointer.swap(pnew_v);
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
            TSystemVectorPointerType pnew_v = TSystemVectorPointerType(new TSystemVectorType(0));
            rSystemVectorPointer.swap(pnew_v);
        }

        TSystemVectorType& v = *rSystemVectorPointer;

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();

        if (v.size() != pBuilderAndSolver->GetEquationSystemSize())
            v.resize(pBuilderAndSolver->GetEquationSystemSize(), true);
    }

    /**
    * @brief This method applies the boundary conditions of the system
    * @param rA The LHS of the system
    * @param rDx The increment of solution
    * @param rb The RHS of the system
    */

    void BuildWithDirichlet(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
        KRATOS_TRY

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();

        pBuilderAndSolver->Build(pScheme, StrategyBaseType::GetModelPart(), rA, rb);
        pBuilderAndSolver->ApplyDirichletConditions(pScheme, StrategyBaseType::GetModelPart(), rA, rDx, rb);

        KRATOS_CATCH( "" )
    }

    /**
     * @brief This method calll the update from the scheme and updates the external loads
     * @param rA The LHS of the system
     * @param rDx The increment of solution
     * @param rb The RHS of the system
     * @param MoveMesh The flag that indicates if the mesh is moved or not
     */

    void UpdateDatabase(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        const bool MoveMesh
        ) override
    {
        KRATOS_TRY

        // We call the base class update database
        BaseType::UpdateDatabase(rA, rDx, rb, MoveMesh);

        // Update External Loads
        UpdateExternalLoads();

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
        // The ratio between the current and previous lambda
        const double lambda_ratio = mLambda/mLambdaOld;

        // Update External Loads
        for(unsigned int i = 0; i < mVariableNames.size(); ++i) {
            ModelPart& r_sub_model_part = *(mSubModelPartList[i]);
            const std::string& r_variable_name = mVariableNames[i];

            auto& nodes_array = r_sub_model_part.Nodes();
            auto& conditions_array = r_sub_model_part.Conditions();

            if( KratosComponents< Variable<double> >::Has( r_variable_name ) ) {
                const Variable<double>& r_variable = KratosComponents<Variable<double>>::Get(r_variable_name);

                // Nodes
                #pragma omp parallel for
                for (int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
                    auto it_node = nodes_array.begin() + i;
                    double& rvalue = it_node->FastGetSolutionStepValue(r_variable);
                    rvalue *= lambda_ratio;
                }
                // Conditions
                #pragma omp parallel for
                for (int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
                    auto it_cond = conditions_array.begin() + i;
                    if (it_cond->Has(r_variable)) {
                        double& rvalue = it_cond->GetValue(r_variable);
                        rvalue *= lambda_ratio;
                    }
                }
            } else if( KratosComponents< Variable<array_1d<double,3> > >::Has(r_variable_name) ) {
                const Variable<array_1d<double,3>>& r_variable = KratosComponents<Variable<array_1d<double,3>>>::Get(r_variable_name);

                // Nodes
                #pragma omp parallel for
                for (int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
                    auto it_node = nodes_array.begin() + i;
                    array_1d<double, 3>& rvalue = it_node->FastGetSolutionStepValue(r_variable);
                    rvalue *= lambda_ratio;
                }
                // Conditions
                #pragma omp parallel for
                for (int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
                    auto it_cond = conditions_array.begin() + i;
                    if (it_cond->Has(r_variable)) {
                        array_1d<double, 3>& rvalue = it_cond->GetValue(r_variable);
                        rvalue *= lambda_ratio;
                    }
                }
            } else {
                KRATOS_ERROR << "One variable of the applied loads has a non supported type. Variable: " << r_variable_name << std::endl;
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

        #pragma omp parallel for reduction(+:reference_dofs_norm)
        for (int i = 0; i < static_cast<int>(rDofSet.size()); ++i) {
            auto it_dof = rDofSet.begin() + i;
            if (it_dof->IsFree()) {
                const double temp = it_dof->GetSolutionStepValue();
                reference_dofs_norm += std::pow(temp, 2);
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
        ModelPart& r_model_part = BaseType::GetModelPart();
        ModelPart& r_external_forces_model_part = r_model_part.HasSubModelPart("EXTERNAL_FORCES_SUBMODELPART") ?  r_model_part.GetSubModelPart("EXTERNAL_FORCES_SUBMODELPART") : r_model_part.CreateSubModelPart("EXTERNAL_FORCES_SUBMODELPART");

        // Set Load SubModelParts and Variable names
        if(mParameters["loads_sub_model_part_list"].size() > 0) {
            mSubModelPartList.resize(mParameters["loads_sub_model_part_list"].size());
            mVariableNames.resize(mParameters["loads_variable_list"].size());

            KRATOS_ERROR_IF(mSubModelPartList.size() != mVariableNames.size()) << "For each SubModelPart there must be a corresponding nodal or conditional Variable" << std::endl;

            for(unsigned int i = 0; i < mVariableNames.size(); ++i) {
                ModelPart& sub_model_part = r_model_part.GetSubModelPart(mParameters["loads_sub_model_part_list"][i].GetString());
                mSubModelPartList[i] = &sub_model_part;
                mVariableNames[i] = mParameters["loads_variable_list"][i].GetString();

                // We add the conditions to our external forces model part
                auto& conditions_array = sub_model_part.Conditions();
                r_external_forces_model_part.AddConditions(conditions_array.begin(), conditions_array.end());
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
