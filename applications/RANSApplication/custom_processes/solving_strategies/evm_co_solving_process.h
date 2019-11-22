//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(EVM_CO_SOLVING_PROCESS_H_INCLUDED)
#define EVM_CO_SOLVING_PROCESS_H_INCLUDED

// System includes
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "rans_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// The base class for all EvmCoSolvingProcess in Kratos.
/** The EvmCoSolvingProcess is the base class for all EvmCoSolvingProcess and defines a simple interface for them.
    Execute method is used to execute the EvmCoSolvingProcess algorithms. While the parameters of this method
  can be very different from one EvmCoSolvingProcess to other there is no way to create enough overridden
  versions of it. For this reason this method takes no argument and all EvmCoSolvingProcess parameters must
  be passed at construction time. The reason is that each constructor can take different set of
  argument without any dependency to other EvmCoSolvingProcess or the base EvmCoSolvingProcess class.
*/
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class EvmCoSolvingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = Process;

    using NodeType = ModelPart::NodeType;

    using NodesContainerType = ModelPart::NodesContainerType;

    using SolvingStrategyType = SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Pointer definition of EvmCoSolvingProcess
    KRATOS_CLASS_POINTER_DEFINITION(EvmCoSolvingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EvmCoSolvingProcess(ModelPart& rModelPart, Parameters rParameters)
        : BaseType(), mrModelPart(rModelPart)
    {
        Parameters default_parameters = Parameters(R"(
        {
            "relative_tolerance"                : 1e-3,
            "absolute_tolerance"                : 1e-5,
            "max_iterations"                    : 10,
            "echo_level"                        : 0,
            "relaxation_factor"                 : 1.0,
            "number_of_parent_solve_iterations" : 0,
            "soft_max_exponent"                 : "inverse_epsilon"
        })");

        if (rParameters.Has("soft_max_exponent"))
        {
            if (rParameters["soft_max_exponent"].IsDouble())
            {
                default_parameters.RemoveValue("soft_max_exponent");
                default_parameters.AddEmptyValue("soft_max_exponent");
                default_parameters["soft_max_exponent"].SetDouble(0.0);
            }
        }

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mEchoLevel = rParameters["echo_level"].GetInt();
        mConvergenceRelativeTolerance = rParameters["relative_tolerance"].GetDouble();
        mConvergenceAbsoluteTolerance = rParameters["absolute_tolerance"].GetDouble();
        mRelaxationFactor = rParameters["relaxation_factor"].GetDouble();
        mMaxIterations = rParameters["max_iterations"].GetInt();
        mSkipIterations = rParameters["number_of_parent_solve_iterations"].GetInt();

        if (rParameters["soft_max_exponent"].IsString())
        {
            if (rParameters["soft_max_exponent"].GetString() ==
                "inverse_epsilon")
            {
                mrModelPart.GetProcessInfo()[RANS_SOFT_MAX_EXPONENT] =
                    1.0 / std::numeric_limits<double>::epsilon();
            }
            else
            {
                KRATOS_ERROR << "Undefined \"soft_max_exponent\"=\"" +
                                    rParameters["soft_max_exponent"].GetString() +
                                    "\". Allowed a double value or "
                                    "\"inverse_epsilon\"\n";
            }
        }
        else if (rParameters["soft_max_exponent"].IsDouble())
        {
            mrModelPart.GetProcessInfo()[RANS_SOFT_MAX_EXPONENT] =
                rParameters["soft_max_exponent"].GetDouble();
        }
        else
        {
            KRATOS_ERROR << "\"soft_max_exponent\" can be passed only as a "
                            "double or a string value.";
        }

        mCurrentParentIteration = 0;
    }

    /// Destructor.
    ~EvmCoSolvingProcess() override
    {
    }

    void ExecuteInitialize() override
    {
        for (Process::Pointer auxiliary_process : mAuxiliaryProcessList)
            auxiliary_process->ExecuteInitialize();
    }

    void ExecuteInitializeSolutionStep() override
    {
        UpdateBeforeSolveEquations();
    }

    /// Execute method is used to execute the ScalarCoSolvingProcess algorithms.
    void Execute() override
    {
        if (mrModelPart.GetProcessInfo()[IS_CO_SOLVING_PROCESS_ACTIVE])
            SolveEquations();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override
    {
        int value = BaseType::Check();

        // KRATOS_CHECK(mrModelPart.HasNodalSolutionStepVariable(KINEMATIC_VISCOSITY));
        // KRATOS_CHECK(mrModelPart.HasNodalSolutionStepVariable(TURBULENT_VISCOSITY));

        for (auto strategy : mrSolvingStrategiesList)
            strategy->Check();

        for (Process::Pointer auxiliary_process : mAuxiliaryProcessList)
            auxiliary_process->Check();

        KRATOS_ERROR_IF(mrSolvingStrategiesList.size() == 0)
            << "No strategies are found for ScalarCoSolvingProcess.";

        return value;
    }

    void AddStrategy(typename SolvingStrategyType::Pointer pStrategy,
                     const Variable<double>& rScalarVariable)
    {
        mrSolvingStrategiesList.push_back(pStrategy);
        mrSolvingVariableNamesList.push_back(rScalarVariable.Name());
    }

    void AddAuxiliaryProcess(Process::Pointer pAuxiliaryProcess)
    {
        mAuxiliaryProcessList.push_back(pAuxiliaryProcess);
    }

    void SetParentSolvingStrategy(typename SolvingStrategyType::Pointer pParentSolvingStrategy)
    {
        mpParentSolvingStrategy = pParentSolvingStrategy;
    }

    bool IsConverged() const
    {
        const Communicator& r_communicator = mrModelPart.GetCommunicator();
        const ModelPart::NodesContainerType& r_nodes =
            r_communicator.LocalMesh().Nodes();
        const int number_of_nodes = r_nodes.size();

        double delta_norm{0.0}, solution_norm{0.0};

#pragma omp parallel for reduction(+ : delta_norm, solution_norm)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const NodeType& r_node = *(r_nodes.begin() + i_node);
            const double new_nu_t = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            delta_norm += std::pow(new_nu_t - mrOldValues[i_node], 2);
            solution_norm += std::pow(new_nu_t, 2);
        }

        // This vector stores norms of the residual
        // index - 0 : increase_norm
        // index - 1 : solution_norm
        // index - 3 : number of nodes
        array_1d<double, 3> residual_norms;
        residual_norms[0] = delta_norm;
        residual_norms[1] = solution_norm;
        residual_norms[2] = static_cast<double>(number_of_nodes);
        array_1d<double, 3> total_residual_norms =
            r_communicator.GetDataCommunicator().SumAll(residual_norms);

        total_residual_norms[0] = std::sqrt(total_residual_norms[0]);
        total_residual_norms[1] = std::sqrt(total_residual_norms[1]);

        const double convergence_relative =
            total_residual_norms[0] /
            (total_residual_norms[1] <= std::numeric_limits<double>::epsilon()
                 ? 1.0
                 : total_residual_norms[1]);
        const double convergence_absolute =
            total_residual_norms[0] / total_residual_norms[2];

        return (convergence_relative < this->mConvergenceRelativeTolerance ||
                convergence_absolute < this->mConvergenceAbsoluteTolerance);
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "EvmCoSolvingProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EvmCoSolvingProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    typename SolvingStrategyType::Pointer mpParentSolvingStrategy;

    std::vector<typename SolvingStrategyType::Pointer> mrSolvingStrategiesList;
    std::vector<std::string> mrSolvingVariableNamesList;
    std::vector<Process::Pointer> mAuxiliaryProcessList;

    double mRelaxationFactor;
    double mConvergenceAbsoluteTolerance;
    double mConvergenceRelativeTolerance;

    int mEchoLevel;
    int mMaxIterations;
    int mSkipIterations;
    int mCurrentParentIteration;

    bool mIsCoSolvingProcessActive;

    Vector mrOldValues;
    Vector mrNewValues;
    Vector mrDeltaValues;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteAuxiliaryProcesses()
    {
        for (Process::Pointer auxiliary_process : mAuxiliaryProcessList)
            auxiliary_process->Execute();
    }

    void UpdateBeforeSolveEquations()
    {
        for (Process::Pointer auxiliary_process : mAuxiliaryProcessList)
            auxiliary_process->ExecuteInitializeSolutionStep();
        this->UpdateConvergenceVariable();
    }

    void UpdateConvergenceVariable()
    {
        const double soft_max_exponent =
            this->mrModelPart.GetProcessInfo()[RANS_SOFT_MAX_EXPONENT];
        Communicator& r_communicator = this->mrModelPart.GetCommunicator();

        ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();

        const int number_of_nodes = r_nodes.size();

        ResizeVector(mrOldValues, number_of_nodes);
        ResizeVector(mrNewValues, number_of_nodes);
        ResizeVector(mrDeltaValues, number_of_nodes);

        RansVariableUtilities::GetNodalVariablesVector(mrOldValues, r_nodes, TURBULENT_VISCOSITY);
        this->ExecuteAuxiliaryProcesses();
        RansVariableUtilities::GetNodalVariablesVector(mrNewValues, r_nodes, TURBULENT_VISCOSITY);
        noalias(mrDeltaValues) = mrNewValues - mrOldValues;

// using soft max for the time being as a test
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            const double nu_t =
                mrOldValues[i_node] + mrDeltaValues[i_node] * this->mRelaxationFactor;
            const double soft_max = RansCalculationUtilities::SoftMax(
                nu_t, std::numeric_limits<double>::epsilon(), soft_max_exponent);
            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = soft_max;
        }

        r_communicator.SynchronizeVariable(TURBULENT_VISCOSITY);

        if (this->mEchoLevel > 0)
        {
            const double min_nu_t = RansVariableUtilities::GetMinimumScalarValue(
                this->mrModelPart, TURBULENT_VISCOSITY);
            const double max_nu_t = RansVariableUtilities::GetMaximumScalarValue(
                this->mrModelPart, TURBULENT_VISCOSITY);
            KRATOS_INFO(this->Info())
                << "TURBULENT_VISCOSITY is bounded between [ " << min_nu_t
                << ", " << max_nu_t << " ].\n";
        }
    }

    void UpdateEffectiveViscosity()
    {
        NodesContainerType& r_nodes = this->mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
            const double nu_t = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);

            r_node.FastGetSolutionStepValue(VISCOSITY) = nu_t + nu;
        }
    }

    void ResizeVector(Vector& rVector, const std::size_t Size)
    {
        if (rVector.size() != Size)
        {
            rVector.resize(Size);
        }
    }

    void SolveEquations()
    {
        ++mCurrentParentIteration;

        if (mCurrentParentIteration > mSkipIterations ||
            mpParentSolvingStrategy->IsConverged())
        {
            mCurrentParentIteration = 0;
            this->UpdateBeforeSolveEquations();

            for (auto p_solving_strategy : this->mrSolvingStrategiesList)
            {
                p_solving_strategy->InitializeSolutionStep();
                p_solving_strategy->Predict();
            }

            bool is_converged = false;
            int iteration = 1;

            const ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();

            int iteration_format_length =
                static_cast<int>(std::log10(this->mMaxIterations)) + 1;

            while (!is_converged && iteration <= this->mMaxIterations)
            {
                for (int i = 0;
                     i < static_cast<int>(this->mrSolvingStrategiesList.size()); ++i)
                {
                    auto p_solving_strategy = this->mrSolvingStrategiesList[i];
                    auto scalar_variable_name = this->mrSolvingVariableNamesList[i];

                    p_solving_strategy->SolveSolutionStep();
                    const unsigned int iterations =
                        r_current_process_info[NL_ITERATION_NUMBER];
                    KRATOS_INFO_IF(this->Info(), this->mEchoLevel > 0)
                        << "Solving " << scalar_variable_name << " used "
                        << iterations << " iterations.\n";
                }

                this->UpdateConvergenceVariable();

                is_converged = this->IsConverged();

                if (this->mEchoLevel > 1 && is_converged)
                {
                    std::stringstream conv_msg;
                    conv_msg << "[Itr.#" << std::setw(iteration_format_length)
                             << iteration << "/" << this->mMaxIterations
                             << "] CONVERGENCE CHECK: TURBULENT_VISCOSITY"
                             << " *** CONVERGENCE IS ACHIEVED ***\n";
                    KRATOS_INFO(this->Info()) << conv_msg.str();
                }

                iteration++;
            }

            this->UpdateEffectiveViscosity();

            KRATOS_INFO_IF(this->Info(), !is_converged && this->mEchoLevel > 2)
                << "\n-------------------------------------------------------"
                << "\n    INFO: Max coupling iterations reached.             "
                << "\n          Please increase coupling max_iterations      "
                << "\n          or decrease coupling                         "
                << "\n          relative_tolerance/absolute tolerance        "
                << "\n-------------------------------------------------------"
                << "\n";

            for (auto p_solving_strategy : this->mrSolvingStrategiesList)
                p_solving_strategy->FinalizeSolutionStep();
        }
        else
        {
            KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
                << "Skipping co-solving process for parent solve to continue, "
                   "since parent solve itertions are less than "
                   "\"number_of_parent_solve_iterations\" [ "
                << mCurrentParentIteration << " <= " << mSkipIterations << " ].\n";
        }
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    EvmCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(
        EvmCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    // EvmCoSolvingProcess(EvmCoSolvingProcess const& rOther);

    ///@}

}; // Class EvmCoSolvingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::istream& operator>>(std::istream& rIStream,
                                EvmCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// output stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::ostream& operator<<(std::ostream& rOStream,
                                const EvmCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // EVM_CO_SOLVING_PROCESS_H_INCLUDED defined
