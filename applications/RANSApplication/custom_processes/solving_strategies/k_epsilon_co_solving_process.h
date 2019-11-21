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

#if !defined(K_EPSILON_CO_SOLVING_PROCESS_H_INCLUDED)
#define K_EPSILON_CO_SOLVING_PROCESS_H_INCLUDED

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

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "processes/find_nodal_neighbours_process.h"
#include "rans_application_variables.h"
#include "scalar_co_solving_process.h"

// debugging
#include "input_output/vtk_output.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// The base class for all KEpsilonCoSolvingProcess in Kratos.
/** The KEpsilonCoSolvingProcess is the base class for all KEpsilonCoSolvingProcess and defines a simple interface for them.
    Execute method is used to execute the KEpsilonCoSolvingProcess algorithms. While the parameters of this method
  can be very different from one KEpsilonCoSolvingProcess to other there is no way to create enough overridden
  versions of it. For this reason this method takes no argument and all KEpsilonCoSolvingProcess parameters must
  be passed at construction time. The reason is that each constructor can take different set of
  argument without any dependency to other KEpsilonCoSolvingProcess or the base KEpsilonCoSolvingProcess class.
*/
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class KEpsilonCoSolvingProcess
    : public ScalarCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = ScalarCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>;

    using NodeType = ModelPart::NodeType;

    using NodesContainerType = ModelPart::NodesContainerType;

    /// Pointer definition of KEpsilonCoSolvingProcess
    KRATOS_CLASS_POINTER_DEFINITION(KEpsilonCoSolvingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    KEpsilonCoSolvingProcess(ModelPart& rModelPart, Parameters rParameters)
        : BaseType(rModelPart, rParameters, TURBULENT_VISCOSITY)
    {
    }

    /// Destructor.
    ~KEpsilonCoSolvingProcess() override
    {
    }

    void ExecuteInitializeSolutionStep() override
    {
        UpdateBeforeSolveEquations();
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

        NodesContainerType& r_nodes = this->mrModelPart.Nodes();
        int number_of_nodes = r_nodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_Y_PLUS, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        }

        return value;
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
        return "KEpsilonCoSolvingProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KEpsilonCoSolvingProcess";
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

    Vector mrOldValues;
    Vector mrNewValues;
    Vector mrDeltaValues;

    double mDeltaNorm;
    double mSolutionNorm;

    ///@}
    ///@name Operations
    ///@{

    void UpdateBeforeSolveEquations() override
    {
        this->ExecuteAuxiliaryProcessesInitializeSolutionStep();
        this->UpdateConvergenceVariable();
    }

    void UpdateAfterSolveEquations() override
    {
        UpdateEffectiveViscosity();
    }

    void UpdateConvergenceVariable() override
    {
        const double soft_max_exponent = this->mrModelPart.GetProcessInfo()[RANS_SOFT_MAX_EXPONENT];
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

        double current_delta_norm{0.0}, current_solution_norm{0.0};
// using soft max for the time being as a test
#pragma omp parallel for reduction(+: current_delta_norm, current_solution_norm)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_nodes.begin() + i_node);
            const double nu_t = mrOldValues[i_node] + mrDeltaValues[i_node] * this->mRelaxationFactor;
            const double soft_max = RansCalculationUtilities::SoftMax(
                nu_t, std::numeric_limits<double>::epsilon(), soft_max_exponent);
            current_delta_norm += std::pow(soft_max - mrOldValues[i_node], 2);
            current_solution_norm += std::pow(soft_max, 2);
            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = soft_max;
        }

        mDeltaNorm = current_delta_norm;
        mSolutionNorm = current_solution_norm;

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

    void CalculateConvergenceNorms(double& rDeltaNorm, double& rSolutionNorm) override
    {
        rDeltaNorm = mDeltaNorm;
        rSolutionNorm = mSolutionNorm;
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

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(
        KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    // KEpsilonCoSolvingProcess(KEpsilonCoSolvingProcess const& rOther);

    ///@}

}; // Class KEpsilonCoSolvingProcess

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
                                KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// output stream function
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
inline std::ostream& operator<<(std::ostream& rOStream,
                                const KEpsilonCoSolvingProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // K_EPSILON_CO_SOLVING_PROCESS_H_INCLUDED defined
