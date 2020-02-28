//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//
//

/* System includes */
#include <cmath>
#include <limits>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"

// Application incldues
#include "rans_application_variables.h"

// Include base h
#include "rans_variable_utilities.h"

namespace Kratos
{
namespace RansVariableUtilities
{
void ClipScalarVariable(unsigned int& rNumberOfNodesBelowMinimum,
                        unsigned int& rNumberOfNodesAboveMaximum,
                        unsigned int& rNumberOfSelectedNodes,
                        const double MinimumValue,
                        const double MaximumValue,
                        const Variable<double>& rVariable,
                        ModelPart& rModelPart)
{
    KRATOS_TRY

    Communicator& r_communicator = rModelPart.GetCommunicator();
    ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();

    const int number_of_nodes = r_nodes.size();

    unsigned int number_of_nodes_below_minimum = 0;
    unsigned int number_of_nodes_above_maximum = 0;
    unsigned int number_of_nodes_selected = 0;

#pragma omp parallel for reduction( +: number_of_nodes_below_minimum, number_of_nodes_above_maximum, number_of_nodes_selected)
    for (int i = 0; i < number_of_nodes; ++i)
    {
        ModelPart::NodeType& r_node = *(r_nodes.begin() + i);
        double& r_value = r_node.FastGetSolutionStepValue(rVariable);

        if (r_value < MinimumValue)
        {
            number_of_nodes_below_minimum++;
            r_value = MinimumValue;
        }
        else if (r_value > MaximumValue)
        {
            number_of_nodes_above_maximum++;
            r_value = MaximumValue;
        }
        number_of_nodes_selected++;
    }

    r_communicator.SynchronizeVariable(rVariable);

    // Stores followings
    // index - 0 : number_of_nodes_below_minimum
    // index - 1 : number_of_nodes_above_maximum
    // index - 2 : number_of_nodes_selected
    std::vector<unsigned int> nodes_count = {number_of_nodes_below_minimum,
                                             number_of_nodes_above_maximum,
                                             number_of_nodes_selected};
    const std::vector<unsigned int>& total_nodes_count =
        r_communicator.GetDataCommunicator().SumAll(nodes_count);

    rNumberOfNodesBelowMinimum = total_nodes_count[0];
    rNumberOfNodesAboveMaximum = total_nodes_count[1];
    rNumberOfSelectedNodes = total_nodes_count[2];

    KRATOS_CATCH("")
}

double GetMinimumScalarValue(const ModelPart& rModelPart, const Variable<double>& rVariable)
{
    KRATOS_TRY

    double min_value = std::numeric_limits<double>::max();

    const Communicator& r_communicator = rModelPart.GetCommunicator();

    const ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();

    const int number_of_nodes = r_nodes.size();

    if (number_of_nodes != 0)
    {
        const int number_of_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector node_partition;
        OpenMPUtils::DivideInPartitions(number_of_nodes, number_of_threads, node_partition);
        Vector min_values(number_of_threads, min_value);

#pragma omp parallel
        {
            const int k = OpenMPUtils::ThisThread();

            auto nodes_begin = r_nodes.begin() + node_partition[k];
            auto nodes_end = r_nodes.begin() + node_partition[k + 1];

            for (auto itNode = nodes_begin; itNode != nodes_end; ++itNode)
            {
                const double value = itNode->FastGetSolutionStepValue(rVariable);
                min_values[k] = std::min(min_values[k], value);
            }
        }

        for (int i = 0; i < number_of_threads; ++i)
        {
            min_value = std::min(min_value, min_values[i]);
        }
    }

    return r_communicator.GetDataCommunicator().MinAll(min_value);

    KRATOS_CATCH("");
}

double GetMaximumScalarValue(const ModelPart& rModelPart, const Variable<double>& rVariable)
{
    KRATOS_TRY

    double max_value = std::numeric_limits<double>::lowest();

    const Communicator& r_communicator = rModelPart.GetCommunicator();

    const ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();

    const int number_of_nodes = r_nodes.size();

    if (number_of_nodes != 0)
    {
        const int number_of_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector node_partition;
        OpenMPUtils::DivideInPartitions(number_of_nodes, number_of_threads, node_partition);
        Vector max_values(number_of_threads, max_value);

#pragma omp parallel
        {
            const int k = OpenMPUtils::ThisThread();

            auto nodes_begin = r_nodes.begin() + node_partition[k];
            auto nodes_end = r_nodes.begin() + node_partition[k + 1];

            for (auto itNode = nodes_begin; itNode != nodes_end; ++itNode)
            {
                const double value = itNode->FastGetSolutionStepValue(rVariable);
                max_values[k] = std::max(max_values[k], value);
            }
        }

        for (int i = 0; i < number_of_threads; ++i)
        {
            max_value = std::max(max_value, max_values[i]);
        }
    }

    return r_communicator.GetDataCommunicator().MaxAll(max_value);

    KRATOS_CATCH("");
}

void GetNodalVariablesVector(Vector& rValues,
                             const ModelPart::NodesContainerType& rNodes,
                             const Variable<double>& rVariable)
{
    const int number_of_nodes = rNodes.size();

    if (static_cast<int>(rValues.size()) != number_of_nodes)
        rValues.resize(number_of_nodes);

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        rValues[i_node] = (rNodes.begin() + i_node)->FastGetSolutionStepValue(rVariable);
}

void GetNodalArray(Vector& rNodalValues, const Element& rElement, const Variable<double>& rVariable)
{
    const Geometry<ModelPart::NodeType>& r_geometry = rElement.GetGeometry();
    std::size_t number_of_nodes = r_geometry.PointsNumber();

    if (rNodalValues.size() != number_of_nodes)
        rNodalValues.resize(number_of_nodes);

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        rNodalValues[i_node] = r_geometry[i_node].FastGetSolutionStepValue(rVariable);
    }
}

void SetNodalVariables(ModelPart::NodesContainerType& rNodes,
                       const Vector& rValues,
                       const Variable<double>& rVariable)
{
    const int number_of_nodes = rNodes.size();

    KRATOS_ERROR_IF(static_cast<int>(rValues.size()) != number_of_nodes)
        << "rValues vector size mismatch with rNodes size in "
           "SetNodalVariables. [ rValues.size = "
        << rValues.size() << ", rNodes.size = " << rNodes.size() << " ]\n";

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        (rNodes.begin() + i_node)->FastGetSolutionStepValue(rVariable) = rValues[i_node];
}

void CopyNodalSolutionStepVariablesList(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart)
{
    KRATOS_TRY

    rDestinationModelPart.GetNodalSolutionStepVariablesList() =
        rOriginModelPart.GetNodalSolutionStepVariablesList();

    KRATOS_CATCH("");
}

void InitializeDuplicatedModelPart(ModelPart& rOriginModelPart, ModelPart& rDuplicatedModelPart)
{
    ModelPart::ElementsContainerType& r_elements = rDuplicatedModelPart.Elements();
    const int number_of_elements = r_elements.size();
#pragma omp parallel for
    for (int i_element = 0; i_element < number_of_elements; ++i_element)
    {
        ModelPart::ElementType& r_element = *(r_elements.begin() + i_element);
        ModelPart::ElementType& r_parent_element =
            rOriginModelPart.GetElement(r_element.Id());
        r_element.SetValue(PARENT_ELEMENT_POINTER, &r_parent_element);
    }
    KRATOS_INFO("DuplicateModelPartInitializer")
        << "Initialized " << rDuplicatedModelPart.Name()
        << " element data from " << rOriginModelPart.Name() << ".\n";

    ModelPart::ConditionsContainerType& r_conditions = rDuplicatedModelPart.Conditions();
    const int number_of_conditions = r_conditions.size();
#pragma omp parallel for
    for (int i_condition = 0; i_condition < number_of_conditions; ++i_condition)
    {
        ModelPart::ConditionType& r_condition = *(r_conditions.begin() + i_condition);
        ModelPart::ConditionType& r_parent_condition =
            rOriginModelPart.GetCondition(r_condition.Id());
        r_condition.SetValue(PARENT_CONDITION_POINTER, &r_parent_condition);
    }

    KRATOS_INFO("DuplicateModelPartInitializer")
        << "Initialized " << rDuplicatedModelPart.Name()
        << " condition data from " << rOriginModelPart.Name() << ".\n";
}

void FixFlaggedDofs(ModelPart& rModelPart,
                    const Variable<double>& rFixingVariable,
                    const Flags& rFlag,
                    const bool CheckValue)
{
    KRATOS_TRY

    const int number_of_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        ModelPart::NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        if (r_node.Is(rFlag) == CheckValue)
        {
            r_node.Fix(rFixingVariable);
        }
    }

    KRATOS_CATCH("");
}

template <typename TDataType>
void CopyFlaggedVariableFromNonHistorical(ModelPart& rModelPart,
                                          const Variable<TDataType>& rVariable,
                                          const Flags& rFlag,
                                          const bool CheckValue)
{
    KRATOS_TRY

    const int number_of_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        ModelPart::NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        if (r_node.Is(rFlag) == CheckValue)
        {
            r_node.FastGetSolutionStepValue(rVariable) = r_node.GetValue(rVariable);
        }
    }

    KRATOS_CATCH("");
}

template <typename TDataType>
void CopyFlaggedVariableToNonHistorical(ModelPart& rModelPart,
                                        const Variable<TDataType>& rVariable,
                                        const Flags& rFlag,
                                        const bool CheckValue)
{
    KRATOS_TRY

    const int number_of_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        ModelPart::NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        if (r_node.Is(rFlag) == CheckValue)
        {
            r_node.SetValue(rVariable, r_node.FastGetSolutionStepValue(rVariable));
        }
    }

    KRATOS_CATCH("");
}

void CalculateMagnitudeSquareFor3DVariable(ModelPart& rModelPart,
                                           const Variable<array_1d<double, 3>>& r3DVariable,
                                           const Variable<double>& rOutputVariable)
{
    const int number_of_nodes = rModelPart.NumberOfNodes();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        ModelPart::NodeType& r_node = *(rModelPart.NodesBegin() + i_node);
        const double magnitude = norm_2(r_node.FastGetSolutionStepValue(r3DVariable));
        r_node.FastGetSolutionStepValue(rOutputVariable) = std::pow(magnitude, 2);
    }
}

template void CopyFlaggedVariableFromNonHistorical<double>(ModelPart&,
                                                           const Variable<double>&,
                                                           const Flags&,
                                                           const bool);
template void CopyFlaggedVariableFromNonHistorical<array_1d<double, 3>>(
    ModelPart&, const Variable<array_1d<double, 3>>&, const Flags&, const bool);

template void CopyFlaggedVariableToNonHistorical<double>(ModelPart&,
                                                         const Variable<double>&,
                                                         const Flags&,
                                                         const bool);
template void CopyFlaggedVariableToNonHistorical<array_1d<double, 3>>(
    ModelPart&, const Variable<array_1d<double, 3>>&, const Flags&, const bool);
// template instantiations

} // namespace RansVariableUtilities

} /* namespace Kratos.*/
