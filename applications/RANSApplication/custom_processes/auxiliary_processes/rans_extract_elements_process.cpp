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

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/find_global_nodal_elemental_neighbours_process.h"
#include "processes/process.h"

#include "custom_utilities/rans_check_utilities.h"

// Include base h
#include "rans_extract_elements_process.h"

namespace Kratos
{
RansExtractElementsProcess::RansExtractElementsProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"                : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"                     : 0,
            "flag_variable_name"             : "PLEASE_SPECIFY_FLAG_VARIABLE_NAME",
            "flag_variable_value"            : true,
            "output_model_part_name"         : "",
            "assign_flag_to_elements"        : true
        })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mModelPartName = mrParameters["model_part_name"].GetString();
    mFlagVariableName = mrParameters["flag_variable_name"].GetString();
    mFlagVariableValue = mrParameters["flag_variable_value"].GetBool();
    mAssignFlagToElements = mrParameters["assign_flag_to_elements"].GetBool();
    mEchoLevel = mrParameters["echo_level"].GetInt();
    mOutputModelPartName = mrParameters["output_model_part_name"].GetString();

    if (mOutputModelPartName == "")
    {
        mOutputModelPartName =
            mModelPartName + ((mFlagVariableValue) ? "_" : "_!") + mFlagVariableName;
    }

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    KRATOS_ERROR_IF(mrModel.GetModelPart(mModelPartName).HasSubModelPart(mOutputModelPartName))
        << mOutputModelPartName << " already exists as a sub-model part in "
        << mModelPartName << ". Please choose a different \"output_model_part_name\".";

    mrModel.GetModelPart(mModelPartName).CreateSubModelPart(mOutputModelPartName);

    KRATOS_CATCH("");
}

int RansExtractElementsProcess::Check()
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH("");
}

void RansExtractElementsProcess::ExecuteInitialize()
{
    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    FindGlobalNodalElementalNeighboursProcess p_global_nodal_elemental_neighbours_process(
        r_model_part.GetCommunicator().GetDataCommunicator(), r_model_part);

    const int number_of_elements = r_model_part.NumberOfElements();
    ModelPart& r_output_model_part = r_model_part.GetSubModelPart(mOutputModelPartName);

    const Flags& r_flag = KratosComponents<Flags>::Get(mFlagVariableName);

    p_global_nodal_elemental_neighbours_process.Execute();
    KRATOS_INFO(this->Info())
        << "Completed nodal elemental neighbours process for " << mModelPartName << ".\n";

    const int number_of_nodes = r_model_part.NumberOfNodes();
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        ModelPart::NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
        if (r_node.Is(r_flag) == mFlagVariableValue)
        {
            GlobalPointersVector<Element>& r_neighbour_elements =
                r_node.GetValue(NEIGHBOUR_ELEMENTS);
            for (auto& r_element : r_neighbour_elements)
            {
                r_element.Set(r_flag, mFlagVariableValue);
            }
        }
    }

// #pragma omp parallel for
//     for (int i_element = 0; i_element < number_of_elements; ++i_element)
//     {
//         ModelPart::ElementType& r_element = *(r_model_part.ElementsBegin() + i_element);
//         if (r_element.Is(r_flag) == mFlagVariableValue)
//         {
//             r_output_model_part.AddElement(r_element.shared_from_this());
//         }
//     }

    KRATOS_INFO(this->Info())
        << "Completed extracting elements from " << mModelPartName << " with "
        << mFlagVariableName << "=" << ((mFlagVariableValue) ? "true" : "false")
        << " to " << mOutputModelPartName << ".\n";
}

std::string RansExtractElementsProcess::Info() const
{
    return std::string("RansExtractElementsProcess");
}

void RansExtractElementsProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansExtractElementsProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
