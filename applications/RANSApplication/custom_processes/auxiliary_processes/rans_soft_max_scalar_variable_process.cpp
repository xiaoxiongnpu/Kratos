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

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"

// Include base h
#include "rans_soft_max_scalar_variable_process.h"

namespace Kratos
{
RansSoftMaxScalarVariableProcess::RansSoftMaxScalarVariableProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"   : "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "echo_level"      : 0,
            "min_value"       : 1e-18,
            "exponent_value"  : 20.0,
            "store_original_value_in_non_historical_data_container" : true
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName = mrParameters["variable_name"].GetString();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mEchoLevel = mrParameters["echo_level"].GetInt();
    mMinValue = mrParameters["min_value"].GetDouble();
    mExponentValue = mrParameters["exponent_value"].GetDouble();
    mStoreOriginalValueInNonHistoricalDataValueContainer =
        mrParameters["store_original_value_in_non_historical_data_container"].GetBool();

    KRATOS_CATCH("");
}

RansSoftMaxScalarVariableProcess::~RansSoftMaxScalarVariableProcess()
{
}

int RansSoftMaxScalarVariableProcess::Check()
{
    KRATOS_TRY

    const Variable<double>& r_scalar_variable =
        KratosComponents<Variable<double>>::Get(mVariableName);

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        mrModel.GetModelPart(mModelPartName), r_scalar_variable);

    return 0;

    KRATOS_CATCH("");
}

void RansSoftMaxScalarVariableProcess::Execute()
{
    KRATOS_TRY

    const Variable<double>& r_scalar_variable =
        KratosComponents<Variable<double>>::Get(mVariableName);

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    const int number_of_nodes = r_model_part.NumberOfNodes();

    if (mStoreOriginalValueInNonHistoricalDataValueContainer)
    {
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            const double scalar_value = r_node.FastGetSolutionStepValue(r_scalar_variable);
            r_node.SetValue(r_scalar_variable, scalar_value);
        }
        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Original solution step values of " << mVariableName << " in " << mModelPartName
            << " is copied to its non-historical data value container.\n";
    }

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
        double& scalar_value = r_node.FastGetSolutionStepValue(r_scalar_variable);
        scalar_value = RansCalculationUtilities::SoftMax(scalar_value, mMinValue, mExponentValue);
    }

    KRATOS_CATCH("");
}

std::string RansSoftMaxScalarVariableProcess::Info() const
{
    return std::string("RansSoftMaxScalarVariableProcess");
}

void RansSoftMaxScalarVariableProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansSoftMaxScalarVariableProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
