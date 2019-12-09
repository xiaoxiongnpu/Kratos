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
#include "processes/process.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_vector_ramping_process.h"

namespace Kratos
{
RansVectorRampingProcess::RansVectorRampingProcess(Model& rModel, Parameters rParameters)
    : Process(), mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"    : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variable_name"      : "VELOCITY",
            "modulus"            : 0.0,
            "constrained"        : true,
            "direction"          : [1.0,0.0,0.0],
            "ramping_iterations" : 20,
            "echo_level"         : 0,
            "local_axes"         : {}
        })");

    if (mrParameters.Has("modulus") && mrParameters["modulus"].IsString())
    {
        default_parameters.RemoveValue("modulus");
        default_parameters.AddEmptyValue("modulus");
        default_parameters["modulus"].SetString("");
        mIsFunction = true;
    }

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    if (mIsFunction)
    {
        mpFunction = PythonGenericFunctionUtility::Pointer(new PythonGenericFunctionUtility(
            mrParameters["modulus"].GetString(), mrParameters["local_axes"]));
    }
    else
    {
        mModulus = mrParameters["modulus"].GetDouble();
    }

    mModelPartName = mrParameters["model_part_name"].GetString();
    mVariableName = mrParameters["variable_name"].GetString();
    mRampingIterations = mrParameters["ramping_iterations"].GetInt();
    mDirection = mrParameters["direction"].GetVector();
    mIsConstrained = mrParameters["constrained"].GetBool();
    mEchoLevel = mrParameters["echo_level"].GetInt();

    KRATOS_ERROR_IF(mDirection.size() != 3)
        << "Direction vector needs to provide values for all three "
           "dimensions.\n";

    const double direction_norm = norm_2(mDirection);
    KRATOS_ERROR_IF(direction_norm == 0.0)
        << "Direction vector should not be zero vector.\n";
    noalias(mDirection) = mDirection * (1.0 / direction_norm);

    KRATOS_ERROR_IF(mRampingIterations <= 0)
        << "Invalid \"ramping_iterations\"=" << mRampingIterations
        << ". Positive number is expected.\n";

    KRATOS_ERROR_IF(!(KratosComponents<Variable<array_1d<double, 3>>>::Has(mVariableName)))
        << "Provided \"variable_name\"=" << mVariableName
        << " is not in the 3D variables list.\n";

    KRATOS_CATCH("");
}

RansVectorRampingProcess::~RansVectorRampingProcess()
{
}

int RansVectorRampingProcess::Check()
{
    KRATOS_TRY

    const Variable<array_1d<double, 3>>& r_vector_variable =
        KratosComponents<Variable<array_1d<double, 3>>>::Get(mVariableName);

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, r_vector_variable);

    if (mIsConstrained)
    {
        const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& r_vector_x =
            KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                mVariableName + "_X");
        const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& r_vector_y =
            KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                mVariableName + "_Y");
        const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& r_vector_z =
            KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                mVariableName + "_Z");

        const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const ModelPart::NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            KRATOS_CHECK_DOF_IN_NODE(r_vector_x, r_node);
            KRATOS_CHECK_DOF_IN_NODE(r_vector_y, r_node);
            KRATOS_CHECK_DOF_IN_NODE(r_vector_z, r_node);
        }
    }

    return 0;

    KRATOS_CATCH("");
}

void RansVectorRampingProcess::ExecuteInitialize()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[IS_RAMPING] = true;
    mCurrentRampingIteration = 1;

    if (mIsConstrained)
    {
        const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& r_vector_x =
            KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                mVariableName + "_X");
        const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& r_vector_y =
            KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                mVariableName + "_Y");
        const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& r_vector_z =
            KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Get(
                mVariableName + "_Z");

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
        const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            ModelPart::NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            r_node.Fix(r_vector_x);
            r_node.Fix(r_vector_y);
            r_node.Fix(r_vector_z);
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Constrained dofs of " << mVariableName << " in "
            << mModelPartName << ".\n";
    }

    KRATOS_CATCH("");
}

void RansVectorRampingProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

    if (r_process_info[IS_RAMPING])
    {
        const Variable<array_1d<double, 3>>& r_vector_variable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mVariableName);

        double weighting_factor = static_cast<double>(mCurrentRampingIteration) /
                                  static_cast<double>(mRampingIterations);

        const double current_time = r_process_info[TIME];

        const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            const double value =
                (mpFunction->CallFunction(r_node.X(), r_node.Y(), r_node.Z(), current_time)) *
                weighting_factor;

            array_1d<double, 3>& r_vector =
                r_node.FastGetSolutionStepValue(r_vector_variable);
            r_vector[0] = mDirection[0] * value;
            r_vector[1] = mDirection[1] * value;
            r_vector[2] = mDirection[2] * value;
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
            << "Applied ramping to " << mVariableName << " in " << mModelPartName
            << " with ramping factor = " << weighting_factor << " [ "
            << mCurrentRampingIteration << "/" << mRampingIterations << " ].\n";

        r_process_info[IS_RAMPING] = (weighting_factor < 1.0);
        ++mCurrentRampingIteration;
    }

    KRATOS_CATCH("");
}

std::string RansVectorRampingProcess::Info() const
{
    return std::string("RansVectorRampingProcess");
}

void RansVectorRampingProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansVectorRampingProcess::PrintData(std::ostream& rOStream) const
{
}
} // namespace Kratos.
