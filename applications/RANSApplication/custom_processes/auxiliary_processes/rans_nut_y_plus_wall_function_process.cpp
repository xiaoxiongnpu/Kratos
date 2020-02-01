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
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_nut_y_plus_wall_function_process.h"

namespace Kratos
{
RansNutYPlusWallFunctionProcess::RansNutYPlusWallFunctionProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "c_mu"            : 0.09,
            "von_karman"      : 0.41,
            "beta"            : 5.2,
            "min_value"       : 1e-18
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mCmu = mrParameters["c_mu"].GetDouble();
    mVonKarman = mrParameters["von_karman"].GetDouble();
    mBeta = mrParameters["beta"].GetDouble();
    mMinValue = mrParameters["min_value"].GetDouble();
    mLimitYPlus =
        RansCalculationUtilities::CalculateLogarithmicYPlusLimit(mVonKarman, mBeta);

    KRATOS_CATCH("");
}

int RansNutYPlusWallFunctionProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, KINEMATIC_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, RANS_Y_PLUS);

    return 0;

    KRATOS_CATCH("");
}

void RansNutYPlusWallFunctionProcess::ExecuteInitialize()
{
    CalculateConditionNeighbourCount();
}

void RansNutYPlusWallFunctionProcess::CalculateConditionNeighbourCount()
{
    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
    VariableUtils().SetNonHistoricalVariableToZero(
        NUMBER_OF_NEIGHBOUR_CONDITIONS, r_model_part.Nodes());

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ConditionType& r_cond = *(r_model_part.ConditionsBegin() + i_cond);
        ConditionGeometryType& r_geometry = r_cond.GetGeometry();
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
        {
            NodeType& r_node = r_geometry[i_node];
            r_node.SetLock();
            r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS) += 1;
            r_node.UnSetLock();
        }
    }

    r_model_part.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_CONDITIONS);
}

void RansNutYPlusWallFunctionProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const double y_plus_limit = r_model_part.GetProcessInfo()[RANS_Y_PLUS_LIMIT];
    const double kappa = mVonKarman;
    const double beta = mBeta;

    VariableUtils().SetHistoricalVariableToZero(TURBULENT_VISCOSITY,
                                                r_model_part.Nodes());

    const int number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for
    for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
    {
        ConditionType& r_cond = *(r_model_part.ConditionsBegin() + i_cond);
        if (r_cond.Is(SLIP))
        {
            ConditionType::GeometryType& r_geometry = r_cond.GetGeometry();
            const double u_tau = r_cond.GetValue(FRICTION_VELOCITY);

            // log region, apply the u_tau which is calculated based on turbulent
            // kinetic energy. gauss point y_plus is estimated assuming the same u_tau
            // derrived from turbulent kinetic energy in the log region equation.
            const std::function<double(double)> log_region_functional =
                [u_tau, kappa, beta](const double velocity_magnitude) -> double {
                return std::exp((velocity_magnitude / u_tau - beta) * kappa);
            };

            // In the linear region, force the velocity to be in the log region
            // lowest y_plus since k - epsilon is only valid in the log region.
            const std::function<double(double)> linear_region_functional =
                [y_plus_limit](const double) -> double { return y_plus_limit; };

            const std::function<double(double)> wall_y_plus_functional =
                (r_cond.Is(MARKER) ? log_region_functional : linear_region_functional);

            for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
            {
                NodeType& r_node = r_geometry[i_node];
                const array_1d<double, 3>& r_velocity =
                    r_node.FastGetSolutionStepValue(VELOCITY);
                const double r_velocity_magnitude = norm_2(r_velocity);
                const double wall_y_plus = wall_y_plus_functional(r_velocity_magnitude);
                const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
                r_node.SetLock();
                r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) +=
                    mVonKarman * wall_y_plus * nu;
                r_node.UnSetLock();
            }
        }
    }
    r_model_part.GetCommunicator().AssembleCurrentData(TURBULENT_VISCOSITY);

    const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
        double& r_nut = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        const double number_of_neighbour_conditions =
            static_cast<double>(r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS));
        r_nut = RansCalculationUtilities::SoftMax(
            r_nut / number_of_neighbour_conditions, mMinValue);
    }

    KRATOS_CATCH("");
}

std::string RansNutYPlusWallFunctionProcess::Info() const
{
    return std::string("RansNutYPlusWallFunctionProcess");
}

void RansNutYPlusWallFunctionProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutYPlusWallFunctionProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
