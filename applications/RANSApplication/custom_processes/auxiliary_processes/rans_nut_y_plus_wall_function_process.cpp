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
            "c1"              : 0.1,
            "c2"              : 0.2,
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
    mC1 = mrParameters["c1"].GetDouble();
    mC2 = mrParameters["c2"].GetDouble();

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
//     VariableUtils().SetHistoricalVariableToZero(TURBULENT_VISCOSITY,
//                                                 r_model_part.Nodes());
//     VariableUtils().SetNonHistoricalVariableToZero(RANS_Y_PLUS, r_model_part.Nodes());

//     const int number_of_conditions = r_model_part.NumberOfConditions();
// #pragma omp parallel for
//     for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
//     {
//         ConditionType& r_cond = *(r_model_part.ConditionsBegin() + i_cond);
//         if (r_cond.Is(SLIP))
//         {
//             const double y_plus = r_cond.GetValue(RANS_Y_PLUS);
//             const double nu = r_cond.GetValue(KINEMATIC_VISCOSITY);
//             ConditionType::GeometryType& r_geometry = r_cond.GetGeometry();
//             if (y_plus > 0.0) // logarithmic region
//             {
//                 for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
//                 {
//                     NodeType& r_node = r_geometry[i_node];
//                     r_node.SetLock();
//                     r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) +=
//                         mVonKarman * y_plus * nu;
//                     r_node.GetValue(RANS_Y_PLUS) += y_plus;
//                     r_node.UnSetLock();
//                 }
//             }
//         }
//     }
//     r_model_part.GetCommunicator().AssembleCurrentData(TURBULENT_VISCOSITY);
//     r_model_part.GetCommunicator().AssembleNonHistoricalData(RANS_Y_PLUS);

//     const int number_of_nodes = r_model_part.NumberOfNodes();
// #pragma omp parallel for
//     for (int i_node = 0; i_node < number_of_nodes; ++i_node)
//     {
//         NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
//         double& r_nut = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
//         const double number_of_neighbour_conditions =
//             static_cast<double>(r_node.GetValue(NUMBER_OF_NEIGHBOUR_CONDITIONS));
//         r_node.GetValue(RANS_Y_PLUS) /= number_of_neighbour_conditions;
//         r_nut = RansCalculationUtilities::SoftMax(
//             r_nut / number_of_neighbour_conditions, mMinValue);
//     }

        const int number_of_nodes = r_model_part.NumberOfNodes();

        unsigned int number_of_modified_nu_t_nodes = 0;

    #pragma omp parallel for reduction(+ : number_of_modified_nu_t_nodes)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
            const double nu = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);

            if (y_plus > mLimitYPlus)
            {
                r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) =
                    std::max(mVonKarman * y_plus * nu, mMinValue);
                ++number_of_modified_nu_t_nodes;
            }
            else
            {
                r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = mMinValue;
            }
        }

        KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
            << "Applied nu_t y_plus wall function to "
            << number_of_modified_nu_t_nodes << " of total "
            << r_model_part.NumberOfNodes() << " nodes in " << mModelPartName << "\n";

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
