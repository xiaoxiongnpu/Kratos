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

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_check_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "rans_omega_y_plus_wall_function_process.h"

namespace Kratos
{
RansOmegaYPlusWallFunctionProcess::RansOmegaYPlusWallFunctionProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "beta_zero"       : 0.072,
            "c_mu"            : 0.09,
            "von_karman"      : 0.41,
            "beta"            : 5.2,
            "is_fixed"        : true
        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mBetaZero = mrParameters["beta_zero"].GetDouble();
    mVonKarman = mrParameters["von_karman"].GetDouble();
    mBeta = mrParameters["beta"].GetDouble();
    mCmu = mrParameters["c_mu"].GetDouble(); // equivalent to beta star
    mIsFixed = mrParameters["is_fixed"].GetBool();
    mLimitYPlus =
        RansCalculationUtilities::CalculateLogarithmicYPlusLimit(mVonKarman, mBeta);

    KRATOS_CATCH("");
}

/// Destructor.
RansOmegaYPlusWallFunctionProcess::~RansOmegaYPlusWallFunctionProcess()
{
}

int RansOmegaYPlusWallFunctionProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, KINEMATIC_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, RANS_Y_PLUS);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_KINETIC_ENERGY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, VELOCITY);

    return 0;

    KRATOS_CATCH("");
}

void RansOmegaYPlusWallFunctionProcess::ExecuteInitialize()
{
    KRATOS_TRY

    if (mIsFixed)
    {
        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
        const int number_of_nodes = r_model_part.NumberOfNodes();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
            r_node.Fix(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        }
    }

    KRATOS_CATCH("");
}

void RansOmegaYPlusWallFunctionProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const int number_of_nodes = r_model_part.NumberOfNodes();
    unsigned int number_of_modified_omega_wall_nodes = 0;
    double c_mu_25 = std::pow(mCmu, 0.25); // 4th rooth of cmu

#pragma omp parallel for reduction(+ : number_of_modified_omega_wall_nodes)
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_model_part.NodesBegin() + i_node);
        // if (r_node.Id() = 1000)
        //     std::cout<< r_node.FastGetSolutionStepValue(VELOCITY)<<std::endl;
        const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
        const double nu_m = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const array_1d<double, 3>& r_velocity = r_node.FastGetSolutionStepValue(VELOCITY);
        double velocity_magnitude = norm_2(r_velocity);

        //         KRATOS_ERROR_IF(y_plus < std::numeric_limits<double>::epsilon())
        // -          << "The input y_plus should be greater than zero.\n";
        // raising flag nor kratos error working on console -> ask suneth
        // if(y_plus < std::numeric_limits<double>::epsilon())
        //        throw std::runtime_error( "The input y_plus should be greater than zero.\n");

        if (y_plus < mLimitYPlus) // if y_plus is less than limit y plus it is in the viscous zones, hence we have to apply wall function
        {
            r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE) =
                (6 * velocity_magnitude * velocity_magnitude) /
                (mBetaZero * std::pow(y_plus, 4) * nu_m);
            number_of_modified_omega_wall_nodes++;
        }
        else // if y plus is greater than the limit y plus, we are not in wall zone so apply the wall function in the log layer
        {
            double u_plus = (1 / mVonKarman) * std::log(y_plus) + mBeta;
            double u_tau_squared = std::pow(velocity_magnitude / (u_plus), 2);
            double denominator_product = c_mu_25 * mVonKarman * y_plus * nu_m;

            r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE) =
                u_tau_squared / denominator_product;
            number_of_modified_omega_wall_nodes++;
        }
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Applied omega y_plus wall function to "
        << number_of_modified_omega_wall_nodes << " of total "
        << r_model_part.NumberOfNodes() << " nodes in " << mModelPartName << "\n";

    KRATOS_CATCH("");
}

std::string RansOmegaYPlusWallFunctionProcess::Info() const
{
    return std::string("RansOmegaYPlusWallFunctionProcess");
}

void RansOmegaYPlusWallFunctionProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansOmegaYPlusWallFunctionProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
