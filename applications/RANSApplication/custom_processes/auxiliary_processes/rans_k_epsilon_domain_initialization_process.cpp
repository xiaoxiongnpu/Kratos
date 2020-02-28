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
#include <cmath>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "rans_application_variables.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/rans_check_utilities.h"

// Include base h
#include "rans_k_epsilon_domain_initialization_process.h"

namespace Kratos
{
RansKEpsilonDomainInitializationProcess::RansKEpsilonDomainInitializationProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name"        : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "turbulent_intensity"    : 0.05,
            "turbulent_mixing_length": 1.0,
            "density"                : 1.2,
            "kinematic_viscosity"    : 1.525e-3,
            "c_mu"                   : 0.09,
            "min_values"             : {
                "turbulent_kinetic_energy"         : 1e-12,
                "turbulent_energy_dissipation_rate": 1e-12,
                "turbulent_viscosity": 1e-12
            }
        })");

    mrParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mTurbulentIntensity = mrParameters["turbulent_intensity"].GetDouble();
    mTurbulentMixingLength = mrParameters["turbulent_mixing_length"].GetDouble();
    mDensity = mrParameters["density"].GetDouble();
    mKinematicViscosity = mrParameters["kinematic_viscosity"].GetDouble();
    mMinTurbulentKineticEnergy =
        mrParameters["min_values"]["turbulent_kinetic_energy"].GetDouble();
    mMinTurbulentEnergyDissipationRate =
        mrParameters["min_values"]["turbulent_energy_dissipation_rate"].GetDouble();
    mMinTurbulentViscosity =
        mrParameters["min_values"]["turbulent_viscosity"].GetDouble();
    mCmu = mrParameters["c_mu"].GetDouble();

    mModelPartName = mrParameters["model_part_name"].GetString();

    KRATOS_ERROR_IF(mTurbulentIntensity < 0.0)
        << "Turbulent intensity needs to be positive in the modelpart "
        << mModelPartName << " [ " << mTurbulentIntensity << " < 0.0 ]\n.";

    KRATOS_ERROR_IF(mTurbulentMixingLength < 0.0)
        << "Turbulent mixing length needs to be positive in the modelpart "
        << mModelPartName << " [ " << mTurbulentMixingLength << " < 0.0 ]\n.";

    KRATOS_ERROR_IF(mDensity < 0.0)
        << "Density needs to be positive in the modelpart " << mModelPartName
        << " [ " << mDensity << " < 0.0 ]\n.";

    KRATOS_ERROR_IF(mKinematicViscosity < 0.0)
        << "Kinematic viscosity needs to be positive in the modelpart "
        << mModelPartName << " [ " << mKinematicViscosity << " < 0.0 ]\n.";

    KRATOS_CATCH("");
}

void RansKEpsilonDomainInitializationProcess::ExecuteInitialize()
{
    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    VariableUtils variable_utilities;
    variable_utilities.SetVariable(DENSITY, mDensity, r_model_part.Nodes());
    variable_utilities.SetVariable(KINEMATIC_VISCOSITY, mKinematicViscosity,
                                   r_model_part.Nodes());
}

void RansKEpsilonDomainInitializationProcess::Execute()
{
    if (!mInitialized)
    {
        mInitialized = true;

        ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

        const double c_mu_75 = std::pow(mCmu, 0.75);

        const int number_of_nodes = r_model_part.NumberOfNodes();
#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            NodeType& r_node = *(r_model_part.NodesBegin() + i_node);

            const array_1d<double, 3>& r_velocity =
                r_node.FastGetSolutionStepValue(VELOCITY);
            const double velocity_magnitude = norm_2(r_velocity);

            double& tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            double& epsilon =
                r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
            double& nu_t = r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            double& nu = r_node.FastGetSolutionStepValue(VISCOSITY);

            tke = std::max(1.5 * std::pow(mTurbulentIntensity * velocity_magnitude, 2),
                           mMinTurbulentKineticEnergy);
            epsilon = std::max(c_mu_75 * std::pow(std::max(tke, 0.0), 1.5) / mTurbulentMixingLength,
                               mMinTurbulentEnergyDissipationRate);
            nu_t = std::max(mCmu * std::pow(tke, 2) / epsilon, mMinTurbulentViscosity);
            nu = nu_t + mKinematicViscosity;
        }

        KRATOS_INFO(this->Info()) << "Initialized " << this->mModelPartName << ".\n";
    }
}

int RansKEpsilonDomainInitializationProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, DENSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, VELOCITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, KINEMATIC_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_KINETIC_ENERGY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE);

    return 0;

    KRATOS_CATCH("");
}

std::string RansKEpsilonDomainInitializationProcess::Info() const
{
    return std::string("RansKEpsilonDomainInitializationProcess");
}

void RansKEpsilonDomainInitializationProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansKEpsilonDomainInitializationProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
