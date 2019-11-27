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

// Application  includes
#include "rans_application_variables.h"
#include "custom_elements/evm_k_omega/evm_k_omega_utilities.h"
#include "custom_utilities/rans_check_utilities.h"

// Include base h
#include "rans_nut_k_omega_high_re_calculation_process.h"
namespace Kratos
{


RansNutKOmegaHighReCalculationProcess::RansNutKOmegaHighReCalculationProcess(
    Model& rModel, Parameters& rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "echo_level"      : 0
    })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();

    KRATOS_CATCH("");
}

/// Destructor.
RansNutKOmegaHighReCalculationProcess::~RansNutKOmegaHighReCalculationProcess()
{
}

int RansNutKOmegaHighReCalculationProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_KINETIC_ENERGY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_VISCOSITY);

    return 0;

    KRATOS_CATCH("");
}

void RansNutKOmegaHighReCalculationProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    NodesContainerType& r_nodes = r_model_part.Nodes();
    const int number_of_nodes = r_nodes.size();

#pragma omp parallel for
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_nodes.begin() + i_node);
        const double tke = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
        const double omega =
            r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        const double nu_t = EvmKomegaModelUtilities::CalculateTurbulentViscosity(
            tke, omega);

        r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = nu_t;
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1)
        << "Calculated nu_t for nodes in " << mModelPartName << "\n";

    KRATOS_CATCH("");
}

/// Turn back information as a string.
std::string RansNutKOmegaHighReCalculationProcess::Info() const
{
    return std::string("RansNutKOmegaHighReCalculationProcess");
}

/// Print information about this object.
void RansNutKOmegaHighReCalculationProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void RansNutKOmegaHighReCalculationProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.

