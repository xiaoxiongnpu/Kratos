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
#include "rans_nut_y_plus_wall_function_omega_process.h"

namespace Kratos
{
RansNutYPlusWallFunctionOmegaProcess::RansNutYPlusWallFunctionOmegaProcess(Model& rModel, Parameters rParameters)
    : mrModel(rModel), mrParameters(rParameters)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
        {
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "echo_level"      : 0,
            "min_value"       : 1e-18,
            "von_karman"      : 0.41,
            "beta_zero"       : 0.072,
            "c_mu"            : 0.09

        })");

    mrParameters.ValidateAndAssignDefaults(default_parameters);

    mEchoLevel = mrParameters["echo_level"].GetInt();
    mModelPartName = mrParameters["model_part_name"].GetString();
    mMinValue = mrParameters["min_value"].GetDouble();
    mCmu = mrParameters["c_mu"].GetDouble();
    mVk = mrParameters["von_karman"].GetDouble();
    mB0 = mrParameters["beta_zero"].GetDouble();

    KRATOS_CATCH("");
}

/// Destructor.
RansNutYPlusWallFunctionOmegaProcess::~RansNutYPlusWallFunctionOmegaProcess()
{
}

int RansNutYPlusWallFunctionOmegaProcess::Check()
{
    KRATOS_TRY

    RansCheckUtilities::CheckIfModelPartExists(mrModel, mModelPartName);

    const ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);
  
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, KINEMATIC_VISCOSITY);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, RANS_Y_PLUS);
    RansCheckUtilities::CheckIfVariableExistsInModelPart(r_model_part, TURBULENT_VISCOSITY);
    return 0;

    KRATOS_CATCH("");
}

void RansNutYPlusWallFunctionOmegaProcess::Execute()
{
    KRATOS_TRY

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const int number_of_nodes = r_model_part.NumberOfNodes();

    unsigned int number_of_modified_nu_t_nodes = 0;
    
    const double c_mu_50= std::pow(mCmu,0.5);
    const double value_1= (36.0/(mB0*mB0));
    const double value_2= 1.0/(mVk*mVk*c_mu_50);

#pragma omp parallel for reduction(+ : number_of_modified_nu_t_nodes)
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        NodeType& r_node = *(r_nodes.begin() + i_node);
        const double nu_k = r_node.FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
        const double y_plus = r_node.FastGetSolutionStepValue(RANS_Y_PLUS);
        double value= std::pow(value_2+(value_1/(y_plus*y_plus*nu_k*nu_k)),0.5);
        const double nu_t = (y_plus * nu_k)/(c_mu_50*value);

        r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) =std::max(nu_t, mMinValue) ;
        ++number_of_modified_nu_t_nodes;

    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0)
        << "Applied nu_t y_plus wall function to "
        << number_of_modified_nu_t_nodes << " of total "
        << r_model_part.NumberOfNodes() << " nodes in " << mModelPartName << "\n";

    KRATOS_CATCH("");
}

std::string RansNutYPlusWallFunctionOmegaProcess::Info() const
{
    return std::string("RansNutYPlusWallFunctionOmegaProcess");
}

void RansNutYPlusWallFunctionOmegaProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void RansNutYPlusWallFunctionOmegaProcess::PrintData(std::ostream& rOStream) const
{
}

} // namespace Kratos.
