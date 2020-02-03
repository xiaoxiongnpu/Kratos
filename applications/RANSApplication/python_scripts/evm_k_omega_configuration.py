from __future__ import print_function, absolute_import, division

import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.RANSApplication.turbulence_eddy_viscosity_model_configuration import TurbulenceEddyViscosityModelConfiguration

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if not CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    msg = "k-omega turbulence model depends on the FluidDynamicsApplication which is not found."
    msg += " Please re-install/compile with FluidDynamicsApplication"
    raise Exception(msg)

if (Kratos.IsDistributedRun()):
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIKEpsilonCoSolvingProcess as k_epsilon_co_solving_process
else:
    from KratosMultiphysics.RANSApplication import KEpsilonCoSolvingProcess as k_epsilon_co_solving_process

class TurbulenceKOmegaConfiguration(
        TurbulenceEddyViscosityModelConfiguration):
    def __init__(self, model, parameters):
        super(TurbulenceKOmegaConfiguration, self).__init__(
            model, parameters)
        self.turbulence_model_process = None

        default_settings = Kratos.Parameters(r'''{
            "scheme_settings": {},
            "echo_level"          :0,
            "turbulent_kinetic_energy_settings":{},
            "turbulent_specific_energy_dissipation_rate_settings":{},
            "is_blended"              : false,
            "constants":
            {
                "sigma_k"                 : 0.5,
                "wall_smoothness_beta"    : 5.2,
                "von_karman"              : 0.41,
                "sigma_omega"             : 0.5,
                "beta_zero"               : 0.072,
                "gamma"                   : 0.52,
                "limit_y_plus_wall"       : 11.06,
                "beta_zero_star"          : 0.09
            },
            "flow_parameters":
            {
                "ramp_up_time"            : 0.5
            },
            "coupling_settings" :{}
         }''')

        parameters["model_settings"].ValidateAndAssignDefaults(
            default_settings)
        parameters["model_settings"]["constants"].ValidateAndAssignDefaults(
            default_settings["constants"])
        parameters["model_settings"][
            "flow_parameters"].ValidateAndAssignDefaults(
                default_settings["flow_parameters"])
        self.model_settings = parameters["model_settings"]
        is_blended = self.model_settings["is_blended"].GetBool()

        self.model_elements_list = ["RansEvmKOmegaK", "RansEvmKOmegaOmega"]
        if is_blended:
            self.model_conditions_list = ["Condition", "RansEvmKOmegaOmegaWallBlended"]
        else:
            self.model_conditions_list = ["Condition", "Condition"]
        
        print("HERE ARE THE CONDITIONS: ", self.model_conditions_list)
        

        self.ramp_up_time = self.model_settings["flow_parameters"][
            "ramp_up_time"].GetDouble()

    def InitializeModelConstants(self):
        # reading constants
        constants = self.model_settings["constants"]
        self.fluid_model_part.ProcessInfo[
            KratosRANS.WALL_SMOOTHNESS_BETA] = constants[
                "wall_smoothness_beta"].GetDouble()
        self.fluid_model_part.ProcessInfo[
            KratosRANS.WALL_VON_KARMAN] = constants["von_karman"].GetDouble()
        self.fluid_model_part.ProcessInfo[
            KratosRANS.TURBULENCE_RANS_BETA_ZERO] = constants["beta_zero"].GetDouble()
        self.fluid_model_part.ProcessInfo[
            KratosRANS.TURBULENCE_RANS_GAMMA] = constants["gamma"].GetDouble()
        self.fluid_model_part.ProcessInfo[
            KratosRANS.TURBULENCE_RANS_BETA_ZERO_STAR] = constants["beta_zero_star"].GetDouble()
        self.fluid_model_part.ProcessInfo[
            KratosRANS.TURBULENCE_RANS_SIGMA_K] = constants["sigma_k"].GetDouble()
        self.fluid_model_part.ProcessInfo[
            KratosRANS.TURBULENCE_RANS_SIGMA_OMEGA] = constants["sigma_omega"].GetDouble()
        self.fluid_model_part.ProcessInfo[
            KratosRANS.TURBULENCE_RANS_Y_PLUS_LIMIT_WALL] = constants["limit_y_plus_wall"].GetDouble()
        # self.fluid_model_part.ProcessInfo[
        #     KratosRANS.TURBULENCE_BLENDING] = constants["is_blended"].GetBool()
        

    def PrepareSolvingStrategy(self):
        scheme_settings = self.model_settings["scheme_settings"]

        # create turbulent kinetic energy strategy
        model_part = self.model_parts_list[0]
        solver_settings = self.model_settings[
            "turbulent_kinetic_energy_settings"]
        scalar_variable = KratosRANS.TURBULENT_KINETIC_ENERGY
        scalar_variable_rate = KratosRANS.TURBULENT_KINETIC_ENERGY_RATE
        relaxed_scalar_variable_rate = KratosRANS.RANS_AUXILIARY_VARIABLE_1
        current_strategy = self.CreateStrategy(
            solver_settings, scheme_settings, model_part, scalar_variable,
            scalar_variable_rate, relaxed_scalar_variable_rate)
        self.strategies_list.append(current_strategy)
        self.GetTurbulenceSolvingProcess().AddStrategy(current_strategy,
                                                       scalar_variable)

        # create turbulent specific energy dissipation rate strategy(omega strategy)
        model_part = self.model_parts_list[1]
        solver_settings = self.model_settings[
            "turbulent_specific_energy_dissipation_rate_settings"]
        scalar_variable = KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE
        scalar_variable_rate = KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2
        relaxed_scalar_variable_rate = KratosRANS.RANS_AUXILIARY_VARIABLE_2
        current_strategy = self.CreateStrategy(
            solver_settings, scheme_settings, model_part, scalar_variable,
            scalar_variable_rate, relaxed_scalar_variable_rate)
        self.strategies_list.append(current_strategy)
        self.GetTurbulenceSolvingProcess().AddStrategy(current_strategy,
                                                       scalar_variable)

        Kratos.Logger.PrintInfo(
            self.__class__.__name__,
            "All turbulence solution strategies are created.")

    def AddVariables(self):
        # adding k-epsilon specific variables
        self.fluid_model_part.AddNodalSolutionStepVariable(
            KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            KratosRANS.TURBULENT_KINETIC_ENERGY_RATE)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            KratosRANS.RANS_AUXILIARY_VARIABLE_1)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            KratosRANS.RANS_AUXILIARY_VARIABLE_2)

        super(TurbulenceKOmegaConfiguration, self).AddVariables()

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_KINETIC_ENERGY,
                                      self.fluid_model_part)
        Kratos.VariableUtils().AddDof(
            KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE,
            self.fluid_model_part)

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "DOFs added successfully.")

    def Initialize(self):
        super(TurbulenceKOmegaConfiguration, self).Initialize()
        self.InitializeModelConstants()

    def InitializeSolutionStep(self):
        if (self.fluid_model_part.ProcessInfo[KratosRANS.
                                              IS_CO_SOLVING_PROCESS_ACTIVE]):
            super(TurbulenceKOmegaConfiguration,
                  self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        if (self.fluid_model_part.ProcessInfo[KratosRANS.
                                              IS_CO_SOLVING_PROCESS_ACTIVE]):
            super(TurbulenceKOmegaConfiguration, self).FinalizeSolutionStep()

        time = self.fluid_model_part.ProcessInfo[Kratos.TIME]
        if (time >= self.ramp_up_time):
            self.fluid_model_part.ProcessInfo[
                KratosRANS.IS_CO_SOLVING_PROCESS_ACTIVE] = True


    def GetTurbulenceSolvingProcess(self):
        if self.turbulence_model_process is None:
            self.turbulence_model_process = k_epsilon_co_solving_process(
                self.fluid_model_part,
                self.model_settings["coupling_settings"])
            Kratos.Logger.PrintInfo(self.__class__.__name__,
                                    "Created turbulence solving process.")

        return self.turbulence_model_process

    def GetFluidVelocityPressureConditionName(self):
        return "RansEvmKEpsilonVmsMonolithicWall"
