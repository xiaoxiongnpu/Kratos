from __future__ import print_function, absolute_import, division

import KratosMultiphysics as Kratos
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.RANSApplication.model_part_factory import CreateDuplicateModelPart

from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics import IntegrationValuesExtrapolationToNodesProcess as extrapolation_process
from KratosMultiphysics import ResidualBasedIncrementalUpdateStaticScheme as scheme
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if (IsDistributedRun()
        and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory as linear_solver_factory
    from KratosMultiphysics.TrilinosApplication import TrilinosResidualCriteria as residual_criteria
    from KratosMultiphysics.TrilinosApplication import TrilinosNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosBlockBuilderAndSolver as block_builder_and_solver
elif (not IsDistributedRun()):
    from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
    from Kratos import ResidualCriteria as residual_criteria
    from Kratos import ResidualBasedNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import BlockBuilderAndSolver as block_builder_and_solver
else:
    raise Exception("Distributed run requires TrilinosApplication")


def CreateSolver(main_model_part, custom_settings):
    return IncompressiblePotentialFlowSolver(main_model_part, custom_settings)


class IncompressiblePotentialFlowSolver(PythonSolver):
    def __init__(self, model, settings):
        self._validate_settings_in_baseclass = True  # To be removed eventually

        super(IncompressiblePotentialFlowSolver,
              self).__init__(model, settings)

        # self.mesh_moving = self.settings["mesh_moving"].GetBool()

        # TODO: Implement stuff for mesh_moving
        self.EpetraCommunicator = None

    def GetDefaultSettings(self):
        return Kratos.Parameters(r'''{
            "velocity_potential_flow_settings":
            {
                "linear_solver_settings": {
                    "solver_type"  : "amgcl"
                },
                "reform_dofs_at_each_step": true,
                "move_mesh_strategy": 0,
                "move_mesh_flag": false,
                "compute_reactions": false
            },
            "pressure_potential_flow_settings":
            {
                "linear_solver_settings": {
                    "solver_type"  : "amgcl"
                },
                "reform_dofs_at_each_step": true,
                "move_mesh_strategy": 0,
                "move_mesh_flag": false,
                "compute_reactions": false
            },
            "echo_level" : 0
        }''')

    def AddVariables(self):
        self.fluid_model_part.AddNodalSolutionStepVariable(
            KratosRANS.VELOCITY_POTENTIAL)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            KratosRANS.PRESSURE_POTENTIAL)

        Kratos.Logger.PrintInfo(
            self.__class__.__name__,
            "Added potential flow initialization variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(KratosRANS.VELOCITY_POTENTIAL,
                                      self.fluid_model_part)
        Kratos.VariableUtils().AddDof(KratosRANS.PRESSURE_POTENTIAL,
                                      self.fluid_model_part)

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Added potential flow initialization dofs.")

    def PrepareModelPart(self):
        self.domain_size = self.fluid_model_part.ProcessInfo[
            Kratos.DOMAIN_SIZE]

        element_suffix = str(
            self.domain_size) + "D" + str(self.domain_size + 1) + "N"
        condition_suffix = str(self.domain_size) + "D" + str(
            self.domain_size) + "N"

        original_condition_name = self.GetFluidVelocityPressureConditionName()

        element_name = "RansIncompressibleVelocityPotentialElement" + element_suffix
        condition_name = "RansIncompressibleVelocityPotentialCondition" + condition_suffix
        self.velocity_model_part = CreateDuplicateModelPart(
            self.fluid_model_part, "IncompressiblePotentialFlow_Velocity",
            element_name, condition_name, original_condition_name)

        element_name = "RansPressurePotentialElement" + element_suffix
        condition_name = "RansIncompressiblePressureCondition" + condition_suffix
        self.pressure_model_part = CreateDuplicateModelPart(
            self.fluid_model_part, "IncompressiblePotentialFlow_Pressure",
            element_name, condition_name, original_condition_name)

        RansVariableUtilities.InitializeDuplicatedModelPart(
            self.fluid_model_part, self.velocity_model_part)
        RansVariableUtilities.InitializeDuplicatedModelPart(
            self.fluid_model_part, self.pressure_model_part)

        Kratos.Logger.PrintInfo(
            self.__class__.__name__,
            "Created potential flow initialization model parts.")

    def Initialize(self):
        # solving for velocity potential
        solver_settings = self.settings["velocity_potential_flow_settings"]
        linear_solver = linear_solver_factory.ConstructSolver(
            solver_settings["linear_solver_settings"])
        builder_and_solver = block_builder_and_solver(linear_solver,
                                                      self.EpetraCommunicator)
        convergence_criteria = residual_criteria(1e-12, 1e-12)
        self.velocity_strategy = newton_raphson_strategy(
            self.velocity_model_part, scheme(), linear_solver,
            convergence_criteria, builder_and_solver, 2,
            solver_settings["compute_reactions"].GetBool(),
            solver_settings["reform_dofs_at_each_step"].GetBool(),
            solver_settings["move_mesh_flag"].GetBool())

        # solving for pressure
        solver_settings = self.settings["pressure_potential_flow_settings"]
        linear_solver = linear_solver_factory.ConstructSolver(
            solver_settings["linear_solver_settings"])
        builder_and_solver = block_builder_and_solver(linear_solver,
                                                      self.EpetraCommunicator)
        convergence_criteria = residual_criteria(1e-12, 1e-12)
        self.pressure_strategy = newton_raphson_strategy(
            self.pressure_model_part, scheme(), linear_solver,
            convergence_criteria, builder_and_solver, 2,
            solver_settings["compute_reactions"].GetBool(),
            solver_settings["reform_dofs_at_each_step"].GetBool(),
            solver_settings["move_mesh_flag"].GetBool())

        self.velocity_strategy.Initialize()
        self.pressure_strategy.Initialize()

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Initialized potential flow solver.")

    def SetCommunicator(self, epetra_communicator):
        self.EpetraCommunicator = epetra_communicator

    def InitializeSolutionStep(self):
        if (not hasattr(self, "is_initialized")):
            self.is_initialized = True

            # solve for velocity potential
            RansVariableUtilities.FixFlaggedDofs(self.velocity_model_part,
                                                 KratosRANS.VELOCITY_POTENTIAL,
                                                 Kratos.OUTLET)
            RansVariableUtilities.FixFlaggedDofs(self.pressure_model_part,
                                                 KratosRANS.PRESSURE_POTENTIAL,
                                                 Kratos.OUTLET)
            RansVariableUtilities.CopyFlaggedVariableToNonHistorical(
                self.velocity_model_part, Kratos.VELOCITY, Kratos.INLET)

            self.velocity_strategy.InitializeSolutionStep()
            self.pressure_strategy.InitializeSolutionStep()

            self.velocity_strategy.Predict()
            self.velocity_strategy.SolveSolutionStep()

            # extrapolate gauss point velocities to nodal velocities
            extrapolation_settings = Kratos.Parameters('''
            {
                "echo_level"                 : 0,
                "area_average"               : true,
                "average_variable"           : "NODAL_AREA",
                "list_of_variables"          : ["VELOCITY"],
                "extrapolate_non_historical" : false
            }''')
            extrapolation_process(self.velocity_model_part,
                                  extrapolation_settings).Execute()

            # take back the original inlet velocities
            RansVariableUtilities.CopyFlaggedVariableFromNonHistorical(
                self.velocity_model_part, Kratos.VELOCITY, Kratos.INLET)

            RansVariableUtilities.CalculateMagnitudeSquareFor3DVariable(
                self.velocity_model_part, Kratos.VELOCITY,
                KratosRANS.VELOCITY_POTENTIAL)

            self.pressure_strategy.Predict()
            self.pressure_strategy.SolveSolutionStep()

            variable_utils = Kratos.VariableUtils()
            variable_utils.CopyModelPartNodalVar(KratosRANS.PRESSURE_POTENTIAL,
                                                 Kratos.PRESSURE,
                                                 self.pressure_model_part,
                                                 self.pressure_model_part, 0)

            self.pressure_strategy.FinalizeSolutionStep()
            self.velocity_strategy.FinalizeSolutionStep()

            self.velocity_strategy.Clear()
            self.pressure_strategy.Clear()

            Kratos.Logger.PrintInfo(
                self.__class__.__name__,
                "Initialized domain with potential flow solution.")

    def GetFluidVelocityPressureConditionName(self):
        return ""
