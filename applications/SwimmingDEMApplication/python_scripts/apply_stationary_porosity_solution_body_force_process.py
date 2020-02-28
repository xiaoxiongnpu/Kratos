# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication

from KratosMultiphysics import Vector
from importlib import import_module
from apply_custom_body_force_process import ApplyCustomBodyForceProcess

import numpy as np

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyStationaryPorositySolutionBodyForceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyStationaryPorositySolutionBodyForceProcess(ApplyCustomBodyForceProcess):
    def __init__(self, model, settings ):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {},
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        super(ApplyCasasSolutionBodyForceProcess, self).__init__(model, settings)

        self.x = np.array([node.X for node in self.model_part.Nodes])
        self.y = np.array([node.Y for node in self.model_part.Nodes])
        self.z = np.array([node.Z for node in self.model_part.Nodes])

        self.center_x = self.settings["benchmark_parameters"]["center_x1"].GetDouble()
        self.center_y = self.settings["benchmark_parameters"]["center_x2"].GetDouble()

    def ExecuteBeforeSolutionLoop(self):
        current_time = 0.0

        value_v  = np.array([self.benchmark.Velocity(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        bf_value = np.array([self.benchmark.BodyForce(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        p_value  = np.array([self.benchmark.Pressure(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])

        self.fluid_fraction_list = np.array([self.benchmark.alpha(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])

        iterator = 0
        for node in self.model_part.Nodes:
            vel_value = Vector(list(value_v[iterator]))
            b_value = Vector(list(bf_value[iterator]))
            press_value = p_value[iterator]
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, vel_value)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE, b_value)
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_VELOCITY, vel_value)
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_PRESSURE, press_value)
            #Fix porosity field due to it's not time dependant
            fluid_fraction = self.fluid_fraction_list[iterator]
            fluid_fraction_rate = self.fluid_fraction_rate_list[iterator]
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION, fluid_fraction)
            node.Fix(KratosMultiphysics.FLUID_FRACTION)
            #Fix Pressure in a point
            if node.X == 1.0 and node.Y == 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, press_value)
                node.Fix(KratosMultiphysics.PRESSURE)
            iterator += 1

    def ExecuteFinalizeSolutionStep(self):
        pass

    def _SetBodyForce(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.bf_value = np.array([self.benchmark.BodyForce(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        self.value_v  = np.array([self.benchmark.Velocity(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])

        iterator = 0
        for node in self.model_part.Nodes:
            #Set BodyForce field
            bodf_value = Vector(list(self.bf_value[iterator]))
            vel_value = Vector(list(self.value_v[iterator]))
            node.SetSolutionStepValue(self.variable, bodf_value)
            #Update BCs
            if node.X == max(self.x) or node.Y == max(self.y) or node.X == min(self.x) or node.Y == min(self.y):
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, vel_value[0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, vel_value[1])
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
            iterator += 1