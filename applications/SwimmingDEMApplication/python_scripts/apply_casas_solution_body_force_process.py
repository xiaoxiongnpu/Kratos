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
    return ApplyCasasSolutionBodyForceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyCasasSolutionBodyForceProcess(ApplyCustomBodyForceProcess):
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
        self.radius   = np.array([np.sqrt((x - self.center_x)**2 + (y - self.center_y)**2) for x, y, z in zip(self.x, self.y, self.z)])

    def ExecuteBeforeSolutionLoop(self):
        current_time = 0.0

        value_v  = np.array([self.benchmark.Velocity(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        bf_value = np.array([self.benchmark.BodyForce(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        p_value  = np.array([self.benchmark.Pressure(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])

        iterator = 0
        for node in self.model_part.Nodes:
            vel_value = Vector(list(value_v[iterator]))
            b_value   = Vector(list(bf_value[iterator]))
            press_value = p_value[iterator]
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, vel_value)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE, b_value)
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_VELOCITY, vel_value)
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_PRESSURE, press_value)
            iterator += 1

    def ExecuteInitializeSolutionStep(self):
        self._SetBodyForceAndPorosityField()

    def ExecuteFinalizeSolutionStep(self):
        pass

    def _SetBodyForceAndPorosityField(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.bf_value = np.array([self.benchmark.BodyForce(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        self.value_v  = np.array([self.benchmark.Velocity(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        self.value_p  = np.array([self.benchmark.Pressure(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])

        self.fluid_fraction_list = np.array([self.benchmark.alpha(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        self.fluid_fraction_rate_list = np.array([self.benchmark.dalphat(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])

        iterator = 0
        for node in self.model_part.Nodes:
            #Set BodyForce field
            bodf_value = Vector(list(self.bf_value[iterator]))
            vel_value  = Vector(list(self.value_v[iterator]))
            press_value = self.value_p[iterator]
            node.SetSolutionStepValue(self.variable, bodf_value)
            #Update BCs
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_PRESSURE, press_value)
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_VELOCITY, vel_value)
            #Set Porosity field
            fluid_fraction = self.fluid_fraction_list[iterator]
            fluid_fraction_rate = self.fluid_fraction_rate_list[iterator]
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION, fluid_fraction)
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION_RATE, fluid_fraction_rate)
            iterator += 1