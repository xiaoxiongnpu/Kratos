# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication
from importlib import import_module
import numpy as np
from KratosMultiphysics import Vector

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCustomBodyForceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyCustomBodyForceProcess(KratosMultiphysics.Process):
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

        self.model_part = model[self.settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())

        self.x = np.array([node.X for node in self.model_part.Nodes])
        self.y = np.array([node.Y for node in self.model_part.Nodes])
        self.z = np.array([node.Z for node in self.model_part.Nodes])

        self.center_x = self.settings["benchmark_parameters"]["center_x1"].GetDouble()
        self.center_y = self.settings["benchmark_parameters"]["center_x2"].GetDouble()

        benchmark_module = import_module(self.settings["benchmark_name"].GetString())
        self.benchmark = benchmark_module.CreateManufacturedSolution(self.settings["benchmark_parameters"])

        self.compute_error = self.settings["compute_nodal_error"].GetBool()

        self.print_output = self.settings["print_convergence_output"].GetBool()
        if self.print_output:
            from custom_body_force.hdf5_output_tool import Hdf5OutputTool
            self.output_process = Hdf5OutputTool(model, self.settings["output_parameters"])

    def ExecuteBeforeSolutionLoop(self):
        current_time = 0.0

        value_v = np.array([self.benchmark.Velocity(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        bf_value = np.array([self.benchmark.BodyForce(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        p_value = np.array([self.benchmark.Pressure(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])

        iterator = 0
        for node in self.model_part.Nodes:
            vel_value = Vector(list(value_v[iterator]))
            b_value = Vector(list(bf_value[iterator]))
            press_value = p_value[iterator]
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, vel_value)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE, b_value)
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_VELOCITY, vel_value)
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_PRESSURE, press_value)
            iterator += 1

    def ExecuteInitializeSolutionStep(self):
        self._SetBodyForce()
        self._SetFluidFractionField()

    def ExecuteFinalizeSolutionStep(self):
        if self.compute_error:
            self._ComputeVelocityError()
            self._ComputePressureError()

    def ExecuteBeforeOutputStep(self):
        self._ComputeVelocityBenchmark()

    def ExecuteFinalize(self):
        if self.compute_error:
            rel_err = self._SumNodalError()
            if self.print_output:
                self.output_process.WriteBodyForceAttributes(self.settings["benchmark_parameters"])
                self.output_process.WriteAverageRelativeError(rel_err)

    def _SetBodyForce(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.bf_value = np.array([self.benchmark.BodyForce(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        self.value_v = np.array([self.benchmark.Velocity(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        self.value_p = np.array([self.benchmark.Pressure(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])

        radius = np.array([np.sqrt((node.X-self.center_x)**2 + (node.Y-self.center_y)**2) for node in self.model_part.Nodes])
            
        iterator = 0
        for node in self.model_part.Nodes:
            bodf_value = Vector(list(self.bf_value[iterator]))
            vel_value = Vector(list(self.value_v[iterator]))
            press_value = self.value_p[iterator]
            node.SetSolutionStepValue(self.variable, bodf_value)
            if node.X == 1.0 and node.Y == 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, press_value)
                node.Fix(KratosMultiphysics.PRESSURE)
            if np.sqrt((node.X-self.center_x)**2 + (node.Y-self.center_y)**2) - min(radius) < 1e-3 or -np.sqrt((node.X-self.center_x)**2 + (node.Y-self.center_y)**2) + max(radius) < 1e-3 or node.X == max(self.x) or node.Y == min(self.y) or node.X - node.Y < 1e-3:
            #if node.X == max(self.x) or node.Y == max(self.y) or node.X == min(self.x) or node.Y == min(self.y):
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, vel_value[0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, vel_value[1])
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)

            iterator += 1

    def _SetFluidFractionField(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        fluid_fraction_list = np.array([self.benchmark.alpha(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])
        fluid_fraction_rate_list = np.array([self.benchmark.dalphat(current_time, x, y, z) for x, y, z in zip(self.x, self.y, self.z)])

        iterator = 0
        for node in self.model_part.Nodes:
            fluid_fraction = fluid_fraction_list[iterator]
            fluid_fraction_rate = fluid_fraction_rate_list[iterator]
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION, fluid_fraction)
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION_RATE, fluid_fraction_rate)
            iterator += 1


    def _ComputeVelocityError(self):
        epsilon = 1e-16

        iterator = 0
        for node in self.model_part.Nodes:
            fem_vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            exact_vel = Vector(list(self.value_v[iterator]))
            fem_vel_modulus = (fem_vel[0]**2 + fem_vel[1]**2 )**0.5
            exact_vel_modulus = (exact_vel[0]**2 + exact_vel[1]**2)**0.5
            error = abs(fem_vel_modulus - exact_vel_modulus)
            error = error / abs(exact_vel_modulus + epsilon)
            node.SetValue(KratosMultiphysics.NODAL_ERROR, error)
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_VELOCITY, exact_vel)
            iterator += 1

    def _ComputePressureError(self):
        epsilon = 1e-16

        iterator = 0
        for node in self.model_part.Nodes:
            fem_press = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
            exact_press = self.value_p[iterator]
            error = abs(fem_press - exact_press)
            error = error / abs(exact_press + epsilon)
            node.SetSolutionStepValue(KratosMultiphysics.SwimmingDEMApplication.EXACT_PRESSURE, exact_press)
            iterator += 1

    def _CopyVelocityAsNonHistorical(self):
        for node in self.model_part.Nodes:
            vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            node.SetValue(KratosMultiphysics.VELOCITY, vel)

    def _ComputeVelocityBenchmark(self):
        pass

    def _SumNodalError(self):
        err_sum = 0.0
        for node in self.model_part.Nodes:
            err_sum = err_sum + node.GetValue(KratosMultiphysics.NODAL_ERROR)
        rel_err = err_sum / self.model_part.Nodes.__len__()
        KratosMultiphysics.Logger.PrintInfo("Benchmark", "The nodal error average is : ", rel_err)
        return rel_err
