# This class can be taken as an abstract template for derivation. It can be used
# for default passive behaviour.

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import numpy as np
from importlib import import_module
from . import recoverer

class ManufacturedSolutionRecoverer(recoverer.DerivativesRecoverer):
    def __init__(self, project_parameters, model_part):
        recoverer.DerivativesRecoverer.__init__(self, project_parameters, model_part)

class ManufacturedFluidFractionSolutionRecoverer(ManufacturedSolutionRecoverer):
    def __init__(self, project_parameters, model_part):
        super(ManufacturedFluidFractionSolutionRecoverer, self).__init__(project_parameters, model_part)
        self.model_part = model_part
        self.fluid_fraction_manufactured_solution = project_parameters["fluid_parameters"]["processes"]["boundary_conditions_process_list"][1]["Parameters"]["benchmark_name"]
        benchmark_module = import_module(self.fluid_fraction_manufactured_solution.GetString())
        self.settings = project_parameters["fluid_parameters"]["processes"]["boundary_conditions_process_list"][1]["Parameters"]["benchmark_parameters"]
        self.fluid_fraction_field = benchmark_module.CreateManufacturedSolution(self.settings)

    def RecoverFluidFractionGradient(self):

        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        x = np.array([node.X for node in self.model_part.Nodes])
        y = np.array([node.Y for node in self.model_part.Nodes])
        z = np.array([node.Z for node in self.model_part.Nodes])
        fluid_fraction_gradient_field = np.array([[self.fluid_fraction_field.alpha1(time, x, y, z), self.fluid_fraction_field.alpha2(time, x, y, z), self.fluid_fraction_field.alpha3(time, x, y, z)] for x, y, z in zip(x, y, z)])
        
        iterator = 0
        for node in self.model_part.Nodes:
            fluid_fraction_gradient_list = fluid_fraction_gradient_field[iterator]
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION_GRADIENT_X, fluid_fraction_gradient_list[0])
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION_GRADIENT_Y, fluid_fraction_gradient_list[1])
            node.SetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION_GRADIENT_Z, fluid_fraction_gradient_list[2])
            iterator += 1