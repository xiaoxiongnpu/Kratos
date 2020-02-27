from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzer_base import AnalyzerBaseClass
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication import structural_response
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics import Parameters, Logger
from KratosMultiphysics.analysis_stage import AnalysisStage
import time as timer
import shutil
import glob, os

if __name__ == "__main__":

    # #Read parameter (VERTICAL LOAD [0, 0, -1])
    # with open("ProjectParameters3.json",'r') as parameter_file:
    #     parameters = KM.Parameters(parameter_file.read())

    # model = KM.Model()
    # simulation = StructuralMechanicsAnalysis(model,parameters)
    # simulation.Run()

    # #Read parameter (HORIZONTAL LOAD [1, 1, 0])
    # with open("ProjectParameters4.json",'r') as parameter_file:
    #     parameters = KM.Parameters(parameter_file.read())

    # model = KM.Model()
    # simulation = StructuralMechanicsAnalysis(model,parameters)
    # simulation.Run()

    #Read parameter (VERTICAL LOAD [0, 0, -1]) Optimized MDPA (plateR)
    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    simulation = StructuralMechanicsAnalysis(model,parameters)
    simulation.Run()
