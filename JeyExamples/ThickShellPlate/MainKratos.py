from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.KratosUnittest import TestCase
import KratosMultiphysics.kratos_utilities as kratos_utilities
import os, csv
"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

if __name__ == "__main__":

    #Read parameter (Structural)
    # with open("ProjectParameters1.json",'r') as parameter_file:
    #     parameters = KM.Parameters(parameter_file.read())

    # model = KM.Model()
    # simulation = StructuralMechanicsAnalysis(model,parameters)
    # simulation.Run()

    # print("###############################")
    # print("-------- Structural Analysis DONE-------------------")
    # print("###############################")
    
    #---------------------------------------------------------------
    print("###############################")
    print("--------LOAD CASE 1-------------------")
    print("###############################")

    # Read parameters (Optimization)
    with open("optimization_parameters1.json",'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    # Create optimizer and perform optimization
    optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model)
    optimizer.Optimize()

    # #---------------------------------------------------------------
    # print("###############################")
    # print("--------LOAD CASE 2-------------------")
    # print("###############################")
    # # Read parameters (Optimization)
    # with open("optimization_parameters2.json",'r') as parameter_file:
    #     parameters = KM.Parameters(parameter_file.read())

    # model = KM.Model()
    # # Create optimizer and perform optimization
    # optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model)
    # optimizer.Optimize()

    # #---------------------------------------------------------------
    # print("###############################")
    # print("--------LOAD CASE 3-------------------")
    # print("###############################")
    # # Read parameters (Optimization)
    # with open("optimization_parameters3.json",'r') as parameter_file:
    #     parameters = KM.Parameters(parameter_file.read())

    # model = KM.Model()
    # # Create optimizer and perform optimization
    # optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model)
    # optimizer.Optimize()

    #---------------------------------------------------------------
    # =====================Multi-Objective-Load 1=================================
    # Read parameters (Optimization)
    # with open("optimization_parameters4.json",'r') as parameter_file:
    #     parameters = KM.Parameters(parameter_file.read())

    # model = KM.Model()
    # # Create optimizer and perform optimization
    # optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model)
    # optimizer.Optimize()