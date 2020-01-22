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

    # =====================Multi-Objective-Load 1=================================
    # Read parameters (Optimization)
    with open("optimization_parameters.json",'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()

    # =======================================================================================================
    # Define external analyzer
    # =======================================================================================================
    class CustomAnalyzer(AnalyzerBaseClass, AnalysisStage):
        # --------------------------------------------------------------------------------------------------
        def __init__( self , model):
            self.model = model

        def DirForOptimization(self, OptiFile, itr , FolderName):

            shutil.copy(OptiFile, str(itr)+ "_ITR" + ".post.bin")            
            dir = "/home/jey/Desktop/CodeWorld/KRATOS/JeyExamples/MultiObjectiveThickShell/ITR_Results/"

            for file in glob.glob("*_ITR.post.bin"):
                dst = dir + "" + file.replace(".post.bin", FolderName)
                os.mkdir(dst)
                print("CUSTOM ANALYZER CALLED !!")
                shutil.move(file, dst)
        # --------------------------------------------------------------------------------------------------
        def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):

            self.jey_response_settings = communicator.jey_response_settings
            
            self.response_settings = []            
            
            for response in self.jey_response_settings:
                if response["identifier"].GetString() == "strain_energy_1":
                    self.response_settings = response["response_settings"]
                    self.identifier = "strain_energy_1"

                if response["identifier"].GetString() == "strain_energy_2":
                    self.response_settings = response["response_settings"]
                    self.identifier = "strain_energy_2"
            
            sr = structural_response.StrainEnergyResponseFunction(self.identifier, self.response_settings, model)
            sr.RunCalculation(calculate_gradient=True)        

            # if self.identifier == "strain_energy_1":
            #     self.DirForOptimization("plate1.post.bin", optimization_iteration, "_ExternalAnalyzer_ID_1")
            # else:
            #     self.DirForOptimization("plate2.post.bin", optimization_iteration, "_ExternalAnalyzer_ID_2")

            # #sr1.Initialize()                   ERROR::ImportModelPart()

            # self.primal_analysis = sr1.primal_analysis
            # self.primal_model_part = sr1.primal_model_part
            # self.response_function_utility = sr1.response_function_utility
            
            # #Structural Response: Initialize()
            # self.Jey_Initialize() #instead of primal_analysis.Initalize()  ERROR :: _ImportModelPart()
            # self.response_function_utility.Initialize()
            
            # #Structural Response: InitializeSolutionStep()
            # #self.Jey_InitializeSolutionStep()
            # sr1.InitializeSolutionStep()

            #Structural Response: Calculate Values
            if communicator.isRequestingValueOf(self.identifier):
                #communicator.reportValue("strain_energy_2", self.__CalculateValue(current_design))
                communicator.reportValue(self.identifier, sr.GetValue())

            #Structural Response: Calculate Gradients
            if communicator.isRequestingGradientOf(self.identifier):
                #communicator.reportGradient("strain_energy_2", self.__CalculateGradient(current_design))
                communicator.reportGradient(self.identifier, sr.GetShapeGradient())
            
            # #Structural Response: FinalizeSolutionStep()
            # #self.Jey_FinalizeSolutionStep()
            # sr1.FinalizeSolutionStep()

            # #Structural Response: Finalize()
            # #self.Jey_Finalize()
            # sr1.Finalize()
        #---------------------------------------------------------------------------
        def Jey_Initialize(self):
            """This function initializes the AnalysisStage
            Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
            This function has to be implemented in deriving classes!
            """
            print("!@$#%^&*^$$:::NOW I'M CALLED, YOUR RUN INITALIZED !")

            self.primal_analysis._GetSolver().PrepareModelPart()
            self.primal_analysis._GetSolver().AddDofs()

            self.primal_analysis.ModifyInitialProperties()
            self.primal_analysis.ModifyInitialGeometry()

            ##here we initialize user-provided processes
            self.primal_analysis._CreateListOfProcesses() # Function MODIFIED from protected to private in analysis_stage.py
            for process in self.primal_analysis._GetListOfProcesses():
                process.ExecuteInitialize()

            self.primal_analysis._GetSolver().Initialize()
            self.primal_analysis.Check()

            self.primal_analysis.ModifyAfterSolverInitialize()

            for process in self.primal_analysis._GetListOfProcesses():
                process.ExecuteBeforeSolutionLoop()
            
            with open(self.response_settings2["primal_settings"].GetString()) as parameters_file:
                self.project_parameters = Parameters(parameters_file.read())

            ## Stepping and time settings
            self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

            #ERROR:: CreateSolver Exceptions
            # if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            #     self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
            # else:
            #     self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()

            #Echo level omitted

            KM.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")
        # --------------------------------------------------------------------------
        def Jey_InitializeSolutionStep(self):
            #primal_analysis -> StructuralMechanicsAnalysis -> NO time() or InitializeSolutionStep()
            # self.primal_analysis.time = self.primal_analysis._GetSolver().AdvanceInTime(self.primal_analysis.time)
            # self.primal_analysis.InitializeSolutionStep()
            print("InitalizeSolutionStep() Not Working!")
        #---------------------------------------------------------------------------
        def __CalculateValue( self, current_design ):
            
            # for node in model["Structure.PointLoad3D_Load"].Nodes:
            #     node.Fix(POINT_LOAD)
            #     node.SetSolutionStepValue(POINT_LOAD, 0, 0)
            #     node.Fix(POINT_LOAD_Y)
            #     node.SetSolutionStepValue(POINT_LOAD_Y, 0, -1)
            #     node.Fix(POINT_LOAD_Z)
            #     node.SetSolutionStepValue(POINT_LOAD_Z, 0, 1)
            #print("JEYARAMANAN :::::",self.primal_analysis)
            #self.primal_analysis.Run()

            Logger.PrintInfo("StrainEnergyResponse", "Starting primal analysis for response", self.identifier2)
            startTime = timer.time()
            self.primal_analysis._GetSolver().Predict()
            self.primal_analysis._GetSolver().SolveSolutionStep()
            Logger.PrintInfo("StrainEnergyResponse", "Time needed for solving the primal analysis",round(timer.time() - startTime,2),"s")
            
            startTime = timer.time()
            value = self.response_function_utility.CalculateValue()
            self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE] = value
            Logger.PrintInfo("StrainEnergyResponse", "Time needed for calculating the response value",round(timer.time() - startTime,2),"s")

            return self.primal_model_part.ProcessInfo[StructuralMechanicsApplication.RESPONSE_VALUE]
        # --------------------------------------------------------------------------
        def __CalculateGradient( self, current_design ):
            
            Logger.PrintInfo("StrainEnergyResponse", "Starting gradient calculation for response", self.identifier2)
            startTime = timer.time()
            self.response_function_utility.CalculateGradient()
            Logger.PrintInfo("StrainEnergyResponse", "Time needed for calculating gradients",round(timer.time() - startTime,2),"s")
            
            gradient = {}
            for node in self.primal_model_part.Nodes:
                gradient[node.Id] = node.GetSolutionStepValue(KM.SHAPE_SENSITIVITY)
            return gradient
        #---------------------------------------------------------------------------
        def Jey_FinalizeSolutionStep(self):
            self.primal_analysis.FinalizeSolutionStep()     #This function is not present structural_mechanics_analaysis.py
            self.primal_analysis.OutputSolutionStep()   
        #---------------------------------------------------------------------------    
        def Jey_Finalize(self):
            self.primal_analysis.Finalize() #This function is not present structural_mechanics_analaysis.py
    # =======================================================================================================

    # Create optimizer and perform optimization
    #optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model, CustomAnalyzer(model))
    
    optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model)
    optimizer.Optimize()