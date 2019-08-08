from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as UnitTest

# Other imports
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.KratosUnittest import isclose as t_isclose

import os

# Check other applications dependency
hdf5_is_available = kratos_utilities.CheckIfApplicationsAvailable("HDF5Application")

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

class PotentialFlowTests(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    def test_Naca0012SmallAdjoint(self):
        if not hdf5_is_available:
            self.skipTest("Missing required application: HDF5Application")
        file_name = "naca0012_small_sensitivities"
        settings_file_name_primal = file_name + "_primal_parameters.json"
        settings_file_name_adjoint = file_name + "_adjoint_parameters.json"
        settings_file_name_adjoint_analytical = file_name + "_adjoint_analytical_parameters.json"
        work_folder = "naca0012_small_adjoint_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name_primal)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.327805503865, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.MOMENT_COEFFICIENT], -0.105810071870, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], 0.3230253050805644, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.32651526722535246, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.DRAG_COEFFICIENT_FAR_FIELD], 0.0036897206842046205, 0.0, 1e-9)
            self._runTest(settings_file_name_adjoint)
            self._runTest(settings_file_name_adjoint_analytical)

            for file_name in os.listdir(os.getcwd()):
                if file_name.endswith(".h5"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

    def test_Naca0012SmallCompressible(self):
        file_name = "naca0012_small_compressible"
        settings_file_name = file_name + "_parameters.json"
        work_folder = "naca0012_small_compressible_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT], 0.4968313580730855, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.MOMENT_COEFFICIENT], -0.1631792300021498, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_JUMP], 0.4876931961465126, 0.0, 1e-9)
            self._check_results(self.main_model_part.ProcessInfo[CPFApp.LIFT_COEFFICIENT_FAR_FIELD], 0.4953997676243705, 0.0, 1e-9)

            for file_name in os.listdir():
                if file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)
    def test_EmbeddedCircleNoWake(self):
        settings_file_name = "embedded_circle_no_wake_parameters.json"
        work_folder = "embedded_test"

        with WorkFolderScope(work_folder):
            self._runTest(settings_file_name)

    def _runTest(self,settings_file_name):
        model = KratosMultiphysics.Model()
        with open(settings_file_name,'r') as settings_file:
            settings = KratosMultiphysics.Parameters(settings_file.read())

        if self.print_output:
            if settings_file_name == "naca0012_small_sensitivities_adjoint_parameters.json":
                settings.AddValue("output_processes", KratosMultiphysics.Parameters(r'''{
                    "gid_output" : [{
                        "python_module" : "gid_output_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"  : "GiDOutputProcess",
                        "help"          : "This process writes postprocessing files for GiD",
                        "Parameters"    : {
                            "model_part_name"        : "MainModelPart",
                            "output_name"            : "naca0012_adjoint",
                            "postprocess_parameters" : {
                                "result_file_configuration" : {
                                    "gidpost_flags"       : {
                                        "GiDPostMode"           : "GiD_PostBinary",
                                        "WriteDeformedMeshFlag" : "WriteDeformed",
                                        "WriteConditionsFlag"   : "WriteConditions",
                                        "MultiFileFlag"         : "SingleFile"
                                    },
                                    "file_label"          : "step",
                                    "output_control_type" : "step",
                                    "output_frequency"    : 1,
                                    "body_output"         : true,
                                    "node_output"         : false,
                                    "skin_output"         : false,
                                    "plane_output"        : [],
                                    "nodal_results"       : ["SHAPE_SENSITIVITY","ADJOINT_VELOCITY_POTENTIAL", "ADJOINT_AUXILIARY_VELOCITY_POTENTIAL"],
                                    "nodal_nonhistorical_results": [],
                                    "elemental_conditional_flags_results": [],
                                    "gauss_point_results" : []
                                },
                                "point_data_configuration"  : []
                            }
                        }
                    }]
                }'''))
            else:
                settings.AddValue("output_processes", KratosMultiphysics.Parameters(r'''{
                    "gid_output" : [{
                        "python_module" : "gid_output_process",
                        "kratos_module" : "KratosMultiphysics",
                        "process_name"  : "GiDOutputProcess",
                        "help"          : "This process writes postprocessing files for GiD",
                        "Parameters"    : {
                            "model_part_name"        : "MainModelPart",
                            "output_name"            : "naca0012",
                            "postprocess_parameters" : {
                                "result_file_configuration" : {
                                    "gidpost_flags"       : {
                                        "GiDPostMode"           : "GiD_PostBinary",
                                        "WriteDeformedMeshFlag" : "WriteDeformed",
                                        "WriteConditionsFlag"   : "WriteConditions",
                                        "MultiFileFlag"         : "SingleFile"
                                    },
                                    "file_label"          : "step",
                                    "output_control_type" : "step",
                                    "output_frequency"    : 1,
                                    "body_output"         : true,
                                    "node_output"         : false,
                                    "skin_output"         : false,
                                    "plane_output"        : [],
                                    "nodal_results"       : ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"],
                                    "nodal_nonhistorical_results": ["TRAILING_EDGE","wAKE_DISTANCE"],
                                    "elemental_conditional_flags_results": ["STRUCTURE"],
                                    "gauss_point_results" : ["PRESSURE_COEFFICIENT","VELOCITY","VELOCITY_LOWER","PRESSURE_LOWER","WAKE","WAKE_ELEMENTAL_DISTANCES","KUTTA"]
                                },
                                "point_data_configuration"  : []
                            }
                        }
                    }]
                }'''))

        potential_flow_analysis = PotentialFlowAnalysis(model, settings)
        potential_flow_analysis.Run()
        self.main_model_part = model.GetModelPart(settings["solver_settings"]["model_part_name"].GetString())

    def _check_results(self, result, reference, rel_tol, abs_tol):
        isclosethis = t_isclose(result, reference, rel_tol, abs_tol)

        full_msg =  "Failed with following parameters:\n"
        full_msg += str(result) + " != " + str(reference) + ", rel_tol = "
        full_msg += str(rel_tol) + ", abs_tol = " + str(abs_tol)

        self.assertTrue(isclosethis, msg=full_msg)

if __name__ == '__main__':
    UnitTest.main()