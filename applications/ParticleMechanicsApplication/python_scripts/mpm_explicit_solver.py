from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Importing the base class
from KratosMultiphysics.ParticleMechanicsApplication.mpm_solver import MPMSolver

def CreateSolver(model, custom_settings):
    return MPMExplicitSolver(model, custom_settings)

class MPMExplicitSolver(PythonSolver):

    def __init__(self, model, custom_settings):
        super(MPMExplicitSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[MPMExplicitSolver]:: ", "Construction is finished.")


    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "scheme_type"   : "forward_euler",
            "stress_update" : "USL"
        }""")
        this_defaults.AddMissingParameters(super(MPMExplicitSolver, cls).GetDefaultSettings())
        return this_defaults

    def AddVariables(self):
        super(MPMExplicitSolver, self).AddVariables()
        self._AddDynamicVariables(self.grid_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[MPMExplicitSolver]:: ", "Variables are all added.")

    ### Protected functions ###

    def _CreateSolutionScheme(self):
        grid_model_part = self.GetGridModelPart()
        domain_size = self._GetDomainSize()
        block_size  = domain_size
        if (self.settings["pressure_dofs"].GetBool()):
            block_size += 1

        # Setting the time integration schemes
        scheme_type = self.settings["scheme_type"].GetString()
        isCentralDifference = False
        StressUpdateOption = 0

        if(scheme_type == "forward_euler"):
            stress_update = self.settings["stress_update"].GetString() #0 = USF, 1 = USL, 2 = MUSL
            if(stress_update == "USF"):
                StressUpdateOption = 0
            elif(stress_update == "USL"):
                StressUpdateOption = 1
            elif(stress_update == "MUSL"):
                StressUpdateOption = 2
            else:
                err_msg = "The requested stress update \"" + stress_update + "\" is not available!\n"
                err_msg += "Available options are: \"USF\", \"USL\",\"MUSL\""
        elif(scheme_type == "central_difference"):
            isCentralDifference = True
        else:
            err_msg = "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"forward_euler\", \"central_difference\""
            raise Exception(err_msg)

        is_dynamic = self._IsDynamic()

        return KratosParticle.MPMExplicitScheme( grid_model_part,
                                                 MaximumDeltaTime,
                                                 DeltaTimeFraction,
                                                 DeltaTimePredictionLevel,
                                                 StressUpdateOption,
                                                 isCentralDifference)


    def _IsDynamic(self):
        return True