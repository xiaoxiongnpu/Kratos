from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

def CreateSolver(model, custom_settings):
    return FluidDEMSolver(model, custom_settings)

class FluidDEMSolver(FluidSolver):

    ## FluidDEMSolver specific methods.

    def _ReplaceElementsAndConditions(self):
        ## Get number of nodes and domain size
        elem_num_nodes = self._GetElementNumNodes()
        cond_num_nodes = self._GetConditionNumNodes()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## If there are no elements and/or conditions, default to triangles/tetra meshes to avoid breaking the ReplaceElementsAndConditionsProcess
        ## This only affects the input name (if there are no elements or conditions to replace, nothing is replaced).
        if elem_num_nodes == 0:
            elem_num_nodes = domain_size + 1
        if cond_num_nodes == 0:
            cond_num_nodes = domain_size

        ## Complete the element name
        if (self.element_name is not None):
            new_elem_name = self.element_name + str(int(domain_size)) + "D" + str(int(elem_num_nodes)) + "N"
        else:
            raise Exception("There is no element name. Define the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if (self.condition_name is not None):
            new_cond_name = self.condition_name + str(int(domain_size)) + "D" + str(int(cond_num_nodes)) + "N"
        else:
            raise Exception("There is no condition name. Define the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        #self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""{}""")
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()