import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return UpdateThermalModelPartProcess(Model, settings["Parameters"])

# All the processes python should be derived from "Process"
class UpdateThermalModelPartProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)
        #print("****************************************")#
        #print("__init__ of UpdateThermalModelPartProcess")
        #print(settings)#
        #print("****************************************")#
        self.origin_model_part = Model[settings["input_model_part"].GetString()]
        #print(input_model_part)#
        self.destination_model_part = Model[settings["output_model_part"].GetString()]
        #print(output_model_part)#


        #self.components_process_list = []

        #if settings["active"][0].GetBool() == True:
        self.params = KratosMultiphysics.Parameters("{}")

        self.params.AddValue("domain_size",settings["domain_size"])
        self.params.AddEmptyValue("two").SetString("EulerianConvDiff2D")
        self.params.AddEmptyValue("three").SetString("EulerianConvDiff3D")
        #self.params.AddValue("two","EulerianConvDiff2D")
        #self.params.AddValue("three","EulerianConvDiff3D")
        #print(self.params)
        #print(1)
        #params.AddValue("2D_reference_element","EulerianConvDiff2D")
        #params.AddValue("3D_reference_element","EulerianConvDiff3D")

    #def ExecuteInitialize(self):
    #    pass
        #for component in self.components_process_list:
        #    component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        #print("start")

        KratosPfemFluid.UpdateThermalModelPartProcess(self.origin_model_part, self.destination_model_part, self.params).ExecuteInitializeSolutionStep()
            #"EulerianConvDiff2D", "EulerianConvDiff3D", self.params).ExecuteInitializeSolutionStep()
        #KratosPfemFluid.UpdateThermalModelPartProcess(input_model_part, output_model_part).ExecuteInitializeSolutionStep()

        # Next steps:
        # 1) call condiff and check id there are inverted elements and conditions
        # 2) call convdiff and rebuild the thermal computing domain
        # these last steps could be problematic as the solver is not defined inside this process
        # anyway, they could be easily moved to the coupled solver

        # TODO: ic c++ there must be a loop over the nodes to update the mesh velocity!