import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")
    return OutputForcesProcess(Model, settings["Parameters"])

# all the processes python processes should be derived from "python_process"

class OutputForcesProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters(r'''{
            "model_part_name": ""
        }''')

        settings.ValidateAndAssignDefaults(default_parameters)

        self.body_model_part = Model[settings["model_part_name"].GetString()]
        self.model_part_name = settings["model_part_name"].GetString()

    def ExecuteFinalizeSolutionStep(self):
        self.Execute()

    def Execute(self):
        KratosMultiphysics.Logger.PrintInfo('Creating Output for ', self.model_part_name)

        LOAD_file = open("LOAD_result_file_for_" + str(self.model_part_name) + ".mdpa", "w+")
        nodes_list = []

        for node in self.body_model_part.Nodes:
            var_x = node.GetSolutionStepValue(SMA.LINE_LOAD_X)
            var_y = node.GetSolutionStepValue(SMA.LINE_LOAD_Y)
            norm_var = node.GetSolutionStepValue(SMA.LINE_LOAD).norm_2()
            nodes_list.append(node.Id)

            LOAD_file.write(str(node.Id) + ' ' + str(var_x) + ' ' + str(var_y) + '  norm=  ' + str(norm_var) + '\n')

        LOAD_file.close()

        sum_load = open("Sum_LOAD_file_for_" + str(self.model_part_name) + ".mdpa", "w+")

        sum_load.write('The Total Load for ' + str(self.model_part_name) + ' is:\n')
        sum_load.write('TOTAL_LOAD_X    TOTAL_LOAD_Y    |LOAD|\n')
        total_load = KratosMultiphysics.VariableUtils().SumHistoricalNodeVectorVariable(SMA.LINE_LOAD, self.body_model_part, 0)
        norm_total_load = total_load.norm_2()
        sum_load.write(str(total_load[0]) + ' ' + str(total_load[1]) + ' ' + str(norm_total_load))
