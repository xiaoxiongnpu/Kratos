import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")
    return OutputResultsProcess(Model, settings["Parameters"])

# all the processes python processes should be derived from "python_process"

class OutputResultsProcess(KratosMultiphysics.Process):
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

        REACTION_file = open("REACTION_result_file_for_" + str(self.model_part_name) + ".mdpa", "w+")
        DISPLACEMENT_file = open("DISPLACEMENT_result_file_for_" + str(self.model_part_name) + ".mdpa", "w+")
        nodes_list = []

        for node in self.body_model_part.Nodes:
            var_x = node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
            var_y = node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y)
            norm_var = node.GetSolutionStepValue(KratosMultiphysics.REACTION).norm_2()
            nodes_list.append(node.Id)

            disp_x = node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_X)
            disp_y = node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y)
            norm_disp = node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT).norm_2()

            REACTION_file.write(str(node.Id) + ' ' + str(var_x) + ' ' + str(var_y) + '  norm=  ' + str(norm_var) + '\n')
            DISPLACEMENT_file.write(str(node.Id) + ' ' + str(disp_x) + ' ' + str(disp_y) + '  norm=  ' + str(norm_disp) + '\n')

        REACTION_file.close()
        DISPLACEMENT_file.close()

        sum_reaction = open("Sum_REACTION_file_for_" + str(self.model_part_name) + ".mdpa", "w+")
        sum_disp = open("Sum_DISPLACEMENT_file_for_" + str(self.model_part_name) + ".mdpa", "w+")

        sum_reaction.write('The Total Reaction for ' + str(self.model_part_name) + ' is:\n')
        sum_reaction.write('TOTAL_REACTION_X    TOTAL_REACTION_Y    |REACTION|\n')
        total_reaction = KratosMultiphysics.VariableUtils().SumHistoricalNodeVectorVariable(KratosMultiphysics.REACTION, self.body_model_part, 0)
        norm_total_reaction = total_reaction.norm_2()
        sum_reaction.write(str(total_reaction[0]) + ' ' + str(total_reaction[1]) + ' ' + str(norm_total_reaction))

        sum_disp.write('The Total Displacement for ' +  str(self.model_part_name) + ' is:\n')
        sum_disp.write('TOTAL_DISPLACEMENT_X    TOTAL_DISPLACEMENT_Y    |DISPLACEMENT|\n')
        total_disp = KratosMultiphysics.VariableUtils().SumHistoricalNodeVectorVariable(KratosMultiphysics.MESH_DISPLACEMENT, self.body_model_part, 0)
        norm_total_disp = total_disp.norm_2()
        sum_disp.write(str(total_disp[0]) + ' ' + str(total_disp[1]) + ' ' + str(norm_total_disp))