from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.MappingApplication import empire_mapper_wrapper
import ctypes as ctp


class TestEmpireMapperWrapperHelpers(KratosUnittest.TestCase):
    def setUp(self):
        self.test_model = KM.Model()
        self.model_part = self.test_model.CreateModelPart("dummy")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.FORCE)

        for i in range(5):
            node = self.model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # coordinates do not matter for this test
            node.SetSolutionStepValue(KM.PRESSURE, NodeScalarHistValue(node.Id))
            node.SetSolutionStepValue(KM.FORCE, NodeVectorHistValue(node.Id))
            node.SetValue(KM.TEMPERATURE, NodeScalarNonHistValue(node.Id))
            node.SetValue(KM.VELOCITY, NodeVectorNonHistValue(node.Id))


    def test_KratosFieldToCArray_scalar_hist(self):
        self.__Check_KratosFieldToCArray(is_scalar=True, historical=True)

    def test_KratosFieldToCArray_scalar_non_hist(self):
        self.__Check_KratosFieldToCArray(is_scalar=True, historical=False)

    def test_KratosFieldToCArray_vector_hist(self):
        self.__Check_KratosFieldToCArray(is_scalar=False)

    def test_KratosFieldToCArray_vector_non_hist(self):
        self.__Check_KratosFieldToCArray(is_scalar=False, historical=False)


    def test_CArrayToKratosField_scalar_hist(self):
        self.__Check_CArrayToKratosField(is_scalar=True)

    def test_CArrayToKratosField_scalar_non_hist(self):
        self.__Check_CArrayToKratosField(is_scalar=True, historical=False)

    def test_CArrayToKratosField_scalar_swap_sign_hist(self):
        self.__Check_CArrayToKratosField(is_scalar=True, swap_sign=True)

    def test_CArrayToKratosField_scalar_add_values_hist(self):
        self.__Check_CArrayToKratosField(is_scalar=True, add_values=True)

    def test_CArrayToKratosField_scalar_swap_sign_add_values_hist(self):
        self.__Check_CArrayToKratosField(is_scalar=True, add_values=True, swap_sign=True)

    def test_CArrayToKratosField_vector_hist(self):
        self.__Check_CArrayToKratosField(is_scalar=False)

    def test_CArrayToKratosField_vector_non_hist(self):
        self.__Check_CArrayToKratosField(is_scalar=False, historical=False)

    def test_CArrayToKratosField_vector_swap_sign_hist(self):
        self.__Check_CArrayToKratosField(is_scalar=False, swap_sign=True)

    def test_CArrayToKratosField_vector_add_values_hist(self):
        self.__Check_CArrayToKratosField(is_scalar=False, add_values=True)

    def test_CArrayToKratosField_vector_swap_sign_add_values_hist(self):
        self.__Check_CArrayToKratosField(is_scalar=False, add_values=True, swap_sign=True)


    def __Check_KratosFieldToCArray(self, is_scalar, historical=True):
        if is_scalar:
            if historical:
                variable = KM.PRESSURE
                fct_ptr = NodeScalarHistValue
            else:
                variable = KM.TEMPERATURE
                fct_ptr = NodeScalarNonHistValue
        else:
            if historical:
                variable = KM.FORCE
                fct_ptr = NodeVectorHistValue
            else:
                variable = KM.VELOCITY
                fct_ptr = NodeVectorNonHistValue

        c_array = empire_mapper_wrapper.KratosFieldToCArray(self.model_part.Nodes, variable, historical)

        # checking the values
        if is_scalar:
            for i, node in enumerate(self.model_part.Nodes):
                self.assertAlmostEqual(c_array[i], fct_ptr(node.Id))
        else:
            for i_node, node in enumerate(self.model_part.Nodes):
                node_val = fct_ptr(node.Id)
                for i in range(3):
                    self.assertAlmostEqual(c_array[i_node*3+i], node_val[i])


    def __Check_CArrayToKratosField(self, is_scalar, historical=True, add_values=False, swap_sign=False):
        if is_scalar:
            if historical:
                variable = KM.PRESSURE
                dim = 1
                fct_ptr = GetSolutionStepValue
            else:
                variable = KM.TEMPERATURE
                dim = 1
                fct_ptr = GetValue
        else:
            if historical:
                variable = KM.FORCE
                dim = 3
                fct_ptr = GetSolutionStepValue
            else:
                variable = KM.VELOCITY
                dim = 3
                fct_ptr = GetValue

        size = self.model_part.NumberOfNodes() * dim
        c_array = (ctp.c_double * size)(0.0)

        for i in range(size):
            c_array[i] = i+0.5

        vals_to_check_against = ArrayDeepCopy(c_array, size) # copy bcs the array might be modified inside "CArrayToKratosField"
        current_values = empire_mapper_wrapper.KratosFieldToCArray(self.model_part.Nodes, variable, historical)

        empire_mapper_wrapper.CArrayToKratosField(c_array, size, self.model_part.Nodes, variable, historical, add_values, swap_sign)

        # checking the values
        if swap_sign:
            for i in range(size):
                vals_to_check_against[i] *= (-1)
        if add_values:
            for i in range(size):
                vals_to_check_against[i] += current_values[i]

        if is_scalar:
            for i, node in enumerate(self.model_part.Nodes):
                self.assertAlmostEqual(vals_to_check_against[i], fct_ptr(node, variable))
        else:
            for i_node, node in enumerate(self.model_part.Nodes):
                node_val = fct_ptr(node, variable)
                for i in range(3):
                    self.assertAlmostEqual(vals_to_check_against[i_node*3+i], node_val[i])


def NodeScalarHistValue(the_id):
    return the_id*1.5

def NodeVectorHistValue(the_id):
    return [the_id*1.1, the_id*1.1+2.5, the_id*2.3]

def NodeScalarNonHistValue(the_id):
    return the_id*2.1+6.2

def NodeVectorNonHistValue(the_id):
    return [the_id*14.7, the_id*19.2-10.5, the_id*303.9]

def GetValue(node, variable):
    return node.GetValue(variable)

def GetSolutionStepValue(node, variable):
    return node.GetSolutionStepValue(variable)

def ArrayDeepCopy(c_array, size):
    c_array_copy = (ctp.c_double * size)(0.0)
    for i in range(size):
        c_array_copy[i] = c_array[i]
    return c_array_copy


if __name__ == '__main__':
    KratosUnittest.main()