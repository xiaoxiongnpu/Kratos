from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
data_comm = KM.DataCommunicator.GetDefault()
import beam_mapper_test_case
from math import cos
import math
import os
import numpy as np

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)

class BladeMappingTests(beam_mapper_test_case.MapperTestCase):
    '''This class contains basic tests for mapping on real geometries
    In this case it is a remodeled NREL Phase VI wind turbine blade
    It also serves as a showcase on how to use the Mapper in FSI
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters):
        structure_mdpa_file_name = "blade_line"
        fluid_mdpa_file_name     = "blade_quad"
        super(BladeMappingTests, cls).setUpModelParts(structure_mdpa_file_name, fluid_mdpa_file_name)

        # TODO ATTENTION: currently the MapperFactory removes some keys, hence those checks have to be done beforehand => improve this!

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

        mapper_parameters_p = mapper_parameters.Clone()
        mapper_parameters_s = mapper_parameters.Clone()

        #mapper_parameters_p.AddEmptyValue("interface_submodel_part_origin").SetString("pressure_side_line")
        mapper_parameters_p.AddEmptyValue("interface_submodel_part_destination").SetString("pressure_side_quad")

        #mapper_parameters_s.AddEmptyValue("interface_submodel_part_origin").SetString("suction_side_line")
        mapper_parameters_s.AddEmptyValue("interface_submodel_part_destination").SetString("suction_side_quad")

        cls.model_part_structure = cls.model_part_origin
        cls.model_part_fluid = cls.model_part_destination

        if data_comm.IsDistributed():
            cls.mapper_pressure_side = KratosMapping.MapperFactory.CreateMPIMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_p)
            cls.mapper_suction_side = KratosMapping.MapperFactory.CreateMPIMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_s)
        else:
            cls.mapper_pressure_side = KratosMapping.MapperFactory.CreateMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_p)
            cls.mapper_suction_side = KratosMapping.MapperFactory.CreateMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_s)

        cls.print_output = False # this can be overridden in derived classes to print the output

    def test_map_displacements(self):
        SetDisplacements(self.model_part_structure)
        if self.print_output:
            beam_mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.DISPLACEMENT, "Blade_" + self.mapper_type + "_Structure_prescr_disp")
            beam_mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.ROTATION, "Blade_" + self.mapper_type + "_Structure_prescr_rot")

        self.mapper_pressure_side.Map(KM.DISPLACEMENT, KM.ROTATION, KM.MESH_DISPLACEMENT)
        self.mapper_suction_side.Map(KM.DISPLACEMENT, KM.ROTATION, KM.MESH_DISPLACEMENT)

        if self.print_output:
            beam_mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.MESH_DISPLACEMENT, "Blade_" + self.mapper_type + "_Fluid_mapped_disp")

        beam_mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_fluid, KM.MESH_DISPLACEMENT, GetFilePath(self.__GetFileName("blade_map_disp")))

    def test_map_forces(self):
        SetReactions(self.model_part_fluid)
        if self.print_output:
            beam_mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.REACTION, "Blade_" + self.mapper_type + "_Fluid_prescr_force")

        # this would be POINT_LOAD in regular StructuralMechanics (using FORCE to avoid the StructuralMechanics import)
        self.mapper_pressure_side.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION)
        self.mapper_suction_side.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION)

        if self.print_output:
            beam_mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.FORCE, "Blade_" + self.mapper_type + "_Structure_mapped_force")
            beam_mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.MOMENT, "Blade_" + self.mapper_type + "_Structure_mapped_moment")

        beam_mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_structure, KM.FORCE, GetFilePath(self.__GetFileName("blade_map_force")))
        beam_mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_structure, KM.MOMENT, GetFilePath(self.__GetFileName("blade_map_moment")))

    def test_map_forces_conservative(self):
        SetReactions(self.model_part_fluid)
        if self.print_output:
            beam_mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.REACTION, "Blade_" + self.mapper_type + "_Fluid_prescr_force")

         # this would be POINT_LOAD in regular StructuralMechanics (using FORCE to avoid the StructuralMechanics import)
        self.mapper_pressure_side.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION)
         # should I fix this???
        #self.__CheckValuesSum(self.model_part_fluid.GetSubModelPart("pressure_side_quad"), self.model_part_structure.GetSubModelPart("pressure_side_line"), KM.REACTION, KM.FORCE)

        # Note: Setting the solution again because in this case some nodes are shared and hence
        # would slightly influence the computation
        SetReactions(self.model_part_fluid)
        self.mapper_suction_side.InverseMap(KM.FORCE, KM.MOMENT, KM.REACTION)

        if self.print_output:
            beam_mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.FORCE, "Blade_" + self.mapper_type + "_Structure_mapped_force_conserv")
            beam_mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.MOMENT, "Blade_" + self.mapper_type + "_Structure_mapped_moment_conserv")

        beam_mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_structure, KM.FORCE, GetFilePath(self.__GetFileName("blade_map_force_conserv")))
        beam_mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_structure, KM.MOMENT, GetFilePath(self.__GetFileName("blade_map_moment_conserv")))

        #self.__CheckValuesSum(self.model_part_fluid.GetSubModelPart("suction_side_quad"), self.model_part_structure.GetSubModelPart("suction_side_line"), KM.REACTION, KM.FORCE)

    def __CheckValuesSum(self, mp1, mp2, var1, var2):
        val_1 = KM.VariableUtils().SumHistoricalNodeVectorVariable(var1, mp1, 0)
        val_2 = KM.VariableUtils().SumHistoricalNodeVectorVariable(var2, mp2, 0)
        # minus because SWAP_SIGN is used
        self.assertAlmostEqual(val_1[0], -val_2[0])
        self.assertAlmostEqual(val_1[1], -val_2[1])
        self.assertAlmostEqual(val_1[2], -val_2[2])

    def __GetFileName(self, file_appendix):
        return os.path.join("result_files", self.mapper_type, file_appendix)

def SetDisplacements(model_part_structure):
    for node in model_part_structure.Nodes:
        #disp_x = 0.0 # edgewise
        #disp_y = 0.008*(node.Z)**3 # flapwise
        #disp_z = 0.0 # spanwise
        #node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Vector([disp_x, disp_y, disp_z]))
    
        lenght_beam = 4.521
        alfa = 0.000020929 
        beta = 0.0
        r = lenght_beam / alfa

        theta_X = - node.Z / r
        theta_Y = 0.0
        theta_Z = (beta * node.Z) / lenght_beam

        e_x = np.array([1.0, 0.0, 0.0])
        e_y = np.array([0.0, 1.0, 0.0])
        e_z = np.array([0.0, 0.0, 1.0])

        Rx = CalculateRotationMatrixWithAngle(e_x, theta_X)
        Ry = CalculateRotationMatrixWithAngle(e_y, theta_Y)
        Rz = CalculateRotationMatrixWithAngle(e_z, theta_Z)

        R_temp = np.dot(Ry, Rz)
        R = np.dot(Rx, R_temp)

        _ROTATION = getRotationVector(R)

        node.SetSolutionStepValue(KM.DISPLACEMENT_X, 0.0)
        node.SetSolutionStepValue(KM.DISPLACEMENT_Y, r - r*math.cos(-theta_X))
        node.SetSolutionStepValue(KM.DISPLACEMENT_Z, r * math.sin(-theta_X) - node.Z )
        node.SetSolutionStepValue(KM.ROTATION_X, _ROTATION[0] )
        node.SetSolutionStepValue(KM.ROTATION_Y, _ROTATION[1] )
        node.SetSolutionStepValue(KM.ROTATION_Z, _ROTATION[2] )

def SetReactions(model_part_fluid):
    for node in model_part_fluid.Nodes:
        react_x = cos(node.X*5)
        react_y = 0.0
        react_z = cos(node.Z*2)
        node.SetSolutionStepValue(KM.REACTION, KM.Vector([react_x, react_y, react_z]))

def CalculateRotationMatrixWithAngle( Axis, Angle ):
    rotationMatrix = np.zeros((3, 3))

    rotationMatrix[0][0] = math.cos( Angle ) + Axis[0]**2 * (1 - math.cos( Angle ))
    rotationMatrix[0][1] = Axis[0] * Axis[1] * (1 - math.cos( Angle )) - Axis[2] * math.sin( Angle )
    rotationMatrix[0][2] = Axis[0] * Axis[2] * (1 - math.cos( Angle )) + Axis[1] * math.sin( Angle )

    rotationMatrix[1][0] = Axis[0] * Axis[1] * (1 - math.cos( Angle )) + Axis[2] * math.sin( Angle )
    rotationMatrix[1][1] = math.cos( Angle ) + Axis[1]**2 * (1 - math.cos( Angle ))
    rotationMatrix[1][2] = Axis[1] * Axis[2] * (1 - math.cos( Angle )) - Axis[0] * math.sin( Angle )

    rotationMatrix[2][0] = Axis[0] * Axis[2] * (1 - math.cos( Angle )) - Axis[1] * math.sin( Angle )
    rotationMatrix[2][1] = Axis[1] * Axis[2] * (1 - math.cos( Angle )) + Axis[0] * math.sin( Angle )
    rotationMatrix[2][2] = math.cos( Angle ) + Axis[2]**2 * (1 - math.cos( Angle ))

    return rotationMatrix

def getRotationVector(rotationMatrix):
    # see Non-linear Modeling and Analysis of Solids and Structures (Steen Krenk 2009) P52
    rotationVector = np.array([0.0, 0.0, 0.0])
    angle = rotationMatrix[0][0] + rotationMatrix[1][1] + rotationMatrix[2][2] - 1.0

    angle = angle/2.0
    if (angle > 1.0):
        angle = 1.0
    elif (angle < -1.0):
        angle = -1.0

    angle = math.acos(angle) # between 0 and pi

    EPS = 1E-6
    M_PI = math.pi
    if (angle < EPS):
        rotationVector[0] = 0.0
        rotationVector[1] = 0.0
        rotationVector[2] = 0.0
        return rotationVector
    elif ((M_PI - angle) < EPS):
        product11 = (rotationMatrix[0][0] + 1.0) / 2.0
        product22 = (rotationMatrix[1][1] + 1.0) / 2.0
        product33 = (rotationMatrix[2][2] + 1.0) / 2.0
        product12 = (rotationMatrix[0][1] + 1.0) / 2.0
        product23 = (rotationMatrix[1][2] + 1.0) / 2.0
        product13 = (rotationMatrix[0][2] + 1.0) / 2.0
        tmp1 = math.sqrt(product11)
        tmp2 = math.sqrt(product22)
        tmp3 = math.sqrt(product33)

        rotationVector[0] = tmp1
        rotationVector[1] = tmp2
        rotationVector[2] = tmp3
        tmp12 = rotationVector[0] * rotationVector[1]
        tmp13 = rotationVector[0] * rotationVector[2]
        tmp23 = rotationVector[1] * rotationVector[2]
        if (math.fabs(tmp12) < EPS or math.fabs(tmp12 - product12) < math.fabs(tmp12 + product12)):
            if (math.fabs(tmp13) < EPS or math.fabs(tmp13 - product13) < math.fabs(tmp13 + product13)):
                if (math.fabs(tmp23) < EPS or math.fabs(tmp23 - product23) < math.fabs(tmp23 + product23)):
                    rotationVector[0] *= M_PI
                    rotationVector[1] *= M_PI
                    rotationVector[2] *= M_PI
                    return rotationVector

        rotationVector[0] =  tmp1
        rotationVector[1] = -tmp2
        rotationVector[2] = -tmp3
        tmp12 = rotationVector[0] * rotationVector[1]
        tmp13 = rotationVector[0] * rotationVector[2]
        tmp23 = rotationVector[1] * rotationVector[2]
        if (math.fabs(tmp12) < EPS or math.fabs(tmp12 - product12) < math.fabs(tmp12 + product12)):
            if (math.fabs(tmp13) < EPS or math.fabs(tmp13 - product13) < math.fabs(tmp13 + product13)):
                if (math.fabs(tmp23) < EPS or math.fabs(tmp23 - product23) < math.fabs(tmp23 + product23)):
                    rotationVector[0] *= M_PI
                    rotationVector[1] *= M_PI
                    rotationVector[2] *= M_PI
                    return rotationVector

        rotationVector[0] = -tmp1
        rotationVector[1] =  tmp2
        rotationVector[2] = -tmp3
        tmp12 = rotationVector[0] * rotationVector[1]
        tmp13 = rotationVector[0] * rotationVector[2]
        tmp23 = rotationVector[1] * rotationVector[2]
        if (math.fabs(tmp12) < EPS or math.fabs(tmp12 - product12) < math.fabs(tmp12 + product12)):
            if (math.fabs(tmp13) < EPS or math.fabs(tmp13 - product13) < math.fabs(tmp13 + product13)):
                if (math.fabs(tmp23) < EPS or math.fabs(tmp23 - product23) < math.fabs(tmp23 + product23)):
                    rotationVector[0] *= M_PI
                    rotationVector[1] *= M_PI
                    rotationVector[2] *= M_PI
                    return rotationVector

        rotationVector[0] = -tmp1
        rotationVector[1] = -tmp2
        rotationVector[2] =  tmp3
        tmp12 = rotationVector[0] * rotationVector[1]
        tmp13 = rotationVector[0] * rotationVector[2]
        tmp23 = rotationVector[1] * rotationVector[2]
        if (math.fabs(tmp12) < EPS or math.fabs(tmp12 - product12) < math.fabs(tmp12 + product12)):
            if (math.fabs(tmp13) < EPS or math.fabs(tmp13 - product13) < math.fabs(tmp13 + product13)):
                if (math.fabs(tmp23) < EPS or math.fabs(tmp23 - product23) < math.fabs(tmp23 + product23)):
                    rotationVector[0] *= M_PI
                    rotationVector[1] *= M_PI
                    rotationVector[2] *= M_PI
                    return rotationVector
        assert(false)
    
    tmp = angle / 2.0 / math.sin(angle)
    rotationVector[0] = -(rotationMatrix[1][2] - rotationMatrix[2][1]) * tmp
    rotationVector[1] =  (rotationMatrix[0][2] - rotationMatrix[2][0]) * tmp
    rotationVector[2] = -(rotationMatrix[0][1] - rotationMatrix[1][0]) * tmp

    return rotationVector
