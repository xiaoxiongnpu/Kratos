from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
data_comm = KM.DataCommunicator.GetDefault()
import mapper_test_case
from math import cos
import os

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)

class BladeMappingTests(mapper_test_case.MapperTestCase):
    '''This class contains basic tests for mapping on real geometries
    In this case it is a remodeled NREL Phase VI wind turbine blade
    It also serves as a showcase on how to use the Mapper in FSI
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters):
        structure_mdpa_file_name = "blade_quad"
        fluid_mdpa_file_name     = "blade_line"
        super(BladeMappingTests, cls).setUpModelParts(structure_mdpa_file_name, fluid_mdpa_file_name)

        # TODO ATTENTION: currently the MapperFactory removes some keys, hence those checks have to be done beforehand => improve this!

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

        mapper_parameters_p = mapper_parameters.Clone()
        mapper_parameters_s = mapper_parameters.Clone()

        mapper_parameters_p.AddEmptyValue("interface_submodel_part_origin").SetString("pressure_side_quad")
        mapper_parameters_p.AddEmptyValue("interface_submodel_part_destination").SetString("pressure_side_line")

        mapper_parameters_s.AddEmptyValue("interface_submodel_part_origin").SetString("suction_side_quad")
        mapper_parameters_s.AddEmptyValue("interface_submodel_part_destination").SetString("suction_side_line")

        cls.model_part_structure = cls.model_part_origin
        cls.model_part_fluid = cls.model_part_destination
        if data_comm.IsDistributed():
            cls.mapper_pressure_side = KratosMapping.MapperFactory.CreateMPIMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_p)
            cls.mapper_suction_side = KratosMapping.MapperFactory.CreateMPIMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_s)
        else:
            cls.mapper_pressure_side = KratosMapping.MapperFactory.CreateMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_p)
            cls.mapper_suction_side = KratosMapping.MapperFactory.CreateMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_s)

        cls.print_output = False # this can be overridden in derived classes to print the output