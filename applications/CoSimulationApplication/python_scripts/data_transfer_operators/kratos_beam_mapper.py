from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as KSM
import KratosMultiphysics.MappingApplication as KratosMapping

def Create(settings):
    return KratosBeamMapper(settings)

class KratosBeamMapper(CoSimulationDataTransferOperator):

    # currently available mapper-flags aka transfer-options
    __mapper_flags_dict = {
        "swap_sign"     : KratosMapping.Mapper.SWAP_SIGN
    }

    def __init__(self, settings):
        super(KratosBeamMapper, self).__init__(settings)
        self.__mappers = {}

        self.second_variable_displacement = KM.KratosGlobals.GetVariable(self.settings["second_variable_displacement"].GetString())
        self.second_variable_force = KM.KratosGlobals.GetVariable(self.settings["second_variable_force"].GetString())

    def TransferData(self, from_solver_data, to_solver_data, transfer_options):
        # TODO check location of data => should coincide with the one for the mapper
        # or throw if it is not in a suitable location (e.g. on the ProcessInfo)

        self._CheckAvailabilityTransferOptions(transfer_options)

        model_part_origin      = from_solver_data.GetModelPart()
        model_part_origin_name = model_part_origin.Name
        variable_origin        = from_solver_data.variable

        model_part_destinatinon      = to_solver_data.GetModelPart()
        model_part_destinatinon_name = model_part_destinatinon.Name
        variable_destination         = to_solver_data.variable

        mapper_flags = self.__GetMapperFlags(transfer_options)

        name_tuple         = (model_part_origin_name, model_part_destinatinon_name)
        inverse_name_tuple = (model_part_destinatinon_name, model_part_origin_name)

        if name_tuple in self.__mappers:
            if not variable_origin == KM.DISPLACEMENT:
                raise Exception("Wrong variable")
            self.__mappers[name_tuple].Map(variable_origin, self.second_variable_displacement, variable_destination, mapper_flags)
        elif inverse_name_tuple in self.__mappers:
            if not variable_origin == KSM.POINT_LOAD:
                raise Exception("Wrong variable")
            self.__mappers[inverse_name_tuple].InverseMap(variable_destination, self.second_variable_force, variable_origin, mapper_flags)
        else:
            if model_part_origin.IsDistributed() or model_part_destinatinon.IsDistributed():
                raise NotImplementedError("Beam mapper does not yet work in MPI")
                mapper_create_fct = KratosMapping.MapperFactory.CreateMPIMapper
            else:
                mapper_create_fct = KratosMapping.MapperFactory.CreateMapper

            beam_mapper_settings = self.settings["mapper_settings"].Clone()
            beam_mapper_settings.AddEmptyValue("mapper_type").SetString("beam_mapper")
            if not variable_origin == KM.DISPLACEMENT:
                raise Exception("Wrong variable")
            self.__mappers[name_tuple] = mapper_create_fct(model_part_origin, model_part_destinatinon, beam_mapper_settings) # Clone is necessary here bcs settings are influenced among mappers otherwise. TODO check in the MapperFactory how to solve this better
            self.__mappers[name_tuple].Map(variable_origin, self.second_variable_displacement, variable_destination, mapper_flags)

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "mapper_settings" : { },
            "second_variable_displacement" : "ROTATION"
            "second_variable_force"        : "POINT_MOMENT"
        }""")
        this_defaults.AddMissingParameters(super(KratosBeamMapper, cls)._GetDefaultSettings())
        return this_defaults

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return cls.__mapper_flags_dict.keys()

    def __GetMapperFlags(self, transfer_options):
        mapper_flags = KM.Flags()
        for flag_name in transfer_options.GetStringArray():
            mapper_flags |= self.__mapper_flags_dict[flag_name]

        return mapper_flags
