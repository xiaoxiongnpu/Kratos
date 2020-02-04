import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication.python_mapper import PythonMapper

import os
import ctypes as ctp

def CreateMapper(model_part_origin, model_part_destination, mapper_settings):
    empire_mapper_type =  mapper_settings["mapper_type"].GetString()[7:] # returns the mapper-type without preceeding "empire_"

    if empire_mapper_type == "nearest_neighbor":
        return EmpireNearestNeighborMapper(model_part_origin, model_part_destination, mapper_settings)

    elif empire_mapper_type == "nearest_element":
        return EmpireNearestElementMapper(model_part_origin, model_part_destination, mapper_settings)

    elif empire_mapper_type == "barycentric":
        return EmpireBarycentricMapper(model_part_origin, model_part_destination, mapper_settings)

    elif empire_mapper_type == "mortar":
        return EmpireMortarMapper(model_part_origin, model_part_destination, mapper_settings)

    else:
        raise Exception('Mapper "{}" not available in empire!'.format(empire_mapper_type))


class EmpireMapperWrapper(PythonMapper):
    """Wrapper for the mappers of Empire with the same Interface as the Mappers in Kratos
    This wrapper requires the development version of Empire which has the mapperlib exposed separately
    """

    mapper_count = 0
    mapper_lib = None

    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        super(EmpireMapperWrapper, self).__init__(model_part_origin, model_part_destination, mapper_settings)

        if model_part_origin.IsDistributed() or model_part_destination.IsDistributed():
            raise Exception('{} does not support mapping with distributed ModelParts!'.format(self._ClassName()))

            # name that is used inside Empire, has to be unique, hence the counter
        self.mapper_name = (self._ClassName()+str(EmpireMapperWrapper.mapper_count)).encode(encoding='UTF-8')

        EmpireMapperWrapper.mapper_count += 1 # required for identification purposes
        self.__inverse_mapper = None

        if EmpireMapperWrapper.mapper_lib == None: # load it only once
            KM.Logger.PrintInfo("EmpireMapperWrapper", "Loading mapper lib ...")
            EmpireMapperWrapper.mapper_lib = self.__LoadMapperLib()

        self.__CreateMeshes()
        self._CreateMapper()
        self.__BuildCouplingMatrices()

    def __del__(self):
        EmpireMapperWrapper.mapper_count -= 1

        EmpireMapperWrapper.mapper_lib.deleteMapper(self.mapper_name)

        if EmpireMapperWrapper.mapper_count == 0: # last mapper was destoyed
            if self.echo_level > 1:
                KM.Logger.PrintInfo('EmpireMapperWrapper', 'Destroying last instance, deleting all meshes & mappers')
            #  delete everything to make sure nothing is left
            EmpireMapperWrapper.mapper_lib.deleteAllMeshes()
            EmpireMapperWrapper.mapper_lib.deleteAllMappers()

    # public methods, same as in "custom_mappers/mapper.h"
    def Map(self, variable_origin, variable_destination, mapper_flags=KM.Flags()):
        self.__CheckMapperExists()
        # if not type() # check variables are matching

        var_dim = GetVariableDimension(variable_origin)

        origin_data_size = len(self.model_part_origin.Nodes)*var_dim
        destination_data_size = len(self.model_part_destination.Nodes)*var_dim

        c_origin_array = KratosFieldToCArray(self.model_part_origin.Nodes, variable_origin, True)
        c_destination_array = (ctp.c_double * destination_data_size)(0.0)

        EmpireMapperWrapper.mapper_lib.doConsistentMapping(
            ctp.c_char_p(self.mapper_name),
            ctp.c_int(var_dim),
            ctp.c_int(origin_data_size),
            c_origin_array,
            ctp.c_int(destination_data_size),
            c_destination_array
            )

        CArrayToKratosField(
            c_destination_array,
            destination_data_size,
            self.model_part_destination.Nodes,
            variable_destination,
            True, False, False)

    def InverseMap(self, variable_origin, variable_destination, mapper_flags=KM.Flags()):
        # TODO check if using transpose => conservative

        if self.__inverse_mapper == None:
            self.__CreateInverseMapper()
        self.__inverse_mapper.Map(variable_destination, variable_origin, mapper_flags) # TODO check this!
        raise NotImplementedError('"InverseMap" was not implemented for "{}"'.format(self._ClassName()))

    def UpdateInterface(self):
        # this requires recreating the mapper completely!
        super(EmpireMapperWrapper, self).UpdateInterface()

    # protected methods
    def _CreateMapper(self):
        raise NotImplementedError('"_CreateMapper" was not implemented for "{}"'.format(self._ClassName()))

    # private methods
    def __LoadMapperLib(self):
        KM.Logger.PrintInfo("EmpireMapperWrapper", "Determining path to mapper lib")
        # first try automatic detection using the environment that is set by Empire => startEMPIRE
        if ('EMPIRE_MAPPER_LIBSO_ON_MACHINE' in os.environ):
            KM.Logger.PrintInfo("EmpireMapperWrapper", "EMPIRE_MAPPER_LIBSO_ON_MACHINE found in environment")
            mapper_lib_path = os.environ['EMPIRE_MAPPER_LIBSO_ON_MACHINE']

        else:
            KM.Logger.PrintInfo("EmpireMapperWrapper", "EMPIRE_MAPPER_LIBSO_ON_MACHINE NOT found in environment, using manually specified path to load the mapper lib")
            mapper_lib_path = self.mapper_settings["mapper_lib"].GetString()
            if mapper_lib_path == "":
                raise Exception('The automatic detection of the mapper lib failed, the path to the mapper lib has to be specified with "mapper_lib"')

        KM.Logger.PrintInfo("EmpireMapperWrapper", "Attempting to load the mapper lib")
        # TODO check if still both are needed! (does the mapperlib link to MPI?)
        try:
            try: # OpenMPI
                loaded_mapper_lib = ctp.CDLL(mapper_lib_path, ctp.RTLD_GLOBAL)
                KM.Logger.PrintInfo('EmpireMapperWrapper', 'Using standard OpenMPI')
            except: # Intel MPI or OpenMPI compiled with "–disable-dlopen" option
                loaded_mapper_lib = ctp.cdll.LoadLibrary(mapper_lib_path)
                KM.Logger.PrintInfo('EmpireMapperWrapper', 'Using Intel MPI or OpenMPI compiled with "–disable-dlopen" option')
        except OSError:
            raise Exception("Mapper lib could not be loaded!")

        KM.Logger.PrintInfo("EmpireMapperWrapper", "Successfully loaded the mapper lib")

        return loaded_mapper_lib

    def __BuildCouplingMatrices(self):
        self.__CheckMapperExists()
        EmpireMapperWrapper.mapper_lib.buildCouplingMatrices(ctp.c_char_p(self.mapper_name))

    def __CreateInverseMapper(self):
        return self.__class__(self.model_part_destination, self.model_part_origin, self.mapper_settings) # TODO check this!

    def __CreateMeshes(self):
        self.mesh_name_origin      = ctp.c_char_p((self.model_part_origin.FullName()).encode(encoding='UTF-8'))
        self.mesh_name_destination = ctp.c_char_p((self.model_part_destination.FullName()).encode(encoding='UTF-8'))

        if not EmpireMapperWrapper.mapper_lib.hasMesh(self.mesh_name_origin):
            self.__MakeEmpireFEMesh(self.mesh_name_origin, self.model_part_origin)

        if not EmpireMapperWrapper.mapper_lib.hasMesh(self.mesh_name_destination):
            self.__MakeEmpireFEMesh(self.mesh_name_destination, self.model_part_destination)

    def __MakeEmpireFEMesh(self, c_mesh_name, model_part):
        # c_mesh_name          = ctp.c_char_p(mesh_name)
        c_num_nodes          = ctp.c_int(len(model_part.Nodes))
        c_num_elems          = ctp.c_int(len(model_part.Conditions))
        c_node_ids           = (ctp.c_int * c_num_nodes.value) (0)
        c_node_coords        = (ctp.c_double * (3 * c_num_nodes.value))(0.0)
        c_num_nodes_per_elem = (ctp.c_int * c_num_elems.value) (0)

        for i_node, node in enumerate(model_part.Nodes):
            c_node_ids[i_node] = node.Id
            c_node_coords[i_node*3]   = node.X
            c_node_coords[i_node*3+1] = node.Y
            c_node_coords[i_node*3+2] = node.Z

        elem_node_ctr = 0
        for elem_ctr, elem in enumerate(model_part.Conditions):
            c_num_nodes_per_elem[elem_ctr] = len(elem.GetNodes())
            elem_node_ctr += c_num_nodes_per_elem[elem_ctr]

        elem_index = 0
        c_elems = (ctp.c_int * elem_node_ctr) (0)
        for elem_ctr, elem in enumerate(model_part.Conditions):
            for elem_node_ctr, elem_node in enumerate(elem.GetNodes()):
                c_elems[elem_index + elem_node_ctr] = elem_node.Id
            elem_index += len(elem.GetNodes())

        EmpireMapperWrapper.mapper_lib.initFEMesh(c_mesh_name, c_num_nodes, c_num_elems, False)
        EmpireMapperWrapper.mapper_lib.setNodesToFEMesh(c_mesh_name, c_node_ids, c_node_coords)
        EmpireMapperWrapper.mapper_lib.setElementsToFEMesh(c_mesh_name, c_num_nodes_per_elem, c_elems)

    def __CheckMapperExists(self):
        if not EmpireMapperWrapper.mapper_lib.hasMapper(ctp.c_char_p(self.mapper_name)):
            raise Exception('Mapper "{}" does not exist!'.format(self.mapper_name))

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "mapper_lib"  : ""
        }""")
        this_defaults.AddMissingParameters(super(EmpireMapperWrapper, cls)._GetDefaultSettings())
        return this_defaults


class EmpireNearestNeighborMapper(EmpireMapperWrapper):
    """Wrapper for the nearest neighbor mapper of Empire"""

    def _CreateMapper(self):
        EmpireMapperWrapper.mapper_lib.initFEMNearestNeighborMapper(
            ctp.c_char_p(self.mapper_name),
            self.mesh_name_origin,
            self.mesh_name_destination
            )


class EmpireNearestElementMapper(EmpireMapperWrapper):
    """Wrapper for the nearest element mapper of Empire"""

    def _CreateMapper(self):
        EmpireMapperWrapper.mapper_lib.initFEMNearestElementMapper(
            ctp.c_char_p(self.mapper_name),
            self.mesh_name_origin,
            self.mesh_name_destination
            )


class EmpireBarycentricMapper(EmpireMapperWrapper):
    """Wrapper for the barycentric mapper of Empire"""

    def _CreateMapper(self):
        EmpireMapperWrapper.mapper_lib.initFEMBarycentricInterpolationMapper(
            ctp.c_char_p(self.mapper_name),
            self.mesh_name_origin,
            self.mesh_name_destination
            )


class EmpireMortarMapper(EmpireMapperWrapper):
    """Wrapper for the mortar mapper of Empire"""

    def _CreateMapper(self):
        dual               = int(self.mapper_settings["dual"].GetBool())
        enforceConsistency = int(self.mapper_settings["enforce_consistency"].GetBool())
        opposite_normals   = int(self.mapper_settings["opposite_normals"].GetBool())

        EmpireMapperWrapper.mapper_lib.initFEMMortarMapper(
            ctp.c_char_p(self.mapper_name),
            self.mesh_name_origin,
            self.mesh_name_destination,
            ctp.c_int(opposite_normals),
            ctp.c_int(dual),
            ctp.c_int(enforceConsistency)
            )

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "dual"                 : false,
            "enforce_consistency"  : false,
            "opposite_normals"     : false
        }""")
        this_defaults.AddMissingParameters(super(EmpireMapperWrapper, cls)._GetDefaultSettings())
        return this_defaults


# Helper functions
def GetVariableDimension(variable):
    var_type = KM.KratosGlobals.GetVariableType(variable.Name())
    if var_type == "Array":
        return 3
    elif var_type in ["Double", "Component"]:
        return 1
    else:
        raise Exception('Wrong variable type: "{}". Only "Array", "Double" and "Component" are allowed'.format(var_type))

def GetValue(node, variable):
    return node.GetValue(variable)

def GetSolutionStepValue(node, variable):
    return node.GetSolutionStepValue(variable)

def SetValue(node, variable, value):
    return node.SetValue(variable, value)

def SetSolutionStepValue(node, variable, value):
    return node.SetSolutionStepValue(variable, 0, value)


def KratosFieldToCArray(nodes, variable, historical):
    dim = GetVariableDimension(variable)
    size = dim * len(nodes)
    c_array = (ctp.c_double * size)(0.0)

    if historical:
        fct_ptr = GetSolutionStepValue
    else:
        fct_ptr = GetValue

    if dim == 1:
        for i_node, node in enumerate(nodes):
            c_array[i_node] = fct_ptr(node, variable)
    else:
        for i_node, node in enumerate(nodes):
            node_value = fct_ptr(node, variable)
            for i_dim in range(dim):
                c_array[i_node*dim + i_dim] = node_value[i_dim]

    return c_array


def CArrayToKratosField(c_array, c_array_size, nodes, variable, historical, add_values, swap_sign):
    dim = GetVariableDimension(variable)
    if c_array_size != dim * len(nodes):
        raise RuntimeError("Wrong size!")

    if historical:
        fct_ptr = SetSolutionStepValue
    else:
        fct_ptr = SetValue

    if swap_sign:
        for i in range(c_array_size):
            c_array[i] *= (-1)

    if add_values:
        current_values = KratosFieldToCArray(nodes, variable, historical)
        for i in range(c_array_size):
            c_array[i] += current_values[i]

    if dim == 1:
        for i_node, node in enumerate(nodes):
            fct_ptr(node, variable, c_array[i_node])
    else:
        for i_node, node in enumerate(nodes):
            values = [c_array[i_node*dim], c_array[i_node*dim+1], c_array[i_node*dim+2]]
            fct_ptr(node, variable, values)
