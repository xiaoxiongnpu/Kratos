from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.KratosUnittest as KratosUnittest

import numpy as np
import math

#mdpa_file_name_beam    = "mdpa_files/beam_mesh_1"
#mdpa_file_name_surface = "mdpa_files/beam_surface_10_elements"

#mdpa_file_name_beam    = "mdpa_files/beam_new"
#mdpa_file_name_surface = "mdpa_files/surface_mesh_Tianyang"

mdpa_file_name_beam    = "mdpa_files/line_10_elements"
mdpa_file_name_surface = "mdpa_files/surface_mesh_Tianyang"

def WriteGiDOutput(model_part):
    from gid_output_process import GiDOutputProcess
    gid_output = GiDOutputProcess(model_part,
        "gid_output_"+model_part.Name,
        KM.Parameters("""
            {
                "result_file_configuration" : {
                    "gidpost_flags": {
                        "GiDPostMode"           : "GiD_PostAscii",
                        "WriteDeformedMeshFlag" : "WriteUndeformed",
                        "WriteConditionsFlag"   : "WriteConditions",
                        "MultiFileFlag"         : "SingleFile"
                    },
                    "nodal_results"       : ["DISPLACEMENT"],
                    "gauss_point_results" : []
                }
            }
            """)
        )

    gid_output.ExecuteInitialize()
    gid_output.ExecuteBeforeSolutionLoop()
    gid_output.ExecuteInitializeSolutionStep()
    gid_output.PrintOutput()
    gid_output.ExecuteFinalizeSolutionStep()
    gid_output.ExecuteFinalize()

def WriteVtkOutputStructure(model_part):
    default_parameters = KM.Parameters("""{
        "file_format"                        : "binary",
        "output_precision"                   : 7,
        "output_control_type"                : "step",
        "output_sub_model_parts"             : false,
        "save_output_files_in_folder"        : false,
        "nodal_solution_step_data_variables" : ["DISPLACEMENT", "ROTATION" , "REACTION", "REACTION_MOMENT"]
    }""")

    vtk_io = KM.VtkOutput(model_part, default_parameters)
    vtk_io.PrintOutput()

def WriteVtkOutputFluid(model_part):
    default_parameters = KM.Parameters("""{
        "file_format"                        : "binary",
        "output_precision"                   : 7,
        "output_control_type"                : "step",
        "output_sub_model_parts"             : false,
        "save_output_files_in_folder"        : false,
        "nodal_solution_step_data_variables" : ["DISPLACEMENT", "REACTION"]
    }""")

    vtk_io = KM.VtkOutput(model_part, default_parameters)
    vtk_io.PrintOutput()

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

class TestBeamMapper(KratosUnittest.TestCase):
    def setUp(self):
        self.current_model = KM.Model()
        self.model_part_beam = self.current_model.CreateModelPart("beam")
        self.model_part_surface = self.current_model.CreateModelPart("surface")

        # list of variables involved in the Mapper-Tests
        self.model_part_beam.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.model_part_beam.AddNodalSolutionStepVariable(KM.ROTATION)
        self.model_part_beam.AddNodalSolutionStepVariable(KM.REACTION)
        self.model_part_beam.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)

        self.model_part_surface.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.model_part_surface.AddNodalSolutionStepVariable(KM.REACTION)

        KM.ModelPartIO(mdpa_file_name_beam).ReadModelPart(self.model_part_beam)
        KM.ModelPartIO(mdpa_file_name_surface).ReadModelPart(self.model_part_surface)
    
    def addDofs(self):
        for node in self.model_part_beam.Nodes:
            node.AddDof(KM.DISPLACEMENT_X, KM.REACTION_X)
            node.AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y)
            node.AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z)
            node.AddDof(KM.ROTATION_X, KM.REACTION_MOMENT_X)
            node.AddDof(KM.ROTATION_Y, KM.REACTION_MOMENT_Y)
            node.AddDof(KM.ROTATION_Z, KM.REACTION_MOMENT_Z)

        for node in self.model_part_surface.Nodes:
            node.AddDof(KM.DISPLACEMENT_X, KM.REACTION_X)
            node.AddDof(KM.DISPLACEMENT_Y, KM.REACTION_Y)
            node.AddDof(KM.DISPLACEMENT_Z, KM.REACTION_Z)

    def test_beam_mapper(self):
        mapper_settings = KM.Parameters("""{
            "mapper_type": "beam_mapper",
            "echo_level" : 3,
            "local_coord_tolerance" : 0.25
        }""")

        self.mapper = KratosMapping.MapperFactory.CreateMapper(self.model_part_beam, self.model_part_surface, mapper_settings)

        for node in self.model_part_beam.Nodes:
            lenght_beam = 100
            alfa = 3.1415 # 20° = 0.3491 rad, 40° = 0.6981, 60° = 1.0472 alfa is the slope of the right end
            beta = alfa
            r = lenght_beam / alfa
            
            theta_X = (beta * node.X) / lenght_beam
            theta_Y = 0.0
            theta_Z = - node.X / r 

            e_x = np.array([1.0, 0.0, 0.0])
            e_y = np.array([0.0, 1.0, 0.0])
            e_z = np.array([0.0, 0.0, 1.0])

            Rx = CalculateRotationMatrixWithAngle(e_x, theta_X)
            Ry = CalculateRotationMatrixWithAngle(e_y, theta_Y)
            Rz = CalculateRotationMatrixWithAngle(e_z, theta_Z)

            R_temp = np.dot(Rx, Ry)
            R = np.dot(Rz, R_temp)

            _rotation = getRotationVector(R)
            
            node.SetSolutionStepValue(KM.DISPLACEMENT_X, r * math.sin(-theta_Z) - node.X)
            node.SetSolutionStepValue(KM.DISPLACEMENT_Y, -r + r*math.cos(-theta_Z))
            node.SetSolutionStepValue(KM.DISPLACEMENT_Z, 0 )
            node.SetSolutionStepValue(KM.ROTATION_X, _rotation[0] )
            node.SetSolutionStepValue(KM.ROTATION_Y, _rotation[1] )
            node.SetSolutionStepValue(KM.ROTATION_Z, _rotation[2] )

        for node in self.model_part_surface.Nodes:
            node.SetSolutionStepValue(KM.REACTION_X, 1.0)
            node.SetSolutionStepValue(KM.REACTION_Y, 0.0)
            node.SetSolutionStepValue(KM.REACTION_Z, 0.0)
            
        self.mapper.Map(KM.DISPLACEMENT, KM.ROTATION, KM.DISPLACEMENT)
        self.mapper.InverseMap(KM.REACTION, KM.REACTION_MOMENT, KM.REACTION) 

        #WriteGiDOutput(self.model_part_beam)
        #WriteGiDOutput(self.model_part_surface)

        WriteVtkOutputStructure(self.model_part_beam)
        WriteVtkOutputFluid(self.model_part_surface)


if __name__ == '__main__':
    KratosUnittest.main()