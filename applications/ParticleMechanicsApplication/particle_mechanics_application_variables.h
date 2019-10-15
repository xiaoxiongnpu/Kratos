//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra, Ilaria Iaconeta
//
//


#if !defined(KRATOS_PARTICLE_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_PARTICLE_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "includes/mat_variables.h"


namespace Kratos
{

    typedef array_1d<double,3> Vector3;
    typedef array_1d<double,6> Vector6;

    // Variables definition

    /* MATERIAL POINT ELEMENTS VARIABLES */
    // Indexing
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, MP_MATERIAL_ID )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, PARTICLES_PER_ELEMENT )

    // Physical
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_MASS )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_DENSITY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_VOLUME )

    // Energy
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_POTENTIAL_ENERGY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_KINETIC_ENERGY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_STRAIN_ENERGY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_TOTAL_ENERGY )

    // Pressure
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_PRESSURE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, PRESSURE_REACTION )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, NODAL_MPRESSURE )

    // Position and kinematics
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_COORD )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_DISPLACEMENT )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_VELOCITY )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MIDDLE_VELOCITY )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MIDDLE_ANGULAR_VELOCITY )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_ACCELERATION )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MP_VOLUME_ACCELERATION )
    KRATOS_DEFINE_APPLICATION_VARIABLE(STRUCTURAL_MECHANICS_APPLICATION, double, NODAL_DISPLACEMENT_DAMPING)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(STRUCTURAL_MECHANICS_APPLICATION, NODAL_ROTATION_DAMPING)

    // Stress Measures
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, Vector, MP_CAUCHY_STRESS_VECTOR )
	KRATOS_DEFINE_APPLICATION_VARIABLE(PARTICLE_MECHANICS_APPLICATION, bool, MUSL_VELOCITY_FIELD_COMPUTED)

    // Strain Measures
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, Vector, MP_ALMANSI_STRAIN_VECTOR )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_DELTA_PLASTIC_STRAIN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_DELTA_PLASTIC_DEVIATORIC_STRAIN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_EQUIVALENT_PLASTIC_STRAIN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN )

    // Constitutive law
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER )
    // CL: Solid
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, RAYLEIGH_ALPHA )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, RAYLEIGH_BETA )
    // CL: Mohr Coulomb
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, COHESION )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, INTERNAL_DILATANCY_ANGLE )
    // CL: Mohr Coulomb Strain Softening
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, INTERNAL_FRICTION_ANGLE_RESIDUAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, COHESION_RESIDUAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, INTERNAL_DILATANCY_ANGLE_RESIDUAL )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, SHAPE_FUNCTION_BETA )

    /* NODAL VARIABLES */
    // Conditions
    // Particle Conditions
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_COORD )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, MPC_CONDITION_ID )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, bool, MPC_IS_NEUMANN )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, MPC_AREA )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_NORMAL )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_DISPLACEMENT )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_VELOCITY )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, MPC_ACCELERATION )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, int, PARTICLES_PER_CONDITION )

    // Essential Boundary Conditions
    KRATOS_DEFINE_APPLICATION_VARIABLE( PARTICLE_MECHANICS_APPLICATION, double, PENALTY_FACTOR )

    // Natural Boundary Conditions
    // Nodal load variables
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PARTICLE_MECHANICS_APPLICATION, POINT_LOAD )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PARTICLE_MECHANICS_APPLICATION, LINE_LOAD )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PARTICLE_MECHANICS_APPLICATION, SURFACE_LOAD )

    // Momentum
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, NODAL_MOMENTUM )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, NODAL_INERTIA )
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(PARTICLE_MECHANICS_APPLICATION, NODAL_INTERNAL_FORCE )

}

#endif // KRATOS_PARTICLE_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED  defined