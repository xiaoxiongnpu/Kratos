#include "meshless_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PRESSURE_ACC );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( VISCOUS_ACC );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( BODYFORCE_ACC );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( BOUNDARY_ACC );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( XSPH_VELOCITY );

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( TEMP_POS );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( TEMP_VEL );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( TEMP_RHS );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( TEMP_DISP );


KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( OLD_VEL );
//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LAGRANGE_DISPLACEMENT );

KRATOS_CREATE_VARIABLE(double, EFFECTIVE_RADIUS);

KRATOS_CREATE_VARIABLE(double, DENSITY_NORM_PARAM);
KRATOS_CREATE_VARIABLE(double, DIV_OF_VEL);

KRATOS_CREATE_VARIABLE(double, OLD_DENSITY);

KRATOS_CREATE_VARIABLE(double, DENS_VARIATION);
KRATOS_CREATE_VARIABLE(double, DENS_DIFF);
KRATOS_CREATE_VARIABLE(double, IS_WET);
KRATOS_CREATE_VARIABLE(double, INI_PRESSURE);

KRATOS_CREATE_VARIABLE(double, OUT_OF_SYSTEM);




KRATOS_CREATE_VARIABLE(double, VER_WALL_LEFT);
KRATOS_CREATE_VARIABLE(double, VER_WALL_RIGHT);
KRATOS_CREATE_VARIABLE(double, HOR_WALL_BOTTOM);



KRATOS_CREATE_VARIABLE(double, DUMMY_NORMALIZE_RHS);
KRATOS_CREATE_VARIABLE(double, DUMMY_APPLY_XSPH);
KRATOS_CREATE_VARIABLE(double, DUMMY_BOUNDARY_PRESSURES);
KRATOS_CREATE_VARIABLE(double, DUMMY_CATCH_FREESURFACE);
KRATOS_CREATE_VARIABLE(double, DUMMY_INTERMEDIATE_RHS);
KRATOS_CREATE_VARIABLE(double, DELTA_TIME_ISPH);

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( BODYFORCE);
//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_DISPLACEMENT);

KRATOS_CREATE_VARIABLE(double, DAMAGE);
KRATOS_CREATE_VARIABLE(double, FRACTURE_ENERGY);

//KRATOS_CREATE_VARIABLE(Vector, CAUCHY_STRESS_VECTOR );
//KRATOS_CREATE_VARIABLE(Vector, PK2_STRESS_VECTOR );


//KRATOS_CREATE_VARIABLE(Matrix, ALMANSI_STRAIN_TENSOR );
//KRATOS_CREATE_VARIABLE(Vector, GREEN_LAGRANGE_STRAIN_VECTOR );
//KRATOS_CREATE_VARIABLE(Vector, ALMANSI_STRAIN_VECTOR );

//KRATOS_CREATE_VARIABLE(double, VON_MISES_STRESS );
KRATOS_CREATE_VARIABLE(Matrix, NODAL_STRESS );
KRATOS_CREATE_VARIABLE(Matrix, NODAL_STRAIN );
KRATOS_CREATE_VARIABLE(double, NODAL_VOLUME);
KRATOS_CREATE_VARIABLE(double, NODAL_DAMAGE);
KRATOS_CREATE_VARIABLE(double,GAUSS_AREA);

KRATOS_CREATE_VARIABLE(Vector, STRESSES);
KRATOS_CREATE_VARIABLE(Vector, PRESTRESS);
KRATOS_CREATE_VARIABLE(Vector, PLASTIC_STRAIN_VECTOR);

KRATOS_CREATE_VARIABLE(ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER );

KRATOS_CREATE_VARIABLE( double, ALPHA );
KRATOS_CREATE_VARIABLE( double, RETRACTION_TIME );


KRATOS_CREATE_VARIABLE(double, YIELD_RATIO);
KRATOS_CREATE_VARIABLE(double, FC);
KRATOS_CREATE_VARIABLE(double, FT);
KRATOS_CREATE_VARIABLE(double, CRUSHING_ENERGY);
KRATOS_CREATE_VARIABLE(double, DILATANCY_ANGLE);
KRATOS_CREATE_VARIABLE(double, COHESION);
KRATOS_CREATE_VARIABLE(Matrix, GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR);
KRATOS_CREATE_VARIABLE(Vector, ALMANSI_PLASTIC_STRAIN);
KRATOS_CREATE_VARIABLE(Vector, ALMANSI_ELASTIC_STRAIN);
KRATOS_CREATE_VARIABLE(Matrix, GRADIENT_DEFORMATION_TENSOR);

KRATOS_CREATE_VARIABLE(double, EQUIVALENT_PLASTIC_STRAIN);
KRATOS_CREATE_VARIABLE(double, YIELD_SURFACE);
KRATOS_CREATE_VARIABLE(double, ISOTROPIC_HARDENING_MODULUS);

KRATOS_CREATE_VARIABLE(double, DAMPING_RATIO);
KRATOS_CREATE_VARIABLE(double, INTERNAL_FRICTION_ANGLE);

KRATOS_CREATE_VARIABLE(Vector, EXTERNAL_FORCES_VECTOR );
KRATOS_CREATE_VARIABLE(Vector, INTERNAL_FORCES_VECTOR );
KRATOS_CREATE_VARIABLE(Vector, CONTACT_FORCES_VECTOR );

KRATOS_CREATE_VARIABLE(Vector, CAUCHY_STRESS_VECTOR );
KRATOS_CREATE_VARIABLE(Vector, PK2_STRESS_VECTOR );

KRATOS_CREATE_VARIABLE(Matrix, ALMANSI_STRAIN_TENSOR );
KRATOS_CREATE_VARIABLE(Vector, GREEN_LAGRANGE_STRAIN_VECTOR );
KRATOS_CREATE_VARIABLE(Vector, ALMANSI_STRAIN_VECTOR );

KRATOS_CREATE_VARIABLE(Matrix, MATERIAL_STIFFNESS_MATRIX );
KRATOS_CREATE_VARIABLE(Matrix, GEOMETRIC_STIFFNESS_MATRIX );

KRATOS_CREATE_VARIABLE(double, VON_MISES_STRESS );

//KRATOS_CREATE_VARIABLE( bool,COMPUTE_RHS_VECTOR );
//KRATOS_CREATE_VARIABLE( bool,COMPUTE_LHS_MATRIX );
//KRATOS_CREATE_VARIABLE( bool,COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
//KRATOS_CREATE_VARIABLE( bool,COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( INTERNAL_FORCE );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EXTERNAL_FORCE );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( CONTACT_FORCE );


KRATOS_CREATE_VARIABLE(  Matrix, CONSTITUTIVE_MATRIX );
KRATOS_CREATE_VARIABLE(  double, NORM_ISOCHORIC_STRESS );
KRATOS_CREATE_VARIABLE(  Matrix, DEFORMATION_GRADIENT );
KRATOS_CREATE_VARIABLE( double, DETERMINANT_F );
KRATOS_CREATE_VARIABLE(  double, NODE_AREA );
KRATOS_CREATE_VARIABLE( double, DomainSize );

KRATOS_CREATE_VARIABLE(double, PLASTIC_STRAIN);
KRATOS_CREATE_VARIABLE( double,NODE_EQUIVALENT_PLASTIC_STRAIN);



}