//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Masó Sotomayor
//


#if !defined(KRATOS_SHALLOW_WAER_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_SHALLOW_WAER_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

namespace Kratos
{
KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, BATHYMETRY)                           // Geometric definition of the problem
KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, RAIN)                                 // Source term

KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, HEIGHT)                               // Main variable
KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, PROJECTED_HEIGHT)                     // Convected variable
KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, DELTA_HEIGHT)                         // Variable to update particles
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SHALLOW_WATER_APPLICATION, PROJECTED_VELOCITY )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SHALLOW_WATER_APPLICATION, DELTA_VELOCITY )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SHALLOW_WATER_APPLICATION, PROJECTED_MOMENTUM )
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SHALLOW_WATER_APPLICATION, DELTA_MOMENTUM )

KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, MEAN_SIZE)                            // Specific variable for PFEM2
KRATOS_DEFINE_APPLICATION_VARIABLE( SHALLOW_WATER_APPLICATION, double, MEAN_VEL_OVER_ELEM_SIZE)              // Specific variable for PFEM2
}

#endif	/* KRATOS_SHALLOW_WAER_APPLICATION_VARIABLES_H_INCLUDED */

