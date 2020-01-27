//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_WALL_ELEMENT_UTILITIES_H_INCLUDED)
#define KRATOS_WALL_ELEMENT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

namespace WallElementUtilities
{
template <unsigned int TDim>
void CalculateRotationMatrix(BoundedMatrix<double, TDim, TDim>& rRotationMatrix,
                             const array_1d<double, 3>& rNormal);

template <unsigned int TDim>
array_1d<double, 3> GetRotatedVector(const array_1d<double, 3>& rVector,
                                     const BoundedMatrix<double, TDim, TDim>& rRotationMatrix);

array_1d<double, 3> GetWallPosition(const ModelPart::ElementType::GeometryType& rGeometry);

array_1d<double, 3> GetGaussPosition(const ModelPart::ElementType::GeometryType& rGeometry,
                                     const Vector& rGaussShapeFunctions);

} // namespace WallElementUtilities

///@}
///@name Kratos Classes
///@{
///@}

} // namespace Kratos

#endif