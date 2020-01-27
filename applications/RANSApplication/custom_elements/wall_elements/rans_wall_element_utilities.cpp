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

// System includes

// External includes

// Project includes

// Application includes

// Include base h
#include "rans_wall_element_utilities.h"

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
template <>
void CalculateRotationMatrix<2>(BoundedMatrix<double, 2, 2>& rRotationMatrix,
                                const array_1d<double, 3>& rNormal)
{
    const double normal_magnitude = norm_2(rNormal);

    rRotationMatrix(0, 0) = rNormal[0] / normal_magnitude;
    rRotationMatrix(0, 1) = rNormal[1] / normal_magnitude;
    rRotationMatrix(1, 0) = -rNormal[1] / normal_magnitude;
    rRotationMatrix(1, 1) = rNormal[0] / normal_magnitude;
}

template <>
void CalculateRotationMatrix<3>(BoundedMatrix<double, 3, 3>& rRotationMatrix,
                                const array_1d<double, 3>& rNormal)
{
    // Get the normal evaluated at the node
    const double normal_magnitude = norm_2(rNormal);

    rRotationMatrix(0, 0) = rNormal[0] / normal_magnitude;
    rRotationMatrix(0, 1) = rNormal[1] / normal_magnitude;
    rRotationMatrix(0, 2) = rNormal[2] / normal_magnitude;
    // Define the new coordinate system, where the first vector is aligned with the normal

    // To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
    array_1d<double, 3> rT1;
    rT1(0) = 1.0;
    rT1(1) = 0.0;
    rT1(2) = 0.0;
    double dot = rRotationMatrix(0, 0); // this->Dot(rN,rT1);

    // It is possible that the normal is aligned with (1,0,0), resulting in
    // norm(rT1) = 0 If this is the case, repeat the procedure using (0,1,0)
    if (fabs(dot) > 0.99)
    {
        rT1(0) = 0.0;
        rT1(1) = 1.0;
        rT1(2) = 0.0;

        dot = rRotationMatrix(0, 1); // this->Dot(rN,rT1);
    }

    // calculate projection and normalize
    rT1[0] -= dot * rRotationMatrix(0, 0);
    rT1[1] -= dot * rRotationMatrix(0, 1);
    rT1[2] -= dot * rRotationMatrix(0, 2);
    noalias(rT1) = rT1 * (1.0 / norm_2(rT1));

    rRotationMatrix(1, 0) = rT1[0];
    rRotationMatrix(1, 1) = rT1[1];
    rRotationMatrix(1, 2) = rT1[2];

    // The third base component is choosen as N x T1, which is normalized by construction
    rRotationMatrix(2, 0) =
        rRotationMatrix(0, 1) * rT1[2] - rRotationMatrix(0, 2) * rT1[1];
    rRotationMatrix(2, 1) =
        rRotationMatrix(0, 2) * rT1[0] - rRotationMatrix(0, 0) * rT1[2];
    rRotationMatrix(2, 2) =
        rRotationMatrix(0, 0) * rT1[1] - rRotationMatrix(0, 1) * rT1[0];
}

template <unsigned int TDim>
array_1d<double, 3> GetRotatedVector(const array_1d<double, 3>& rVector,
                                     const BoundedMatrix<double, TDim, TDim>& rRotationMatrix)
{
    array_1d<double, 3> value = ZeroVector(3);
    for (unsigned int i = 0; i < TDim; ++i)
    {
        for (unsigned int j = 0; j < TDim; ++j)
        {
            value[i] += rVector[j] * rRotationMatrix(i, j);
        }
    }
    return value;
}

array_1d<double, 3> GetWallPosition(const ModelPart::ElementType::GeometryType& rGeometry)
{
    array_1d<double, 3> wall_position = ZeroVector(3);

    const int number_of_nodes = rGeometry.PointsNumber();
    int number_of_wall_nodes = 0;
    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const ModelPart::NodeType& r_node = rGeometry[i_node];
        if (r_node.Is(STRUCTURE))
        {
            noalias(wall_position) += r_node.Coordinates();
            ++number_of_wall_nodes;
        }
    }

    noalias(wall_position) =
        wall_position * (1.0 / static_cast<double>(number_of_wall_nodes));

    return wall_position;
}

array_1d<double, 3> GetGaussPosition(const ModelPart::ElementType::GeometryType& rGeometry,
                                     const Vector& rGaussShapeFunctions)
{
    array_1d<double, 3> gauss_position = ZeroVector(3);
    const int number_of_nodes = rGeometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        noalias(gauss_position) =
            rGeometry[i_node].Coordinates() * rGaussShapeFunctions[i_node];
    }

    return gauss_position;
}

// template instantiations
template array_1d<double, 3> GetRotatedVector<2>(const array_1d<double, 3>&,
                                                 const BoundedMatrix<double, 2, 2>&);

template array_1d<double, 3> GetRotatedVector<3>(const array_1d<double, 3>&,
                                                 const BoundedMatrix<double, 3, 3>&);

} // namespace WallElementUtilities

///@}
///@name Kratos Classes
///@{
///@}

} // namespace Kratos
