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
#include "rans_application.h"
#include "geometries/line_2d_2.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "rans_application_variables.h"

namespace Kratos
{
KratosRANSApplication::KratosRANSApplication()
    : KratosApplication("RANSApplication"),
      mRansEvmKEpsilonLowReK2D(0,
                               Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                                   Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonLowReK3D(0,
                               Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                                   Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonLowReEpsilon2D(
          0,
          Element::GeometryType::Pointer(
              new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonLowReEpsilon3D(
          0,
          Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
              Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonK2D(0,
                          Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                              Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonK3D(0,
                          Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                              Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonEpsilon2D(0,
                                Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                                    Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonEpsilon3D(0,
                                Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                                    Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKEpsilonEpsilonWall2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonEpsilonWall3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKEpsilonVmsMonolithicWall2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmKEpsilonVmsMonolithicWall3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
    // k-omega elements
      mRansEvmOmega2D(0,
                      Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                          Element::GeometryType::PointsArrayType(3)))),
      mRansEvmOmega3D(0,
                      Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                          Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaK2D(0,
                        Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(
                            Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaK3D(0,
                        Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(
                            Element::GeometryType::PointsArrayType(4)))),
      mRansEvmKOmegaOmegaWall2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmKOmegaOmegaWall3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
      mRansEvmKOmegaOmegaWallBlended2D2N(
          0,
          Element::GeometryType::Pointer(
              new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
      mRansEvmKOmegaOmegaWallBlended3D3N(
          0,
          Element::GeometryType::Pointer(
              new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3))))                        
{
}

void KratosRANSApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosRANSApplication..." << std::endl;

    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE)
    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY_RATE)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE_2)
    KRATOS_REGISTER_VARIABLE(IS_CO_SOLVING_PROCESS_ACTIVE)
    KRATOS_REGISTER_VARIABLE(RANS_Y_PLUS)
    KRATOS_REGISTER_VARIABLE(RANS_AUXILIARY_VARIABLE_1)
    KRATOS_REGISTER_VARIABLE(RANS_AUXILIARY_VARIABLE_2)
    KRATOS_REGISTER_VARIABLE(WALL_SMOOTHNESS_BETA)
    KRATOS_REGISTER_VARIABLE(WALL_VON_KARMAN)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C_MU)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C1)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_C2)
    KRATOS_REGISTER_VARIABLE(TURBULENT_KINETIC_ENERGY_SIGMA)
    KRATOS_REGISTER_VARIABLE(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA)
    KRATOS_REGISTER_VARIABLE(NUMBER_OF_NEIGHBOUR_CONDITIONS)

    //for the k-omega element
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_SIGMA_K)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_SIGMA_OMEGA)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_GAMMA)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_BETA_ZERO_STAR)
    KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_BETA_ZERO)
    KRATOS_REGISTER_VARIABLE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE)
	KRATOS_REGISTER_VARIABLE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2)
	KRATOS_REGISTER_VARIABLE(TURBULENCE_RANS_Y_PLUS_LIMIT_WALL)
	// KRATOS_REGISTER_VARIABLE(TURBULENCE_BLENDING)

    // Register Elements
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReK2D3N", mRansEvmKEpsilonLowReK2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReK3D4N", mRansEvmKEpsilonLowReK3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReEpsilon2D3N", mRansEvmKEpsilonLowReEpsilon2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonLowReEpsilon3D4N", mRansEvmKEpsilonLowReEpsilon3D);

    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonK2D3N", mRansEvmKEpsilonK2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonK3D4N", mRansEvmKEpsilonK3D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilon2D3N", mRansEvmKEpsilonEpsilon2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKEpsilonEpsilon3D4N", mRansEvmKEpsilonEpsilon3D);

    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonWall2D2N", mRansEvmKEpsilonEpsilonWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonEpsilonWall3D3N", mRansEvmKEpsilonEpsilonWall3D3N);

    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonVmsMonolithicWall2D2N",
                              mRansEvmKEpsilonVmsMonolithicWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKEpsilonVmsMonolithicWall3D3N",
                              mRansEvmKEpsilonVmsMonolithicWall3D3N);

    // k-omega elements
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmega2D3N", mRansEvmOmega2D);
    KRATOS_REGISTER_ELEMENT("RansEvmKOmegaOmega3D4N", mRansEvmOmega3D);
	KRATOS_REGISTER_ELEMENT("RansEvmKOmegaK2D3N", mRansEvmKOmegaK2D);
	KRATOS_REGISTER_ELEMENT("RansEvmKOmegaK3D4N", mRansEvmKOmegaK3D);
    // k-omega conditions
    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaWall2D2N", mRansEvmKOmegaOmegaWall2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaWall3D3N", mRansEvmKOmegaOmegaWall3D3N);
    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaWallBlended2D2N", mRansEvmKOmegaOmegaWallBlended2D2N);
    KRATOS_REGISTER_CONDITION("RansEvmKOmegaOmegaWallBlended3D3N", mRansEvmKOmegaOmegaWallBlended3D3N);
}
} // namespace Kratos.
