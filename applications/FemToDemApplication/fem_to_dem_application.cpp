//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"

#include "fem_to_dem_application.h"
#include "fem_to_dem_application_variables.h"


namespace Kratos {

KratosFemToDemApplication::KratosFemToDemApplication(): KratosApplication("FemToDemApplication"),
mSmallStrainModifiedMohrCoulombFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainModifiedMohrCoulombFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainRankineFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainRankineFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainSimoJuFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainSimoJuFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainDruckerPragerFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainDruckerPragerFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainVonMisesFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainVonMisesFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainTrescaFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainTrescaFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainMohrCoulombFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainMohrCoulombFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mLargeDisplacementModifiedMohrCoulombFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mLargeDisplacementModifiedMohrCoulombFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mLargeDisplacementRankineFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mLargeDisplacementRankineFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mLargeDisplacementSimoJuFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mLargeDisplacementSimoJuFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mLargeDisplacementDruckerPragerFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mLargeDisplacementDruckerPragerFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mLargeDisplacementVonMisesFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mLargeDisplacementVonMisesFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mLargeDisplacementTrescaFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mLargeDisplacementTrescaFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mLargeDisplacementMohrCoulombFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mLargeDisplacementMohrCoulombFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4))))
{}

void KratosFemToDemApplication::Register() 
{
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
	
	//REGISTER VARIABLES FEM2DEM
	KRATOS_REGISTER_VARIABLE(VOLUME_COUNTED)
	KRATOS_REGISTER_VARIABLE(ERASED_VOLUME)
	KRATOS_REGISTER_VARIABLE(COHESION)
	KRATOS_REGISTER_VARIABLE(RECOMPUTE_NEIGHBOURS)
	KRATOS_REGISTER_VARIABLE(GENERATE_DEM)
	KRATOS_REGISTER_VARIABLE(DISPLACEMENT_INCREMENT)
	KRATOS_REGISTER_VARIABLE(DAMAGE_ELEMENT)
	KRATOS_REGISTER_VARIABLE(TIME_UNIT_CONVERTER)
	KRATOS_REGISTER_VARIABLE(STRESS_VECTOR)
	KRATOS_REGISTER_VARIABLE(YIELD_STRESS_C)
	KRATOS_REGISTER_VARIABLE(YIELD_STRESS_T)
	KRATOS_REGISTER_VARIABLE(FRAC_ENERGY_T)
	KRATOS_REGISTER_VARIABLE(FRAC_ENERGY_C)
	KRATOS_REGISTER_VARIABLE(INTERNAL_PRESSURE_ITERATION)
	KRATOS_REGISTER_VARIABLE(STRESS_VECTOR_INTEGRATED)
	KRATOS_REGISTER_VARIABLE(THRESHOLD)
	KRATOS_REGISTER_VARIABLE(SMOOTHED_STRESS_VECTOR)
	KRATOS_REGISTER_VARIABLE(YIELD_SURFACE)
	KRATOS_REGISTER_VARIABLE(STRAIN_VECTOR)
	KRATOS_REGISTER_VARIABLE(SMOOTHING)
	KRATOS_REGISTER_VARIABLE(IS_DAMAGED)
	KRATOS_REGISTER_VARIABLE(TANGENT_CONSTITUTIVE_TENSOR)
	KRATOS_REGISTER_VARIABLE(RECONSTRUCT_PRESSURE_LOAD)
	KRATOS_REGISTER_VARIABLE(IS_DYNAMIC)
	KRATOS_REGISTER_VARIABLE(STRESS_THRESHOLD)
	KRATOS_REGISTER_VARIABLE(INTEGRATION_COEFFICIENT)
	KRATOS_REGISTER_VARIABLE(MAPPING_PROCEDURE)
	KRATOS_REGISTER_VARIABLE(INITIAL_THRESHOLD)
	KRATOS_REGISTER_VARIABLE(IS_DEM)
	KRATOS_REGISTER_VARIABLE(DEM_RADIUS)
	KRATOS_REGISTER_VARIABLE(DEM_GENERATED)
	KRATOS_REGISTER_VARIABLE(INACTIVE_NODE)
	KRATOS_REGISTER_VARIABLE(NUMBER_OF_ACTIVE_ELEMENTS)
	KRATOS_REGISTER_VARIABLE(NODAL_FORCE_APPLIED)
	KRATOS_REGISTER_VARIABLE(NODAL_FORCE_X)
	KRATOS_REGISTER_VARIABLE(NODAL_FORCE_Y)
	KRATOS_REGISTER_VARIABLE(NODAL_FORCE_Z)
	KRATOS_REGISTER_VARIABLE(NODAL_STRESS_VECTOR)
	KRATOS_REGISTER_VARIABLE(EQUIVALENT_NODAL_STRESS)
	KRATOS_REGISTER_VARIABLE(PRESSURE_EXPANDED)
	KRATOS_REGISTER_VARIABLE(IS_SKIN)
	KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EQUIVALENT_NODAL_STRESS_GRADIENT)
	KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUXILIAR_GRADIENT)

	KRATOS_REGISTER_VARIABLE(STRAIN_TENSOR);
	KRATOS_REGISTER_VARIABLE(STRESS_TENSOR);
	KRATOS_REGISTER_VARIABLE(STRESS_TENSOR_INTEGRATED);
	
	// Composite
	KRATOS_REGISTER_VARIABLE(CONCRETE_STRESS_TENSOR);
	KRATOS_REGISTER_VARIABLE(STEEL_STRESS_TENSOR);
	KRATOS_REGISTER_VARIABLE(CONCRETE_STRESS_VECTOR);
	KRATOS_REGISTER_VARIABLE(STEEL_STRESS_VECTOR);
	KRATOS_REGISTER_VARIABLE(YOUNG_MODULUS_STEEL);
	KRATOS_REGISTER_VARIABLE(DENSITY_STEEL);
	KRATOS_REGISTER_VARIABLE(POISSON_RATIO_STEEL);
	KRATOS_REGISTER_VARIABLE(STEEL_VOLUMETRIC_PART);
	KRATOS_REGISTER_VARIABLE(CONCRETE_STRESS_TENSOR_INTEGRATED);
	
	KRATOS_REGISTER_VARIABLE(YIELD_STRESS_C_STEEL);
	KRATOS_REGISTER_VARIABLE(YIELD_STRESS_T_STEEL);
	KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_STEEL);
	KRATOS_REGISTER_VARIABLE(PLASTIC_DISSIPATION_CAPAP);
	KRATOS_REGISTER_VARIABLE(EQUIVALENT_STRESS_VM);
	KRATOS_REGISTER_VARIABLE(NODAL_DAMAGE);
	KRATOS_REGISTER_VARIABLE(IS_TAKEN);
	KRATOS_REGISTER_VARIABLE(PRESSURE_ID);

	// Hardening variables plasticity
	KRATOS_REGISTER_VARIABLE(HARDENING_LAW);
	KRATOS_REGISTER_VARIABLE(MAXIMUM_STRESS);
	KRATOS_REGISTER_VARIABLE(MAXIMUM_STRESS_POSITION);	
	
	KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD);
	KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD);
	
	//Register element
	KRATOS_REGISTER_ELEMENT("SmallStrainModifiedMohrCoulombFemDemElement2D", mSmallStrainModifiedMohrCoulombFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("SmallStrainModifiedMohrCoulombFemDemElement3D", mSmallStrainModifiedMohrCoulombFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("SmallStrainRankineFemDemElement2D", mSmallStrainRankineFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("SmallStrainRankineFemDemElement3D", mSmallStrainRankineFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("SmallStrainSimoJuFemDemElement2D", mSmallStrainSimoJuFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("SmallStrainSimoJuFemDemElement3D", mSmallStrainSimoJuFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("SmallStrainDruckerPragerFemDemElement2D", mSmallStrainDruckerPragerFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("SmallStrainDruckerPragerFemDemElement3D", mSmallStrainDruckerPragerFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("SmallStrainVonMisesFemDemElement2D", mSmallStrainVonMisesFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("SmallStrainVonMisesFemDemElement3D", mSmallStrainVonMisesFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("SmallStrainTrescaFemDemElement2D", mSmallStrainTrescaFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("SmallStrainTrescaFemDemElement3D", mSmallStrainTrescaFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("SmallStrainMohrCoulombFemDemElement3D", mSmallStrainMohrCoulombFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("SmallStrainMohrCoulombFemDemElement2D", mSmallStrainMohrCoulombFemDemElement2D)

	KRATOS_REGISTER_ELEMENT("LargeDisplacementModifiedMohrCoulombFemDemElement2D", mLargeDisplacementModifiedMohrCoulombFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementModifiedMohrCoulombFemDemElement3D", mLargeDisplacementModifiedMohrCoulombFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementRankineFemDemElement2D", mLargeDisplacementRankineFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementRankineFemDemElement3D", mLargeDisplacementRankineFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementSimoJuFemDemElement2D", mLargeDisplacementSimoJuFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementSimoJuFemDemElement3D", mLargeDisplacementSimoJuFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementDruckerPragerFemDemElement2D", mLargeDisplacementDruckerPragerFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementDruckerPragerFemDemElement3D", mLargeDisplacementDruckerPragerFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementVonMisesFemDemElement2D", mLargeDisplacementVonMisesFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementVonMisesFemDemElement3D", mLargeDisplacementVonMisesFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementTrescaFemDemElement2D", mLargeDisplacementTrescaFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementTrescaFemDemElement3D", mLargeDisplacementTrescaFemDemElement3D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementMohrCoulombFemDemElement2D", mLargeDisplacementMohrCoulombFemDemElement2D)
	KRATOS_REGISTER_ELEMENT("LargeDisplacementMohrCoulombFemDemElement3D", mLargeDisplacementMohrCoulombFemDemElement3D)
	
	//Register Constitutive Laws
	Serializer::Register("ZarateLaw", mZarateLaw);
	Serializer::Register("FemDemElasticLaw", mFemDemElasticLaw);

	KRATOS_REGISTER_CONSTITUTIVE_LAW("FemDemElasticLaw", mFemDemElasticLaw)

	

}

}  // namespace Kratos.