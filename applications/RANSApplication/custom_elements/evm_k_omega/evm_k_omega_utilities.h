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

#if !defined(KRATOS_EVM_K_OMEGA_UTILITIES_H_INCLUDED)
#define KRATOS_EVM_K_OMEGA_UTILITIES_H_INCLUDED //ALL Epsilons renamed for omega

// System includes

// Project includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/ublas_interface.h"

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

namespace EvmKomegaModelUtilities
{
	//change calculation in this one:
double CalculateTurbulentViscosity(const double turbulent_kinetic_energy,
		const double turbulent_specific_energy_dissipation_rate);


// these are changing parameters from the k-e model:
//double CalculateFmu(const double y_plus);
//
//double CalculateF2(const double turbulent_kinetic_energy,
//                   const double kinematic_viscosity,
//                   const double turbulent_energy_dissipation_rate);

//we need to calculate our own parameters, which are: Beta and Beta_Star --> parameters of the CDR equations.
template <unsigned int TDim>
double CalculateBeta(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
					const double TurbulentSpecificDissipationRate,
					const double BetaZero,
					const double BetaZeroStar); //define the input parameters

double CalculateBetaStar(const double Omega,
							const array_1d<double, 3> &OmegaGradient,
							const array_1d<double, 3> &TkeGradient,
							const double BetaZeroStar); //define the input parameters



template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                           const double turbulent_kinematic_viscosity,
                           const double turbulent_kinetic_energy);

double CalculateTheta(const double turbulent_kinematic_viscosity);


} // namespace EvmKomegaModelUtilities

///@}

} // namespace Kratos

#endif