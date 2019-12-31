//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah (https://github.com/sdharmin)
//                   Bence Rochlitz (https://github.com/bencerochlitz)
//
//  Supervised by:   Jordi Cotela (https://github.com/jcotela)
//                   Suneth Warnakulasuriya (https://github.com/sunethwarna)

#include "evm_k_omega_utilities.h"
#include <cmath>
#include <iostream>
#include <limits>

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
double CalculateTurbulentViscosity(const double turbulent_kinetic_energy,
	const double turbulent_specific_energy_dissipation_rate)
{
	return turbulent_kinetic_energy / turbulent_specific_energy_dissipation_rate;
}

//we need to calculate our own parameters, which are: Beta and Beta_Star --> parameters of the CDR equations.
template <unsigned int TDim>
double CalculateBeta(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
	const double TurbulentSpecificDissipationRate,
	const double BetaZero,
	const double BetaZeroStar) //omega element
{
	BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
	noalias(symmetric_velocity_gradient) = 0.5 * (rVelocityGradient + trans(rVelocityGradient));

	BoundedMatrix<double, TDim, TDim> antisymmetric_velocity_gradient;
	noalias(antisymmetric_velocity_gradient) = 0.5 * (rVelocityGradient - trans(rVelocityGradient));
	
	// x_omega = abs(Omega_ij * Omega_jk * S_ki / (BetaZeroStar * omega)^3)
	BoundedMatrix<double, TDim, TDim> antisymmetric_velocity_gradient_product;
	noalias(antisymmetric_velocity_gradient_product) = prod(antisymmetric_velocity_gradient, trans(antisymmetric_velocity_gradient));

	double x_omega = 0.0;
	for (unsigned int i = 0; i < TDim; ++i)
		for (unsigned int j = 0; j < TDim; ++j)
			x_omega += antisymmetric_velocity_gradient_product(i, j) * symmetric_velocity_gradient(i, j);
	x_omega = std::abs(x_omega / std::pow(BetaZeroStar * TurbulentSpecificDissipationRate, 3));
	
	double f_beta = (1 + 70 * x_omega) / (1 + 80 * x_omega);

	return BetaZero * f_beta;
}

double CalculateBetaStar(const double Omega,
						const array_1d<double, 3> &OmegaGradient,
						const array_1d<double, 3> &TkeGradient,
						const double BetaZeroStar) //--> ask Suneth
{
	double X_k = (TkeGradient[0] * OmegaGradient[0] + 
					TkeGradient[1] * OmegaGradient[1] +
					TkeGradient[2] * OmegaGradient[2])
					/ std::pow(Omega, 3);

	double f_beta_star;
	if (X_k <= 0)
	{
		f_beta_star = 1.0;
	}
	else
	{
		f_beta_star = (1 + 680 * X_k * X_k) / (1 + 400 * X_k * X_k);
	}

	return BetaZeroStar * f_beta_star;
}



template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                           const double turbulent_kinematic_viscosity,
                           const double turbulent_kinetic_energy)
{
    const double velocity_divergence = 
        RansCalculationUtilities::CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);
    
    BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
    noalias(symmetric_velocity_gradient) = rVelocityGradient + trans(rVelocityGradient);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;

	noalias(reynolds_stress_tensor) =		
		// turbulent_kinematic_viscosity *
		// (symmetric_velocity_gradient - (2.0 / 3.0) * velocity_divergence * identity);
		turbulent_kinematic_viscosity *
		(symmetric_velocity_gradient - (2.0 / 3.0) * velocity_divergence * identity - (2.0 / 3.0) * turbulent_kinetic_energy * identity);


	double source = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            {
				source += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);
			}
    return source;
} // calculates "P" tensor for the production term - elementwise multiplication

double CalculateTheta(const double turbulent_kinematic_viscosity)
{
    return 1.0 / turbulent_kinematic_viscosity;
}


template double CalculateSourceTerm<2>(const BoundedMatrix<double, 2, 2>&, const double, const double);
template double CalculateSourceTerm<3>(const BoundedMatrix<double, 3, 3>&, const double, const double);

template double CalculateBeta<3>(const BoundedMatrix<double, 3, 3>&,
	const double,
	const double,
	const double); //omega element

template double CalculateBeta<2>(const BoundedMatrix<double, 2, 2>&,
	const double,
	const double,
	const double);
} // namespace EvmKOmegaModelUtilities

///@}

} // namespace Kratos
