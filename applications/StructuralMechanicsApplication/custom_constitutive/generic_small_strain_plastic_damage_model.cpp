// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Sergio Jimenez
//  
//

// System includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_constitutive/generic_small_strain_plastic_damage_model.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress Damage
    Vector& r_integrated_stress_vector = rValues.GetStressVector();
    const double characteristic_length = ConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
    array_1d<double, VoigtSize> auxiliar_integrated_stress_vector = r_integrated_stress_vector;
    Matrix& r_tangent_tensor = rValues.GetConstitutiveMatrix(); // todo modify after integration
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if( r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
        this->CalculateValue(rValues, STRAIN, r_strain_vector);
    }

    // Elastic Matrix
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);
    }

    // We compute the stress
    if(r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }

        // Converged values
        double threshold_plasticity = mThresholdPlasticity;
        double threshold_damage = mThresholdDamage;
        double damage = mDamage;
        double plastic_dissipation = mPlasticDissipation;
        Vector plastic_strain = mPlasticStrain;

        // Stress Predictor S = (1-d)C:(E-Ep)
        array_1d<double, VoigtSize> predictive_stress_vector = (1.0 - damage) * prod(r_constitutive_matrix, r_strain_vector - plastic_strain);

        // Initialize Plastic Parameters
        double uniaxial_stress_plasticity = 0.0, plastic_denominator = 0.0, uniaxial_stress_damage = 0.0;
        BoundedArrayType f_flux = ZeroVector(VoigtSize); // DF/DS
        BoundedArrayType g_flux = ZeroVector(VoigtSize); // DG/DS
        BoundedArrayType plastic_strain_increment = ZeroVector(VoigtSize);

        // Compute the plastic parameters
        const double plasticity_indicator = TPlasticityIntegratorType::CalculatePlasticParameters(
                                        predictive_stress_vector, r_strain_vector, uniaxial_stress_plasticity,
                                        threshold_plasticity, plastic_denominator, f_flux, g_flux,
                                        plastic_dissipation, plastic_strain_increment,
                                        r_constitutive_matrix, rValues, characteristic_length,
                                        plastic_strain);

        // Compute Damage Parameters
        TDamageIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress_damage, rValues);
        const double damage_indicator = uniaxial_stress_damage - threshold_damage;

        // Verification threshold for the plastic-damage process
        if (plasticity_indicator >= std::abs(1.0e-4 * threshold_plasticity) && damage_indicator >= std::abs(1.0e-4 * threshold_damage)) {
            array_1d<double, VoigtSize> damage_yield_flux = ZeroVector(VoigtSize);
            TDamageIntegratorType::CalculateYieldSurfaceDerivative(predictive_stress_vector, damage_yield_flux, rValues);

            const double scalar_prod_dam_yield_sigma_eff = inner_prod(damage_yield_flux, predictive_stress_vector / (1.0 - damage));
            const double scalar_prod_plast_yield_sigma_eff = inner_prod(f_flux, predictive_stress_vector / (1.0 - damage));

            double innerprod_dam_yield_elastic_tensor = 0.0, HKG = 0.0, fact1 = 0.0;
            Vector hcapa = ZeroVector(VoigtSize);
            for (IndexType i = 0; i < VoigtSize; ++i) {
                for (IndexType j = 0; j < VoigtSize; ++j) {
                    innerprod_dam_yield_elastic_tensor += damage_yield_flux[j] * r_constitutive_matrix(i,j);
                }
                hcapa[i] = predictive_stress_vector[i] / uniaxial_stress_plasticity;
                HKG +=  hcapa[i] * g_flux[i];
                fact1 += innerprod_dam_yield_elastic_tensor * g_flux[i]; // ?????
            }
            fact1 *= (1.0 - damage);
            const double factorA = scalar_prod_plast_yield_sigma_eff;
            const double factorB = 1 / plastic_denominator; // ?? ABETA
            //const double factorC = scalar_prod_dam_yield_sigma_eff + hard_damage_slope;

        }


    }

} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    ProcessInfo dummy_process_info;
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, dummy_process_info);

    // // We call the integrator
    double initial_threshold_plast, initial_threshold_damage;
    TPlasticityIntegratorType::GetInitialUniaxialThreshold(aux_param, initial_threshold_plast);
    TDamageIntegratorType::GetInitialUniaxialThreshold(aux_param, initial_threshold_damage);
    mThresholdPlasticity = initial_threshold_plast;
    mThresholdDamage = initial_threshold_damage;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{











}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
bool GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
bool GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        return true;
    } else {
        return BaseType::Has(rThisVariable);
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
bool GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::Has(const Variable<Matrix>& rThisVariable)
{
    return BaseType::Has(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        mPlasticDissipation = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        mPlasticStrain = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
double& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == PLASTIC_DISSIPATION) {
        rValue = mPlasticDissipation;
    } else {
        BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
Vector& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_VECTOR) {
        rValue = mPlasticStrain;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
Matrix& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == PLASTIC_STRAIN_TENSOR) {
        rValue = MathUtils<double>::StrainVectorToTensor(mPlasticStrain);
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
double& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // todo
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
Vector& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
Matrix& GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (this->Has(rThisVariable)) {
        return this->GetValue(rThisVariable, rValue);
    } else {
        return BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TPlasticityIntegratorType, class TDamageIntegratorType>
int GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    const int check_integrator_plasticity = TPlasticityIntegratorType::Check(rMaterialProperties);
	const int check_integrator_damage = TDamageIntegratorType::Check(rMaterialProperties);
    KRATOS_ERROR_IF_NOT(VoigtSize == this->GetStrainSize()) << "You are combining not compatible constitutive laws" << std::endl;
    if ((check_base + check_integrator_plasticity + check_integrator_damage) > 0) return 1;
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void CalculateDamageParameters(
    array_1d<double, 6>& rPredictiveStressVector,
    Vector& rStrainVector,
    double& rUniaxialStress,
    double& rThreshold,
    double& rDamageDissipation,
    const Matrix& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues,
    const double CharacteristicLength,  
    array_1d<double, 6>& rFflux,
    const Vector& rPlasticStrain,
    const double Damage,
    const double DamageIncrement
)
{
    array_1d<double, VoigtSize> deviator = ZeroVector(6);
    array_1d<double, VoigtSize> h_capa = ZeroVector(6);
    double J2, tensile_indicator_factor, compression_indicator_factor, slope, hardening_parameter, suma;

    YieldSurfaceType::CalculateEquivalentStress( rPredictiveStressVector, rStrainVector, rUniaxialStress, rValues);
    const double I1 = rPredictiveStressVector[0] + rPredictiveStressVector[1] + rPredictiveStressVector[2];
    ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
    YieldSurfaceType::CalculateYieldSurfaceDerivative(rPredictiveStressVector, deviator, J2, rFflux, rValues);
    CalculateIndicatorsFactors(rPredictiveStressVector, tensile_indicator_factor, compression_indicator_factor, suma);
    // TPlasticityIntegratorType::CalculatePlasticDissipation(rPredictiveStressVector, tensile_indicator_factor,compression_indicator_factor, rPlasticStrainIncrement,rPlasticDissipation, h_capa, rValues, CharacteristicLength);
    // TPlasticityIntegratorType::CalculateHardeningParameter(rFflux, slope, h_capa, hardening_parameter);

    auto& r_matProps = rValues.GetMaterialProperties();
    const bool has_symmetric_yield_stress = r_matProps.Has(YIELD_STRESS);
    const double yield_compression = has_symmetric_yield_stress ? r_matProps[YIELD_STRESS] : r_matProps[YIELD_STRESS_COMPRESSION];
    const double yield_tension = has_symmetric_yield_stress ? r_matProps[YIELD_STRESS] : r_matProps[YIELD_STRESS_TENSION];
    const double yield_ratio = yield_compression / yield_tension;

    const double fracture_energy_tension = r_matProps[FRACTURE_ENERGY] / CharacteristicLength;
    const double fracture_energy_compression = fracture_energy_tension * std::pow(yield_ratio, 2.0);

    const double constant1 = tensile_indicator_factor * std::abs(rUniaxialStress / yield_ratio) / (fracture_energy_tension * suma);
    const double constant2 = compression_indicator_factor * std::abs(rUniaxialStress) / (fracture_energy_compression * suma);
    const double constant = constant1 + constant2;

    // Free Energy Undamaged
    const double free_energy_undamaged = 0.5 * (inner_prod(rStrainVector - rPlasticStrain), rPredictiveStressVector / (1.0 - Damage));
    const double hcapd = constant * free_energy_undamaged;
    const double damage_dissipation_increment = hcapd * DamageIncrement;
    if (damage_dissipation_increment > 1.0 || damage_dissipation_increment <tolerance) damage_dissipation_increment = 0.0;
    rDamageDissipation += damage_dissipation_increment;
    if (rDamageDissipation > 1.0) rDamageDissipation = 0.99999;
    Vector slopes(2), thresholds(2);

}


/***********************************************************************************/
/***********************************************************************************/
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
void CalculateIndicatorsFactors(
    const array_1d<double, 6>& rPredictiveStressVector,
    double& rTensileIndicatorFactor,
    double& rCompressionIndicatorFactor,
    const double SumPrincipalStresses
)
{
        // We do an initial check
        if (norm_2(rPredictiveStressVector) < 1.0e-8) {
            rTensileIndicatorFactor = 0.0;
            rCompressionIndicatorFactor = 0.0;
            return;
        }

        // We proceed as usual
        array_1d<double, Dimension> principal_stresses = ZeroVector(Dimension);
        ConstitutiveLawUtilities<VoigtSize>::CalculatePrincipalStresses(principal_stresses, rPredictiveStressVector);

        double suma = 0.0, sumb = 0.0, sumc = 0.0;
        double aux_sa;

        for (IndexType i = 0; i < Dimension; ++i) {
            aux_sa = std::abs(principal_stresses[i]);
            suma += aux_sa;
            sumb += 0.5 * (principal_stresses[i] + aux_sa);
            sumc += 0.5 * (-principal_stresses[i] + aux_sa);
        }
        SumPrincipalStresses = suma;

        if (std::abs(suma) > tolerance) {
            rTensileIndicatorFactor = sumb / suma;
            rCompressionIndicatorFactor = sumc / suma;
        } else {
            rTensileIndicatorFactor = sumb;
            rCompressionIndicatorFactor = sumc;
        }

        // Final check
        if ((std::abs(rTensileIndicatorFactor) + std::abs(rCompressionIndicatorFactor)) < tolerance) {
            rTensileIndicatorFactor = 0.0;
            rCompressionIndicatorFactor = 0.0;
            return;
        }
}



/***********************************************************************************/
/***********************************************************************************/
template class GenericSmallStrainPlasticDamageModel<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;


} // namespace Kratos
