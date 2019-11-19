// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Thomas Helfer
//                   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/mat_variables.h"
#include "custom_advanced_constitutive/mgis_constitutive_law.h"

namespace Kratos {

  /******************************CONSTRUCTOR******************************************/
  /***********************************************************************************/

  MGISConstitutiveLaw::MGISConstitutiveLaw(const Kratos::shared_ptr<Behaviour>& b)
      : ConstitutiveLaw(), behaviour(b), data(*(this->behaviour)) {
    if ((this->behaviour->btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) &&
        (this->behaviour->kinematic == Behaviour::SMALLSTRAINKINEMATIC)) {
      this->strain_measure = StrainMeasure_Infinitesimal;
      this->stress_measure = StressMeasure_Cauchy;
    } else if ((this->behaviour->btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) &&
               (this->behaviour->kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY)) {
      this->strain_measure = StrainMeasure_Deformation_Gradient;
      this->stress_measure = StressMeasure_Cauchy;
    } else {
      KRATOS_ERROR << "unsupported behaviour type or behaviour kinematic\n";
    }
  }

  /******************************COPY CONSTRUCTOR*************************************/
  /***********************************************************************************/

  MGISConstitutiveLaw::MGISConstitutiveLaw(const MGISConstitutiveLaw& rOther) = default;

  /********************************CLONE**********************************************/
  /***********************************************************************************/

  ConstitutiveLaw::Pointer MGISConstitutiveLaw::Clone() const {
    return Kratos::make_shared<MGISConstitutiveLaw>(*this);
  }

  /*******************************DESTRUCTOR******************************************/
  /***********************************************************************************/

  MGISConstitutiveLaw::~MGISConstitutiveLaw() = default;

  /***********************************************************************************/
  /***********************************************************************************/

  MGISConstitutiveLaw::SizeType MGISConstitutiveLaw::WorkingSpaceDimension() {
    return mgis::behaviour::getSpaceDimension(this->behaviour->hypothesis);
  }  // end of MGISConstitutiveLaw::WorkingSpaceDimension

  MGISConstitutiveLaw::SizeType MGISConstitutiveLaw::GetStrainSize() {
    if (this->strain_measure == StrainMeasure_Deformation_Gradient) {
      return mgis::behaviour::getTensorSize(this->behaviour->hypothesis);
    } else if (this->strain_measure == StrainMeasure_Infinitesimal) {
      return mgis::behaviour::getStensorSize(this->behaviour->hypothesis);
    }
    KRATOS_ERROR << "unsupported behaviour\n";
    return 0;
  }  // end of MGISConstitutiveLaw::GetStrainSize

  MGISConstitutiveLaw::StrainMeasure MGISConstitutiveLaw::GetStrainMeasure() {
    return this->strain_measure;
  }

  MGISConstitutiveLaw::StressMeasure MGISConstitutiveLaw::GetStressMeasure() {
    return this->stress_measure;
  }

  void MGISConstitutiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) {
    KRATOS_TRY;
    //     // b.- Get Values to compute the constitutive law:
    //     Flags& r_constitutive_law_options = rValues.GetOptions();
    //
    //     Vector& r_strain_vector = rValues.GetStrainVector();
    //
    //     // NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN
    //     // MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    //     if (r_constitutive_law_options.IsNot(
    //             ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //       CalculateCauchyGreenStrain(rValues, r_strain_vector);
    //     }
    //
    //     if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
    //       Vector& r_stress_vector = rValues.GetStressVector();
    //       CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
    //     }
    //
    //     if (r_constitutive_law_options.Is(
    //             ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
    //       Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    //       CalculateElasticMatrix(r_constitutive_matrix, rValues);
    //     }
    //
    KRATOS_CATCH("");
  }

  /***********************************************************************************/
  /***********************************************************************************/

  // NOTE: Since we are in the hypothesis of small strains we can use the same
  // function for everything

  void MGISConstitutiveLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) {}

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::CalculateMaterialResponseKirchhoff(
      ConstitutiveLaw::Parameters& rValues) {
    //    CalculateMaterialResponsePK2(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) {
    //    CalculateMaterialResponsePK2(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) {
    //     // Small deformation so we can call the Cauchy method
    //     InitializeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) {
    //     // Small deformation so we can call the Cauchy method
    //     InitializeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) {
    // TODO: Add if necessary
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::InitializeMaterialResponseKirchhoff(
      ConstitutiveLaw::Parameters& rValues) {
    // Small deformation so we can call the Cauchy method
    //    InitializeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) {
    // Small deformation so we can call the Cauchy method
    //    FinalizeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) {
    // Small deformation so we can call the Cauchy method
    //    FinalizeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) {
    // TODO: Add if necessary
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::FinalizeMaterialResponseKirchhoff(
      ConstitutiveLaw::Parameters& rValues) {
    // Small deformation so we can call the Cauchy method
    //    FinalizeMaterialResponseCauchy(rValues);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  double& MGISConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                              const Variable<double>& rThisVariable,
                                              double& rValue) {
    //     Vector& r_strain_vector = rParameterValues.GetStrainVector();
    //     Vector& r_stress_vector = rParameterValues.GetStressVector();
    //
    //     if (rThisVariable == STRAIN_ENERGY) {
    //       this->CalculateCauchyGreenStrain(rParameterValues,
    //       r_strain_vector);
    //       this->CalculatePK2Stress(r_strain_vector, r_stress_vector,
    //                                rParameterValues);
    //
    //       rValue = 0.5 * inner_prod(r_strain_vector,
    //                                 r_stress_vector);  // Strain energy =
    //                                 0.5*E:C:E
    //     }
    //
    return (rValue);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  Vector& MGISConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                              const Variable<Vector>& rThisVariable,
                                              Vector& rValue) {
    if (rThisVariable == STRAIN || rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {
      // ???
      //      this->CalculateCauchyGreenStrain(rParameterValues, rValue);
    } else if (rThisVariable == STRESSES || rThisVariable == CAUCHY_STRESS_VECTOR ||
               rThisVariable == KIRCHHOFF_STRESS_VECTOR || rThisVariable == PK2_STRESS_VECTOR) {
      // Get Values to compute the constitutive law:
      Flags& r_flags = rParameterValues.GetOptions();

      // Previous flags saved
      const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
      const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

      r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
      r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

      // We compute the stress
      MGISConstitutiveLaw::CalculateMaterialResponseCauchy(rParameterValues);
      rValue = rParameterValues.GetStressVector();

      // Previous flags restored
      r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
      r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    }

    return (rValue);
  }

  /***********************************************************************************/
  /***********************************************************************************/

  Matrix& MGISConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                              const Variable<Matrix>& rThisVariable,
                                              Matrix& rValue) {
    //     if (rThisVariable == CONSTITUTIVE_MATRIX ||
    //         rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
    //         rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
    //       this->CalculateElasticMatrix(rValue, rParameterValues);
    //     }

    return (rValue);
  }

  //*************************CONSTITUTIVE LAW GENERAL FEATURES
  //*************************
  /***********************************************************************************/

  void MGISConstitutiveLaw::GetLawFeatures(Features& rFeatures) {
    // Set the type of law
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    if (this->behaviour->symmetry == Behaviour::ISOTROPIC) {
      rFeatures.mOptions.Set(ISOTROPIC);
    } else if (this->behaviour->symmetry == Behaviour::ORTHOTROPIC) {
      rFeatures.mOptions.Set(ANISOTROPIC);
    }
    // Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(this->GetStrainMeasure());
    // Set the strain size
    rFeatures.mStrainSize = this->GetStrainSize();
    // Set the spacedimension
    rFeatures.mSpaceDimension = this->WorkingSpaceDimension();
  }

  /***********************************************************************************/
  /***********************************************************************************/

  int MGISConstitutiveLaw::Check(const Properties& rMaterialProperties,
                                 const GeometryType& rElementGeometry,
                                 const ProcessInfo& rCurrentProcessInfo) {
    return 0;
  }

}  // Namespace Kratos
