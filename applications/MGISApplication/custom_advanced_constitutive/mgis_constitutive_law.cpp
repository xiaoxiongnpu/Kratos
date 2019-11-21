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
#include "MGIS/Behaviour/Integrate.hxx"

// Project includes
#include "includes/checks.h"
#include "includes/mat_variables.h"
#include "custom_advanced_constitutive/mgis_constitutive_law.h"

namespace Kratos {

  namespace MGIS {

    template <class TMatrix>
    static void ConvertDeformationGradientToMGIS(double* const F,
                                                 const TMatrix& Fk,
                                                 const mgis::behaviour::Hypothesis h) {
      F[0] = Fk(0, 0);
      F[1] = Fk(1, 1);
      F[2] = Fk(2, 2);
      if ((h == mgis::behaviour::Hypothesis::AXISYMMETRICAL) ||
          (h == mgis::behaviour::Hypothesis::PLANESTRESS) ||
          (h == mgis::behaviour::Hypothesis::PLANESTRAIN) ||
          (h == mgis::behaviour::Hypothesis::GENERALISEDPLANESTRAIN)) {
        F[3] = Fk(0, 1);
        F[4] = Fk(1, 0);
      }
      if (h == mgis::behaviour::Hypothesis::TRIDIMENSIONAL) {
        F[3] = Fk(0, 1);
        F[4] = Fk(1, 0);
        F[5] = Fk(0, 2);
        F[6] = Fk(2, 0);
        F[7] = Fk(2, 1);
        F[8] = Fk(1, 2);
      }
    }  // end of ConvertDeformationGradientToMGIS

    template <class TVector>
    static void ConvertStrainToMGIS(double* const e, const TVector& ek) {
      // \note strain ordering conventions used by Kratos:
      // \f$ [ e11, e22, e33, 2*e12, 2*e23, 2*e13 ] \f$ for 3D case and
      // \f$ [ e11, e22, e33, 2*e12 ] \f$ for 2D case.
      // \f$ [ e11, e22, 2*e12 ] \f$ for 2D case.
      constexpr const double isqrt2 = 0.70710678118654752440;
      if (ek.size() == 6) {
        e[0] = ek[0];
        e[1] = ek[1];
        e[2] = ek[2];
        e[3] = ek[3] * isqrt2;
        e[4] = ek[5] * isqrt2;
        e[5] = ek[4] * isqrt2;
      } else if (ek.size() == 4) {
        e[0] = ek[0];
        e[1] = ek[1];
        e[2] = ek[2];
        e[3] = ek[3] * isqrt2;
      } else if (ek.size() == 3) {
        e[0] = ek[0];
        e[1] = ek[1];
        e[2] = double{};
        e[3] = ek[2] * isqrt2;
      } else {
        KRATOS_ERROR << "ConvertStrainToMGIS: unsupported tensor size\n";
      }
    }  // end of ConvertStrainToMGIS

    template <class TVector>
    static void ConvertStressToKratos(TVector& sk, double* const s) {
      // \note stress ordering conventions used by Kratos:
      // \f$ [ s11, s22, s33, s12, s23, s13 ] \f$ for 3D case and
      // \f$ [ s11, s22, s33, s12 ] \f$ for 2D case.
      // \f$ [ s11, s22, s12 ] \f$ for 2D case.
      constexpr const double isqrt2 = 0.70710678118654752440;
      if (sk.size() == 6) {
        sk[0] = s[0];
        sk[1] = s[1];
        sk[2] = s[2];
        sk[3] = s[3] * isqrt2;
        sk[4] = s[5] * isqrt2;
        sk[5] = s[4] * isqrt2;
      } else if (sk.size() == 4) {
        sk[0] = s[0];
        sk[1] = s[1];
        sk[2] = s[2];
        sk[3] = s[3] * isqrt2;
      } else if (sk.size() == 3) {
        sk[0] = s[0];
        sk[1] = s[1];
        sk[2] = s[3] * isqrt2;
      } else {
        KRATOS_ERROR << "ConvertStressToKratos: unsupported tensor size\n";
      }
    }  // end of ConvertStressToKratos

    template <class Matrix>
    static void ConvertTangentOperatorToKratos(Matrix& Kk,
                                               double* const K,
                                               const mgis::behaviour::Hypothesis h) {
      constexpr const double cste = 1.41421356237309504880;
      if (h == mgis::behaviour::Hypothesis::TRIDIMENSIONAL) {
        for (unsigned short i = 0; i != 6u; ++i) {
          for (unsigned short j = 0; j != 6u; ++j) {
            Kk(i, j) = K[i * 6 + j];
          }
        }
        Kk(0, 3) *= cste;
        Kk(1, 3) *= cste;
        Kk(2, 3) *= cste;
        Kk(0, 4) *= cste;
        Kk(1, 4) *= cste;
        Kk(2, 4) *= cste;
        Kk(0, 5) *= cste;
        Kk(1, 5) *= cste;
        Kk(2, 5) *= cste;
        Kk(3, 0) *= cste;
        Kk(3, 1) *= cste;
        Kk(3, 2) *= cste;
        Kk(4, 0) *= cste;
        Kk(4, 1) *= cste;
        Kk(4, 2) *= cste;
        Kk(5, 0) *= cste;
        Kk(5, 1) *= cste;
        Kk(5, 2) *= cste;
        Kk(3, 3) *= 2;
        Kk(3, 4) *= 2;
        Kk(3, 5) *= 2;
        Kk(4, 3) *= 2;
        Kk(4, 4) *= 2;
        Kk(4, 5) *= 2;
        Kk(5, 3) *= 2;
        Kk(5, 4) *= 2;
        Kk(5, 5) *= 2;
        // 
        for (unsigned short i = 0; i != 6u; ++i) {
          std::swap(Kk(4, i), Kk(5, i));
        }
        for (unsigned short i = 0; i != 6u; ++i) {
          std::swap(Kk(i, 4), Kk(i, 5));
        }
      } else {
        KRATOS_ERROR << "ConvertTangentOperatorToKratos: unsupported hypothesis";
      }
    }  // end of ConvertTangentOperatorToKratos
  }

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
  }  // end of MGISConstitutiveLaw::MGISConstitutiveLaw

  /******************************COPY CONSTRUCTOR*************************************/
  /***********************************************************************************/

  MGISConstitutiveLaw::MGISConstitutiveLaw(const MGISConstitutiveLaw& rOther) = default;

  /********************************CLONE**********************************************/
  /***********************************************************************************/

  ConstitutiveLaw::Pointer MGISConstitutiveLaw::Clone() const {
    return Kratos::make_shared<MGISConstitutiveLaw>(*this);
  }  // end of MGISConstitutiveLaw::Clone

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

  void MGISConstitutiveLaw::Integrate(ConstitutiveLaw::Parameters& rValues) {
    const auto& opts = rValues.GetOptions();
    const auto& pi = rValues.GetProcessInfo();
    // getting the gradients
    if (this->strain_measure == StrainMeasure_Infinitesimal) {
      if (opts.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        KRATOS_ERROR << "MGISConstitutiveLaw::Integrate: the strain tensor must be provided by the "
                        "element\n";
      }
      const auto& ek = rValues.GetStrainVector();
      Kratos::MGIS::ConvertStrainToMGIS(this->data.s1.gradients.data(), ek);
    } else if (this->strain_measure == StrainMeasure_Deformation_Gradient) {
      const auto& Fk = rValues.GetDeformationGradientF();
      Kratos::MGIS::ConvertDeformationGradientToMGIS(this->data.s1.gradients.data(), Fk,
                                                     this->behaviour->hypothesis);
    } else {
      KRATOS_ERROR << "MGISConstitutiveLaw::Integrate: unimplemented yet\n";
    }
    const auto k = opts.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if (opts.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
      this->data.K[0] = k ? 4 : 0;
    } else {
      if (!k) {
        KRATOS_ERROR << "MGISConstitutiveLaw::Integrate: don't know what to do, neither the stress "
                        "nor the tangent operator is requested\n";
      }
      // return the elastic operator
      this->data.K[0] = -1;
    }
    data.dt = pi[DELTA_TIME];
    auto data_view = mgis::behaviour::make_view(this->data);
    if (mgis::behaviour::integrate(data_view, *(this->behaviour)) != 0) {
      KRATOS_ERROR << "MGISConstitutiveLaw::Integrate: behaviour integration failed\n";
    }
  }  // end of MGISConstitutiveLaw::Integrate

  void MGISConstitutiveLaw::Update(ConstitutiveLaw::Parameters& rValues) {
    // copy the state at the end of the time step on the state at the beginning of the time step
    mgis::behaviour::update(data);
  }  // end of MGISConstitutiveLaw::Update

  void MGISConstitutiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) {
    KRATOS_TRY;
    const auto& opts = rValues.GetOptions();
    this->Integrate(rValues);
    if (opts.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
      // exporting the stress
      auto& pk2 = rValues.GetStressVector();
      Kratos::MGIS::ConvertStressToKratos(pk2, this->data.s1.thermodynamic_forces.data());
    }
    if (opts.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
      auto& K = rValues.GetConstitutiveMatrix();
      Kratos::MGIS::ConvertTangentOperatorToKratos(K, this->data.K.data(),
                                                   this->behaviour->hypothesis);
    }
    KRATOS_CATCH("");
  }

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) {
    KRATOS_ERROR << "PK1 stress measure is not supported yet\n";
  }  // end of MGISConstitutiveLaw::CalculateMaterialResponsePK1

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::CalculateMaterialResponseKirchhoff(
      ConstitutiveLaw::Parameters& rValues) {
    KRATOS_TRY;
    if (this->strain_measure == StrainMeasure_Deformation_Gradient) {
      KRATOS_ERROR << "Kirchhoff stress measure is not supported yet\n";
    }
    this->Integrate(rValues);
    KRATOS_CATCH("");
  }  // end of MGISConstitutiveLaw::CalculateMaterialResponseKirchhoff

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) {
    KRATOS_TRY;
    if (this->strain_measure == StrainMeasure_Deformation_Gradient) {
      KRATOS_ERROR << "Kirchhoff stress measure is not supported yet\n";
    }
    this->Integrate(rValues);
    KRATOS_CATCH("");
  }  // end of MGISConstitutiveLaw::CalculateMaterialResponseCauchy

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) {
  }  // end of MGISConstitutiveLaw::InitializeMaterialResponsePK1

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) {
  }  // end of MGISConstitutiveLaw::InitializeMaterialResponsePK2

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::InitializeMaterialResponseKirchhoff(
      ConstitutiveLaw::Parameters& rValues) {
  }  // end of  MGISConstitutiveLaw::InitializeMaterialResponseKirchhoff

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) {
  }  // end of  MGISConstitutiveLaw::InitializeMaterialResponseCauchy

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) {
    this->Update(rValues);
  }  // end of MGISConstitutiveLaw::FinalizeMaterialResponsePK1

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) {
    this->Update(rValues);
  }  // end of MGISConstitutiveLaw::FinalizeMaterialResponsePK2

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::FinalizeMaterialResponseKirchhoff(
      ConstitutiveLaw::Parameters& rValues) {
    this->Update(rValues);
  }  // end of MGISConstitutiveLaw::FinalizeMaterialResponseKirchhoff

  /***********************************************************************************/
  /***********************************************************************************/

  void MGISConstitutiveLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) {
    this->Update(rValues);
  }  // end of MGISConstitutiveLaw::FinalizeMaterialResponseCauchy

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
    //     if (rThisVariable == STRAIN || rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
    //         rThisVariable == ALMANSI_STRAIN_VECTOR) {
    //       // ???
    //       //      this->CalculateCauchyGreenStrain(rParameterValues, rValue);
    //     } else if (rThisVariable == STRESSES || rThisVariable == CAUCHY_STRESS_VECTOR ||
    //                rThisVariable == KIRCHHOFF_STRESS_VECTOR || rThisVariable ==
    //                PK2_STRESS_VECTOR) {
    //       // Get Values to compute the constitutive law:
    //       auto& r_flags = rParameterValues.GetOptions();
    //
    //       // Previous flags saved
    //       const auto flag_const_tensor =
    //       r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    //       const auto flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    //
    //       r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    //       r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    //
    //       // We compute the stress
    //       MGISConstitutiveLaw::CalculateMaterialResponseCauchy(rParameterValues);
    //       rValue = rParameterValues.GetStressVector();
    //
    //       // Previous flags restored
    //       r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    //       r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    //     }

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
