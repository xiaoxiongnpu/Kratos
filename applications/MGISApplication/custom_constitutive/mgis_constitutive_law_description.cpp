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

// External includes

// Project includes
#include "kratos_to_tfel_glossary_mapping.h"
#include "custom_constitutive/mgis_constitutive_law_description.h"

namespace Kratos {

  namespace MGIS {

    MGISConstitutiveLawDescription::MGISConstitutiveLawDescription(
        const mgis::behaviour::FiniteStrainBehaviourOptions opts,
        const std::string& l,
        const std::string& b,
        const mgis::behaviour::Hypothesis h)
        : mgis::behaviour::Behaviour(mgis::behaviour::load(opts, l, b, h)) {
      this->init();
    }  // end of MGISConstitutiveLawDescription::MGISConstitutiveLawDescription

    MGISConstitutiveLawDescription::MGISConstitutiveLawDescription(
        const std::string& l, const std::string& b, const mgis::behaviour::Hypothesis h)
        : mgis::behaviour::Behaviour(mgis::behaviour::load(l, b, h)) {
      this->init();
    }  // end of MGISConstitutiveLawDescription::MGISConstitutiveLawDescription

    void MGISConstitutiveLawDescription::init() {
      // treating internal state variables
      auto pos = mgis::size_type{};
      for (const auto& iv : this->isvs) {
        const auto& n = Kratos::MGIS::GetKratosVariableName(iv.name);
        if (n.empty()) {
          pos += mgis::behaviour::getVariableSize(iv, this->hypothesis);
          continue;
        }
        if (iv.type == mgis::behaviour::Variable::SCALAR) {
          const auto& kv = Kratos::MGIS::GetKratosVariable<double>(iv.name);
          this->SetInternalStateVariablePosition<double>(kv, pos);
        } else if (iv.type == mgis::behaviour::Variable::STENSOR) {
          const auto& kv = Kratos::MGIS::GetKratosVariable<Kratos::Vector>(iv.name);
          this->SetInternalStateVariablePosition<Kratos::Vector>(kv, pos);
        } else {
          KRATOS_ERROR << "unsupported variable type for variable '" << iv.name << "'\n";
        }
        pos += mgis::behaviour::getVariableSize(iv, this->hypothesis);
      }
    }  // end of MGISConstitutiveLawDescription::init

    void MGISConstitutiveLawDescription::throwAlreadyRegistredVariable(const std::string& n) const {
      KRATOS_ERROR << "MGISConstitutiveLawDescription::SetInternalStateVariablePosition: "
                   << "variable '" << n << "' is already registred\n";
    }  // end of MGISConstitutiveLawDescription::throwAlreadyRegistredVariable

  }  // // end of namespace MGIS

}  // end of namespace Kratos
