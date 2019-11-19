// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Thomas Helfer
//  Collaborator:    Riccardo Rossi, Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"

// Project includes
#include "custom_advanced_constitutive/mgis_constitutive_law_factory.h"

namespace Kratos
{
ConstitutiveLaw::Pointer MGISConstitutiveLawFactory::Create(Kratos::Parameters NewParameters) const
{
  using mgis::behaviour::Behaviour;
  const auto hypothesis = NewParameters.Has("hypothesis")
                              ? mgis::behaviour::fromString(NewParameters["hypothesis"].GetString())
                              : mgis::behaviour::Hypothesis::TRIDIMENSIONAL;
  const auto& library = NewParameters["library"].GetString();
  const auto& name = NewParameters["behaviour"].GetString();
  // mgis behaviour
  auto b = Kratos::shared_ptr<Behaviour>{};
  // finite strain behaviour options, if any
  if (NewParameters.Has("stress_measure") || NewParameters.Has("tangent_operator")) {
    using FiniteStrainBehaviourOptions = mgis::behaviour::FiniteStrainBehaviourOptions;
    auto opts = FiniteStrainBehaviourOptions{};
    if (NewParameters.Has("stress_measure")) {
      const auto sm = NewParameters["stress_measure"].GetString();
      if (sm == "Cauchy") {
        opts.stress_measure = FiniteStrainBehaviourOptions::StressMeasure::CAUCHY;
      } else if ((sm == "PK2") || (sm == "first_Piola_Kirchhoff_stress")) {
        opts.stress_measure = FiniteStrainBehaviourOptions::StressMeasure::PK1;
      } else if ((sm == "PK2") || (sm == "second_Piola_Kirchhoff_stress")) {
        opts.stress_measure = FiniteStrainBehaviourOptions::StressMeasure::PK2;
      } else {
        KRATOS_ERROR << "MGISConstitutiveLawFactory::Create: invalid stress measure '" << sm
                     << "'\n";
      }
    }
    if (NewParameters.Has("tangent_operator")) {
      const auto to = NewParameters["tangent_operator"].GetString();
      if (to == "DSIG_DF") {
        opts.tangent_operator = FiniteStrainBehaviourOptions::TangentOperator::DSIG_DF;
      } else if (to == "DPK1_DF") {
        opts.tangent_operator = FiniteStrainBehaviourOptions::TangentOperator::DPK1_DF;
      } else if (to == "DS_DEGL") {
        opts.tangent_operator = FiniteStrainBehaviourOptions::TangentOperator::DS_DEGL;
      } else {
        KRATOS_ERROR << "MGISConstitutiveLawFactory::Create: invalid tangent operator '" << to
                     << "'\n";
      }
    }
    // loading the behaviour
    b = make_shared<Behaviour>(mgis::behaviour::load(opts, library, name, hypothesis));
  } else {
    // loading the behaviour
    b = make_shared<Behaviour>(mgis::behaviour::load(library, name, hypothesis));
  }
  //   return KratosComponents<ConstitutiveLaw>::Get(name).Clone();
}

} // namespace Kratos
