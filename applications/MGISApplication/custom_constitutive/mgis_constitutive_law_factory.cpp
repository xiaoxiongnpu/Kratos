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
#include "MGIS/LibrariesManager.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"

// Project includes
#include "custom_constitutive/mgis_constitutive_law.h"
#include "custom_constitutive/mgis_constitutive_law_factory.h"

namespace Kratos {

  ConstitutiveLaw::Pointer MGISConstitutiveLawFactory::Create(
      Kratos::Parameters NewParameters) const {
    using mgis::behaviour::Behaviour;
    using mgis::behaviour::FiniteStrainBehaviourOptions;
    const auto hypothesis =
        NewParameters.Has("hypothesis")
            ? mgis::behaviour::fromString(NewParameters["hypothesis"].GetString())
            : mgis::behaviour::Hypothesis::TRIDIMENSIONAL;
    const auto& library = NewParameters["library"].GetString();
    const auto& name = NewParameters["behaviour"].GetString();
    // mgis behaviour
    auto b = Kratos::shared_ptr<Behaviour>{};
    const auto btype = [&library, &name] {
      auto& lm = mgis::LibrariesManager::get();
      /* - 0 : general behaviour
       * - 1 : strain based behaviour *
       * - 2 : standard finite strain behaviour *
       * - 3 : cohesive zone model */
      switch (lm.getBehaviourType(library, name)) {
        case 0:
          return Behaviour::GENERALBEHAVIOUR;
        case 1:
          return Behaviour::STANDARDSTRAINBASEDBEHAVIOUR;
        case 2:
          return Behaviour::STANDARDFINITESTRAINBEHAVIOUR;
        case 3:
          return Behaviour::COHESIVEZONEMODEL;
      }
      KRATOS_ERROR << "unsupported behaviour type\n";
    }();
    if (btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      auto opts = FiniteStrainBehaviourOptions{};
      opts.stress_measure = FiniteStrainBehaviourOptions::StressMeasure::PK2;
      opts.tangent_operator = FiniteStrainBehaviourOptions::TangentOperator::DS_DEGL;
      b = Kratos::make_shared<Behaviour>(mgis::behaviour::load(opts, library, name, hypothesis));
    } else {
      b = Kratos::make_shared<Behaviour>(mgis::behaviour::load(library, name, hypothesis));
    }
    return Kratos::make_shared<MGISConstitutiveLaw>(b);
  }  // end of MGISConstitutiveLawFactory::Create

}  // end of  namespace Kratos
