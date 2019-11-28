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

#if !defined(KRATOS_MGIS_CONSTITUTIVE_LAW_DESCRIPTION_H_INCLUDED)
#define KRATOS_MGIS_CONSTITUTIVE_LAW_DESCRIPTION_H_INCLUDED

// System includes
#include <map>

// External includes
#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos {

  namespace MGIS {

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

    ///@}
    ///@name Kratos Classes
    ///@{

    /**
     * @class MGISPositionContainer
     * @ingroup MGISApplication
     * @brief A class mapping a `Kratos`' variable to its position
     */
    template <typename Type>
    struct MGISPositionContainer
        : public std::map<typename Variable<Type>::KeyType, mgis::size_type> {
    };  // end of struct MGISPositionContainer

    /**
     * @class MGISConstitutiveLawDescription
     * @ingroup MGISApplication
     * @brief This class gathers information to ease the interoperation of a `MGIS` constitutive
     * law
     * in `Kratos`.
     * @author Thomas Helfer
     * @author Riccardo Rossi
     * @author Vicente Mataix Ferrandiz
     */
    class KRATOS_API(MGIS_APPLICATION) MGISConstitutiveLawDescription
        : protected MGISPositionContainer<double>,
          public mgis::behaviour::Behaviour {
     public:
      ///@name Type Definitions
      ///@{

      ///@}
      ///@name Lyfe Cycle
      ///@{

      MGISConstitutiveLawDescription(const mgis::behaviour::FiniteStrainBehaviourOptions,
                                     const std::string&,
                                     const std::string&,
                                     const mgis::behaviour::Hypothesis);

      MGISConstitutiveLawDescription(const std::string&,
                                     const std::string&,
                                     const mgis::behaviour::Hypothesis);

      ///@}
      ///@name Public member Variables
      ///@{

      ///@}
      ///@name Operators
      ///@{

      ///@}
      ///@name Operations
      ///@{

      template <typename Type>
      std::pair<bool, mgis::size_type> GetInternalStateVariablePosition(
          const Variable<Type>& v) const {
        const auto& c = static_cast<const MGISPositionContainer<double>&>(*this);
        auto p = c.find(v.Key());
        if (p == c.end()) {
          return {false, mgis::size_type{}};
        }
        return {true, p->second};
      }  // end of GetInternalStateVariablePosition

      ///@}

     protected:
      ///@name Protected static Member Variables
      ///@{

      ///@}
      ///@name Protected member Variables
      ///@{

      ///@}
      ///@name Protected Operators
      ///@{

      ///@}
      ///@name Protected Operations
      ///@{

      ///@}

     private:
      ///@name Static Member Variables
      ///@{

      ///@}
      ///@name Member Variables
      ///@{

      ///@}
      ///@name Private Operators
      ///@{

      ///@}
      ///@name Private Operations
      ///@{

      void init();

      template <typename Type>
      void SetInternalStateVariablePosition(const Variable<Type>& v, const mgis::size_type p) {
        auto& c = static_cast<MGISPositionContainer<double>&>(*this);
        if (!c.insert({v.Key(), p}).second) {
          throwAlreadyRegistredVariable(v.Name());
        }
      }  // end of SetInternalStateVariablePosition

      void throwAlreadyRegistredVariable(const std::string&) const;

      ///@}
      ///@name Private  Access
      ///@{
      ///@}

      ///@}
      ///@name Serialization
      ///@{
      ///@}

    };  // end of class MGISConstitutiveLawDescription

  }  // end of MGIS

}  // namespace Kratos.

#endif  // KRATOS_MGIS_CONSTITUTIVE_LAW_DESCRIPTION_H_INCLUDED  defined
