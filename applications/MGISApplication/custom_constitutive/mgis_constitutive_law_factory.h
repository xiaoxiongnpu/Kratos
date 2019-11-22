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

#if !defined(KRATOS_MGIS_CONSTITUTIVE_LAW_FACTORY_H_INCLUDED)
#define KRATOS_MGIS_CONSTITUTIVE_LAW_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

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

///@}
///@name Kratos Classes
///@{
/**
 * @class MGISConstitutiveLawFactory
 * @ingroup StructuralMechanicsApplication
 * @brief Dummy class to register, only implements the `create` method
 * @author Thomas Helfer
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) MGISConstitutiveLawFactory
    : public ConstitutiveLaw
{
  public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(MGISConstitutiveLawFactory);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    MGISConstitutiveLawFactory()
    {
    }

    /**
    * Destructor.
    */
    ~MGISConstitutiveLawFactory() override
    {
    }

    /**
     * @brief Creates a new constitutive law pointer
     * @param NewParameters The configuration parameters of the new constitutive law
     * @return a Pointer to the new constitutive law
     */
    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

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
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
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

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif
