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

#if !defined(KRATOS_MGIS_CONSTITUTIVE_LAW_H_INCLUDED)
#define KRATOS_MGIS_CONSTITUTIVE_LAW_H_INCLUDED

// System includes

// External includes
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"

// Project includes
#include "includes/constitutive_law.h"

// Application includes
#include "mgis_constitutive_law_description.h"

namespace Kratos {
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
   * @class MGISConstitutiveLaw
   * @ingroup MGISApplication
   * @brief This class defines a constitutive
   * law imported through MGIS
   * @details This class derives from the base constitutive law
   * @author Thomas Helfer
   * @author Riccardo Rossi
   * @author Vicente Mataix Ferrandiz
   */
  class KRATOS_API(MGIS_APPLICATION) MGISConstitutiveLaw : public ConstitutiveLaw {
   public:
    ///@name Type Definitions
    ///@{

    /// The process info type definition
    typedef ProcessInfo ProcessInfoType;

    /// The base class ConstitutiveLaw type definition
    typedef ConstitutiveLaw BaseType;

    /// The size type definition
    typedef std::size_t SizeType;

    /// A simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /// A simple alias
    using Behaviour = mgis::behaviour::Behaviour;
    /// A simple alias
    using BehaviourData = mgis::behaviour::BehaviourData;

    /// Counted pointer of MGISConstitutiveLaw
    KRATOS_CLASS_POINTER_DEFINITION(MGISConstitutiveLaw);

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief constructor.
     * @param[in] b: MGIS behaviour description
     */
    MGISConstitutiveLaw(const Kratos::shared_ptr<Kratos::MGIS::MGISConstitutiveLawDescription>&);

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    MGISConstitutiveLaw(const MGISConstitutiveLaw& rOther);

    /**
     * @brief Destructor.
     */
    ~MGISConstitutiveLaw() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called once to check compatibility
     * with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
    * @brief Dimension of the law:
    */
    SizeType WorkingSpaceDimension() override;

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override;

    /**
     * @brief Returns the expected strain measure of this constitutive law (by
     * default Green-Lagrange)
     * @return the expected strain measure
     */
    StrainMeasure GetStrainMeasure() override;

    /**
     * @brief Returns the stress measure of this constitutive law (by default
     * 2st Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
    StressMeasure GetStressMeasure() override;

    /**
     * @brief Computes the material response:
     * @details PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response:
     * @details PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response:
     * @details Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response:
     * @details Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 1st Piola-Kirchhoff
     * stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff
     * stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
      * @brief Updates the material response:
      * @details Cauchy stresses and Internal Variables
      * @param rValues The internal values of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Updates the material response:
     * @details Cauchy stresses and Internal Variables
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Updates the material response:
     * @details Cauchy stresses and Internal Variables
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief It calculates the value of a specified variable (double case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<double>& rThisVariable,
                           double& rValue) override;

    /**
     * @brief It calculates the value of a specified variable (Vector case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<Vector>& rThisVariable,
                           Vector& rValue) override;

    /**
     * @brief It calculates the value of a specified variable (Matrix case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                           const Variable<Matrix>& rThisVariable,
                           Matrix& rValue) override;

    /**
     * @brief This function provides the place to perform checks on the
     * completeness of the input.
     * @details It is designed to be called only once (or anyway, not often)
     * typically at the beginning of the calculations, so to verify that nothing
     * is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(const Properties& rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo& rCurrentProcessInfo) override;
    ///@}

   protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    Kratos::shared_ptr<Kratos::MGIS::MGISConstitutiveLawDescription> behaviour;
    BehaviourData data;
    StrainMeasure strain_measure;
    StressMeasure stress_measure;

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

    /**
     * @brief integrate the behaviour over the time step.
     * @param rValues The internal values of the law
     * @see   Parameters
     * @note  This method does not update the stress and
     * the tangent operator on the Kratos side. This must
     * be done in the calling methods
     * (`CalculateMaterialResponseCauchy` for example)
     * which are responsible for converting the stress and
     * the tangent operator returned by `MGIS` (which is
     * respectively the second Piola-Kirchhoff stress and
     * its derivative with respect to the Green-Lagrange
     * strain.
     */
    void Integrate(ConstitutiveLaw::Parameters& rValues);

    /**
     * @brief update the material state.
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void Update(ConstitutiveLaw::Parameters& rValues);

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }
    ///@}

  };  // Class MGISConstitutiveLaw
}  // namespace Kratos.

#endif  // KRATOS_MGIS_CONSTITUTIVE_LAW_H_INCLUDED  defined
