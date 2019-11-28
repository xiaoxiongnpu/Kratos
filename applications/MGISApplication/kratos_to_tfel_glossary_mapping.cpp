//   __  __  _____ _____  _____                     _ _           _   _
//  |  \/  |/ ____|_   _|/ ____|  /\               | (_)         | | (_)
//  | \  / | |  __  | | | (___   /  \   _ __  _ __ | |_  ___ __ _| |_ _  ___  _ __
//  | |\/| | | |_ | | |  \___ \ / /\ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//  | |  | | |__| |_| |_ ____) / ____ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//  |_|  |_|\_____|_____|_____/_/    \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                     | |   | |
//                                     |_|   |_|
//
//  License: BSD License
//   license: MGISApplication/license.txt
//
//  Main authors:  Thomas Helfer, Vicente Mataix Ferrandiz
//

// System includes
#include <map>

// External includes

// Project includes
#include "includes/mat_variables.h"
// Application includes
#include "structural_mechanics_application_variables.h"
#include "kratos_to_tfel_glossary_mapping.h"

namespace Kratos {

  namespace MGIS {

    /**
     * @brief class in charge of mapping a TFEL glossary names to Kratos' variable name.
     */
    struct KratosTFELMapping
        : private std::map<std::string, std::string>,
          private std::map<std::string, Kratos::Variable<double>>,
          private std::map<std::string, Kratos::Variable<Kratos::Vector>>,
          private std::map<Kratos::Variable<Kratos::Vector>::KeyType, StensorStorageConvention> {
      //! \return the unique instance of this class
      static const KratosTFELMapping& Get() {
        static const KratosTFELMapping i;
        return i;
      }  // end of get

      const std::map<std::string, std::string>& GetNamesMapping() const {
        return static_cast<const std::map<std::string, std::string>&>(*this);
      }  // end of GetNamesMapping

      template <typename Type>
      const std::map<std::string, Kratos::Variable<Type>>& GetVariablesMapping() const {
        return static_cast<const std::map<std::string, Kratos::Variable<Type>>&>(*this);
      }  // end of GetVariablesMapping

      StensorStorageConvention GetStensorStorageConvention(const Kratos::Variable<Kratos::Vector>& v) const {
        const auto& m = static_cast<
            const std::map<Kratos::Variable<double>::KeyType, StensorStorageConvention>&>(*this);
        const auto p = m.find(v.Key());
        if (p == m.end()) {
          KRATOS_ERROR << "no storage convention associated with variable '" << v.Name() << "'\n";
        }
        return p->second;
      } // end of GetStensorStorageConvention

     private:
      //! \brief default constructor
      KratosTFELMapping() {
        // external state variables
        this->Add("Temperature", Kratos::TEMPERATURE);
        // material properties
        this->Add("YoungModulus", Kratos::YOUNG_MODULUS);
        this->Add("YoungModulus1", Kratos::YOUNG_MODULUS_X);
        this->Add("YoungModulus2", Kratos::YOUNG_MODULUS_Y);
        this->Add("YoungModulus3", Kratos::YOUNG_MODULUS_Z);
        this->Add("PoissonRatio", Kratos::POISSON_RATIO);
        this->Add("PoissonRatio12", Kratos::POISSON_RATIO_XY);
        this->Add("PoissonRatio13", Kratos::POISSON_RATIO_XZ);
        this->Add("PoissonRatio23", Kratos::POISSON_RATIO_YZ);
        this->Add("ShearModulus", Kratos::SHEAR_MODULUS);
        this->Add("ShearModulus12", Kratos::SHEAR_MODULUS_XY);
        this->Add("ShearModulus13", Kratos::SHEAR_MODULUS_XZ);
        this->Add("ShearModulus23", Kratos::SHEAR_MODULUS_YZ);
        // internal state variables
        this->Add("Damage", Kratos::DAMAGE);
        this->Add("EquivalentPlasticStrain", Kratos::EQUIVALENT_PLASTIC_STRAIN);
        this->Add("PlasticStrain", Kratos::PLASTIC_STRAIN, StensorStorageConvention::STRAIN);
      }  // end of KratosTFELMapping

      /**
       * @brief add an new mapping
       * @tparam Type: type of the Kratos' variable
       * @param[in] k: TFEL glossary name
       * @param[in] k: Kratos' variable
       */
      template <typename Type>
      void Add(const std::string& k, const Variable<Type>& v) {
        auto& names_mapping = static_cast<std::map<std::string, std::string>&>(*this);
        auto& variables_mapping = static_cast<std::map<std::string, Variable<Type>>&>(*this);
        names_mapping.insert({k, v.Name()});
        variables_mapping.insert({k, v});
      }  // end of add

      /**
       * @brief add an new mapping
       * @tparam Type: type of the Kratos' variable
       * @param[in] k: TFEL glossary name
       * @param[in] k: Kratos' variable
       * @param[in] c: storage conventions
       */
      template <typename Type>
      void Add(const std::string& k, const Variable<Type>& v, const StensorStorageConvention c) {
        auto& m = static_cast<
            std::map<Kratos::Variable<Kratos::Vector>::KeyType, StensorStorageConvention>&>(*this);
        this->Add(k, v);
        if (!m.insert({v.Key(), c}).second) {
          KRATOS_ERROR << "storage convention has already been defined for variable '" << v.Name()
                       << "'";
        }
      }  // end of add

      KratosTFELMapping(KratosTFELMapping&&) = delete;
      KratosTFELMapping(const KratosTFELMapping&) = delete;
      KratosTFELMapping& operator=(KratosTFELMapping&&) = delete;
      KratosTFELMapping& operator=(const KratosTFELMapping&) = delete;
    };  // end of struct KratosTFELMapping

    std::string GetKratosVariableName(const std::string& g) {
      const auto& m = KratosTFELMapping::Get().GetNamesMapping();
      auto p = m.find(g);
      if (p != m.end()) {
        return p->second;
      }
      return "";
    }  // end of getKratosVariableName

    template <typename Type>
    static const Variable<Type>& GetKratosVariableImpl(const std::string& g) {
      const auto& m = KratosTFELMapping::Get().GetVariablesMapping<Type>();
      auto p = m.find(g);
      if (p == m.end()) {
        KRATOS_ERROR << "Kratos::MGIS::GetKratosVariableImpl: "
                     << "no variable associated for glossary name '" << g << "'";
      }
      return p->second;
    }  // end of getKratosVariableName

    template <>
    const Variable<double>& GetKratosVariable<double>(const std::string& g) {
      return GetKratosVariableImpl<double>(g);
    }  // end of GetKratosVariable<double>

    template <>
    const Variable<Kratos::Vector>& GetKratosVariable<Kratos::Vector>(const std::string& g) {
      return GetKratosVariableImpl<Kratos::Vector>(g);
    }  // end of GetKratosVariable<Kratos::Vector>

    StensorStorageConvention GetStensorStorageConvention(const Kratos::Variable<Kratos::Vector>& v) {
      const auto& m = KratosTFELMapping::Get();
      return m.GetStensorStorageConvention(v);
    } // end of GetStensorStorageConvention

  }  // end of namespace MGIS

}  // end of namespace Kratos

