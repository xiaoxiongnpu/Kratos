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

// Application includes
#include "kratos_to_tfel_glossary_mapping.h"

namespace Kratos {

  /**
   * @brief class in charge of mappy a TFEL glossary names to Kratos' variable name.
   */
  struct KratosToTFELMapping : std::map<std::string, std::string> {
    //! \return the unique instance of this class
    static const KratosToTFELMapping& get() {
      static const KratosToTFELMapping i;
      return i;
    }  // end of get
   private:
    //! \brief default constructor
    KratosToTFELMapping() {
      // external state variables
      add("Temperature", TEMPERATURE);
      // material properties
      add("YoungModulus", YOUNG_MODULUS);
      add("YoungModulus1", YOUNG_MODULUS_X);
      add("YoungModulus2", YOUNG_MODULUS_Y);
      add("YoungModulus3", YOUNG_MODULUS_Z);
      add("PoissonRatio", POISSON_RATIO);
      add("PoissonRatio12", POISSON_RATIO_XY);
      add("PoissonRatio13", POISSON_RATIO_XZ);
      add("PoissonRatio23", POISSON_RATIO_YZ);
      add("ShearModulus", SHEAR_MODULUS);
      add("ShearModulus12", SHEAR_MODULUS_XY);
      add("ShearModulus13", SHEAR_MODULUS_XZ);
      add("ShearModulus23", SHEAR_MODULUS_YZ);
      // internal state variables
      add("Damage", DAMAGE);
      add("EquivalentPlasticStrain", EQUIVALENT_PLASTIC_STRAIN);
      add("PlasticStrain", PLASTIC_STRAIN);
    }  // end of KratosToTFELMapping

    /**
     * @brief add an new mapping
     * @tparam Type: type of the Kratos' variable
     * @param[in] k: TFEL glossary name
     * @param[in] k: Kratos' variable
     */
    template <typename Type>
    void add(const std::string& k, const Variable<Type>& v) {
      this->insert({k, v.Name()});
    }  // end of add

    KratosToTFELMapping(KratosToTFELMapping&&) = delete;
    KratosToTFELMapping(const KratosToTFELMapping&) = delete;
    KratosToTFELMapping& operator=(KratosToTFELMapping&&) = delete;
    KratosToTFELMapping& operator=(const KratosToTFELMapping&) = delete;
  };  // end of struct KratosToTFELMapping

  std::string getKratosVariableNameFromTFELGlossaryName(const std::string& g) {
    const auto& m = KratosToTFELMapping::get();
    auto p = m.find(g);
    if (p != m.end()) {
      return p->second;
    }
    return "";
  }  // end of getKratosVariableNameFromTFELGlossaryName

}  // end of namespace Kratos