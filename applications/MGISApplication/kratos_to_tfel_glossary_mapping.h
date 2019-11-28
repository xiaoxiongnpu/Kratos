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

#if !defined(KRATOS_TO_TFEL_GLOSSARY_MAPPING_H_INCLUDED)
#define KRATOS_TO_TFEL_GLOSSARY_MAPPING_H_INCLUDED

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

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
    enum struct StensorStorageConvention { STRAIN, STRESS };
    ///@}
    ///@name  Functions
    ///@{

    /**
     * @brief this function maps TFEL glossary names to Kratos' variables names
     * @param[in] g: glossary name
     * @return Kratos' variable name if found, an empty string otherwise.
     */
    KRATOS_API(MGIS_APPLICATION)
    std::string GetKratosVariableName(const std::string&);

    /**
     * @brief This function maps TFEL glossary names to Kratos' variables
     * @param[in] g: glossary name
     * @return Kratos' variable
     */
    template <typename Type>
    const Kratos::Variable<Type>& GetKratosVariable(const std::string&);

    /**
     * @brief Specialisation of the `getKratosVariable` for `double`
     * @param[in] g: glossary name
     * @return Kratos' variable
     */
    template<>
    KRATOS_API(MGIS_APPLICATION)
    const Kratos::Variable<double>& GetKratosVariable<double>(const std::string&);

    /**
     * @brief Specialisation of the `getKratosVariable` for `Vector`
     * @param[in] g: glossary name
     * @return Kratos' variable
     */
    template <>
    KRATOS_API(MGIS_APPLICATION)
    const Kratos::Variable<Kratos::Vector>& GetKratosVariable<Kratos::Vector>(const std::string&);

    /**
     * @brief Get the storage convention associated with a `Kratos`' variable
     * @param[in] Kratos' variable
     */
    KRATOS_API(MGIS_APPLICATION)
    StensorStorageConvention GetStensorStorageConvention(const Kratos::Variable<Kratos::Vector>&);

    ///@}
    ///@name Kratos Classes
    ///@{

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}

  }  // end of namespace MGIS

}  // namespace Kratos.

#endif /* KRATOS_TO_TFEL_GLOSSARY_MAPPING_H_INCLUDED */
