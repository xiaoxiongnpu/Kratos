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
#include <iostream>

// Project includes
#include "includes/define.h"

// Application includes
#include "mgis_application_variables.h"

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

  /**
   * @brief this function maps TFEL glossary names to Kratos' variables names
   * @param[in] g: glossary name
   * @return Kratos' variable name if found, an empty string otherwise.
   */
  KRATOS_API(MGIS_APPLICATION)
  std::string getKratosVariableNameFromTFELGlossaryName(const std::string&);

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

}  // namespace Kratos.

#endif /* KRATOS_TO_TFEL_GLOSSARY_MAPPING_H_INCLUDED */
