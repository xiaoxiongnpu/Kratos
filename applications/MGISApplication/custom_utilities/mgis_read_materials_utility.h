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
//  Main authors:  Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MGIS_APPLICATION_H_INCLUDED )
#define  KRATOS_MGIS_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "utilities/read_materials_utility.h"

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
 * @class MGISReadMaterialsUtility
 * @ingroup MGISApplication
 * @brief Process to read constitutive law and material properties from a json file
 * @details This process reads constitutive law and material properties from a json file
 * and assign them to elements and conditions. Adaptaed for MGIS variables
 * The definition includes the creation of subproperties
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(MGIS_APPLICATION) MGISReadMaterialsUtility 
    : public ReadMaterialsUtility
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef ReadMaterialsUtility BaseType;
    
    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of ReadMaterialProcess
    KRATOS_CLASS_POINTER_DEFINITION(MGISReadMaterialsUtility);

    ///@}
    ///@name Life Cycle
    ///@{
    
    /**
     * @brief Default constructor
     * @param rModel The model containing the problem to solve
     */
    MGISReadMaterialsUtility(Model& rModel) : BaseType(rModel) {};
    
    /**
     * @brief Constructor reading directly from file, via parameters
     * @param Params The configuration parameters telling where the configuration file can be found
     * @param rModel The model containing the problem to solve
     */
    MGISReadMaterialsUtility(
        Parameters Params,
        Model& rModel
        ) : BaseType(Params, rModel) 
    {
        
    };
    
    /**
     * @brief Constructor reading directly from file, via text
     * @param Params The string telling where the configuration file can be found
     * @param rModel The model containing the problem to solve
     */
    MGISReadMaterialsUtility(
        const std::string& rParametersName,
        Model& rModel
        ) : BaseType(rParametersName, rModel) 
    {
        
    };
    
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

    /// Turn back information as a string.
    std::string Info() const  {
        return "MGISReadMaterialsUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const  {
        rOStream << "MGISReadMaterialsUtility";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const  {
    }

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

    /**
     * @brief This method assigns the variables to a property from configuration parameters
     * @param MaterialData The parameters containing all the configurations of the materials
     * @param rProperty The reference to the property for which the materials are to be assigned
     */
    void AssignVariablesToProperty(
        const Parameters MaterialData,
        Properties& rProperty
        ) override;
        
    /**
     * @brief This method creates an auxiliar Parameters when reading properties in order to avoid error, so these non-registered properties can be processed later
     * @param VariablesParameters The original variable parameters
     * @param PropertyId The current property Id (for a warning)
     * @return The variables filtered if required
     */
    Parameters FilterVariables(
        const Parameters VariablesParameters,
        const IndexType PropertyId = 0
        ) override;

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

}; // Class MGISReadMaterialsUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_MGIS_READ_MATERIALS_H_INCLUDED defined
