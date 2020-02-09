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

// System includes

// External includes

// Project includes
#include "custom_utilities/mgis_read_materials_utility.h"

namespace Kratos {
namespace {

template <class TValueType>
void CheckIfOverwritingValue(const Properties& rProps,
                             const Variable<TValueType>& rVariable,
                             const TValueType& rValue)
{
    KRATOS_WARNING_IF("ReadMaterialsUtility", rProps.Has(rVariable)) << "The properties ID: "
        << rProps.Id() << " already has " << rVariable.Name() << ".\nOverwriting "
        << rProps[rVariable] << " with " << rValue << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template <class TValueType>
void CheckIfOverwritingTable(const Properties& rProps,
                             const Variable<TValueType>& rInputVariable,
                             const Variable<TValueType>& rOutputVariable)
{
    KRATOS_WARNING_IF("ReadMaterialsUtility", rProps.HasTable(rInputVariable, rOutputVariable))
        << "The properties ID: " << rProps.Id() << " already has a table for "
        << rInputVariable.Name() << " and " << rOutputVariable.Name()
        << ".\nIt is overwritten." << std::endl;
}

}

/***********************************************************************************/
/***********************************************************************************/

void MGISReadMaterialsUtility::AssignVariablesToProperty(
    const Parameters MaterialData,
    Properties& rProperty
    )
{
    KRATOS_TRY;
 
    ReadMaterialsUtility::AssignVariablesToProperty(MaterialData, rProperty);
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Parameters MGISReadMaterialsUtility::FilterVariables(
    const Parameters VariablesParameters,
    const IndexType PropertyId
    )
{
    KRATOS_TRY;

    Parameters variables_considered;
    for (auto iter = VariablesParameters.begin(); iter != VariablesParameters.end(); ++iter) {
        const Parameters value = VariablesParameters.GetValue(iter.name());
        std::string variable_name = iter.name();
        TrimComponentName(variable_name);
        // We create an axiliar parameters for adding the variables
        if (KratosComponents<Variable<double> >::Has(variable_name)) {
            variables_considered.AddValue(variable_name, value);
        } else if(KratosComponents<Variable<bool> >::Has(variable_name)) {
            variables_considered.AddValue(variable_name, value);
        } else if(KratosComponents<Variable<int> >::Has(variable_name)) {
            variables_considered.AddValue(variable_name, value);
        } else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name)) {
            variables_considered.AddValue(variable_name, value);
        } else if(KratosComponents<Variable<array_1d<double, 6> > >::Has(variable_name)) {
            variables_considered.AddValue(variable_name, value);
        } else if(KratosComponents<Variable<Vector > >::Has(variable_name)) {
            variables_considered.AddValue(variable_name, value);
        } else if(KratosComponents<Variable<Matrix> >::Has(variable_name)) {
            variables_considered.AddValue(variable_name, value);
        } else if(KratosComponents<Variable<std::string> >::Has(variable_name)) {
            variables_considered.AddValue(variable_name, value);
        } else {
            KRATOS_WARNING("Read materials") << "The variable property: " << variable_name << " for material ID: " << PropertyId << " is not a registered variable. It will be skipped, and assumed that this value will be processed later" << std::endl;
        }
    }
    return variables_considered;

    KRATOS_CATCH("");
}

}  // namespace Kratos.
