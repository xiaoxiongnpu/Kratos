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

void MGISReadMaterialsUtility::AssingVariablesToProperty(
    const Parameters MaterialData,
    Properties& rProperty
    )
{
    KRATOS_TRY;
 
    // Add / override the values of material parameters in the p_properties
    if (MaterialData.Has("Variables")) {
        Parameters variables = MaterialData["Variables"];
        for (auto iter = variables.begin(); iter != variables.end(); ++iter) {
            const Parameters value = variables.GetValue(iter.name());

            std::string variable_name = iter.name();
            TrimComponentName(variable_name);

            // We don't just copy the values, we do some tyransformation depending of the destination variable
            if (KratosComponents<Variable<double> >::Has(variable_name)) {
                const Variable<double>& r_variable = KratosComponents<Variable<double>>().Get(variable_name);
                CheckIfOverwritingValue(rProperty, r_variable, value.GetDouble());
                rProperty.SetValue(r_variable, value.GetDouble());
            } else if(KratosComponents<Variable<bool> >::Has(variable_name)) {
                const Variable<bool>& r_variable = KratosComponents<Variable<bool>>().Get(variable_name);
                CheckIfOverwritingValue(rProperty, r_variable, value.GetBool());
                rProperty.SetValue(r_variable, value.GetBool());
            } else if(KratosComponents<Variable<int> >::Has(variable_name)) {
                const Variable<int>& r_variable = KratosComponents<Variable<int>>().Get(variable_name);
                CheckIfOverwritingValue(rProperty, r_variable, value.GetInt());
                rProperty.SetValue(r_variable, value.GetInt());
            } else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name)) {
                const Variable<array_1d<double, 3>>& r_variable = KratosComponents<Variable<array_1d<double, 3>>>().Get(variable_name);
                array_1d<double, 3> temp = ZeroVector(3);
                const Vector& r_value_variable = value.GetVector();
                KRATOS_ERROR_IF(r_value_variable.size() != 3) << "The vector of variable " << variable_name << " has size " << r_value_variable.size() << " and it is supposed to be 3" << std::endl;
                for (IndexType index = 0; index < 3; ++index)
                    temp[index] = r_value_variable[index];
                CheckIfOverwritingValue(rProperty, r_variable, temp);
                rProperty.SetValue(r_variable, temp);
            } else if(KratosComponents<Variable<array_1d<double, 6> > >::Has(variable_name)) {
                const Variable<array_1d<double, 6>>& r_variable = KratosComponents<Variable<array_1d<double, 6>>>().Get(variable_name);
                array_1d<double, 6> temp(6, 0.0);
                const Vector& r_value_variable = value.GetVector();
                KRATOS_ERROR_IF(r_value_variable.size() != 6) << "The vector of variable " << variable_name << " has size " << r_value_variable.size() << " and it is supposed to be 6" << std::endl;
                for (IndexType index = 0; index < 6; ++index)
                    temp[index] = r_value_variable[index];
                CheckIfOverwritingValue(rProperty, r_variable, temp);
                rProperty.SetValue(r_variable, temp);
            } else if(KratosComponents<Variable<Vector > >::Has(variable_name)) {
                const Variable<Vector>& r_variable = KratosComponents<Variable<Vector>>().Get(variable_name);
                CheckIfOverwritingValue(rProperty, r_variable, value.GetVector());
                rProperty.SetValue(r_variable, value.GetVector());
            } else if(KratosComponents<Variable<Matrix> >::Has(variable_name)) {
                const Variable<Matrix>& r_variable = KratosComponents<Variable<Matrix>>().Get(variable_name);
                CheckIfOverwritingValue(rProperty, r_variable, value.GetMatrix());
                rProperty.SetValue(r_variable, value.GetMatrix());
            } else if(KratosComponents<Variable<std::string> >::Has(variable_name)) {
                const Variable<std::string>& r_variable = KratosComponents<Variable<std::string>>().Get(variable_name);
                CheckIfOverwritingValue(rProperty, r_variable, value.GetString());
                rProperty.SetValue(r_variable, value.GetString());
            } else {
                // TODO: Here the MGIS varaibles are been read
            }
        }
    } else {
        KRATOS_INFO("Read materials") << "No variables defined for material ID: " << rProperty.Id() << std::endl;
    }
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

}  // namespace Kratos.
