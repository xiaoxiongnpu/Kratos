//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_parameters.h"
#include "constitutive/constitutive_law_with_imposed_deformation.h"

namespace Kratos
{
    namespace Testing
    {
        /**
        * Checks the correct work of the Has methods
        */
        KRATOS_TEST_CASE_IN_SUITE(CreateMethodForImposedDeformation, KratosCoreFastSuite)
        {
            const auto& r_law = KratosComponents<ConstitutiveLaw>::Get("ConstitutiveLaw");

            KRATOS_CHECK_IS_FALSE(r_law.Is(ConstitutiveLaw::INTERNAL_IMPOSED_DEFORMATION));

            const auto& r_law_with_imposed_deformation = KratosComponents<ConstitutiveLaw>::Get("ConstitutiveLawWithImposedDeformation");

            KRATOS_CHECK(r_law_with_imposed_deformation.Is(ConstitutiveLaw::INTERNAL_IMPOSED_DEFORMATION));

//             Parameters this_parameters = Parameters(R"(
//             {
//                 "name" : "ConstitutiveLawWithImposedDeformation"
//             })" );
//             auto p_law = r_law.Create(this_parameters);
//
//             KRATOS_CHECK(p_law->Is(ConstitutiveLaw::INTERNAL_IMPOSED_DEFORMATION));
        }

    } // namespace Testing
}  // namespace Kratos.
