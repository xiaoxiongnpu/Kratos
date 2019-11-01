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
#include "includes/imposed_deformation.h"
#include "constitutive/constitutive_law_with_imposed_deformation.h"

namespace Kratos
{
    namespace Testing
    {
        /**
        * Checks the correct work of the create methods
        */
        KRATOS_TEST_CASE_IN_SUITE(InitializeMaterialForConstitutiveLawWithImposedDeformation, KratosCoreFastSuite)
        {
            const auto& r_law = KratosComponents<ConstitutiveLaw>::Get("ConstitutiveLaw");

            KRATOS_CHECK_IS_FALSE(r_law.Is(ConstitutiveLaw::INTERNAL_IMPOSED_DEFORMATION));

            auto& r_law_with_imposed_deformation = const_cast<ConstitutiveLaw&>(KratosComponents<ConstitutiveLaw>::Get("ConstitutiveLawWithImposedDeformation"));

            KRATOS_CHECK(r_law_with_imposed_deformation.Is(ConstitutiveLaw::INTERNAL_IMPOSED_DEFORMATION));

            // Check initialize material works
            ImposedDeformation::Pointer p_imposed_deformation = Kratos::make_shared<ImposedDeformation>();

            Properties this_prop;
            this_prop.SetValue(IMPOSED_DEFORMATION, p_imposed_deformation);
            Geometry<Node<3>> geometry;
            Vector vector;

            r_law_with_imposed_deformation.InitializeMaterial(this_prop, geometry, vector);

            ConstitutiveLaw::Parameters dummy_parameters;
            KRATOS_CHECK_IS_FALSE(r_law_with_imposed_deformation.GetImposedDeformation(dummy_parameters) == NULL);
        }

    } // namespace Testing
}  // namespace Kratos.
