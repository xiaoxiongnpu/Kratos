//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Michael Breitenberger
//                   Riccardo Rossi
//

// System includes

// External includes
#include "custom_conditions/penalty_coupling_condition.h"

// Project includes

namespace Kratos
{

    void PenaltyCouplingCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        const auto& r_geometry_master = GetGeometry();
        const auto& r_geometry_slave = r_geometry_master.GetSlaveGeometry();

        const int number_of_nodes_master = r_geometry_master.size();
        const int number_of_nodes_slave = r_geometry_slave.size();

        const int mat_size = 3 * (number_of_nodes_master + number_of_nodes_slave);

        const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry_master.IntegrationPoints(integration_method);

        const double Penalty = GetProperties()[PENALTY_FACTOR];

        for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
        {
            Matrix N_master = r_geometry_master.ShapeFunctionsValues(integration_method);
            Matrix N_slave = r_geometry_slave.ShapeFunctionsValues(integration_method);

            //FOR DISPLACEMENTS
            Matrix Hcomplete = ZeroMatrix(3, mat_size);
            for (unsigned int i = 0; i < number_of_nodes_master; i++)
            {
                int index = 3 * i;
                if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                    H(0, index) = N_master[i];
                if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                    H(1, index + 1) = N_master[i];
                if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                    H(2, index + 2) = N_master[i];
            }

            for (unsigned int i = 0; i < number_of_nodes_slave; i++)
            {
                int index = 3 * i + number_of_nodes_master;
                if (Is(IgaFlags::FIX_DISPLACEMENT_X))
                    H(0, index) = N_slave[i];
                if (Is(IgaFlags::FIX_DISPLACEMENT_Y))
                    H(1, index + 1) = N_slave[i];
                if (Is(IgaFlags::FIX_DISPLACEMENT_Z))
                    H(2, index + 2) = N_slave[i];
            }

            Vector TDisplacements(mat_size);
            for (unsigned int i = 0; i < number_of_nodes_master; i++)
            {
                const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
                int index = 3 * i;
                TDisplacements[index] = disp[0];
                TDisplacements[index + 1] = disp[1];
                TDisplacements[index + 2] = disp[2];
            }
            for (unsigned int i = 0; i < number_of_nodes_slave; i++)
            {
                const array_1d<double, 3> disp = GetGeometry().GetSlaveGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
                int index = 3 * i + number_of_nodes_master;
                TDisplacements[index] = disp[0];
                TDisplacements[index + 1] = disp[1];
                TDisplacements[index + 2] = disp[2];
            }

            const double integration_weight = integration_points[point_number].Weight();
            const double determinat_jacobian = GetGeometry().DeterminatJacobian(point_number);

            noalias(rLeftHandSideMatrix) += prod(trans(Hcomplete), Hcomplete)
                * integration_weight * determinat_jacobian * Penalty;
            noalias(rRightHandSideVector) -= prod(prod(trans(Hcomplete), Hcomplete), TDisplacements)
                * integration_weight * determinat_jacobian * Penalty;
        }

        KRATOS_CATCH("")
    }

    void PenaltyCouplingCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int number_of_nodes_master = GetGeometry().size();
        const int number_of_nodes_slave = GetGeometry().GetSlaveGeometry().size();

        if (rResult.size() != 3 * (number_of_nodes_master + number_of_nodes_slave))
            rResult.resize(3 * (number_of_nodes_master + number_of_nodes_slave), false);

        for (unsigned int i = 0; i < number_of_nodes_master; ++i) {
            const unsigned int index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }

        for (unsigned int i = 0; i < number_of_nodes_slave; ++i) {
            const unsigned int index = i * 3 + number_of_nodes_master;
            rResult[index]     = GetGeometry().GetSlaveGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry().GetSlaveGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry().GetSlaveGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }

        KRATOS_CATCH("")
    }

    void PenaltyCouplingCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int number_of_nodes_master = GetGeometry().size();
        const int number_of_nodes_slave = GetGeometry().GetSlaveGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * (number_of_nodes_master + number_of_nodes_slave));

        for (unsigned int i = 0; i < number_of_nodes_master; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        for (unsigned int i = 0; i < number_of_nodes_slave; ++i) {
            rElementalDofList.push_back(GetGeometry().GetSlaveGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry().GetSlaveGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry().GetSlaveGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos


