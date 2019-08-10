//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "factories/standard_builder_and_solver_factory.h"
#include "spaces/ublas_space.h"

// Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_constraints_elementwise.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h"

namespace Kratos
{
    void RegisterBuilderAndSolvers()
    {
        typedef TUblasSparseSpace<double> SpaceType;
        typedef TUblasDenseSpace<double> LocalSpaceType;
        typedef LinearSolver<SpaceType, LocalSpaceType> LinearSolverType;

        typedef BuilderAndSolver<SpaceType,  LocalSpaceType, LinearSolverType> BuilderAndSolverType;
        typedef ResidualBasedEliminationBuilderAndSolver<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedEliminationBuilderAndSolverType;
        typedef ResidualBasedEliminationBuilderAndSolverWithConstraints<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedEliminationBuilderAndSolverWithConstraintsType;
        typedef ResidualBasedBlockBuilderAndSolver<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverType;
        typedef ResidualBasedBlockBuilderAndSolverWithConstraintsElementWise<SpaceType,  LocalSpaceType, LinearSolverType> ResidualBasedBlockBuilderAndSolverWithConstraintsElementWiseType;

        //NOTE: here we must create persisting objects for the builder and solvers
        static auto ResidualBasedEliminationBuilderAndSolverFactory = StandardBuilderAndSolverFactory<BuilderAndSolverType, LinearSolverType, ResidualBasedEliminationBuilderAndSolverType>();
        static auto ResidualBasedEliminationBuilderAndSolverWithConstraintsFactory = StandardBuilderAndSolverFactory<BuilderAndSolverType, LinearSolverType, ResidualBasedEliminationBuilderAndSolverWithConstraintsType>();
        static auto ResidualBasedBlockBuilderAndSolverFactory = StandardBuilderAndSolverFactory<BuilderAndSolverType, LinearSolverType, ResidualBasedBlockBuilderAndSolverType>();
        static auto ResidualBasedBlockBuilderAndSolverWithConstraintsElementWiseFactory = StandardBuilderAndSolverFactory<BuilderAndSolverType, LinearSolverType, ResidualBasedBlockBuilderAndSolverWithConstraintsElementWiseType>();

        // Registration of convergence solvers
        KRATOS_REGISTER_BUILDER_AND_SOLVER("elimination_builder_and_solver", ResidualBasedEliminationBuilderAndSolverFactory);
        KRATOS_REGISTER_BUILDER_AND_SOLVER("elimination_builder_and_solver_with_constraints", ResidualBasedEliminationBuilderAndSolverWithConstraintsFactory);
        KRATOS_REGISTER_BUILDER_AND_SOLVER("block_builder_and_solver", ResidualBasedBlockBuilderAndSolverFactory);
        KRATOS_REGISTER_BUILDER_AND_SOLVER("block_builder_and_solver_with_constraints_element_wise", ResidualBasedBlockBuilderAndSolverWithConstraintsElementWiseFactory);
    };
} // Namespace Kratos
