//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//


// System includes


// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

// Strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/custom_strategies/linear_mor_matrix_output_strategy.hpp"
#include "custom_strategies/custom_strategies/mor_offline_strategy.hpp"
#include "custom_strategies/custom_strategies/mor_online_strategy.hpp"

// Builders and solvers
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

// Linear solvers
#include "linear_solvers/linear_solver.h"


namespace Kratos {
namespace Python {

void  AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef LinearSolverType::Pointer LinearSolverPointer;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef BuilderAndSolverType::Pointer BuilderAndSolverPointer;

    // Custom strategy types
    typedef LinearMorMatrixOutputStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > LinearMorMatrixOutputStrategyType;
    typedef MorOnlineStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > MorOnlineStrategyType;
    typedef MorOfflineStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > MorOfflineStrategyType;

    // Custom builder and solver types
    typedef SystemMatrixBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > SystemMatrixBuilderAndSolverType;


    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    py::class_< LinearMorMatrixOutputStrategyType, typename LinearMorMatrixOutputStrategyType::Pointer, BaseSolvingStrategyType >(m,"LinearMorMatrixOutputStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer, bool >())
        ;

    py::class_< MorOnlineStrategyType, typename MorOnlineStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorOnlineStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer, BuilderAndSolverPointer, BaseSolvingStrategyType::Pointer, vector<double>, bool >())
        ;

    py::class_< MorOfflineStrategyType, typename MorOfflineStrategyType::Pointer, BaseSolvingStrategyType >(m,"MorOfflineStrategy")
        .def(py::init < ModelPart&, BaseSchemeType::Pointer, LinearSolverPointer, vector<double>, bool >())
        ;

    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

    py::class_< SystemMatrixBuilderAndSolverType, typename SystemMatrixBuilderAndSolverType::Pointer, BuilderAndSolverType >(m, "SystemMatrixBuilderAndSolver")
        .def(py::init < LinearSolverPointer >())
        ;

}

} // namespace Python.
} // Namespace Kratos
