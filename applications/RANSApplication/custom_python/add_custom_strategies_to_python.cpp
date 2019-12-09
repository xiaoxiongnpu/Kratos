//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/define_python.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

// strategies
#include "custom_strategies/generic_residual_based_bossak_velocity_scalar_scheme.h"
#include "custom_strategies/generic_residualbased_simple_steady_scalar_scheme.h"
#include "custom_strategies/rans_residualbased_newton_raphson_strategy.h"

// convergence criterians
#include "custom_strategies/generic_convergence_criteria.h"

namespace Kratos
{
namespace Python
{
void AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using BaseSchemeType = Scheme<SparseSpaceType, LocalSpaceType>;
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    using ConvergenceCriteriaType = ConvergenceCriteria<SparseSpaceType, LocalSpaceType>;

    // Convergence criteria
    py::class_<GenericConvergenceCriteria<SparseSpaceType, LocalSpaceType>,
               typename GenericConvergenceCriteria<SparseSpaceType, LocalSpaceType>::Pointer, ConvergenceCriteriaType>(
        m, "GenericScalarConvergenceCriteria")
        .def(py::init<double, double>());

    py::class_<GenericResidualBasedBossakVelocityScalarScheme<SparseSpaceType, LocalSpaceType>,
               typename GenericResidualBasedBossakVelocityScalarScheme<SparseSpaceType, LocalSpaceType>::Pointer, BaseSchemeType>(
        m, "GenericResidualBasedBossakVelocityDynamicScalarScheme")
        .def(py::init<const double, const double, const Variable<double>&,
                      const Variable<double>&, const Variable<double>&>());

    py::class_<GenericResidualBasedSimpleSteadyScalarScheme<SparseSpaceType, LocalSpaceType>,
               typename GenericResidualBasedSimpleSteadyScalarScheme<SparseSpaceType, LocalSpaceType>::Pointer, BaseSchemeType>(
        m, "GenericResidualBasedSimpleSteadyScalarScheme")
        .def(py::init<const double>());

    using BaseSolvingStrategyType =
        SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using BuilderAndSolverType =
        BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    using RansResidualBasedNewtonRaphsonStrategyType =
        RansResidualBasedNewtonRaphsonStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>;

    py::class_<RansResidualBasedNewtonRaphsonStrategyType, typename RansResidualBasedNewtonRaphsonStrategyType::Pointer, BaseSolvingStrategyType>(
        m, "RansResidualBasedNewtonRaphsonStrategy")
        .def(py::init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer,
                      ConvergenceCriteriaType::Pointer, int, bool, bool, bool>())
        .def(py::init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
                      BuilderAndSolverType::Pointer, int, bool, bool, bool>())
        .def(py::init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer,
                      ConvergenceCriteriaType::Pointer, Parameters>())
        .def(py::init<ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer,
                      ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, Parameters>())
        .def("SetMaxIterationNumber", &RansResidualBasedNewtonRaphsonStrategyType::SetMaxIterationNumber)
        .def("GetMaxIterationNumber", &RansResidualBasedNewtonRaphsonStrategyType::GetMaxIterationNumber)
        .def("SetKeepSystemConstantDuringIterations",
             &RansResidualBasedNewtonRaphsonStrategyType::SetKeepSystemConstantDuringIterations)
        .def("GetKeepSystemConstantDuringIterations",
             &RansResidualBasedNewtonRaphsonStrategyType::GetKeepSystemConstantDuringIterations)
        .def("SetInitializePerformedFlag",
             &RansResidualBasedNewtonRaphsonStrategyType::SetInitializePerformedFlag)
        .def("GetInitializePerformedFlag",
             &RansResidualBasedNewtonRaphsonStrategyType::GetInitializePerformedFlag)
        .def("GetSystemMatrix", &RansResidualBasedNewtonRaphsonStrategyType::GetSystemMatrix,
             py::return_value_policy::reference_internal)
        .def("GetSystemVector", &RansResidualBasedNewtonRaphsonStrategyType::GetSystemVector,
             py::return_value_policy::reference_internal)
        .def("GetSolutionVector", &RansResidualBasedNewtonRaphsonStrategyType::GetSolutionVector,
             py::return_value_policy::reference_internal);
}

} // namespace Python.
} // Namespace Kratos
