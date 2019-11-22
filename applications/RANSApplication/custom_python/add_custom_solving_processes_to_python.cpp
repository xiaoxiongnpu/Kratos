// System includes

#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif                                    // KRATOS_USE_AMATRIX

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

#include "custom_python/add_custom_solving_processes_to_python.h"

// Application includes
#include "custom_processes/solving_strategies/evm_co_solving_process.h"

namespace Kratos
{
namespace Python
{
void AddCustomSolvingProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;

    // Adding solving strategies
    using EvmCoSolvingProcessType =
        EvmCoSolvingProcess<SparseSpaceType, LocalSpaceType, LinearSolverType>;
    py::class_<EvmCoSolvingProcessType, EvmCoSolvingProcessType::Pointer, Process>(
        m, "EvmCoSolvingProcess")
        .def(py::init<ModelPart&, Parameters>())
        .def("AddStrategy", &EvmCoSolvingProcessType::AddStrategy)
        .def("AddAuxiliaryProcess", &EvmCoSolvingProcessType::AddAuxiliaryProcess)
        .def("SetParentSolvingStrategy", &EvmCoSolvingProcessType::SetParentSolvingStrategy);
}

} // namespace Python.
} // Namespace Kratos
