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

#include "Epetra_FEVector.h"
#include "trilinos_space.h"

#include "add_trilinos_solving_processes_to_python.h"

// Application includes
#include "custom_processes/solving_strategies/evm_co_solving_process.h"

namespace Kratos
{
namespace Python
{
void AddTrilinosSolvingProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using MPISparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using MPILinearSolverType = LinearSolver<MPISparseSpaceType, LocalSpaceType>;

    using MPIEvmCoSolvingProcessType =
        EvmCoSolvingProcess<MPISparseSpaceType, LocalSpaceType, MPILinearSolverType>;
    py::class_<MPIEvmCoSolvingProcessType, MPIEvmCoSolvingProcessType::Pointer, Process>(
        m, "MPIEvmCoSolvingProcess")
        .def(py::init<ModelPart&, Parameters>())
        .def("AddStrategy", &MPIEvmCoSolvingProcessType::AddStrategy)
        .def("AddAuxiliaryProcess", &MPIEvmCoSolvingProcessType::AddAuxiliaryProcess)
        .def("SetParentSolvingStrategy", &MPIEvmCoSolvingProcessType::SetParentSolvingStrategy);
}

} // namespace Python.
} // Namespace Kratos
