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
#include "custom_python/add_custom_utilities_to_python.h"

#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"

namespace Kratos
{
namespace Python
{
void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def_submodule("RansVariableUtilities")
        .def("ClipScalarVariable", &RansVariableUtilities::ClipScalarVariable)
        .def("GetMinimumScalarValue", &RansVariableUtilities::GetMinimumScalarValue)
        .def("GetMaximumScalarValue", &RansVariableUtilities::GetMaximumScalarValue)
        .def("CopyNodalSolutionStepVariablesList",
             &RansVariableUtilities::CopyNodalSolutionStepVariablesList)
        .def("InitializeDuplicatedModelPart", &RansVariableUtilities::InitializeDuplicatedModelPart)
        .def("CalculateMagnitudeSquareFor3DVariable", &RansVariableUtilities::CalculateMagnitudeSquareFor3DVariable)
        .def("FixFlaggedDofs", &RansVariableUtilities::FixFlaggedDofs, py::arg("model_part"), py::arg("variable"), py::arg("check_flag"), py::arg("check_value") = true)
        .def("CopyFlaggedVariableToNonHistorical", &RansVariableUtilities::CopyFlaggedVariableToNonHistorical<double>, py::arg("model_part"), py::arg("variable"), py::arg("check_flag"), py::arg("check_value") = true)
        .def("CopyFlaggedVariableToNonHistorical", &RansVariableUtilities::CopyFlaggedVariableToNonHistorical<array_1d<double, 3>>, py::arg("model_part"), py::arg("variable"), py::arg("check_flag"), py::arg("check_value") = true)
        .def("CopyFlaggedVariableFromNonHistorical", &RansVariableUtilities::CopyFlaggedVariableFromNonHistorical<double>, py::arg("model_part"), py::arg("variable"), py::arg("check_flag"), py::arg("check_value") = true)
        .def("CopyFlaggedVariableFromNonHistorical", &RansVariableUtilities::CopyFlaggedVariableFromNonHistorical<array_1d<double, 3>>, py::arg("model_part"), py::arg("variable"), py::arg("check_flag"), py::arg("check_value") = true)
        ;

    m.def_submodule("RansCalculationUtilities")
        .def("CalculateLogarithmicYPlusLimit",
             &RansCalculationUtilities::CalculateLogarithmicYPlusLimit,
             py::arg("kappa") = 0.41, py::arg("beta") = 5.2,
             py::arg("max_iterations") = 20, py::arg("tolerance") = 1e-6);
}

} // namespace Python.
} // Namespace Kratos
