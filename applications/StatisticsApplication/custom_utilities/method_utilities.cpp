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
#include <cmath>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/ublas_interface.h"

// Application includes

// Include base h
#include "method_utilities.h"

namespace Kratos
{
namespace MethodsUtilities
{
template <typename TDataType>
TDataType RaiseToPower(const TDataType& rData, const double Power)
{
    return std::pow(rData, Power);
}

template <>
array_1d<double, 3> RaiseToPower(const array_1d<double, 3>& rData, const double Power)
{
    array_1d<double, 3> output;
    for (int i = 0; i < 3; ++i)
    {
        output[i] = RaiseToPower(rData[i], Power);
    }
    return output;
}

template <>
Vector RaiseToPower(const Vector& rData, const double Power)
{
    const int n = rData.size();
    Vector output(n);
    for (int i = 0; i < n; ++i)
    {
        output[i] = RaiseToPower(rData[i], Power);
    }
    return output;
}

template <>
Matrix RaiseToPower(const Matrix& rData, const double Power)
{
    const int n1 = rData.size1();
    const int n2 = rData.size2();
    Matrix output(n1, n2);
    for (int i = 0; i < n1; ++i)
    {
        for (int j = 0; j < n2; ++j)
        {
            output(i, j) = RaiseToPower(rData(i, j), Power);
        }
    }
    return output;
}

template <typename TDataType>
void DataTypeSizeInitializer(TDataType& rData, const TDataType& rReferenceData)
{
    // do nothing in the case of double and int or array_1d<double, 3>
}

template <>
void DataTypeSizeInitializer(Vector& rData, const Vector& rReferenceData)
{
    KRATOS_TRY

    const std::size_t size = rReferenceData.size();
    KRATOS_ERROR_IF(size == 0) << "Vector size is zero.\n";

    if (rData.size() != size)
        rData.resize(size, false);

    for (std::size_t i = 0; i < size; ++i)
    {
        rData[i] = 0.0;
    }

    KRATOS_CATCH("");
}

template <>
void DataTypeSizeInitializer(Matrix& rData, const Matrix& rReferenceData)
{
    KRATOS_TRY

    const std::size_t size1 = rReferenceData.size1();
    const std::size_t size2 = rReferenceData.size2();

    KRATOS_ERROR_IF(size1 == 0 || size2 == 0) << "Matrix size is zero.\n";

    if ((rData.size1() != size1) || (rData.size2() != size2))
        rData.resize(size1, size2, false);

    for (std::size_t i = 0; i < size1; ++i)
    {
        for (std::size_t j = 0; j < size2; ++j)
            rData(i, j) = 0.0;
    }

    KRATOS_CATCH("");
}

template <typename TDataType>
void DataTypeSizeChecker(const TDataType& rData, const TDataType& rReferenceData)
{
    // do nothing in the case of double and int or array_1d<double, 3>
}

template <>
void DataTypeSizeChecker(const Vector& rData, const Vector& rReferenceData)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rData.size() != rReferenceData.size())
        << "Data size and reference data vector size mismatch. [ "
        << rData.size() << " != " << rReferenceData.size() << " ].\n";

    KRATOS_CATCH("");
}

template <>
void DataTypeSizeChecker(const Matrix& rData, const Matrix& rReferenceData)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rData.size1() != rReferenceData.size1())
        << "Data size and reference data matrix size1 mismatch. [ "
        << rData.size1() << " != " << rReferenceData.size1() << " ].\n";
    KRATOS_ERROR_IF(rData.size2() != rReferenceData.size2())
        << "Data size and reference data matrix size2 mismatch. [ "
        << rData.size2() << " != " << rReferenceData.size2() << " ].\n";

    KRATOS_CATCH("");
}

template <>
NodesContainerType& GetDataContainer(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Nodes();
}

template <>
ElementsContainerType& GetDataContainer(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Elements();
}

template <>
ConditionsContainerType& GetDataContainer(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Conditions();
}

template <>
const NodesContainerType& GetDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Nodes();
}

template <>
const ElementsContainerType& GetDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Elements();
}

template <>
const ConditionsContainerType& GetDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Conditions();
}

template <>
const std::function<double(const double&)> GetNormMethod(const Variable<double>& rVariable,
                                                         const std::string& rNormType)
{
    KRATOS_TRY

    if (rNormType == "value")
    {
        return [](const double& rValue) -> double { return rValue; };
    }
    else if (rNormType == "magnitude")
    {
        return [](const double& rValue) -> double { return std::abs(rValue); };
    }
    else
    {
        KRATOS_ERROR << "Unknown norm type for double variable "
                     << rVariable.Name() << ". [ NormType = " << rNormType << " ]\n"
                     << "   Allowed norm types are:\n"
                     << "        magnitude\n"
                     << "        value\n";
    }

    return [](const double&) -> double { return 0.0; };

    KRATOS_CATCH("");
}

template <>
const std::function<double(const array_1d<double, 3>&)> GetNormMethod(
    const Variable<array_1d<double, 3>>& rVariable, const std::string& rNormType)
{
    KRATOS_TRY

    if (rNormType == "magnitude")
    {
        return [](const array_1d<double, 3>& rValue) -> double {
            return norm_2(rValue);
        };
    }
    else if (rNormType == "component_x")
    {
        return
            [](const array_1d<double, 3>& rValue) -> double { return rValue[0]; };
    }
    else if (rNormType == "component_y")
    {
        return
            [](const array_1d<double, 3>& rValue) -> double { return rValue[1]; };
    }
    else if (rNormType == "component_z")
    {
        return
            [](const array_1d<double, 3>& rValue) -> double { return rValue[2]; };
    }
    else
    {
        KRATOS_ERROR << "Unknown norm type for 3d variable " << rVariable.Name()
                     << ". [ NormType = " << rNormType << " ]\n"
                     << "   Allowed norm types are:\n"
                     << "        magnitude\n"
                     << "        component_x\n"
                     << "        component_y\n"
                     << "        component_z\n";
    }

    return [](const array_1d<double, 3>&) -> double { return 0.0; };

    KRATOS_CATCH("");
}

// method template instantiations

template double RaiseToPower(const double&, const double);
template int RaiseToPower(const int&, const double);

template void DataTypeSizeInitializer(double&, const double&);
template void DataTypeSizeInitializer(int&, const int&);
template void DataTypeSizeInitializer(array_1d<double, 3>&, const array_1d<double, 3>&);

template void DataTypeSizeChecker(const double&, const double&);
template void DataTypeSizeChecker(const int&, const int&);
template void DataTypeSizeChecker(const array_1d<double, 3>&, const array_1d<double, 3>&);

// class template instantiations
template class NonHistoricalDataValueRetrievalFunctor<NodeType>;
template class HistoricalDataValueRetrievalFunctor<NodeType>;
template class NonHistoricalDataValueRetrievalFunctor<ConditionType>;
template class NonHistoricalDataValueRetrievalFunctor<ElementType>;

//

} // namespace MethodsUtilities

} // namespace Kratos.