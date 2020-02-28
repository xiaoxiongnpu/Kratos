//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "incompressible_pressure_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
IncompressiblePressureCondition<TDim, TNumNodes>& IncompressiblePressureCondition<TDim, TNumNodes>::operator=(
    IncompressiblePressureCondition<TDim, TNumNodes> const& rOther)
{
    Condition::operator=(rOther);

    return *this;
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer IncompressiblePressureCondition<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<IncompressiblePressureCondition>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer IncompressiblePressureCondition<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<IncompressiblePressureCondition>(NewId, pGeom, pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer IncompressiblePressureCondition<TDim, TNumNodes>::Clone(
    IndexType NewId, NodesArrayType const& rThisNodes) const
{
    Condition::Pointer pNewCondition = Create(
        NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template <unsigned int TDim, unsigned int TNumNodes>
int IncompressiblePressureCondition<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = BaseType::Check(rCurrentProcessInfo);

    const GeometryType& r_geometry = this->GetGeometry();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const NodeType& r_node = r_geometry[i_node];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE_POTENTIAL, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE, r_node);
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE_POTENTIAL, r_node);
    }

    return Check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::Initialize()
{
    KRATOS_TRY;

    this->CalculateNormal();

    KRATOS_ERROR_IF(!(this->GetValue(PARENT_CONDITION_POINTER)))
        << "Parent condition not found for duplicated " << this->Info() << "\n";

    const ConditionType& r_parent_condition = *(this->GetValue(PARENT_CONDITION_POINTER));

    KRATOS_ERROR_IF(r_parent_condition.GetValue(NEIGHBOUR_ELEMENTS).size() == 0)
        << r_parent_condition.Info() << " cannot find parent element\n";

    mNormalLength = RansCalculationUtilities::CalculateWallHeight(
        r_parent_condition, this->GetValue(NORMAL));

    KRATOS_CATCH("");
}

template <>
void IncompressiblePressureCondition<2, 2>::CalculateNormal()
{
    KRATOS_TRY

    Geometry<Node<3>>& pGeometry = this->GetGeometry();

    array_1d<double, 3> normal;

    normal[0] = pGeometry[1].Y() - pGeometry[0].Y();
    normal[1] = -(pGeometry[1].X() - pGeometry[0].X());
    normal[2] = 0.00;

    const double normal_magnitude = norm_2(normal);
    KRATOS_ERROR_IF(normal_magnitude == 0.0)
        << "Normal magnitude is zero in element " << this->Info();

    noalias(normal) = normal * (1 / normal_magnitude);
    this->SetValue(NORMAL, normal);

    KRATOS_CATCH("");
}

template <>
void IncompressiblePressureCondition<3, 3>::CalculateNormal()
{
    Geometry<Node<3>>& pGeometry = this->GetGeometry();

    array_1d<double, 3> v1, v2, normal;

    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(normal, v1, v2);

    const double normal_magnitude = norm_2(normal);
    KRATOS_ERROR_IF(normal_magnitude == 0.0)
        << "Normal magnitude is zero in element " << this->Info();

    noalias(normal) = normal * (1 / normal_magnitude);
    this->SetValue(NORMAL, normal);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (unsigned int i = 0; i < TNumNodes; ++i)
        rResult[i] = Condition::GetGeometry()[i].GetDof(PRESSURE_POTENTIAL).EquationId();
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::GetDofList(DofsVectorType& ConditionDofList,
                                                                  ProcessInfo& CurrentProcessInfo)
{
    if (ConditionDofList.size() != TNumNodes)
        ConditionDofList.resize(TNumNodes);

    for (unsigned int i = 0; i < TNumNodes; ++i)
        ConditionDofList[i] = Condition::GetGeometry()[i].pGetDof(PRESSURE_POTENTIAL);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::GetValuesVector(VectorType& rValues, int Step)
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] =
            rGeom[iNode].FastGetSolutionStepValue(PRESSURE_POTENTIAL, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod IncompressiblePressureCondition<TDim, TNumNodes>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // Calculate RHS
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Calculate LHS
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes, false);

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

    const GeometryType& r_geometry = this->GetGeometry();
    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    RansCalculationUtilities::CalculateConditionGeometryData(
        r_geometry, this->GetIntegrationMethod(), gauss_weights, shape_functions);
    const IndexType num_gauss_points = gauss_weights.size();

    const array_1d<double, 3>& r_normal = this->GetValue(NORMAL);

    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const Vector gauss_shape_functions = row(shape_functions, g);

        const array_1d<double, 3>& r_body_force = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, BODY_FORCE, gauss_shape_functions);
        const double density = RansCalculationUtilities::EvaluateInPoint(
            r_geometry, DENSITY, gauss_shape_functions);

        const double value =
            density * inner_prod(r_body_force, r_normal) * gauss_weights[g];

        for (IndexType a = 0; a < TNumNodes; ++a)
        {
            rRightHandSideVector[a] += value * gauss_shape_functions[a];
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string IncompressiblePressureCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "IncompressiblePressureCondition" << TDim << "D"
           << " ID: " << this->Id();
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "IncompressiblePressureCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void IncompressiblePressureCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// template instantiations

template class IncompressiblePressureCondition<2, 2>;
template class IncompressiblePressureCondition<3, 3>;

} // namespace Kratos.
