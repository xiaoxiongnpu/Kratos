//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
#include "log_law_element.h"

namespace Kratos
{
///@name Specialized implementation of LogLawElement for functions that depend on TDim
///@{

/**
 * @see LogLawElement::EquationIdVector
 */
template <>
void LogLawElement<2>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(3), LocalSize(9);
    unsigned int LocalIndex = 0;

    unsigned int vpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(VELOCITY_X, vpos).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(VELOCITY_Y, vpos + 1).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(PRESSURE, ppos).EquationId();
    }
}

/**
 * @see LogLawElement::EquationIdVector
 */
template <>
void LogLawElement<3>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(4), LocalSize(16);
    unsigned int LocalIndex = 0;
    unsigned int vpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(VELOCITY_X, vpos).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(VELOCITY_Y, vpos + 1).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(VELOCITY_Z, vpos + 2).EquationId();
        rResult[LocalIndex++] =
            this->GetGeometry()[iNode].GetDof(PRESSURE, ppos).EquationId();
    }
}

/**
 * @see LogLawElement::GetDofList
 */
template <>
void LogLawElement<2>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(3), LocalSize(9);
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/**
 * @see LogLawElement::GetDofList
 */
template <>
void LogLawElement<3>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(4), LocalSize(16);
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/**
 * @see LogLawElement::GetFirstDerivativesVector
 */
template <>
void LogLawElement<2>::GetFirstDerivativesVector(Vector& Values, int Step)
{
    const unsigned int NumNodes(3), LocalSize(9);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        array_1d<double, 3>& rVelocity =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
        Values[LocalIndex++] = rVelocity[0];
        Values[LocalIndex++] = rVelocity[1];
        Values[LocalIndex++] =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
    }
}

/**
 * @see LogLawElement::GetFirstDerivativesVector
 */
template <>
void LogLawElement<3>::GetFirstDerivativesVector(Vector& Values, int Step)
{
    const unsigned int NumNodes(4), LocalSize(16);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        array_1d<double, 3>& rVelocity =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
        Values[LocalIndex++] = rVelocity[0];
        Values[LocalIndex++] = rVelocity[1];
        Values[LocalIndex++] = rVelocity[2];
        Values[LocalIndex++] =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
    }
}

/**
 * @see LogLawElement::GetSecondDerivativesVector
 */
template <>
void LogLawElement<2>::GetSecondDerivativesVector(Vector& Values, int Step)
{
    const unsigned int NumNodes(3), LocalSize(9);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        array_1d<double, 3>& rAcceleration =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
        Values[LocalIndex++] = rAcceleration[0];
        Values[LocalIndex++] = rAcceleration[1];
        Values[LocalIndex++] = 0.0; // Pressure Dof
    }
}

/**
 * @see LogLawElement::GetSecondDerivativesVector
 */
template <>
void LogLawElement<3>::GetSecondDerivativesVector(Vector& Values, int Step)
{
    const unsigned int NumNodes(4), LocalSize(16);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        array_1d<double, 3>& rAcceleration =
            this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
        Values[LocalIndex++] = rAcceleration[0];
        Values[LocalIndex++] = rAcceleration[1];
        Values[LocalIndex++] = rAcceleration[2];
        Values[LocalIndex++] = 0.0; // Pressure Dof
    }
}

///@} // Specialized implementations
template class LogLawElement<2>;
template class LogLawElement<3>;

} // namespace Kratos
