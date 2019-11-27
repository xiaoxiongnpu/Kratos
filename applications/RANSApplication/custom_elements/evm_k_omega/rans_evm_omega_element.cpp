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
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/properties.h"

#include "custom_elements/evm_k_omega/evm_k_omega_utilities.h"

#include "includes/cfd_variables.h"
#include "rans_application_variables.h"
// Include Base h
#include "rans_evm_omega_element.h" //renamed

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * Constructor.
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmOmegaElement<TDim, TNumNodes>::RansEvmOmegaElement(IndexType NewId)
    : StabilizedConvectionDiffusionReaction<TDim, TNumNodes, RansEvmOmegaElementData>(NewId)
{
}

/**
 * Constructor using an array of nodes
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmOmegaElement<TDim, TNumNodes>::RansEvmOmegaElement(IndexType NewId,
                                                              const NodesArrayType& ThisNodes)
    : StabilizedConvectionDiffusionReaction<TDim, TNumNodes, RansEvmOmegaElementData>(
          NewId, ThisNodes)
{
}

/**
 * Constructor using Geometry
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmOmegaElement<TDim, TNumNodes>::RansEvmOmegaElement(IndexType NewId,
                                                              GeometryType::Pointer pGeometry)
    : StabilizedConvectionDiffusionReaction<TDim, TNumNodes, RansEvmOmegaElementData>(
          NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmOmegaElement<TDim, TNumNodes>::RansEvmOmegaElement(IndexType NewId,
                                                              GeometryType::Pointer pGeometry,
                                                              PropertiesType::Pointer pProperties)
    : StabilizedConvectionDiffusionReaction<TDim, TNumNodes, RansEvmOmegaElementData>(
          NewId, pGeometry, pProperties)
{
}

/**
 * Copy Constructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmOmegaElement<TDim, TNumNodes>::RansEvmOmegaElement(
	RansEvmOmegaElement<TDim, TNumNodes> const& rOther)
    : StabilizedConvectionDiffusionReaction<TDim, TNumNodes, RansEvmOmegaElementData>(rOther)
{
}

/**
 * Destructor
 */
template <unsigned int TDim, unsigned int TNumNodes>
RansEvmOmegaElement<TDim, TNumNodes>::~RansEvmOmegaElement()
{
}

///@}
///@name Operators
///@{

///@}
///@name Operations
///@{

/**
 * ELEMENTS inherited from this class have to implement next
 * Create and Clone methods: MANDATORY
 */

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer RansEvmOmegaElement<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmOmegaElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer RansEvmOmegaElement<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmOmegaElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer RansEvmOmegaElement<TDim, TNumNodes>::Clone(IndexType NewId,
                                                               NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<RansEvmOmegaElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmOmegaElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                              ProcessInfo& CurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (unsigned int i = 0; i < TNumNodes; i++)
        rResult[i] =
            Element::GetGeometry()[i].GetDof(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE).EquationId(); // TURBULENT_ENERGY_DISSIPATION_RATE --> omega (turbulent specific dissipation rate)
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmOmegaElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                        ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != TNumNodes)
        rElementalDofList.resize(TNumNodes);

    for (unsigned int i = 0; i < TNumNodes; i++) // node indexes of the elements
        rElementalDofList[i] =
            Element::GetGeometry()[i].pGetDof(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE); // get the omega DoFs from the elements DoF list
	//(this also includes u, p, k, omega) but we only need the position of the omega DoFs to build up the omega DoF list of the element [w1, w2, w3] --> in the nodes
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmOmegaElement<TDim, TNumNodes>::GetValuesVector(VectorType& rValues, int Step)
{
    this->GetFirstDerivativesVector(rValues, Step);
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmOmegaElement<TDim, TNumNodes>::GetFirstDerivativesVector(VectorType& rValues, int Step) // it is called in the function above
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] = rGeom[iNode].FastGetSolutionStepValue(
			TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, Step); //omega
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmOmegaElement<TDim, TNumNodes>::GetSecondDerivativesVector(VectorType& rValues, int Step) // this is called by a scheme, not the element
{
    if (rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    GeometryType& rGeom = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rValues[LocalIndex++] = rGeom[iNode].FastGetSolutionStepValue(
			TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, Step); // TURBULENT_ENERGY_DISSIPATION_RATE_2 is epsilon dot (derivative of eps.) --> same for omega
    } //register : TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod RansEvmOmegaElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
	// definition of how many Gauss points we use -->GI_GAUSS_2 --> 2 points
	return GeometryData::GI_GAUSS_2;
}

/**
 * This method provides the place to perform checks on the completeness of the
 * input and the compatibility with the problem options as well as the
 * contitutive laws selected It is designed to be called only once (or anyway,
 * not often) typically at the beginning of the calculations, so to verify that
 * nothing is missing from the input or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
template <unsigned int TDim, unsigned int TNumNodes>
int RansEvmOmegaElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Check(rCurrentProcessInfo);

    for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
    {
        NodeType& r_node = this->GetGeometry()[iNode];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_2, r_node); //RANS_AUXILIARY_VARIABLE_2 is only for transient simulations

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node); // we don't need information about "k" here --> only in k element
    } // TURBULENT_ENERGY_DISSIPATION_RATE should be changed for omega

    return 0;

    KRATOS_CATCH("");
}

///@}
///@name Access
///@{

///@}
///@name Inquiry
///@{

///@}
///@name Input and output
///@{

/// Turn back information as a string.

// from here we don't need to change anything (only debugging)

template <unsigned int TDim, unsigned int TNumNodes>
std::string RansEvmOmegaElement<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "RansEvmOmegaElement #" << Element::Id();
    return buffer.str();
}

/// Print information about this object.

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmOmegaElement<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "RansEvmOmegaElement #" << Element::Id();
}

/// Print object's data.

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmOmegaElement<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
    Element::pGetGeometry()->PrintData(rOStream);
}

// until here

///@}
///@name Friends
///@{

///@}

///@name Protected static Member Variables
///@{

///@}
///@name Protected member Variables
///@{

///@}
///@name Protected Operators
///@{

///@}
///@name Protected Operations
///@{

///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{

///@}
///@name Protected LifeCycle
///@{

///@}

///@name Static Member Variables
///@{

///@}
///@name Member Variables
///@{

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{
template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& RansEvmOmegaElement<TDim, TNumNodes>::GetPrimalVariable() const
{
	return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& RansEvmOmegaElement<TDim, TNumNodes>::GetPrimalRelaxedRateVariable() const
{
	return RANS_AUXILIARY_VARIABLE_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmOmegaElement<TDim, TNumNodes>::CalculateElementData(
    RansEvmOmegaElementData& rData,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo,
    const int Step) const
{
	//specific for omega element
	rData.SigmaOmega = rCurrentProcessInfo[TURBULENCE_RANS_SIGMA_OMEGA];
	rData.Gamma = rCurrentProcessInfo[TURBULENCE_RANS_GAMMA];
	rData.BetaZero = rCurrentProcessInfo[TURBULENCE_RANS_BETA_ZERO];
	rData.BetaZeroStar = rCurrentProcessInfo[TURBULENCE_RANS_BETA_ZERO_STAR];

	//both k and omega elements
    const double nu = this->EvaluateInPoint(KINEMATIC_VISCOSITY, rShapeFunctions);
    const double nu_t = this->EvaluateInPoint(TURBULENT_VISCOSITY, rShapeFunctions); //figure out where this is calculated
    const double tke = this->EvaluateInPoint(TURBULENT_KINETIC_ENERGY, rShapeFunctions);
	const double omega = this->EvaluateInPoint(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, rShapeFunctions); //implement TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE !!
    const double theta = EvmKomegaModelUtilities::CalculateTheta(nu_t); //-->theta

	BoundedMatrix<double, TDim, TDim> velocity_gradient_matrix; //--> ask Suneth how this works
	this->CalculateGradient(velocity_gradient_matrix, VELOCITY, rShapeFunctionDerivatives);
	const double beta = EvmKomegaModelUtilities::CalculateBeta<TDim>(velocity_gradient_matrix,
																omega,
																rData.BetaZero,
																rData.BetaZeroStar);
	rData.Beta = beta;

    rData.ShapeFunctionDerivatives = rShapeFunctionDerivatives;
    rData.TurbulentKinematicViscosity = nu_t;
    rData.KinematicViscosity = nu;
    rData.TurbulentKineticEnergy = tke;
	rData.TurbulentSpecificDissipationRate = omega;
	rData.Theta = theta;
    rData.VelocityDivergence =
        this->GetDivergenceOperator(VELOCITY, rShapeFunctionDerivatives);
}

template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmOmegaElement<TDim, TNumNodes>::CalculateEffectiveKinematicViscosity(
									const RansEvmOmegaElementData& rData,
									const Vector& rShapeFunctions,
									const Matrix& rShapeFunctionDerivatives,
									const ProcessInfo& rCurrentProcessInfo,
									const int Step) const
{
    return rData.KinematicViscosity + rData.SigmaOmega * rData.TurbulentKinematicViscosity;
}



//this is done
template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmOmegaElement<TDim, TNumNodes>::CalculateReactionTerm(
    const RansEvmOmegaElementData& rData, const Vector& rShapeFunctions,
	const Matrix& rShapeFunctionDerivatives,
	const ProcessInfo& rCurrentProcessInfo,
	const int Step) const
{
    return std::max(rData.Beta * rData.TurbulentSpecificDissipationRate + (2.0 / 3.0) *rData.Gamma* rData.VelocityDivergence, 0.0);
}

//this is done
template <unsigned int TDim, unsigned int TNumNodes>
double RansEvmOmegaElement<TDim, TNumNodes>::CalculateSourceTerm(
    const RansEvmOmegaElementData& rData, const Vector& rShapeFunctions,
	const Matrix& rShapeFunctionDerivatives,
	const ProcessInfo& rCurrentProcessInfo,
	const int Step) const
{
    double production = 0.0;

    BoundedMatrix<double, TDim, TDim> velocity_gradient_matrix;
    this->CalculateGradient(velocity_gradient_matrix, VELOCITY, rShapeFunctionDerivatives);

    production = EvmKomegaModelUtilities::CalculateSourceTerm<TDim>(
        velocity_gradient_matrix, rData.TurbulentKinematicViscosity, rData.TurbulentKineticEnergy);

    production *= (rData.Gamma * rData.Theta);
    return production;
}

///@}
///@name Serialization
///@{

// from here on only change names --> RansEvmEpsilonElement

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmOmegaElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}

template <unsigned int TDim, unsigned int TNumNodes>
void RansEvmOmegaElement<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

    // List
    // To be completed with the class member list
}

///@}
///@name Private  Access
///@{

///@}
///@name Private Inquiry
///@{

///@}
///@name Un accessible methods
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function

template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                RansEvmOmegaElement<TDim, TNumNodes>& rThis);

/// output stream function

template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEvmOmegaElement<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

// Class template instantiation

template class RansEvmOmegaElement<2, 3>;
template class RansEvmOmegaElement<3, 4>;
template class RansEvmOmegaElement<2, 4>;
template class RansEvmOmegaElement<3, 8>;



} // namespace Kratos.
