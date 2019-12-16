//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes


// External includes


// Project includes
#include "includes/variables.h"
#include "includes/checks.h"
#include "mass_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"


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
MassElement::MassElement(IndexType NewId)
    : Element(NewId)
{
}

/**
 * Constructor using an array of nodes
 */
MassElement::MassElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
{
}

/**
 * Constructor using Geometry
 */
MassElement::MassElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
MassElement::MassElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

/**
 * Copy Constructor
 */
MassElement::MassElement(MassElement const& rOther)
    : Element(rOther)
{
}

/**
 * Destructor
 */
MassElement::~MassElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
MassElement & MassElement::operator=(MassElement const& rOther)
{
    BaseType::operator=(rOther);
    Flags::operator =(rOther);
    // mpProperties = rOther.mpProperties;
    return *this;
}

///@}
///@name Operations
///@{

Element::Pointer MassElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MassElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

Element::Pointer MassElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MassElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

Element::Pointer MassElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MassElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

void MassElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType local_size = GetGeometry().PointsNumber() * 3;

    if (rResult.size() != local_size) {
        rResult.resize(local_size, false);
    }

    SizeType local_index = 0;
    const SizeType d_pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

    for (const auto& r_node : GetGeometry()) {
        rResult[local_index++] = r_node.GetDof(DISPLACEMENT_X, d_pos).EquationId();
        rResult[local_index++] = r_node.GetDof(DISPLACEMENT_Y, d_pos + 1).EquationId();
        rResult[local_index++] = r_node.GetDof(DISPLACEMENT_Z, d_pos + 2).EquationId();
    }

    KRATOS_CATCH("")
}

void MassElement::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType local_size = GetGeometry().PointsNumber() * 3;

    if (rElementalDofList.size() != local_size) {
        rElementalDofList.resize(local_size);
    }

    SizeType local_index = 0;

    for (const auto& r_node : GetGeometry()){
        rElementalDofList[local_index++] = r_node.pGetDof(DISPLACEMENT_X);
        rElementalDofList[local_index++] = r_node.pGetDof(DISPLACEMENT_Y);
        rElementalDofList[local_index++] = r_node.pGetDof(DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
}


void MassElement::GenericGetValuesVector(Vector& rValues, int Step, const ArrayVariableType& rVariable)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType local_size = number_of_nodes * 3;

    if (rValues.size() != local_size) {
        rValues.resize(local_size, false);
    }

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        int index = i * 3;
        const auto& r_vals = GetGeometry()[i].FastGetSolutionStepValue(rVariable, Step);

        rValues[index] =     r_vals[0];
        rValues[index + 1] = r_vals[1];
        rValues[index + 2] = r_vals[2];
    }

    KRATOS_CATCH("")
}

void MassElement::GetValuesVector(Vector& rValues, int Step)
{
    GenericGetValuesVector(rValues, Step, DISPLACEMENT);
}

void MassElement::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    GenericGetValuesVector(rValues, Step, VELOCITY);
}

void MassElement::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    GenericGetValuesVector(rValues, Step, ACCELERATION);
}

void MassElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
}

void MassElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const SizeType local_size = GetGeometry().PointsNumber() * 3;
    if (rLeftHandSideMatrix.size1() != local_size) {
        rLeftHandSideMatrix.resize(local_size, local_size, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(local_size, local_size);
}

void MassElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
}

void MassElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const auto& r_geom = GetGeometry();

    // lumped mass matrix
    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType local_size = number_of_nodes * 3;

    if (rMassMatrix.size1() != local_size) {
        rMassMatrix.resize(local_size, local_size, false);
    }

    noalias(rMassMatrix) = ZeroMatrix(local_size, local_size);

    const SizeType local_dim = GetGeometry().LocalSpaceDimension();

    double total_mass;

    if (local_dim == 1) { // line
        total_mass = GetProperties()[CROSS_AREA] * StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);

    } else if (local_dim == 2) { // tri / quad
        total_mass = GetProperties()[THICKNESS]* r_geom.Area(); // TODO Area is calculated o the current configuration, fix this!

    } else {
        KRATOS_ERROR << "Wrong local space dimension!" << std::endl;
    }

    Vector lump_fact = ZeroVector(number_of_nodes);
    r_geom.LumpingFactors(lump_fact);

    total_mass *= StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);

    for (SizeType i=0; i<number_of_nodes; ++i) {
        const double temp = lump_fact[i] * total_mass;

        for (SizeType j=0; j<3; ++j) {
            const SizeType index = i * 3 + j;
            rMassMatrix(index, index) = temp;
        }
    }



    KRATOS_CATCH("")
}

void MassElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        GetGeometry().PointsNumber() * 3);
}

int MassElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1) << "Element found with Id 0 or negative" << std::endl;

    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const SizeType local_dim = GetGeometry().WorkingSpaceDimension();

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : GetGeometry()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node);
    }

    // verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(CROSS_AREA);
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS);
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)

    if (GetProperties().Has(DENSITY) == false) {
        KRATOS_ERROR << "DENSITY not provided for this element #" << Id() << std::endl;
    }

    // dimension specific checks
    if (dim == 1) { // line
        KRATOS_ERROR_IF(number_of_nodes != 2) << "Wrong number of nodes!" << std::endl;

        KRATOS_ERROR_IF(StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this) <= 0) << "On element #" << Id() << "; Length cannot be less than or equal to 0" << std::endl;

        KRATOS_ERROR_IF(GetProperties().Has(CROSS_AREA) == false || GetProperties()[CROSS_AREA] <= numerical_limit) << "CROSS_AREA not provided for this element #" << Id() << std::endl;
    } else if (local_dim == 2) { // tri / quad
        KRATOS_ERROR_IF(number_of_nodes != 3 && number_of_nodes != 4) << "Wrong number of nodes!" << std::endl;

        KRATOS_ERROR_IF(GetGeometry().Area() <= 0) << "On element #" << Id() << "; Area cannot be less than or equal to 0" << std::endl;

        KRATOS_ERROR_IF(GetProperties().Has(THICKNESS) == false || GetProperties()[THICKNESS] <= numerical_limit) << "THICKNESS not provided for this element #" << Id() << std::endl;

    } else {
        KRATOS_ERROR << "Wrong local space dimension, can only be 2 (line) or 3 (surface)!" << std::endl;
    }

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

std::string MassElement::Info() const {
    std::stringstream buffer;
    buffer << "MassElement #" << Id();
    return buffer.str();
}

/// Print information about this object.

void MassElement::PrintInfo(std::ostream& rOStream) const {
    rOStream << "MassElement #" << Id();
}

/// Print object's data.

void MassElement::PrintData(std::ostream& rOStream) const {
    pGetGeometry()->PrintData(rOStream);
}

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

///@}
///@name Serialization
///@{

void MassElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

void MassElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );

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
inline std::istream & operator >> (std::istream& rIStream, MassElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const MassElement& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.
