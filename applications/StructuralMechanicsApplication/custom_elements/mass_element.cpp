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

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer MassElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MassElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer MassElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MassElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer MassElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<MassElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
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

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
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

/**
 * this is called during the assembling process in order
 * to calculate all elemental contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void MassElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void MassElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void MassElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental mass matrix
 * @param rMassMatrix: the elemental mass matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void MassElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const auto& r_geom = GetGeometry();

    // lumped mass matrix
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = number_of_nodes * 3;

    if (rMassMatrix.size1() != mat_size) {
        rMassMatrix.resize(mat_size, mat_size, false);
    }

    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);

    const SizeType local_dim = GetGeometry().LocalSpaceDimension();

    if (local_dim == 2) { // line
        // Clear matrix
        if (rMassVector.size() != msLocalSize) {
            rMassVector.resize(msLocalSize, false);
        }

        const double A = GetProperties()[CROSS_AREA];
        const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        const double rho = GetProperties()[DENSITY];

        const double total_mass = A * L * rho;

        for (int i = 0; i < msNumberOfNodes; ++i) {
            for (int j = 0; j < msDimension; ++j) {
                int index = i * msDimension + j;

                rMassVector[index] = total_mass * 0.50;
            }
        }

        // Compute lumped mass matrix
        VectorType temp_vector(msLocalSize);
        CalculateLumpedMassVector(temp_vector);

        // Clear matrix
        if (rMassMatrix.size1() != msLocalSize || rMassMatrix.size2() != msLocalSize) {
            rMassMatrix.resize(msLocalSize, msLocalSize, false);
        }
        rMassMatrix = ZeroMatrix(msLocalSize, msLocalSize);

        // Fill the matrix
        for (IndexType i = 0; i < msLocalSize; ++i) {
            rMassMatrix(i, i) = temp_vector[i];
        }

    } else if (local_dim == 3 || local_dim == 4) { // tri / quad
        const double total_mass = r_geom.Area() * GetProperties()[THICKNESS] * StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);

        Vector lump_fact = ZeroVector(number_of_nodes);
        r_geom.LumpingFactors(lump_fact);

        for (SizeType i=0; i<number_of_nodes; ++i) {
            const double temp = lump_fact[i] * total_mass;

            for (SizeType j=0; j<3; ++j) {
                const SizeType index = i * 3 + j;
                rMassMatrix(index, index) = temp;
            }
        }
    } else {
        KRATOS_ERROR << "Wrong local space dimension!" << std::endl;
    }



    KRATOS_CATCH("")
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental damping matrix
 * @param rDampingMatrix: the elemental damping matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void MassElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rDampingMatrix.size1() != 0)
        rDampingMatrix.resize(0, 0, false);
}

/**
 * This method provides the place to perform checks on the completeness of the input
 * and the compatibility with the problem options as well as the contitutive laws selected
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
int MassElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1) <<"MassElement found with Id 0 or negative" << std::endl;

    // KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0) << "On MassElement #" << this->Id() <<  "; Area cannot be less than or equal to 0" << std::endl;

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
