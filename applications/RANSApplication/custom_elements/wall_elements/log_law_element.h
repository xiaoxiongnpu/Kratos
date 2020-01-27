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

#if !defined(KRATOS_LOG_LAW_ELEMENT_H_INCLUDED)
#define KRATOS_LOG_LAW_ELEMENT_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/vms.h"
#include "custom_elements/wall_elements/rans_wall_element_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

template <unsigned int TDim, unsigned int TNumNodes = TDim + 1>
class LogLawElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LogLawElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LogLawElement);

    /// base type: an IndexedObject that automatically has a unique number
    typedef IndexedObject BaseType;

    /// definition of node type (default is: Node<3>)
    typedef Node<3> NodeType;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    typedef Properties PropertiesType;

    /// definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;

    /// definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef array_1d<double, TNumNodes> ShapeFunctionsType;
    typedef BoundedMatrix<double, TNumNodes, TDim> ShapeFunctionDerivativesType;

    using FluidElement = VMS<TDim, TNumNodes>;

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    LogLawElement(IndexType NewId = 0) : Element(NewId)
    {
    }

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    LogLawElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    LogLawElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    LogLawElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    /// Destructor.
    ~LogLawElement() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new LogLawElement element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LogLawElement<TDim, TNumNodes>>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LogLawElement<TDim, TNumNodes>>(NewId, pGeom, pProperties);
    }

    /// Provides local contributions from body forces and OSS projection terms
    /**
     * This is called during the assembly process and provides the terms of the
     * system that are either constant or computed explicitly (from the 'old'
     * iteration variables). In this case this means the body force terms and the
     * OSS projections, that are treated explicitly.
     * @param rLeftHandSideMatrix the elemental left hand side matrix. Not used here, required for compatibility purposes only.
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int LocalSize = (TDim + 1) * TNumNodes;

        // Check sizes and initialize matrix
        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

        // Calculate RHS
        this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    /// Returns a zero matrix of appropiate size (provided for compatibility with scheme)
    /**
     * @param rLeftHandSideMatrix Local matrix, will be filled with zeros
     * @param rCurrentProcessInfo Process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int LocalSize = (TDim + 1) * TNumNodes;

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    }

    /// Provides local contributions from body forces and projections to the RHS
    /**
     * This is called during the assembly process and provides the RHS terms of the
     * system that are either constant or computed explicitly (from the 'old'
     * iteration variables). In this case this means the body force terms and the
     * OSS projections, that are treated explicitly.
     * @param rRightHandSideVector Will be filled with the elemental right hand side
     * @param rCurrentProcessInfo ProcessInfo instance from the ModelPart. It is
     * expected to contain values for OSS_SWITCH, DYNAMIC_TAU and DELTA_TIME
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        FluidElement fluid_element(this->Id(), this->pGetGeometry());
        fluid_element.CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    /// Computes local contributions to the mass matrix
    /**
     * Provides the local contributions to the mass matrix, which is defined
     * here as the matrix associated to velocity derivatives. Note that the mass
     * matrix implemented here is lumped.
     * @param rMassMatrix Will be filled with the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
    }

    /// Computes the local contribution associated to 'new' velocity and pressure values
    /**
     * Provides local contributions to the system associated to the velocity and
     * pressure terms (convection, diffusion, pressure gradient/velocity divergence
     * and stabilization).
     * @param rDampingMatrix Will be filled with the velocity-proportional "damping" matrix
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo) override
    {
        const double kappa = 0.41;
        const double beta = 5.2;

        constexpr unsigned int LocalSize = (TDim + 1) * TNumNodes;
        constexpr unsigned int BlockSize = (TDim + 1);

        FluidElement fluid_element(this->Id(), this->pGetGeometry());
        fluid_element.CalculateLocalVelocityContribution(
            rDampingMatrix, rRightHandSideVector, rCurrentProcessInfo);

        // Now get the original right hand side
        BoundedVector<double, LocalSize> U(0.0);
        int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
        {
            array_1d<double, 3>& rVel =
                this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
            {
                U[LocalIndex] = rVel[d];
                ++LocalIndex;
            }
            U[LocalIndex] =
                this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
            ++LocalIndex;
        }
        const BoundedVector<double, LocalSize>& LHS_Vector = prod(rDampingMatrix, U);
        noalias(rRightHandSideVector) += LHS_Vector;

        Vector gauss_weights;
        Matrix shape_functions;
        GeometryData::ShapeFunctionsGradientsType shape_function_gradients;
        RansCalculationUtilities::CalculateGeometryData(
            this->GetGeometry(), fluid_element.GetIntegrationMethod(),
            gauss_weights, shape_functions, shape_function_gradients);

        const Vector& gauss_shape_functions = row(shape_functions, 0);
        const Matrix& gauss_shape_function_derivatives = shape_function_gradients[0];
        const double gauss_weight = gauss_weights[0];

        array_1d<double, 3> normal = RansCalculationUtilities::EvaluateInPoint(
            this->GetGeometry(), NORMAL, gauss_shape_functions);
        noalias(normal) = normal * (1.0 / norm_2(normal));

        BoundedMatrix<double, TDim, TDim> rotation_matrix;
        WallElementUtilities::CalculateRotationMatrix<TDim>(rotation_matrix, normal);

        // calculating velocity gradient
        const array_1d<double, 3> velocity = RansCalculationUtilities::EvaluateInPoint(
            this->GetGeometry(), VELOCITY, gauss_shape_functions);
        BoundedMatrix<double, TDim, TDim> velocity_gradient;
        RansCalculationUtilities::CalculateGradient<TDim>(
            velocity_gradient, this->GetGeometry(), VELOCITY, gauss_shape_function_derivatives);
        BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
        noalias(symmetric_velocity_gradient) =
            0.5 * (velocity_gradient + trans(velocity_gradient));
        const double kinematic_viscosity = RansCalculationUtilities::EvaluateInPoint(
            this->GetGeometry(), VISCOSITY, gauss_shape_functions);

        // calculating skin friction stress
        array_1d<double, 3> skin_friction = ZeroVector(3);
        for (unsigned int i = 0; i < TDim; ++i)
        {
            for (unsigned int j = 0; j < TDim; ++j)
            {
                skin_friction[i] += symmetric_velocity_gradient(i, j) * normal[j];
            }
        }
        const double temp = inner_prod(skin_friction, normal);
        noalias(skin_friction) = skin_friction - normal * temp;
        noalias(skin_friction) = skin_friction * kinematic_viscosity;
        const double u_tau = std::sqrt(norm_2(skin_friction));
        const double y =
            norm_2(WallElementUtilities::GetWallPosition(this->GetGeometry()) -
                   WallElementUtilities::GetGaussPosition(
                       this->GetGeometry(), gauss_shape_functions));

        const double y_plus = u_tau * y / kinematic_viscosity;
        const double y_plus_limit =
            RansCalculationUtilities::CalculateLogarithmicYPlusLimit(kappa, beta);

        const double u_star = (y_plus > y_plus_limit)
                                  ? u_tau * (std::log(y_plus) / kappa + beta)
                                  : u_tau * y_plus;

        const array_1d<double, 3>& rotated_velocity =
            WallElementUtilities::GetRotatedVector<TDim>(velocity, rotation_matrix);

        const double rotated_velocity_tangential_magnitude =
            std::sqrt(rotated_velocity[1] * rotated_velocity[1] +
                      rotated_velocity[2] * rotated_velocity[2]);

        BoundedVector<double, LocalSize> rotated_right_hand_side_vector =
            ZeroVector(LocalSize);
        for (unsigned int a = 0; a < TNumNodes; ++a)
        {
            // surface normal direction case
            double value = 0.0;
            for (unsigned int m = 0; m < TDim; ++m)
            {
                value += rRightHandSideVector[a * BlockSize + m] * rotation_matrix(0, m);
            }
            rotated_right_hand_side_vector[a * BlockSize] +=
                value + gauss_weight * gauss_shape_functions[0] * rotated_velocity[0];

            // surface tangential directions
            if (rotated_velocity_tangential_magnitude > 0.0)
            {
                for (unsigned int i = 1; i < TDim; ++i)
                {
                    // add the boundary layer velocity contribution
                    rotated_right_hand_side_vector[a * BlockSize + i] +=
                        gauss_weight * gauss_shape_functions[a] * u_star *
                        rotated_velocity[i] / rotated_velocity_tangential_magnitude;

                    for (unsigned int m = 0; m < TDim; ++m)
                    {
                        rotated_right_hand_side_vector[a * BlockSize + i] +=
                            rotation_matrix(i, m) * LHS_Vector[a * BlockSize + m];
                    }
                }
            }

            // Adding the additional term for damping matrix
            for (unsigned int b = 0; b < TNumNodes; ++b)
            {
                for (unsigned int i = 0; i < TDim; ++i)
                {
                    for (unsigned int j = 0; j < TDim; ++j)
                    {
                        rDampingMatrix(a * BlockSize + i, b * BlockSize + j) +=
                            gauss_weight * gauss_shape_functions[a] *
                            gauss_shape_functions[b] * rotation_matrix(i, j);
                    }
                }
            }
        }

        noalias(rRightHandSideVector) =
            rotated_right_hand_side_vector - prod(rDampingMatrix, U);
    }

    // The following methods have different implementations depending on TDim
    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& Values, int Step = 0) override;

    /// Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& Values, int Step = 0) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Elemental Data
    ///@{

    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "LogLawElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "LogLawElement" << TDim << "D";
    }

    //        /// Print object's data.
    //        virtual void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    LogLawElement& operator=(LogLawElement const& rOther);

    /// Copy constructor.
    LogLawElement(LogLawElement const& rOther);

    ///@}

}; // Class LogLawElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream, LogLawElement<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const LogLawElement<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_LOG_LAW_ELEMENT_H_INCLUDED  defined
