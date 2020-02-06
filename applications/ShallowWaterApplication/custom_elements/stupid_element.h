//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_STUPID_ELEMENT_H_INCLUDED
#define KRATOS_STUPID_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "shallow_water_application_variables.h"


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

class StupidElement : public Element
{
public:

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of StupidElement
    KRATOS_CLASS_POINTER_DEFINITION(StupidElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    StupidElement(IndexType NewId = 0)
    : Element(NewId)
    {}

    /**
     * Constructor using an array of nodes
     */
    StupidElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
    {}

    /**
     * Constructor using Geometry
     */
    StupidElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    /**
     * Constructor using Properties
     */
    StupidElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /**
     * Destructor
     */
    ~StupidElement(){};

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
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<StupidElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<StupidElement>(NewId, pGeom, pProperties);
    }

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        Element::Pointer p_new_elem = Kratos::make_intrusive<StupidElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
        p_new_elem->SetData(this->GetData());
        p_new_elem->Set(Flags(*this));
        return p_new_elem;
        KRATOS_CATCH("")
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override
    {
        const size_t element_size = TNumNodes*3;
        if(rResult.size() != element_size)
            rResult.resize(element_size, false); // False says not to preserve existing storage!!

        GeometryType& geom = GetGeometry();

        int counter=0;
        for (size_t i = 0; i < TNumNodes; i++)
        {
            rResult[counter++] = geom[i].GetDof(VELOCITY_X).EquationId();
            rResult[counter++] = geom[i].GetDof(VELOCITY_Y).EquationId();
            rResult[counter++] = geom[i].GetDof(HEIGHT).EquationId();
        }
    }

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override
    {
        const size_t element_size = TNumNodes*3;
        if(rElementalDofList.size() != element_size)
            rElementalDofList.resize(element_size);

        GeometryType& geom = GetGeometry();

        int counter=0;
        for (size_t i = 0; i < TNumNodes; i++)
        {
            rElementalDofList[counter++] = geom[i].pGetDof(VELOCITY_X);
            rElementalDofList[counter++] = geom[i].pGetDof(VELOCITY_Y);
            rElementalDofList[counter++] = geom[i].pGetDof(HEIGHT);
        }
    }

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override
    {}

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override
    {
        size_t element_size = TNumNodes*3;
    
        if(rRightHandSideVector.size() != element_size)
            rRightHandSideVector.resize(element_size, false);

        for (size_t i = 0; i < element_size; ++i)
        {
            rRightHandSideVector[i] = 1.0;
        }
    }

    void AddExplicitContribution(ProcessInfo& rCurrentProcessInfo) override
    {    
        VectorType local_vector;
        this->CalculateRightHandSide(local_vector, rCurrentProcessInfo);

        auto& r_geom = this->GetGeometry();
        size_t count = 0;

        // Adding the explicit contribution to the RHS variable
        for (size_t i = 0; i < TNumNodes; ++i)
        {
            r_geom[i].SetLock();
            r_geom[i].FastGetSolutionStepValue(RHS_VELOCITY_X) += local_vector[count++];
            r_geom[i].FastGetSolutionStepValue(RHS_VELOCITY_Y) += local_vector[count++];
            r_geom[i].FastGetSolutionStepValue(RHS_HEIGHT)     += local_vector[count++];
            r_geom[i].UnSetLock();
        }
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
    std::string Info() const override
    {
        return "StupidElement";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    static constexpr size_t TNumNodes = 3;

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
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
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
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class StupidElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_STUPID_ELEMENT_H_INCLUDED  defined
