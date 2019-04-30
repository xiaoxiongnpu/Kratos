//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//

#if !defined(KRATOS_GEOMETRICAL_OBJECT_H_INCLUDED )
#define  KRATOS_GEOMETRICAL_OBJECT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/flags.h"
#include "geometries/geometry.h"

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
 * @class GeometricalObject
 * @ingroup KratosCore
 * @brief This defines the geometrical object, base definition of the element and condition entities
 * @details Derives from IndexedObject, so it has an ID, and from Flags
 * @author Pooyan Dadvand
*/
class GeometricalObject : public IndexedObject, public Flags, public std::intrusive_base<GeometricalObject>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometricalObject
    typedef Kratos::intrusive_ptr<GeometricalObject> Pointer;
    typedef Kratos::intrusive_weak_ptr<GeometricalObject> WeakPointer;
    typedef Kratos::unique_ptr<GeometricalObject> UniquePointer;

    /// Definition of the node type
    typedef Node <3> NodeType;

    /// The geometry type definition
    typedef Geometry<NodeType> GeometryType;

    /// Defines the index type
    typedef std::size_t IndexType;

    /// Defines the result type
    typedef std::size_t result_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit GeometricalObject(IndexType NewId = 0)
        : IndexedObject(NewId),
          Flags(),
          mpGeometry()
    {}

    /// Default constructor.
    GeometricalObject(IndexType NewId, GeometryType::Pointer pGeometry)
        : IndexedObject(NewId),
          Flags(),
          mpGeometry(pGeometry)
    {}

    /// Destructor.
    ~GeometricalObject() override {}

    /// Copy constructor.
    GeometricalObject(GeometricalObject const& rOther)
        : IndexedObject(rOther.Id()),
          Flags(rOther),
          mpGeometry(rOther.mpGeometry)
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GeometricalObject& operator=(GeometricalObject const& rOther)
    {
        IndexedObject::operator=(rOther);
        Flags::operator =(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Returns the pointer to the geometry
     * @return The pointer of the geometry
     */
    GeometryType::Pointer pGetGeometry()
    {
        return mpGeometry;
    }

    /**
     * @brief Returns the pointer to the geometry (const version)
     * @return The pointer of the geometry
     */
    const GeometryType::Pointer pGetGeometry() const
    {
        return mpGeometry;
    }

    /**
     * @brief Returns the reference of the geometry
     * @return The reference of the geometry
     */
    GeometryType& GetGeometry()
    {
        return *mpGeometry;
    }

    /**
     * @brief Returns the reference of the geometry (const version)
     * @return The reference of the geometry
     */
    GeometryType const& GetGeometry() const
    {
        return *mpGeometry;
    }

    /**
     * @brief Returns the flags of the object
     * @return The  flags of the object
     */
    Flags& GetFlags()
    {
        return *this;
    }

    /**
     * @brief Returns the flags of the object (const version)
     * @return The  flags of the object
     */
    Flags const& GetFlags() const
    {
        return *this;
    }

    /**
     * @brief Sets the flags of the object
     * @param rThisFlags The flags to be set
     */
    void SetFlags(Flags const& rThisFlags)
    {
        Flags::operator=(rThisFlags);
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Checks if two GeometricalObject have the same type
     * @return True if the objects are the same type, false otherwise
     */
    inline static bool HasSameType(const GeometricalObject& rLHS, const GeometricalObject& rRHS) {
        return (typeid(rLHS) == typeid(rRHS));
    }

    /**
     * @brief Checks if two GeometricalObject have the same type (pointer version)
     * @return True if the objects are the same type, false otherwise
     */
    inline static bool HasSameType(const GeometricalObject * rLHS, const GeometricalObject* rRHS) {
        return GeometricalObject::HasSameType(*rLHS, *rRHS);
    }

    /**
     * @brief Checks if two GeometricalObject have the same geometry type
     * @return True if the geometries are the same type, false otherwise
     */
    inline static bool HasSameGeometryType(const GeometricalObject& rLHS, const GeometricalObject& rRHS) {
        return (rLHS.GetGeometry().GetGeometryType() == rRHS.GetGeometry().GetGeometryType());
    }

    /**
     * @brief Checks if two GeometricalObject have the same geometry type (pointer version)
     * @return True if the geometries are the same type, false otherwise
     */
    inline static bool HasSameGeometryType(const GeometricalObject* rLHS, const GeometricalObject* rRHS) {
        return GeometricalObject::HasSameGeometryType(*rLHS, *rRHS);
    }

    /**
     * @brief Checks if two GeometricalObject are the same
     * @return True if the object is the same, false otherwise
     */
    inline static bool IsSame(const GeometricalObject& rLHS, const GeometricalObject& rRHS) {
        return GeometricalObject::HasSameType(rLHS, rRHS) && GeometricalObject::HasSameGeometryType(rLHS, rRHS);
    }

    /**
     * @brief Checks if two GeometricalObject are the same (pointer version)
     * @return True if the object is the same, false otherwise
     */
    inline static bool IsSame(const GeometricalObject* rLHS, const GeometricalObject* rRHS) {
        return GeometricalObject::HasSameType(*rLHS, *rRHS) && GeometricalObject::HasSameGeometryType(*rLHS, *rRHS);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Geometrical object # "
               << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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

    GeometryType::Pointer mpGeometry; /// Pointer to the entity geometry

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        rSerializer.save("Geometry",mpGeometry);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
        rSerializer.load("Geometry",mpGeometry);
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

}; // Class GeometricalObject

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GeometricalObject& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GeometricalObject& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GEOMETRICAL_OBJECT_H_INCLUDED  defined


