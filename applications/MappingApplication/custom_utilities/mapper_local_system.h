//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_MAPPER_LOCAL_SYSTEM_H_INCLUDED )
#define  KRATOS_MAPPER_LOCAL_SYSTEM_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "mapper_interface_info.h"

namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// This is the "Condition" of the mappers
/** This class assembles the local system for the mappers, using the Information that
 * was provided by the "MapperInterfaceInfo"
*/
class MapperLocalSystem
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperLocalSystem
    KRATOS_CLASS_POINTER_DEFINITION(MapperLocalSystem);

    typedef Kratos::shared_ptr<MapperInterfaceInfo> MapperInterfaceInfoPointerType;
    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemUniquePointer;

    typedef typename MapperInterfaceInfo::CoordinatesArrayType CoordinatesArrayType;

    typedef Matrix MatrixType;
    typedef Vector VectorType;
    typedef std::vector<int> EquationIdVectorType; // int bcs of mpi

    typedef InterfaceObject::NodePointerType NodePointerType;
    typedef InterfaceObject::GeometryType GeometryType; 
    typedef InterfaceObject::GeometryPointerType GeometryPointerType;

    ///@}
    ///@name  Enum's
    ///@{

    enum class PairingStatus
    {
        NoInterfaceInfo,
        Approximation,
        InterfaceInfoFound
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~MapperLocalSystem() = default;

    ///@}
    ///@name Operations
    ///@{

    void EquationIdVectors(EquationIdVectorType& rOriginIds,
                           EquationIdVectorType& rDestinationIds)
    {
        if (!mIsComputed) {
            CalculateAll(mLocalMappingMatrix, mOriginIds, mDestinationIds, mPairingStatus);
            mIsComputed = true;
        }

        rOriginIds      = mOriginIds;
        rDestinationIds = mDestinationIds;
    }

    void CalculateLocalSystem(MatrixType& rLocalMappingMatrix,
                              EquationIdVectorType& rOriginIds,
                              EquationIdVectorType& rDestinationIds) const
    {
        if (mIsComputed) {
            std::cout << "EquationIdVectors was computed already" << std::endl;
            // This will be called if the EquationIdVectors have been querried before
            // i.e. matrix-based mapping
            rLocalMappingMatrix = mLocalMappingMatrix;
            rOriginIds      = mOriginIds;
            rDestinationIds = mDestinationIds;
        }
        else {            
            std::cout << "EquationIdVectors was NOT computed already" << std::endl;
            // This will be called if the EquationIdVectors have NOT been querried before
            // i.e. matrix-free mapping
            CalculateAll(rLocalMappingMatrix, rOriginIds, rDestinationIds, mPairingStatus);
        }
    }

    /**
    * @brief Resizing the output if no InterfaceInfo is available
    * This function resizes the system vectors to zero and also sets that no valid
    * Information from the other side could be found to compute the local system
    * @param rLocalMappingMatrix The vector conatining the mapping weights
    * @param rOriginIds The vector containing the ids on the origin
    * @param rDestinationIds The vector containing the ids on the destination
    * @param rPairingStatus The pairingstatus of the MapperLocalSystem
    * @see CalculateAll
    * @author Philipp Bucher
    */
    void ResizeToZero(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const
    {
        rPairingStatus = MapperLocalSystem::PairingStatus::NoInterfaceInfo;

        rLocalMappingMatrix.resize(0, 0, false);
        rOriginIds.resize(0);
        rDestinationIds.resize(0);
    }

    virtual CoordinatesArrayType& Coordinates() const = 0;


    void AddInterfaceInfo(MapperInterfaceInfoPointerType pInterfaceInfo) // TODO pass by const ref?
    {
        mInterfaceInfos.push_back(pInterfaceInfo);
    }

    bool HasInterfaceInfo() const
    {
        return mInterfaceInfos.size() > 0;
    }


    virtual void Clear()
    {
        mInterfaceInfos.clear();
        mLocalMappingMatrix.clear();
        mOriginIds.clear();
        mDestinationIds.clear();
        mIsComputed = false;
    }

    PairingStatus GetPairingStatus() const
    {
        return mPairingStatus;
    }

    ///@}
    ///@name Input and output
    ///@{

    virtual std::string PairingInfo(const int EchoLevel) const = 0;

    /// Turn back information as a string.
    virtual std::string Info() const {return "MapperLocalSystem";}

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

    virtual void CalculateRotationMatrixInterfaceInfos(MatrixType& _rotationMatrix_B,
                                                       VectorType& _translationVector_G_B,
                                                       VectorType& _linearShapeValues,
                                                       VectorType& _hermitianShapeValues,
                                                       VectorType& _hermitanDerShapeValues,
                                                       GeometryType& r_geom){ KRATOS_ERROR << "Base class function called! This function is used only for the Beam Mapper. " << std::endl;}

protected:
    ///@name Protected Life Cycle
    ///@{

    MapperLocalSystem() = default; // only accessbíble by derived classes

    ///@}
    ///@name Protected member Variables
    ///@{

    std::vector<MapperInterfaceInfoPointerType> mInterfaceInfos;

    bool mIsComputed = false;

    MatrixType mLocalMappingMatrix;
    EquationIdVectorType mOriginIds;
    EquationIdVectorType mDestinationIds;

    mutable PairingStatus mPairingStatus = PairingStatus::NoInterfaceInfo;

    ///@}
    ///@name Protected Operations
    ///@{

    // This function calculates the components necessary for the mapping
    // Note that it is "const", therefore it can NOT modify its members
    // Whether members are to be saved is decided in other functions of this class
    virtual void CalculateAll(MatrixType& rLocalMappingMatrix,
                              EquationIdVectorType& rOriginIds,
                              EquationIdVectorType& rDestinationIds,
                              MapperLocalSystem::PairingStatus& rPairingStatus) const = 0;

    ///@}

}; // Class MapperLocalSystem

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_LOCAL_SYSTEM_H_INCLUDED  defined
