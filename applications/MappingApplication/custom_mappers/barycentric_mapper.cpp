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

// System includes

// External includes

// Project includes
#include "barycentric_mapper.h"
#include "mapping_application_variables.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{

typedef Node<3> NodeType;
typedef Geometry<NodeType> GeometryType;

namespace { // internal functions only used in this cpp

array_1d<double, 3> CreateArrayFromVector(const std::vector<double>& rVector,
                                          const std::size_t StartPosition)
{
    array_1d<double,3> the_array;
    the_array[0] = rVector[StartPosition];
    the_array[1] = rVector[StartPosition+1];
    the_array[2] = rVector[StartPosition+2];
    return the_array;
}

void ComputeNeighborDistances(const array_1d<double,3>& rCoords,
                              const std::vector<double>& rNeighborCoords,
                              std::vector<double>& rDistances)
{
    for (std::size_t i=0; i<rDistances.size(); ++i) {
        const auto neighbor_coords = CreateArrayFromVector(rNeighborCoords, i*3);
        rDistances[i] = MapperUtilities::ComputeDistance(rCoords, neighbor_coords);
    }
}

void InsertIfCloser(const array_1d<double,3>& rRefCoords,
                    const int CandidateEquationId,
                    const array_1d<double,3>& rCandidateCoords,
                    std::vector<int>& rNeighborIds,
                    std::vector<double>& rNeighborCoods)
{
    const std::size_t num_interpolation_nodes = rNeighborIds.size();

    std::vector<double> neighbor_distances(num_interpolation_nodes, std::numeric_limits<double>::max());

    // compute the distances to the currently closest nodes to the candidate to check for coinciding nodes
    ComputeNeighborDistances(rCandidateCoords, rNeighborCoods, neighbor_distances);
    for (const double dist : neighbor_distances) {
        if (dist < 1e-12) return; // the candidate is coinciding with an already used point, hence we don't use it
    }

    // compute the distances to the currently closest nodes based on the coordinates
    ComputeNeighborDistances(rRefCoords, rNeighborCoods, neighbor_distances);

    // compute distance to candidate
    const double candidate_distance = MapperUtilities::ComputeDistance(rRefCoords, rCandidateCoords);

    // check if the candidate is closer than the previously found nodes and save it if it is closer
    for (std::size_t i=0; i<num_interpolation_nodes; ++i) {
        if (candidate_distance < neighbor_distances[i]) {
            rNeighborIds.insert(rNeighborIds.begin()+i, CandidateEquationId);
            rNeighborCoods.insert(rNeighborCoods.begin()+(i*3), std::begin(rCandidateCoords), std::end(rCandidateCoords));
            break;
        }
    }
    // resize is required because insert increases the size
    if (rNeighborIds.size() != num_interpolation_nodes) {
        rNeighborIds.resize(num_interpolation_nodes);
    }
    if (rNeighborCoods.size() != 3*num_interpolation_nodes) {
        rNeighborCoods.resize(3*num_interpolation_nodes);
    }
}

bool BarycentricInterpolateInEntity(const array_1d<double,3>& rRefCoords,
                                    const std::vector<double>& rCoordinates,
                                    Vector& rShapeFunctionValues,
                                    std::vector<int>& rEquationIds,
                                    ProjectionUtilities::PairingIndex& rPairingIndex)
{
    // Check how many "proper" results were found
    const std::size_t num_interpolation_nodes = rEquationIds.size() - std::count(rEquationIds.begin(), rEquationIds.end(), -1);

    const bool is_full_projection = num_interpolation_nodes == rEquationIds.size();

    KRATOS_DEBUG_ERROR_IF(num_interpolation_nodes < 2 || num_interpolation_nodes > 4) << "Wrong number of interpolation nodes" << std::endl;
    KRATOS_DEBUG_ERROR_IF(rCoordinates.size() < num_interpolation_nodes*3) << "Not enough coords" << std::endl;
    KRATOS_DEBUG_ERROR_IF(rCoordinates.size()%3 != 0) << "Coords have wrong size" << std::endl;

    GeometryType::PointsArrayType geom_points;
    for (std::size_t i=0; i<num_interpolation_nodes; ++i) {
        geom_points.push_back(Kratos::make_intrusive<NodeType>(0, rCoordinates[i*3], rCoordinates[i*3+1], rCoordinates[i*3+2]));
        geom_points[i].SetValue(INTERFACE_EQUATION_ID, rEquationIds[i]);
    }

    Kratos::unique_ptr<GeometryType> p_geom;
    if      (num_interpolation_nodes == 2) p_geom = Kratos::make_unique<Line2D2<NodeType>>(geom_points);
    else if (num_interpolation_nodes == 3) p_geom = Kratos::make_unique<Triangle3D3<NodeType>>(geom_points);
    else if (num_interpolation_nodes == 4) p_geom = Kratos::make_unique<Tetrahedra3D4<NodeType>>(geom_points);

    double dummy_dist;

    return is_full_projection && ProjectionUtilities::ComputeProjection(*p_geom, Point(rRefCoords), 0.25, rShapeFunctionValues, rEquationIds, dummy_dist, rPairingIndex, true);
}

}

void BarycentricInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject,
                                                   const double NeighborDistance)
{
    SetLocalSearchWasSuccessful();

    const auto p_node = rInterfaceObject.pGetBaseNode();
    InsertIfCloser(
        Coordinates(),
        p_node->GetValue(INTERFACE_EQUATION_ID),
        p_node->Coordinates(),
        mNodeIds,
        mNeighborCoordinates);
}

void BarycentricLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    if (mInterfaceInfos.size() < 1) {
        ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
        return;
    }

    std::vector<int> node_ids;
    std::vector<double> neighbor_coods;

    // allocate final vectors, using the max possible size (in case of a volume interpolation)
    std::vector<int> final_node_ids(4);
    std::vector<double> final_neighbor_coords(12);
    std::fill(final_node_ids.begin(), final_node_ids.end(), -1);
    std::fill(final_neighbor_coords.begin(), final_neighbor_coords.end(), std::numeric_limits<double>::max());

    for (std::size_t i=0; i<mInterfaceInfos.size(); ++i) {
        mInterfaceInfos[i]->GetValue(node_ids, MapperInterfaceInfo::InfoType::Dummy);
        mInterfaceInfos[i]->GetValue(neighbor_coods, MapperInterfaceInfo::InfoType::Dummy);

        const std::size_t num_interpolation_nodes = node_ids.size();

        if (final_node_ids.size() != num_interpolation_nodes) {
            final_node_ids.resize(num_interpolation_nodes);
        }

        for (std::size_t j=0; j<num_interpolation_nodes; ++j) {
            InsertIfCloser(
                Coordinates(),
                node_ids[j],
                CreateArrayFromVector(neighbor_coods, j*3),
                final_node_ids,
                final_neighbor_coords);
        }
    }

    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    if (rDestinationIds.size() != 1) rDestinationIds.resize(1);
    rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);

    const std::size_t num_interpolation_nodes = final_node_ids.size();

    KRATOS_DEBUG_ERROR_IF(num_interpolation_nodes < 2 || num_interpolation_nodes > 4) << "Wrong number of interpolation nodes" << std::endl;

    KRATOS_ERROR_IF(final_node_ids[0] == -1) << "Not even an approximation is found, this should not happen!" << std::endl; // TODO should this be an error?

    if (final_node_ids[1] == -1) { // only one node was found => using nearest-neighbor
        rPairingStatus = MapperLocalSystem::PairingStatus::Approximation;
        mPairingIndex = ProjectionUtilities::PairingIndex::Closest_Point;
        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != 1) rLocalMappingMatrix.resize(1, 1, false);
        rLocalMappingMatrix(0,0) = 1.0;
        if (rOriginIds.size() != 1) rOriginIds.resize(1);
        rOriginIds[0] = final_node_ids[0];
    } else {
        Vector shape_function_values;
        const bool is_full_projection = BarycentricInterpolateInEntity(Coordinates(), final_neighbor_coords, shape_function_values, final_node_ids, mPairingIndex);
        rPairingStatus = is_full_projection ? MapperLocalSystem::PairingStatus::InterfaceInfoFound : MapperLocalSystem::PairingStatus::Approximation;

        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != shape_function_values.size()) {
            rLocalMappingMatrix.resize(1, shape_function_values.size(), false);
        }
        for (std::size_t i=0; i<shape_function_values.size(); ++i) {
            rLocalMappingMatrix(0,i) = shape_function_values[i];
        }

        if (rOriginIds.size() != final_node_ids.size()) rOriginIds.resize(final_node_ids.size());
        for (std::size_t i=0; i<final_node_ids.size(); ++i) {
            rOriginIds[i] = final_node_ids[i];
        }
    }
}

std::string BarycentricLocalSystem::PairingInfo(const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    std::stringstream buffer;
    buffer << "NearestElementLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) {// TODO leave here?
        buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
        if (mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) {
            mpNode->SetValue(PAIRING_STATUS, (int)mPairingIndex);
        }
    }
    return buffer.str();
}

}  // namespace Kratos.
