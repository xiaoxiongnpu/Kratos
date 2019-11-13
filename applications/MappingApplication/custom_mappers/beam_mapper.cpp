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
#include "utilities/math_utils.h"

// Project includes
#include "beam_mapper.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapping_matrix_utilities.h"
#include "mapping_application_variables.h"
#include <math.h>
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_searching/interface_communicator_mpi.h"
#endif

namespace Kratos
{

// ************ BeamMapperInterfaceInfo function definitions

void BeamMapperInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject,
                                                      const double NeighborDistance)
{
    SaveSearchResult(rInterfaceObject,false);
}

void BeamMapperInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject,
                                                                  const double NeighborDistance)
{
    KRATOS_ERROR << "ProcessSearchResultForApproximation is not implemented for Beam Mapping yet." << std::endl;
    //SaveSearchResult(rInterfaceObject, true);
}

void BeamMapperInterfaceInfo::SaveSearchResult(const InterfaceObject& rInterfaceObject,
                                               const bool ComputeApproximation)
{
    const auto p_geom = rInterfaceObject.pGetBaseGeometry();
    std::cout << "my first result " << std::endl;
    rInterfaceObject.PrintInfo(std::cout);
    std::cout << "\n";
    double proj_dist;

    const Point point_to_proj(this->Coordinates());
    Point projection_point;

    mCoordinates = point_to_proj;

    Vector linear_shape_function_values;
    Vector hermitian_shape_function_values;
    Vector hermitian_shape_function_derivatives_values;
    std::vector<int> eq_ids;

    ProjectionUtilities::PairingIndex pairing_index;

    const auto geom_family = p_geom->GetGeometryFamily();
    KRATOS_ERROR_IF(geom_family != GeometryData::Kratos_Linear) << "Invalid geometry of the Origin! The geometry should be a beam!" << std::endl;

    //const bool is_full_projection = ProjectionUtilities::ComputeProjection(*p_geom, point_to_proj, mLocalCoordTol, shape_function_values, eq_ids, proj_dist, pairing_index, ComputeApproximation);

    // Calculating and storing the shape function values
    // Linear shape functions
    pairing_index = ProjectionUtilities::ProjectOnLine(*p_geom, point_to_proj, mLocalCoordTol, linear_shape_function_values, eq_ids, proj_dist, ComputeApproximation);
    // Hermitian shape functions
    ProjectionUtilities::ProjectOnLineHermitian(*p_geom, point_to_proj, mLocalCoordTol, hermitian_shape_function_values, hermitian_shape_function_derivatives_values, proj_dist, projection_point);
    const bool is_full_projection = (pairing_index == ProjectionUtilities::PairingIndex::Line_Inside);
    
    std::cout << "rShapeFunctionValues : " << linear_shape_function_values << std::endl;
    std::cout << "rHermitianShapeFunctionValues : " << hermitian_shape_function_values << std::endl;
    std::cout << "rHermitianShapeFunctionValuesDerivatives : " << hermitian_shape_function_derivatives_values << std::endl;
    

    if (is_full_projection) {
        SetLocalSearchWasSuccessful();
    } else {
        if (!ComputeApproximation) {
            return;
        } else {
            SetIsApproximation();
        }
    }

    const std::size_t num_values_linear = linear_shape_function_values.size();
    const std::size_t num_values_hermitian = hermitian_shape_function_values.size();
    const std::size_t num_values_hermitian_der = hermitian_shape_function_derivatives_values.size();

    //KRATOS_ERROR_IF_NOT(num_values == eq_ids.size()) << "Number of equation-ids is not the same as the number of ShapeFunction values, something went wrong!" << std::endl;

    if (pairing_index > mPairingIndex || (pairing_index == mPairingIndex && proj_dist < mClosestProjectionDistance)) {
        mPairingIndex = pairing_index;
        mClosestProjectionDistance = proj_dist;
        mNodeIds = eq_ids;

        if (mLinearShapeFunctionValues.size() != num_values_linear) mLinearShapeFunctionValues.resize(num_values_linear);
        for (std::size_t i=0; i<num_values_linear; ++i) {
            mLinearShapeFunctionValues[i] = linear_shape_function_values[i];
        }

        if (mHermitianShapeFunctionValues.size() != num_values_linear) mHermitianShapeFunctionValues.resize(num_values_hermitian);
        for (std::size_t i=0; i<num_values_hermitian; ++i) {
            mHermitianShapeFunctionValues[i] = hermitian_shape_function_values[i];
        }

        if (mHermitianShapeFunctionValuesDerivatives.size() != num_values_hermitian_der) mHermitianShapeFunctionValuesDerivatives.resize(num_values_hermitian_der);
        for (std::size_t i=0; i<num_values_hermitian_der; ++i) {
            mHermitianShapeFunctionValuesDerivatives[i] = hermitian_shape_function_derivatives_values[i];
        }

        mProjectionOfPoint = projection_point;
        //std::cout << "the projected point defined in the GCS is : " << mProjectionOfPoint << std::endl; 
        
        //mpInterfaceObject = &rInterfaceObject; // save a pointer that points to the beam element
        rInterfaceObject.PrintInfo(std::cout);
        std::cout << "\n-----" << std::endl;
        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());
        mpInterfaceObject->PrintInfo(std::cout);
        const auto g = mpInterfaceObject->pGetBaseGeometry();
        std::cout << "\n coordinate 1  " << (*g)[0].Coordinates() << std::endl;
        std::cout << "\n coordinate 2  " << (*g)[1].Coordinates() << std::endl;
        //std::cout << "interfaceObject geometry was extracted successfully" << std::endl;
        std::cout << "\n------------- END OF SEARCH -----------------" << std::endl;
    }
    
}

void BeamMapperInterfaceInfo::ComputeRotationMatrix()
{
    array_1d<double, 3> axisX;
    array_1d<double, 3> axisY;
    array_1d<double, 3> axisZ;
    
    const auto p_geom = mpInterfaceObject->pGetBaseGeometry();

    auto temp_v = (*p_geom)[1].Coordinates() - (*p_geom)[0].Coordinates();
    double lengthX = sqrt(temp_v[0]*temp_v[0] + temp_v[1]*temp_v[1] + temp_v[2]*temp_v[2]);
    axisX = (temp_v / lengthX); 
        
    double lengthY = sqrt(temp_v[0]*temp_v[0] + temp_v[1]*temp_v[1]);
    axisY[0] = -temp_v[1] / lengthY;
    axisY[1] =  temp_v[0] / lengthY;
    axisY[2] =  0;

    axisZ[0] = axisX[1]*axisY[2] - axisX[2]*axisY[1]; 
    axisZ[1] = axisX[2]*axisY[0] - axisX[0]*axisY[2];
    axisZ[2] = axisX[0]*axisY[1] - axisX[1]*axisY[0];

    //std::cout << "The unitary vector of axis x is : " << axisX << std::endl;
    //std::cout << "The unitary vector of axis y is : " << axisY << std::endl;
    //std::cout << "The unitary vector of axis z is : " << axisZ << std::endl;

    MatrixType _RotationMatrix( 3, 3, 0.0 );

    for(std::size_t j = 0; j < 3; j++)
    {
        _RotationMatrix( j, 0 ) = axisX[j];
        _RotationMatrix( j, 1 ) = axisY[j];
        _RotationMatrix( j, 2 ) = axisZ[j];
    }

    mRotationMatrixOfBeam = _RotationMatrix;

    std::cout << "The mRotationMatrixOfBeam of this BeamMapper Interface Info is" << mRotationMatrixOfBeam << std::endl;

}

// ************ BeamMapperLocalSystem function definitions
        
void BeamMapperLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{ std::cout << "CalculateAll is running here" << std::endl;
    if (mInterfaceInfos.size() > 0) {
        double distance;
        double min_distance = std::numeric_limits<double>::max();
        int found_idx = -1;
        for (IndexType i=0; i<mInterfaceInfos.size(); ++i) {
            // the approximations will be processed in the next step if necessary
            if (!mInterfaceInfos[i]->GetIsApproximation()) {
                mInterfaceInfos[i]->GetValue(distance, MapperInterfaceInfo::InfoType::Dummy);
                if (distance < min_distance) {
                    min_distance = distance;
                    found_idx = static_cast<int>(i); // TODO explicit conversion needed?
                    rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;
                }
            }
        }

        if (found_idx == -1) { // no valid projection exists => using an approximation
            int int_pairing_index;
            ProjectionUtilities::PairingIndex pairing_index;
            for (IndexType i=0; i<mInterfaceInfos.size(); ++i) {
                // now the approximations are being checked
                if (mInterfaceInfos[i]->GetIsApproximation()) {
                    mInterfaceInfos[i]->GetValue(int_pairing_index, MapperInterfaceInfo::InfoType::Dummy);
                    pairing_index = (ProjectionUtilities::PairingIndex)int_pairing_index;
                    mInterfaceInfos[i]->GetValue(distance, MapperInterfaceInfo::InfoType::Dummy);

                    if (pairing_index > mPairingIndex || (pairing_index == mPairingIndex && distance < min_distance)) {
                        mPairingIndex = pairing_index;
                        min_distance = distance;
                        found_idx = static_cast<int>(i); // TODO explicit conversion needed?
                        rPairingStatus = MapperLocalSystem::PairingStatus::Approximation;
                    }
                }
            }
        }

        KRATOS_ERROR_IF(found_idx == -1) << "Not even an approximation is found, this should not happen!"
            << std::endl; // TODO should this be an error?

        KRATOS_ERROR_IF(mPairingIndex == ProjectionUtilities::PairingIndex::Unspecified && mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) << "Not even an approximation is found (enum), this should not happen! " << found_idx << std::endl; // TODO should this be an error?

        std::vector<double> sf_values;

        mInterfaceInfos[found_idx]->GetValue(sf_values, MapperInterfaceInfo::InfoType::Dummy);

        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != sf_values.size()) {
            rLocalMappingMatrix.resize(1, sf_values.size(), false);
        }
        for (IndexType i=0; i<sf_values.size(); ++i) {
            rLocalMappingMatrix(0,i) = sf_values[i];
        }

        mInterfaceInfos[found_idx]->GetValue(rOriginIds, MapperInterfaceInfo::InfoType::Dummy);
        std::cout << "wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww " << std::endl;
        std::cout << "found_idx is : " << found_idx << std::endl;
        std::cout << "My OriginIDs in this interface Info are: " << rOriginIds << std::endl;

        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

        if (rDestinationIds.size() != 1) rDestinationIds.resize(1);
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);

        std::cout << "My DestinationIDs in this interface Info are: " << rDestinationIds << std::endl;
    }
    else ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
}

std::string BeamMapperLocalSystem::PairingInfo(const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    std::stringstream buffer;
    buffer << "BeamMapperLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) { // TODO leave here?
        buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
        if (mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) {
            mpNode->SetValue(PAIRING_STATUS, 0);
        } 
    }
    return buffer.str();
}

// ************ BeamMapper function definitions

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::ValidateInput()
{
    MapperUtilities::CheckInterfaceModelParts(0);

    Parameters mapper_default_settings(GetMapperDefaultSettings());
    mMapperSettings.ValidateAndAssignDefaults(mapper_default_settings);

    if (mMapperSettings["search_radius"].GetDouble() < 0.0) {
        const double search_radius = MapperUtilities::ComputeSearchRadius(
                                        mrModelPartOrigin,
                                        mrModelPartDestination,
                                        mMapperSettings["echo_level"].GetInt());
        mMapperSettings["search_radius"].SetDouble(search_radius);
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::Initialize()
{
    InitializeInterfaceCommunicator();
    InitializeInterface();
}

template<>
void BeamMapper<MapperDefinitions::SparseSpaceType,
    MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicator>(mrModelPartOrigin,
                                                                        mMapperLocalSystems,
                                                                        mMapperSettings);
}

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template<>
void BeamMapper<MapperDefinitions::MPISparseSpaceType,
    MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicatorMPI>(mrModelPartOrigin,
                                                                           mMapperLocalSystems,
                                                                           mMapperSettings);
}
#endif

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::InitializeInformationBeams(const Variable< array_1d<double, 3> >& rOriginVariablesDisplacements,
                                                                       const Variable< array_1d<double, 3> >& rOriginVariablesRotations)
{
    for( auto& r_local_sys : mMapperLocalSystems )
    {
        // Calculates rotation matrices
        if( r_local_sys->HasInterfaceInfo())
        {
            MatrixType _rotationMatrix_B(3, 3);
            VectorType _linearShapeValues(2);
            VectorType _hermitianShapeValues(4);
            VectorType _hermitanDerShapeValues(4);
            GeometryType r_geom;

            r_local_sys->CalculateRotationMatrixInterfaceInfos(_rotationMatrix_B,
                                                               _linearShapeValues,
                                                               _hermitianShapeValues,
                                                               _hermitanDerShapeValues,
                                                               r_geom);
            
            //std::cout << "1 node " << r_geom[0].Coordinates() << std::endl;
            //std::cout << "2 node " << r_geom[1].Coordinates() << std::endl;
            
            const std::vector<std::string> var_comps{"_X", "_Y", "_Z"};
            VectorType displacementNode1_G(3); //Expresses in global coordinates
            VectorType displacementNode2_G(3); //Expresses in global coordinates
            VectorType rotationNode1_G(3); //Expresses in global coordinates
            VectorType rotationNode2_G(3); //Expresses in global coordinates

            VectorType displacementNode1_B(3); //Expresses in beam coordinates
            VectorType displacementNode2_B(3); //Expresses in beam coordinates
            VectorType rotationNode1_B(3); //Expresses in beam coordinates
            VectorType rotationNode2_B(3); //Expresses in beam coordinates

            size_t k = 0;

            for (const auto& var_ext : var_comps)
            {
                const auto& var_origin_disp = KratosComponents<ComponentVariableType>::Get(rOriginVariablesDisplacements.Name() + var_ext);
                displacementNode1_G(k) = r_geom[0].FastGetSolutionStepValue(var_origin_disp);
                displacementNode2_G(k) = r_geom[1].FastGetSolutionStepValue(var_origin_disp);


                const auto& var_origin_rot = KratosComponents<ComponentVariableType>::Get(rOriginVariablesRotations.Name() + var_ext);
                rotationNode1_G(k) = r_geom[0].FastGetSolutionStepValue(var_origin_rot);
                rotationNode2_G(k) = r_geom[1].FastGetSolutionStepValue(var_origin_rot);
                k++;
            }
            //std::cout << "displacement of node 1 is" << displacementNode1 << std::endl;
            //std::cout << "displacement of node 2 is" << displacementNode2 << std::endl;
            //std::cout << "rotation of node 1 is" << rotationNode1 << std::endl;
            //std::cout << "rotation of node 2 is" << rotationNode2 << std::endl; 

            MatrixType _rotationMatrixInverse( 3, 3 );
            double determinant;
            MathUtils<double>::InvertMatrix3(_rotationMatrix_B, _rotationMatrixInverse, determinant );
            
            TDenseSpace::Mult( _rotationMatrixInverse, displacementNode1_G, displacementNode1_B );
            TDenseSpace::Mult( _rotationMatrixInverse, rotationNode1_G, rotationNode1_B );
            TDenseSpace::Mult( _rotationMatrixInverse, displacementNode2_G, displacementNode2_B );
            TDenseSpace::Mult( _rotationMatrixInverse, rotationNode2_G, rotationNode2_B );
            
            std::cout << "rotated displacement in node 1 is : " << displacementNode1_B << std::endl;
            std::cout << "rotated rotation in node 1 is : " << rotationNode1_B << std::endl;
            std::cout << "rotated displacement in node 2 is : " << displacementNode2_B << std::endl;
            std::cout << "rotated rotation in node 2 is : " << rotationNode2_B << std::endl;

            //double _VX, _VY, _VZ, _ThetaX, _ThetaY, _ThetaZ;

            // Initializing matrix of shape functions
            MatrixType _ShapeFunctionsMatrix(6, 12, 0.0);

            _ShapeFunctionsMatrix(0 , 0) = _linearShapeValues(0);
            _ShapeFunctionsMatrix(0 , 6) = _linearShapeValues(1);
            _ShapeFunctionsMatrix(1 , 1) = _hermitianShapeValues(0);
            _ShapeFunctionsMatrix(1 , 5) = _hermitianShapeValues(1);
            _ShapeFunctionsMatrix(1 , 7) = _hermitianShapeValues(2);
            _ShapeFunctionsMatrix(1 , 11) = _hermitianShapeValues(3);
            _ShapeFunctionsMatrix(2 , 2) = _hermitianShapeValues(0);
            _ShapeFunctionsMatrix(2 , 4) = -_hermitianShapeValues(1);
            _ShapeFunctionsMatrix(2 , 8) = _hermitianShapeValues(2);
            _ShapeFunctionsMatrix(2 , 10) = -_hermitianShapeValues(3);
            
            _ShapeFunctionsMatrix(3 , 3) = _linearShapeValues(0);
            _ShapeFunctionsMatrix(3 , 9) = _linearShapeValues(1);
            _ShapeFunctionsMatrix(4 , 2) = -_hermitanDerShapeValues(0);
            _ShapeFunctionsMatrix(4 , 4) = _hermitanDerShapeValues(1);
            _ShapeFunctionsMatrix(4 , 8) = -_hermitanDerShapeValues(2);
            _ShapeFunctionsMatrix(4 , 10) = -_hermitanDerShapeValues(3);
            _ShapeFunctionsMatrix(5 , 1) = _hermitanDerShapeValues(0);
            _ShapeFunctionsMatrix(5 , 5) = _hermitanDerShapeValues(1);
            _ShapeFunctionsMatrix(5 , 7) = _hermitanDerShapeValues(2);
            _ShapeFunctionsMatrix(5 , 11) = _hermitanDerShapeValues(3);
            
            VectorType _DOF_Vector(12);
            for (size_t i = 0; i < 3; i++){
                _DOF_Vector(i) = displacementNode1_B(i);
                _DOF_Vector(i + 3) = rotationNode1_B(i);
                _DOF_Vector(i + 6) = displacementNode2_B(i);
                _DOF_Vector(i + 9) = rotationNode2_B(i);
            }

            std::cout << "_ShapeFunctionsMatrix is : " << _ShapeFunctionsMatrix << std::endl;
            std::cout << "_DOF_Vector is : " << _DOF_Vector << std::endl;
            
            VectorType _DisplacementsRotationsP(6);
            TDenseSpace::Mult(_ShapeFunctionsMatrix, _DOF_Vector, _DisplacementsRotationsP);

            std::cout << "Interpolated displacements and rotations are : " << _DisplacementsRotationsP << std::endl;
            std::cout << "-------------------------- END OF INTERPOLATION ---------------------------" << std::endl; 
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    // Here we can see that the local systems are done with the Destination Model Part 
    CreateMapperLocalSystems(mrModelPartDestination.GetCommunicator(), mMapperLocalSystems);

    // Lets find the information for the local systems CHANGE NAME OF THIS FUNCTION LATER
    BuildMappingMatrix(MappingOptions);

}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::BuildMappingMatrix(Kratos::Flags MappingOptions)
{
    MapperUtilities::AssignInterfaceEquationIds(mrModelPartOrigin.GetCommunicator());
    MapperUtilities::AssignInterfaceEquationIds(mrModelPartDestination.GetCommunicator());

    KRATOS_ERROR_IF_NOT(mpIntefaceCommunicator) << "mpInterfaceCommunicator is a nullptr" << std::endl;
    const MapperInterfaceInfoUniquePointerType p_ref_interface_info = GetMapperInterfaceInfo();
    std::cout << "Calculating shape function values" << std::endl;
    mpIntefaceCommunicator->ExchangeInterfaceData(mrModelPartDestination.GetCommunicator(),
                                                  MappingOptions,
                                                  p_ref_interface_info);
    //const int echo_level = mMapperSettings["echo_level"].GetInt();
    
    //KRATOS_WATCH("BEFORE Building MMatrix")
    //MappingMatrixUtilities::BuildMappingMatrix<TSparseSpace, TDenseSpace>(
    //    mpMappingMatrix,
    //    mpInterfaceVectorContainerOrigin->pGetVector(),
    //    mpInterfaceVectorContainerDestination->pGetVector(),
    //    mpInterfaceVectorContainerOrigin->GetModelPart(),
    //    mpInterfaceVectorContainerDestination->GetModelPart(),
    //    mMapperLocalSystems,
    //    echo_level);
}

// Extern template instantiation

template class BeamMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;



}  // namespace Kratos.
