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
// See Master-Thesis E. G. Loera Villeda
// "  "

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
    
    //std::cout << "rShapeFunctionValues : " << linear_shape_function_values << std::endl;
    //std::cout << "rHermitianShapeFunctionValues : " << hermitian_shape_function_values << std::endl;
    //std::cout << "rHermitianShapeFunctionValuesDerivatives : " << hermitian_shape_function_derivatives_values << std::endl;
    

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

        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());

        //std::cout << "\n coordinate 1  " << (*g)[0].Coordinates() << std::endl;
        //std::cout << "\n coordinate 2  " << (*g)[1].Coordinates() << std::endl;
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

    //std::cout << "The mRotationMatrixOfBeam of this BeamMapper Interface Info is" << mRotationMatrixOfBeam << std::endl;

}

// ************ BeamMapperLocalSystem function definitions
        
void BeamMapperLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{ //std::cout << "CalculateAll is running here" << std::endl;
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
        //std::cout << "wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww " << std::endl;
        //std::cout << "found_idx is : " << found_idx << std::endl;
        //std::cout << "My OriginIDs in this interface Info are: " << rOriginIds << std::endl;

        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

        if (rDestinationIds.size() != 1) rDestinationIds.resize(1);
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);

        //std::cout << "My DestinationIDs in this interface Info are: " << rDestinationIds << std::endl;
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
                                                                       const Variable< array_1d<double, 3> >& rOriginVariablesRotations,
                                                                       const Variable< array_1d<double, 3> >& rDestinationVariableDisplacement)
{   size_t i = 0;
    for( auto& r_local_sys : mMapperLocalSystems )
    {   
        std::cout << "----------- LOCAL SYSTEM " << i << std::endl;
        i++;

        if( r_local_sys->HasInterfaceInfo())
        {
            MatrixType _rotationMatrix_G_B(3, 3);
            VectorType _translationVector_B_P(3);
            VectorType _linearShapeValues(2);
            VectorType _hermitianShapeValues(4);
            VectorType _hermitanDerShapeValues(4);
            GeometryType _r_geom;
            NodePointerType _pNode;

            r_local_sys->CalculateRotationMatrixInterfaceInfos(_rotationMatrix_G_B,
                                                               _translationVector_B_P,
                                                               _linearShapeValues,
                                                               _hermitianShapeValues,
                                                               _hermitanDerShapeValues,
                                                               _r_geom,
                                                               _pNode);

            std::cout << "The coordinates of the surface mesh node are : " << _pNode->Coordinates() << std::endl;            
            std::cout << "The _rotationMatrix_G_B = " << _rotationMatrix_G_B << std::endl;
            KRATOS_ERROR_IF_NOT(_pNode) << "Node is a nullptr"<< std::endl;

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
                displacementNode1_G(k) = _r_geom[0].FastGetSolutionStepValue(var_origin_disp);
                displacementNode2_G(k) = _r_geom[1].FastGetSolutionStepValue(var_origin_disp);


                const auto& var_origin_rot = KratosComponents<ComponentVariableType>::Get(rOriginVariablesRotations.Name() + var_ext);
                rotationNode1_G(k) = _r_geom[0].FastGetSolutionStepValue(var_origin_rot);
                rotationNode2_G(k) = _r_geom[1].FastGetSolutionStepValue(var_origin_rot);
                k++;
            }
            std::cout << "displacement of node 1 is" << displacementNode1_G << std::endl;
            std::cout << "displacement of node 2 is" << displacementNode2_G << std::endl;
            std::cout << "rotation of node 1 is" << rotationNode1_G << std::endl;
            std::cout << "rotation of node 2 is" << rotationNode2_G << std::endl; 

            MatrixType _rotationMatrix_B_G( 3, 3 );
            double determinant;
            MathUtils<double>::InvertMatrix3(_rotationMatrix_G_B, _rotationMatrix_B_G, determinant );
            
            TDenseSpace::Mult( _rotationMatrix_B_G, displacementNode1_G, displacementNode1_B );
            TDenseSpace::Mult( _rotationMatrix_B_G, rotationNode1_G, rotationNode1_B );
            TDenseSpace::Mult( _rotationMatrix_B_G, displacementNode2_G, displacementNode2_B );
            TDenseSpace::Mult( _rotationMatrix_B_G, rotationNode2_G, rotationNode2_B );
            
            //std::cout << "rotated displacement in node 1 is : " << displacementNode1_B << std::endl;
            //std::cout << "rotated rotation in node 1 is : " << rotationNode1_B << std::endl;
            //std::cout << "rotated displacement in node 2 is : " << displacementNode2_B << std::endl;
            //std::cout << "rotated rotation in node 2 is : " << rotationNode2_B << std::endl;

            // Initializing matrix of shape functions
            MatrixType _ShapeFunctionsMatrix(6, 12, 0.0);

            VectorType beamVector = _r_geom[1].Coordinates() - _r_geom[0].Coordinates();
            double length_beamVector = norm_2(beamVector);
            _ShapeFunctionsMatrix(0 , 0) = _linearShapeValues(0);
            _ShapeFunctionsMatrix(0 , 6) = _linearShapeValues(1);
            _ShapeFunctionsMatrix(1 , 1) = _hermitianShapeValues(0);
            _ShapeFunctionsMatrix(1 , 5) = _hermitianShapeValues(1) * length_beamVector;
            _ShapeFunctionsMatrix(1 , 7) = _hermitianShapeValues(2);
            _ShapeFunctionsMatrix(1 , 11) = _hermitianShapeValues(3) * length_beamVector;
            _ShapeFunctionsMatrix(2 , 2) = _hermitianShapeValues(0);
            _ShapeFunctionsMatrix(2 , 4) = -_hermitianShapeValues(1) * length_beamVector;
            _ShapeFunctionsMatrix(2 , 8) = _hermitianShapeValues(2);
            _ShapeFunctionsMatrix(2 , 10) = -_hermitianShapeValues(3) * length_beamVector;
            
            _ShapeFunctionsMatrix(3 , 3) = _linearShapeValues(0);
            _ShapeFunctionsMatrix(3 , 9) = _linearShapeValues(1);
            _ShapeFunctionsMatrix(4 , 2) = -_hermitanDerShapeValues(0) / length_beamVector;
            _ShapeFunctionsMatrix(4 , 4) = _hermitanDerShapeValues(1);
            _ShapeFunctionsMatrix(4 , 8) = -_hermitanDerShapeValues(2) / length_beamVector;
            _ShapeFunctionsMatrix(4 , 10) = -_hermitanDerShapeValues(3);
            _ShapeFunctionsMatrix(5 , 1) = _hermitanDerShapeValues(0) / length_beamVector;
            _ShapeFunctionsMatrix(5 , 5) = _hermitanDerShapeValues(1);
            _ShapeFunctionsMatrix(5 , 7) = _hermitanDerShapeValues(2) / length_beamVector;
            _ShapeFunctionsMatrix(5 , 11) = _hermitanDerShapeValues(3);
            
            VectorType _DOF_Vector(12);
            for (size_t i = 0; i < 3; i++){
                _DOF_Vector(i) = displacementNode1_B(i);
                _DOF_Vector(i + 3) = rotationNode1_B(i);
                _DOF_Vector(i + 6) = displacementNode2_B(i);
                _DOF_Vector(i + 9) = rotationNode2_B(i);
            }

            //std::cout << "_ShapeFunctionsMatrix is : " << _ShapeFunctionsMatrix << std::endl;
            //std::cout << "_DOF_Vector is : " << _DOF_Vector << std::endl;
            
            VectorType _DisplacementsRotationsP(6);
            TDenseSpace::Mult(_ShapeFunctionsMatrix, _DOF_Vector, _DisplacementsRotationsP);

            //std::cout << "Interpolated displacements and rotations are : " << _DisplacementsRotationsP << std::endl;
            //std::cout << "-------------------------- END OF INTERPOLATION ---------------------------" << std::endl; 

            VectorType axisX(3, 0.0);
            VectorType axisY(3, 0.0);
            VectorType axisZ(3, 0.0);
            axisX(0) = 1.0;
            axisY(1) = 1.0;
            axisZ(2) = 1.0; 

            MatrixType rotation_X(3, 3);
            MatrixType rotation_Y(3, 3);
            MatrixType rotation_Z(3, 3);

            CalculateRotationMatrixWithAngle( axisX, _DisplacementsRotationsP(3), rotation_X);
            CalculateRotationMatrixWithAngle( axisY, _DisplacementsRotationsP(4), rotation_Y);
            CalculateRotationMatrixWithAngle( axisZ, _DisplacementsRotationsP(5), rotation_Z);

            MatrixType RotationMatrix_P(3, 3);
            MatrixType tmpRotation(3, 3);
            tmpRotation = prod(rotation_X ,rotation_Z);
            RotationMatrix_P = prod(tmpRotation, rotation_Y);

            std::cout << "R_P is : " << RotationMatrix_P << std::endl;

            VectorType TranslationVector_P(3);
            TranslationVector_P(0) = _DisplacementsRotationsP(0);
            TranslationVector_P(1) = _DisplacementsRotationsP(1);
            TranslationVector_P(2) = _DisplacementsRotationsP(2);

            std::cout << "t_P is : " << TranslationVector_P << std::endl;

            // Rigid body operation 
            // phi_G( X ) = R_G_B * R_P * (R_G_B ^ T) * X - R_G_B * R_P * (R_G_B ^ T) * t_B_P + R_G_B * t_B + t_B_P

            MatrixType MatrixProduct(3, 3); 
            MathUtils<double>::BDBtProductOperation(MatrixProduct, RotationMatrix_P, _rotationMatrix_G_B);
            
            VectorType firstProduct(3), secondProduct(3), thirdProduct(3), phi_G(3), displacement(3), _X(3);

            _X( 0 ) = (r_local_sys->Coordinates())[ 0 ];
            _X( 1 ) = (r_local_sys->Coordinates())[ 1 ];
            _X( 2 ) = (r_local_sys->Coordinates())[ 2 ];
            TDenseSpace::Mult(MatrixProduct, _X, firstProduct );
            TDenseSpace::Mult(MatrixProduct, _translationVector_B_P, secondProduct );
            TDenseSpace::Mult(_rotationMatrix_G_B, TranslationVector_P, thirdProduct );

            phi_G( 0 ) = firstProduct( 0 ) - secondProduct( 0 ) + thirdProduct( 0 ) + _translationVector_B_P( 0 );
            phi_G( 1 ) = firstProduct( 1 ) - secondProduct( 1 ) + thirdProduct( 1 ) + _translationVector_B_P( 1 );
            phi_G( 2 ) = firstProduct( 2 ) - secondProduct( 2 ) + thirdProduct( 2 ) + _translationVector_B_P( 2 ); 

            displacement( 0 ) = phi_G( 0 ) - _X( 0 );
            displacement( 1 ) = phi_G( 1 ) - _X( 1 );
            displacement( 2 ) = phi_G( 2 ) - _X( 2 );

            std::cout << "phi_G is : " << phi_G << std::endl;

            std::cout << "_X is : " << _X << std::endl;

            std::cout << "The final displacement is : " << displacement << std::endl;

            size_t c = 0;
            for( const auto& var_ext : var_comps )
            {
                const auto& var_destination_disp = KratosComponents<ComponentVariableType>::Get(rDestinationVariableDisplacement.Name() + var_ext);
                _pNode->FastGetSolutionStepValue(var_destination_disp) = displacement( c );
                c++;
            }
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::InitializeInformationBeamsCorotation(const Variable< array_1d<double, 3> >& rOriginVariablesDisplacements,
                                                                                 const Variable< array_1d<double, 3> >& rOriginVariablesRotations,
                                                                                 const Variable< array_1d<double, 3> >& rDestinationVariableDisplacement)
{   size_t i = 0;
    for( auto& r_local_sys : mMapperLocalSystems )
    {   
        std::cout << "----------- LOCAL SYSTEM " << i << std::endl;
        i++;

        if( r_local_sys->HasInterfaceInfo())
        {
            MatrixType _rotationMatrix_G_B(3, 3);
            VectorType _translationVector_B_P(3);
            VectorType _linearShapeValues(2);
            VectorType _hermitianShapeValues(4);
            VectorType _hermitanDerShapeValues(4);
            GeometryType _r_geom;
            NodePointerType _pNode;

            r_local_sys->CalculateRotationMatrixInterfaceInfos(_rotationMatrix_G_B,
                                                               _translationVector_B_P,
                                                               _linearShapeValues,
                                                               _hermitianShapeValues,
                                                               _hermitanDerShapeValues,
                                                               _r_geom,
                                                               _pNode);

            std::cout << "The coordinates of the surface mesh node are : " << _pNode->Coordinates() << std::endl;            
            std::cout << "The _rotationMatrix_G_B = " << _rotationMatrix_G_B << std::endl;
            
            KRATOS_ERROR_IF_NOT(_pNode) << "Node is a nullptr"<< std::endl;

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
                displacementNode1_G(k) = _r_geom[0].FastGetSolutionStepValue(var_origin_disp);
                displacementNode2_G(k) = _r_geom[1].FastGetSolutionStepValue(var_origin_disp);


                const auto& var_origin_rot = KratosComponents<ComponentVariableType>::Get(rOriginVariablesRotations.Name() + var_ext);
                rotationNode1_G(k) = _r_geom[0].FastGetSolutionStepValue(var_origin_rot);
                rotationNode2_G(k) = _r_geom[1].FastGetSolutionStepValue(var_origin_rot);
                k++;
            }
            std::cout << "displacement of node 1 is" << displacementNode1_G << std::endl;
            std::cout << "displacement of node 2 is" << displacementNode2_G << std::endl;
            std::cout << "rotation of node 1 is" << rotationNode1_G << std::endl;
            std::cout << "rotation of node 2 is" << rotationNode2_G << std::endl;  

            MatrixType _rotationMatrix_B_G( 3, 3 );
            double determinant;
            MathUtils<double>::InvertMatrix3(_rotationMatrix_G_B, _rotationMatrix_B_G, determinant );
            
            // Transforming the displacements to the BCS
            TDenseSpace::Mult( _rotationMatrix_B_G, displacementNode1_G, displacementNode1_B );
            TDenseSpace::Mult( _rotationMatrix_B_G, displacementNode2_G, displacementNode2_B );

            // Transforming the nodal rotations to the BCS
            double angle1 = norm_2(rotationNode1_G);
            VectorType n1(3, 0.0);
            n1 = rotationNode1_G;
            if (angle1 != 0){
                n1 /= angle1;
            }
            MatrixType Rotation_G_1(3, 3);
            CalculateRotationMatrixWithAngle(n1, angle1, Rotation_G_1);
            MatrixType Rotation_B_1(3, 3);
            MathUtils<double>::BtDBProductOperation(Rotation_B_1, Rotation_G_1,_rotationMatrix_G_B);
            //std::cout << "angle1 = " << angle1 << std::endl;
            //std::cout << "n1 = " << n1 << std::endl;
            std::cout << "Rotation_G_1 = " << Rotation_G_1 << std::endl;
            std::cout << "Rotation_B_1 = " << Rotation_B_1 << std::endl;
            
            double angle2 = norm_2(rotationNode2_G);
            VectorType n2(3 , 0.0);
            n2 = rotationNode2_G;
            if(angle2 != 0){
                n2 /= angle2;
            }
            MatrixType Rotation_G_2(3, 3);
            CalculateRotationMatrixWithAngle(n2, angle2, Rotation_G_2);
            MatrixType Rotation_B_2(3, 3);
            MathUtils<double>::BtDBProductOperation(Rotation_B_2, Rotation_G_2,_rotationMatrix_G_B);
            //std::cout << "angle2 = " << angle2 << std::endl;
            //std::cout << "n2 = " << n2 << std::endl;
            std::cout << "Rotation_G_2 = " << Rotation_G_2 << std::endl;
            std::cout << "Rotation_B_2 = " << Rotation_B_2 << std::endl;

            // Calculating R_d  and t_d for Phi_d 
            VectorType e_x_d_G(3), e_x_d_B(3);
            e_x_d_G = _r_geom[1].Coordinates() + displacementNode2_G - (_r_geom[0].Coordinates() + displacementNode1_G); // this vector is described in global system 
            //std::cout << "_r_geom[0].Coordinates() = " << _r_geom[0] << std::endl;
            //std::cout << "_r_geom[1].Coordinates() = " << _r_geom[1] << std::endl;
            //std::cout << "displacementNode1_G =  " << displacementNode1_G << std::endl;
            //std::cout << "displacementNode2_G =  " << displacementNode2_G << std::endl;
            e_x_d_G /= norm_2(e_x_d_G);
            //std::cout << "e_x_d_G = " << e_x_d_G << std::endl;
            //std::cout << "this takes me to the beam C.S. _rotationMatrix_B_G = " << _rotationMatrix_B_G << std::endl; 
            //std::cout << " _rotationMatrix_G_B = " << _rotationMatrix_G_B << std::endl; 
            TDenseSpace::Mult(_rotationMatrix_B_G, e_x_d_G, e_x_d_B); // transforming e_x_d to the beam coordinate system
            e_x_d_B /= norm_2(e_x_d_B);
            //std::cout << "e_x_d_B = " << e_x_d_B << std::endl;
            
            VectorType e_x(3, 0.0), n_d(3, 0.0);
            MatrixType R_d(3, 3), _I(3, 3, 0.0), _R(3, 3, 0.0), I_2nT(3, 3, 0.0);
            e_x(0) = 1.0;
            _I(0,0) = 1.0;
            _I(1,1) = 1.0;
            _I(2,2) = 1.0;
            _R(0, 0) = -1.0;
            _R(1, 1) = 1.0;
            _R(2, 2) = 1.0;
            n_d =  e_x + e_x_d_B;
            n_d /= norm_2(n_d);
            std::cout << "n_d = " << n_d << std::endl;
            std::cout << "n_d * n_d_T = " << MathUtils<double>::TensorProduct3(n_d, n_d) << std::endl;
            I_2nT = _I - 2 * MathUtils<double>::TensorProduct3(n_d, n_d);
            R_d = prod(I_2nT, _R);
            std::cout << "Rd is = " << R_d << std::endl;
           
            VectorType t_d(3);
            t_d(0) = _linearShapeValues(0) * displacementNode1_B(0) + _linearShapeValues(1) * displacementNode2_B(0);
            t_d(1) = _linearShapeValues(0) * displacementNode1_B(1) + _linearShapeValues(1) * displacementNode2_B(1);
            t_d(2) = _linearShapeValues(0) * displacementNode1_B(2) + _linearShapeValues(1) * displacementNode2_B(2);
            std::cout << "td is = " << t_d << std::endl;

            // First, transform the nodal rotations (expresed in the Beam Coordinate System) to the Co-rotational Coordinate System
            // Calculating Co-rotate basis vectors
            VectorType corotateXAxis(3), corotateYAxis(3), corotateZAxis(3);
            corotateXAxis = e_x_d_B; // Expresed in local system
            corotateYAxis(0) = -2.0 * n_d(0) * n_d(1);
            corotateYAxis(1) = 1.0 - 2.0 * n_d(1) * n_d(1); // diagonal
            corotateYAxis(2) = -2.0 * n_d(2) * n_d(1);
            corotateZAxis(0) = -2.0 * n_d(0) * n_d(2);
            corotateZAxis(1) = -2.0 * n_d(1) * n_d(2);
            corotateZAxis(2) = 1.0 - 2.0 * n_d(2) * n_d(2); // diagonal
            
            MatrixType Rotation_B_COR(3, 3), Rotation_COR_B(3, 3);
            Rotation_B_COR(0, 0) = corotateXAxis(0);
            Rotation_B_COR(1, 0) = corotateXAxis(1);
            Rotation_B_COR(2, 0) = corotateXAxis(2);
            Rotation_B_COR(0, 1) = corotateYAxis(0);
            Rotation_B_COR(1, 1) = corotateYAxis(1);
            Rotation_B_COR(2, 1) = corotateYAxis(2);
            Rotation_B_COR(0, 2) = corotateZAxis(0);
            Rotation_B_COR(1, 2) = corotateZAxis(1);
            Rotation_B_COR(2, 2) = corotateZAxis(2);
            
            //std::cout << "The other version of R_d is = " << Rotation_B_COR << std::endl;
            //MathUtils<double>::InvertMatrix3(Rotation_B_COR, Rotation_COR_B, determinant);
            //MatrixType Rotation_COR_1(3, 3), Rotation_COR_2(3, 3);
            //Rotation_COR_1 = prod(Rotation_COR_B, Rotation_B_1);
            //Rotation_COR_2 = prod(Rotation_COR_B, Rotation_B_2);

            // Calculating R_s
            MatrixType Rl1_Rs(3, 3), Rl2_Rs(3, 3), R_d_T(3, 3);
            MathUtils<double>::InvertMatrix3(R_d, R_d_T, determinant);
            Rl1_Rs = prod(R_d_T, Rotation_B_1); // Rl1 * Rs = R_d_T * Phi_1, displacements cancel out
            Rl2_Rs = prod(R_d_T, Rotation_B_2); // Rl2 * Rs = R_d_T * Phi_2, displacements cancel out
            std::cout << "Rl1_Rs = " << Rl1_Rs << std::endl;
            std::cout << "Rl2_Rs = " << Rl2_Rs << std::endl; 
            VectorType v_Rl1_Rs(3), v_Rl2_Rs(3);
            getRotationVector(Rl1_Rs, v_Rl1_Rs);
            getRotationVector(Rl2_Rs, v_Rl2_Rs);
            std::cout << "v_Rl1_Rs = " << v_Rl1_Rs << std::endl;
            std::cout << "v_Rl2_Rs = " << v_Rl2_Rs << std::endl;

            double theta_s = (v_Rl1_Rs(0) + v_Rl2_Rs(0)) / 2;
            MatrixType R_s(3, 3);
            VectorType e_x_s(3, 0.0);
            e_x_s(0) = 1.0;
            //std::cout << "theta_s = " << theta_s << std::endl;
            CalculateRotationMatrixWithAngle(e_x_s, theta_s, R_s);
            std::cout << "R_s = " << R_s << std::endl;

            // Calculating rotation vector l on node 1 and 2 of the beam
            MatrixType Rl1(3, 3), Rl2(3, 3), R_s_T(3, 3);
            VectorType v_Rl1(3), v_Rl2(3);
            MathUtils<double>::InvertMatrix3(R_s, R_s_T, determinant);
            //std::cout << " R_s_T  = " << R_s_T << std::endl;
            //std::cout << " Rs_Rl1  = " << Rs_Rl1 << std::endl;
            //std::cout << " Rs_Rl2  = " << Rs_Rl2 << std::endl;
            Rl1 = prod(Rl1_Rs, R_s_T);
            Rl2 = prod(Rl2_Rs, R_s_T);
            std::cout << " Rl1  = " << Rl1 << std::endl;
            std::cout << " Rl2  = " << Rl2 << std::endl;
            getRotationVector(Rl1, v_Rl1);
            getRotationVector(Rl2, v_Rl2); 
            std::cout << " v_Rl1  = " << v_Rl1 << std::endl;
            std::cout << " v_Rl2  = " << v_Rl2 << std::endl;

            MatrixType _ShapeFunctionsMatrix(6, 12, 0.0);
            VectorType corBeamVector(3, 0.0);
            corBeamVector = _r_geom[1].Coordinates() + displacementNode2_G - (_r_geom[0].Coordinates() + displacementNode1_G); // this vector is described in global system 
            double length_beamVector = norm_2(corBeamVector);            

            _ShapeFunctionsMatrix(0 , 0) = _linearShapeValues(0);
            _ShapeFunctionsMatrix(0 , 6) = _linearShapeValues(1);
            _ShapeFunctionsMatrix(1 , 1) = _hermitianShapeValues(0);
            _ShapeFunctionsMatrix(1 , 5) = _hermitianShapeValues(1) * length_beamVector;
            _ShapeFunctionsMatrix(1 , 7) = _hermitianShapeValues(2);
            _ShapeFunctionsMatrix(1 , 11) = _hermitianShapeValues(3) * length_beamVector;
            _ShapeFunctionsMatrix(2 , 2) = _hermitianShapeValues(0);
            _ShapeFunctionsMatrix(2 , 4) = -_hermitianShapeValues(1) * length_beamVector;
            _ShapeFunctionsMatrix(2 , 8) = _hermitianShapeValues(2);
            _ShapeFunctionsMatrix(2 , 10) = -_hermitianShapeValues(3) * length_beamVector;
            
            _ShapeFunctionsMatrix(3 , 3) = _linearShapeValues(0);
            _ShapeFunctionsMatrix(3 , 9) = _linearShapeValues(1);
            _ShapeFunctionsMatrix(4 , 2) = -_hermitanDerShapeValues(0) / length_beamVector;
            _ShapeFunctionsMatrix(4 , 4) = _hermitanDerShapeValues(1);
            _ShapeFunctionsMatrix(4 , 8) = -_hermitanDerShapeValues(2) / length_beamVector;
            _ShapeFunctionsMatrix(4 , 10) = -_hermitanDerShapeValues(3);
            _ShapeFunctionsMatrix(5 , 1) = _hermitanDerShapeValues(0) / length_beamVector;
            _ShapeFunctionsMatrix(5 , 5) = _hermitanDerShapeValues(1);
            _ShapeFunctionsMatrix(5 , 7) = _hermitanDerShapeValues(2) / length_beamVector;
            _ShapeFunctionsMatrix(5 , 11) = _hermitanDerShapeValues(3); 

            std::cout << "Shape functions matrix = " << _ShapeFunctionsMatrix << std::endl;
            VectorType _DOF_Vector_l(12);
            for (size_t i = 0; i < 3; i++){
                _DOF_Vector_l(i) = 0.0;
                _DOF_Vector_l(i + 3) = v_Rl1(i);
                _DOF_Vector_l(i + 6) = 0.0;
                _DOF_Vector_l(i + 9) = v_Rl2(i);
            }
            std::cout << "Vector DOF_l is : " << _DOF_Vector_l << std::endl;

            VectorType _DisplacementsRotations_L(6);
            TDenseSpace::Mult(_ShapeFunctionsMatrix, _DOF_Vector_l, _DisplacementsRotations_L);

            std::cout << "displacements and rotations for phi l are : _DisplacementsRotations_L = " << _DisplacementsRotations_L << std::endl;

            // Constructing phi_l
            VectorType axis_l(3), t_l(3);
            MatrixType R_l(3, 3), R_l_temp(3, 3), R_x(3, 3), R_y(3, 3), R_z(3, 3);
            t_l(0) = _DisplacementsRotations_L(0); // t_l
            t_l(1) = _DisplacementsRotations_L(1); // t_l
            t_l(2) = _DisplacementsRotations_L(2); // t_l
            axis_l(0) = _DisplacementsRotations_L(3);
            axis_l(1) = _DisplacementsRotations_L(4);
            axis_l(2) = _DisplacementsRotations_L(5); 
            double angle_l = sqrt( pow(_DisplacementsRotations_L(3), 2) + pow(_DisplacementsRotations_L(4), 2) + pow(_DisplacementsRotations_L(5), 2) );
            axis_l = axis_l/angle_l;

            CalculateRotationMatrixWithAngle(axis_l, angle_l, R_l); // R_l
            std::cout << "t_l = " << t_l << std::endl;
            // Calculating R_P and t_P
            MatrixType R_P(3, 3), R_P_temp(3, 3);
            std::cout << "------ Before multiplying -------- " << std::endl;
            std::cout << "R_d = " << R_d << std::endl;
            std::cout << "R_s = " << R_s << std::endl;
            std::cout << "R_l = " << R_l << std::endl;
            R_P_temp = prod(R_d, R_l);
            R_P = prod(R_P_temp, R_s);

            VectorType t_P(3), t_P_temp(3, 0.0);
            TDenseSpace::Mult(R_d, t_l, t_P_temp);
            t_P = t_P_temp + t_d;

            std::cout << "R_P is : " << R_P << std::endl;
            std::cout << "t_P is : " << t_P << std::endl;

            // Rigid body operation 
            // phi_G( X ) = R_G_B * R_P * (R_G_B ^ T) * X - R_G_B * R_P * (R_G_B ^ T) * t_B_P + R_G_B * t_B + t_B_P

            MatrixType MatrixProduct(3, 3); 
            MathUtils<double>::BDBtProductOperation(MatrixProduct, R_P, _rotationMatrix_G_B);
            
            VectorType firstProduct(3), secondProduct(3), thirdProduct(3), phi_G(3), displacement(3), _X(3);

            _X( 0 ) = (r_local_sys->Coordinates())[ 0 ];
            _X( 1 ) = (r_local_sys->Coordinates())[ 1 ];
            _X( 2 ) = (r_local_sys->Coordinates())[ 2 ];
            TDenseSpace::Mult(MatrixProduct, _X, firstProduct );
            TDenseSpace::Mult(MatrixProduct, _translationVector_B_P, secondProduct );
            TDenseSpace::Mult(_rotationMatrix_G_B, t_P, thirdProduct );

            phi_G( 0 ) = firstProduct( 0 ) - secondProduct( 0 ) + thirdProduct( 0 ) + _translationVector_B_P( 0 );
            phi_G( 1 ) = firstProduct( 1 ) - secondProduct( 1 ) + thirdProduct( 1 ) + _translationVector_B_P( 1 );
            phi_G( 2 ) = firstProduct( 2 ) - secondProduct( 2 ) + thirdProduct( 2 ) + _translationVector_B_P( 2 ); 

            displacement( 0 ) = phi_G( 0 ) - _X( 0 );
            displacement( 1 ) = phi_G( 1 ) - _X( 1 );
            displacement( 2 ) = phi_G( 2 ) - _X( 2 );

            std::cout << "phi_G is : " << phi_G << std::endl;

            std::cout << "_X is : " << _X << std::endl;

            std::cout << "The final displacement is : " << displacement << std::endl;

            size_t c = 0;
            for( const auto& var_ext : var_comps )
            {
                const auto& var_destination_disp = KratosComponents<ComponentVariableType>::Get(rDestinationVariableDisplacement.Name() + var_ext);
                _pNode->FastGetSolutionStepValue(var_destination_disp) = displacement( c );
                c++;
            }
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::CalculateRotationMatrixWithAngle( VectorType& rAxis, double& rAngle , MatrixType& rRotationMatrix)
{
    rRotationMatrix(0, 0) = cos( rAngle ) + pow(rAxis(0), 2) * (1 - cos( rAngle ));
    rRotationMatrix(0, 1) = rAxis(0) * rAxis(1) * (1 - cos( rAngle )) - rAxis(2)*sin(rAngle);
    rRotationMatrix(0, 2) = rAxis(0)*rAxis(2)*(1-cos(rAngle)) + rAxis(1)*sin(rAngle);

    rRotationMatrix(1, 0) = rAxis(0)*rAxis(1)*(1-cos(rAngle)) + rAxis(2)*sin(rAngle);
    rRotationMatrix(1, 1) = cos(rAngle) + pow(rAxis(1), 2)*(1 - cos(rAngle));
    rRotationMatrix(1, 2) = rAxis(1)*rAxis(2)*(1-cos(rAngle)) - rAxis(0)*sin(rAngle);

    rRotationMatrix(2, 0) = rAxis(0)*rAxis(2)*(1-cos(rAngle)) - rAxis(1)*sin(rAngle);
    rRotationMatrix(2, 1) = rAxis(1)*rAxis(2)*(1-cos(rAngle)) + rAxis(0)*sin(rAngle);
    rRotationMatrix(2, 2) = cos(rAngle) + pow(rAxis(2),2)*(1-cos(rAngle));
}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::getRotationVector(const MatrixType& rotationMatrix, VectorType& rotationVector) {
    // see Non-linear Modeling and Analysis of Solids and Structures (Steen Krenk 2009) P52
    double angle = rotationMatrix(0, 0) + rotationMatrix(1, 1) + rotationMatrix(2, 2) - 1.0;
    angle /= 2.0;
    if (angle > 1.0)
        angle = 1.0;
    else if (angle < -1.0)
        angle = -1.0;

    angle = acos( angle ); // between 0 and pi

    const double EPS = 1E-6;
    const double PI = 3.1416;
    if (angle < EPS) {
        rotationVector(0) = 0.0;
        rotationVector(1) = 0.0;
        rotationVector(2) = 0.0;

        return;
    } else if ((PI - angle) < EPS) {
        const double product11 = (rotationMatrix(0, 0) + 1.0) / 2.0;
        const double product22 = (rotationMatrix(1, 1) + 1.0) / 2.0;
        const double product33 = (rotationMatrix(2, 2) + 1.0) / 2.0;
        const double product12 = (rotationMatrix(0, 1) + 1.0) / 2.0;
        const double product23 = (rotationMatrix(1, 2) + 1.0) / 2.0;
        const double product13 = (rotationMatrix(0, 2) + 1.0) / 2.0;
        const double tmp1 = sqrt(product11);
        const double tmp2 = sqrt(product22);
        const double tmp3 = sqrt(product33);

        { // case 1 +++:
            rotationVector(0) = tmp1;
            rotationVector(1) = tmp2;
            rotationVector(2) = tmp3;
            const double tmp12 = rotationVector(0) * rotationVector(1);
            const double tmp13 = rotationVector(0) * rotationVector(2);
            const double tmp23 = rotationVector(1) * rotationVector(2);
            if (fabs(tmp12) < EPS || fabs(tmp12 - product12) < fabs(tmp12 + product12))
                if (fabs(tmp13) < EPS || fabs(tmp13 - product13) < fabs(tmp13 + product13))
                    if (fabs(tmp23) < EPS || fabs(tmp23 - product23) < fabs(tmp23 + product23)) {
                        rotationVector(0) *= PI;
                        rotationVector(1) *= PI;
                        rotationVector(2) *= PI;
                        return;
                    }
        }
        { // case 2 +--:
            rotationVector(0) = tmp1;
            rotationVector(1) = -tmp2;
            rotationVector(2) = -tmp3;
            const double tmp12 = rotationVector(0) * rotationVector(1);
            const double tmp13 = rotationVector(0) * rotationVector(2);
            const double tmp23 = rotationVector(1) * rotationVector(2);
            if (fabs(tmp12) < EPS || fabs(tmp12 - product12) < fabs(tmp12 + product12))
                if (fabs(tmp13) < EPS || fabs(tmp13 - product13) < fabs(tmp13 + product13))
                    if (fabs(tmp23) < EPS || fabs(tmp23 - product23) < fabs(tmp23 + product23)) {
                        rotationVector(0) *= PI;
                        rotationVector(1) *= PI;
                        rotationVector(2) *= PI;
                        return;
                    }
        }
        { // case 3 -+-:
            rotationVector(0) = -tmp1;
            rotationVector(1) = tmp2;
            rotationVector(2) = -tmp3;
            const double tmp12 = rotationVector(0) * rotationVector(1);
            const double tmp13 = rotationVector(0) * rotationVector(2);
            const double tmp23 = rotationVector(1) * rotationVector(2);
            if (fabs(tmp12) < EPS || fabs(tmp12 - product12) < fabs(tmp12 + product12))
                if (fabs(tmp13) < EPS || fabs(tmp13 - product13) < fabs(tmp13 + product13))
                    if (fabs(tmp23) < EPS || fabs(tmp23 - product23) < fabs(tmp23 + product23)) {
                        rotationVector(0) *= PI;
                        rotationVector(1) *= PI;
                        rotationVector(2) *= PI;
                        return;
                    }
        }
        { // case 4 --+:
            rotationVector(0) = -tmp1;
            rotationVector(1) = -tmp2;
            rotationVector(2) = tmp3;
            const double tmp12 = rotationVector(0) * rotationVector(1);
            const double tmp13 = rotationVector(0) * rotationVector(2);
            const double tmp23 = rotationVector(1) * rotationVector(2);
            if (fabs(tmp12) < EPS || fabs(tmp12 - product12) < fabs(tmp12 + product12))
                if (fabs(tmp13) < EPS || fabs(tmp13 - product13) < fabs(tmp13 + product13))
                    if (fabs(tmp23) < EPS || fabs(tmp23 - product23) < fabs(tmp23 + product23)) {
                        rotationVector(0) *= PI;
                        rotationVector(1) *= PI;
                        rotationVector(2) *= PI;
                        return;
                    }
        }
        assert(false);
    }

    double tmp = angle / 2.0 / sin(angle);
    rotationVector(0) = -(rotationMatrix(1, 2) - rotationMatrix(2, 1)) * tmp;
    rotationVector(1) =  (rotationMatrix(0, 2) - rotationMatrix(2, 0)) * tmp;
    rotationVector(2) = -(rotationMatrix(0, 1) - rotationMatrix(1, 0)) * tmp;
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
    //std::cout << "Calculating shape function values" << std::endl;
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
