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

    // Calculating and storing the shape function values
    // Linear shape functions
    pairing_index = ProjectionUtilities::ProjectOnLine(*p_geom, point_to_proj, mLocalCoordTol, linear_shape_function_values, eq_ids, proj_dist, ComputeApproximation);
    // Hermitian shape functions
    ProjectionUtilities::ProjectOnLineHermitian(*p_geom, point_to_proj, mLocalCoordTol, hermitian_shape_function_values, hermitian_shape_function_derivatives_values, proj_dist, projection_point);
    const bool is_full_projection = (pairing_index == ProjectionUtilities::PairingIndex::Line_Inside);
    

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
        
        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());

        std::cout << "\n------------- END OF SEARCH -----------------" << std::endl;
    }
    
}

void BeamMapperInterfaceInfo::ComputeRotationMatrix()
{
    std::vector<double> axisX;
    std::vector<double> axisY;
    std::vector<double> axisZ;

    axisX.resize(3);
    axisY.resize(3);
    axisZ.resize(3);
    
    const auto p_geom = mpInterfaceObject->pGetBaseGeometry();

    auto temp_v = (*p_geom)[1].Coordinates() - (*p_geom)[0].Coordinates();
    double lengthX = sqrt(temp_v[0]*temp_v[0] + temp_v[1]*temp_v[1] + temp_v[2]*temp_v[2]);
    
    KRATOS_ERROR_IF(lengthX < 0.000001) << "Lenght of the beam is 0.0" << std::endl;
    
    axisX[0] = temp_v[0] / lengthX;
    axisX[1] = temp_v[1] / lengthX;
    axisX[2] = temp_v[2] / lengthX;   

    if (axisX[0] == 1.0 && axisX[1] == 0.0 && axisX[2] == 0.0 ){
        axisY[0] = 0.0;
        axisY[1] = 1.0;
        axisY[2] = 0.0;
        axisZ[0] = 0.0;
        axisZ[1] = 0.0;
        axisZ[2] = 1.0;
    }
    else if (axisX[0] == 0.0 && axisX[1] == 1.0 && axisX[2] == 0.0 ){
        axisY[0] = 0.0;
        axisY[1] = 0.0;
        axisY[2] = 1.0;
        axisZ[0] = 1.0;
        axisZ[1] = 0.0;
        axisZ[2] = 0.0;
    }
    else if (axisX[0] == 0.0 && axisX[1] == 0.0 && axisX[2] == 1.0 ){
        axisY[0] = 0.0;
        axisY[1] = 1.0;
        axisY[2] = 0.0;
        axisZ[0] = 1.0;
        axisZ[1] = 0.0;
        axisZ[2] = 0.0;
    }
    else if (axisX[0] != 0.0 && axisX[1] != 0.0 && axisX[2] == 0.0 ){
        axisY[0] = -axisX[1];
        axisY[1] =  axisX[0];
        axisY[2] =  0.0;
        axisZ[0] = axisX[1]*axisY[2] - axisX[2]*axisY[1]; 
        axisZ[1] = axisX[2]*axisY[0] - axisX[0]*axisY[2];
        axisZ[2] = axisX[0]*axisY[1] - axisX[1]*axisY[0];
    }
    else if (axisX[0] != 0.0 && axisX[1] == 0.0 && axisX[2] != 0.0 ){
        axisY[0] = -axisX[2];
        axisY[1] =  0;
        axisY[2] =  axisX[0];
        axisZ[0] = axisX[1]*axisY[2] - axisX[2]*axisY[1]; 
        axisZ[1] = axisX[2]*axisY[0] - axisX[0]*axisY[2];
        axisZ[2] = axisX[0]*axisY[1] - axisX[1]*axisY[0];
    }
    else if (axisX[0] == 0.0 && axisX[1] != 0.0 && axisX[2] != 0.0){
        axisY[0] =  0;
        axisY[1] = -axisX[2];
        axisY[2] =  axisX[1];
        axisZ[0] = axisX[1]*axisY[2] - axisX[2]*axisY[1]; 
        axisZ[1] = axisX[2]*axisY[0] - axisX[0]*axisY[2];
        axisZ[2] = axisX[0]*axisY[1] - axisX[1]*axisY[0];
    }
    else{
        axisY[0] = 1;
        axisY[1] = 1;
        axisY[2] = (-axisX[0] - axisX[1]) / axisX[2];
        double lenghtY = sqrt(axisY[0]*axisY[0] + axisY[1]*axisY[1] + axisY[2]*axisY[2]);
        axisY[0] = axisY[0]/lenghtY;
        axisY[1] = axisY[1]/lenghtY;
        axisY[2] = axisY[2]/lenghtY;
        
        axisZ[0] = axisX[1]*axisY[2] - axisX[2]*axisY[1]; 
        axisZ[1] = axisX[2]*axisY[0] - axisX[0]*axisY[2];
        axisZ[2] = axisX[0]*axisY[1] - axisX[1]*axisY[0];
    }

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
{ 
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

        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

        if (rDestinationIds.size() != 1) rDestinationIds.resize(1);
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
   
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
        std::cout << " ------------------------------- "<< std::endl;
        std::cout << "----------- LOCAL SYSTEM " << i  << "------------"<< std::endl;
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
    {   std::cout << " ---------------------------------- "<< std::endl;
        std::cout << " ----------- LOCAL SYSTEM " << i << "---------"<< std::endl;
        std::cout << " ---------------------------------- "<< std::endl;
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
            std::cout << "rotation of node 1 is" << rotationNode1_G << std::endl;
            std::cout << "displacement of node 2 is" << displacementNode2_G << std::endl;
            std::cout << "rotation of node 2 is" << rotationNode2_G << std::endl;

            // Processing previously the nodal rotations // THIS IS JUST FOR THE TESTS
            MatrixType Rx(3, 3), Ry(3, 3), Rz(3, 3), R(3, 3), R_temp(3, 3);
            VectorType e_1(3, 0.0), e_2(3, 0.0), e_3(3, 0.0);
            e_1(0) = 1.0;
            e_2(1) = 1.0;
            e_3(2) = 1.0;
            // For rotations in node 1
            CalculateRotationMatrixWithAngle(e_1, rotationNode1_G(0), Rx);
            CalculateRotationMatrixWithAngle(e_2, rotationNode1_G(1), Ry);
            CalculateRotationMatrixWithAngle(e_3, rotationNode1_G(2), Rz);
            //R_temp = prod(Rx, Ry); // this works for beam
            //R = prod(Rz, R_temp);
            R_temp = prod(Ry, Rz); // this works for wing
            R = prod(Rx, R_temp);
            getRotationVector(R, rotationNode1_G);
            double _angle1 = norm_2(rotationNode1_G);
            if (_angle1 != 0.0){
                rotationNode1_G /= _angle1;
            }

            // For rotations in node 2
            CalculateRotationMatrixWithAngle(e_1, rotationNode2_G(0), Rx);
            CalculateRotationMatrixWithAngle(e_2, rotationNode2_G(1), Ry);
            CalculateRotationMatrixWithAngle(e_3, rotationNode2_G(2), Rz);
            //R_temp = prod(Rx, Ry); // this works for beam
            //R = prod(Rz, R_temp);
            R_temp = prod(Ry, Rz); // this works for wing
            R = prod(Rx, R_temp);
            std::cout << "RzRxRy in 2 is" << R << std::endl;
            VectorType temp_rotationNode2_G(3);
            temp_rotationNode2_G = rotationNode2_G;
            getRotationVector(R, rotationNode2_G);
            std::cout << "rotation vector of RzRxRy in 2 is : " << rotationNode2_G << std::endl;
            double _angle2 = norm_2(rotationNode2_G);
            if (_angle2 != 0.0){
                rotationNode2_G /= _angle2;
            }

            MatrixType _rotationMatrix_B_G( 3, 3 );
            double determinant;
            //std::cout << "_rotationMatrix_G_B : " << _rotationMatrix_G_B << std::endl;
            MathUtils<double>::InvertMatrix3(_rotationMatrix_G_B, _rotationMatrix_B_G, determinant );
            //std::cout << "_rotationMatrix_B_G : " << _rotationMatrix_B_G << std::endl;
            //std::cout << "displacementNode1_G : " << displacementNode1_G << std::endl;
            //std::cout << "displacementNode2_G : " << displacementNode2_G << std::endl;
            // Transforming the displacements to the BCS
            TDenseSpace::Mult( _rotationMatrix_B_G, displacementNode1_G, displacementNode1_B );
            TDenseSpace::Mult( _rotationMatrix_B_G, displacementNode2_G, displacementNode2_B );
            

            // other option to calculate rotation vectors
            MatrixType Rot2(3,3), Rot2fin_temp(3,3), Rot2fin(3,3);
            double angle_unitary2 = norm_2(temp_rotationNode2_G);
            std::cout << "angle_unitary2 = " << angle_unitary2 << std::endl;
            //if(angle_unitary2 > 3.1415){
            //    angle_unitary2 = angle_unitary2 - 3.1415;
            //}
            //std::cout << "angle_unitary2 = " << angle_unitary2 << std::endl;
            VectorType unitary_rotation2(3);
            unitary_rotation2 = rotationNode2_G / norm_2(rotationNode2_G);
            std::cout << "unitary_rotation2 = " << unitary_rotation2 << std::endl;
            CalculateRotationMatrixWithAngle(unitary_rotation2, angle_unitary2, Rot2);
            Rot2fin_temp = prod(_rotationMatrix_G_B, Rot2);
            Rot2fin = prod(Rot2fin_temp, _rotationMatrix_B_G);
            std::cout << "otra opcion para calcular RzRxRy es : " << Rot2fin << std::endl;
            /////

            // Transforming the nodal rotations to the BCS
            MatrixType Rotation_G_1(3, 3);
            CalculateRotationMatrixWithAngle(rotationNode1_G, _angle1, Rotation_G_1);
            MatrixType Rotation_B_1(3, 3);
            MathUtils<double>::BtDBProductOperation(Rotation_B_1, Rotation_G_1,_rotationMatrix_G_B);
            std::cout << "Rotation_G_1 = " << Rotation_G_1 << std::endl;
            std::cout << "Rotation_B_1 = " << Rotation_B_1 << std::endl;
            
            MatrixType Rotation_G_2(3, 3);
            CalculateRotationMatrixWithAngle(rotationNode2_G, _angle2, Rotation_G_2);
            MatrixType Rotation_B_2(3, 3);
            MathUtils<double>::BtDBProductOperation(Rotation_B_2, Rotation_G_2,_rotationMatrix_G_B);
            std::cout << "Rotation_G_2 = " << Rotation_G_2 << std::endl;
            std::cout << "Rotation_B_2 = " << Rotation_B_2 << std::endl;

            // Calculating R_d  and t_d for Phi_d 
            VectorType e_x_d_G(3), e_x_d_B(3);
            e_x_d_G = _r_geom[1].Coordinates() + displacementNode2_G - (_r_geom[0].Coordinates() + displacementNode1_G); // this vector is described in global system 
            e_x_d_G /= norm_2(e_x_d_G);
            TDenseSpace::Mult(_rotationMatrix_B_G, e_x_d_G, e_x_d_B); // transforming e_x_d to the beam coordinate system
            e_x_d_B /= norm_2(e_x_d_B);
            
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
            I_2nT = _I - 2 * MathUtils<double>::TensorProduct3(n_d, n_d);
            R_d = prod(I_2nT, _R);
            std::cout << "Rd is = " << R_d << std::endl;
           
            VectorType t_d(3);
            t_d(0) = _linearShapeValues(0) * displacementNode1_B(0) + _linearShapeValues(1) * displacementNode2_B(0);
            t_d(1) = _linearShapeValues(0) * displacementNode1_B(1) + _linearShapeValues(1) * displacementNode2_B(1);
            t_d(2) = _linearShapeValues(0) * displacementNode1_B(2) + _linearShapeValues(1) * displacementNode2_B(2);
            std::cout << "td is = " << t_d << std::endl;

            // Calculating R_s
            MatrixType Rs_Rl1(3, 3), Rs_Rl2(3, 3), R_d_T(3, 3);
            MathUtils<double>::InvertMatrix3(R_d, R_d_T, determinant);
            Rs_Rl1 = prod(R_d_T, Rotation_B_1); 
            Rs_Rl2 = prod(R_d_T, Rotation_B_2);
            VectorType v_Rs_Rl1(3), v_Rs_Rl2(3);
            getRotationVector(Rs_Rl1, v_Rs_Rl1);
            getRotationVector(Rs_Rl2, v_Rs_Rl2);


            double theta_s = (v_Rs_Rl1(0) + v_Rs_Rl2(0)) / 2;
            MatrixType R_s(3, 3);
            VectorType e_x_s(3, 0.0);
            e_x_s(0) = 1.0;
            CalculateRotationMatrixWithAngle(e_x_s, theta_s, R_s);
            std::cout << "R_s = " << R_s << std::endl;

            // Calculating rotation vector l on node 1 and 2 of the beam
            MatrixType Rl1(3, 3), Rl2(3, 3), R_s_T(3, 3);
            VectorType v_Rl1(3), v_Rl2(3);
            MathUtils<double>::InvertMatrix3(R_s, R_s_T, determinant);
            Rl1 = prod(R_s_T, Rs_Rl1);
            Rl2 = prod(R_s_T, Rs_Rl2);
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
            _ShapeFunctionsMatrix(4 , 10) = _hermitanDerShapeValues(3);
            _ShapeFunctionsMatrix(5 , 1) = _hermitanDerShapeValues(0) / length_beamVector;
            _ShapeFunctionsMatrix(5 , 5) = _hermitanDerShapeValues(1);
            _ShapeFunctionsMatrix(5 , 7) = _hermitanDerShapeValues(2) / length_beamVector;
            _ShapeFunctionsMatrix(5 , 11) = _hermitanDerShapeValues(3); 

            VectorType _DOF_Vector_l(12);
            for (size_t i = 0; i < 3; i++){
                _DOF_Vector_l(i) = 0.0;
                _DOF_Vector_l(i + 3) = v_Rl1(i);
                _DOF_Vector_l(i + 6) = 0.0;
                _DOF_Vector_l(i + 9) = v_Rl2(i);
            }

            VectorType _DisplacementsRotations_L(6);
            TDenseSpace::Mult(_ShapeFunctionsMatrix, _DOF_Vector_l, _DisplacementsRotations_L);

            // Constructing phi_l
            VectorType axis_l(3), t_l(3);
            MatrixType R_l(3, 3), R_l_temp(3, 3), R_x(3, 3), R_y(3, 3), R_z(3, 3);
            t_l(0) = _DisplacementsRotations_L(0); // t_l
            t_l(1) = _DisplacementsRotations_L(1); // t_l
            t_l(2) = _DisplacementsRotations_L(2); // t_l
            axis_l(0) = _DisplacementsRotations_L(3);
            axis_l(1) = _DisplacementsRotations_L(4);
            axis_l(2) = _DisplacementsRotations_L(5); 

            std::cout << "axis_l : " << axis_l << std::endl;

            //axis_l = axis_l/angle_l;
            CalculateRotationMatrixWithAngle(e_1, axis_l(0), Rx);
            CalculateRotationMatrixWithAngle(e_2, axis_l(1), Ry);
            CalculateRotationMatrixWithAngle(e_3, axis_l(2), Rz);
            //R_temp = prod(Rx, Ry); // this works for beam
            //R_l = prod(Rz, R_temp);
            R_temp = prod(Ry, Rz); // this works for wing
            R_l = prod(Rx, R_temp);

            std::cout << "t_d = " << t_d << std::endl;
            std::cout << "t_s = none " << std::endl;
            std::cout << "t_l = " << t_l << std::endl;

            // Calculating R_P and t_P
            MatrixType R_P(3, 3), R_P_temp(3, 3);
            std::cout << "R_d = " << R_d << std::endl;
            std::cout << "R_s = " << R_s << std::endl;
            std::cout << "R_l = " << R_l << std::endl;
            R_P_temp = prod(R_d, R_s);
            R_P = prod(R_P_temp, R_l);

            VectorType t_P(3), t_P_temp(3, 0.0);
            TDenseSpace::Mult(R_P_temp, t_l, t_P_temp);
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

    angle = acos(angle); // between 0 and pi

    const double EPS = 1E-6;
    if (angle < EPS) {
        rotationVector(0) = 0.0;
        rotationVector(1) = 0.0;
        rotationVector(2) = 0.0;

        return;
    } else if ((M_PI - angle) < EPS) {
        const double product11 = (rotationMatrix(0,0) + 1.0) / 2.0;
        const double product22 = (rotationMatrix(1,1) + 1.0) / 2.0;
        const double product33 = (rotationMatrix(2,2) + 1.0) / 2.0;
        const double product12 = (rotationMatrix(0,1) + 1.0) / 2.0;
        const double product23 = (rotationMatrix(1,2) + 1.0) / 2.0;
        const double product13 = (rotationMatrix(0,2) + 1.0) / 2.0;
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
                        rotationVector(0) *= M_PI;
                        rotationVector(1) *= M_PI;
                        rotationVector(2) *= M_PI;
                        return;
                    }
        }
        { // case 2 +--:
            rotationVector(0) = tmp1;
            rotationVector(1) = -tmp2;
            rotationVector(2) = -tmp3;
            const double tmp12 = rotationVector[0] * rotationVector[1];
            const double tmp13 = rotationVector[0] * rotationVector[2];
            const double tmp23 = rotationVector[1] * rotationVector[2];
            if (fabs(tmp12) < EPS || fabs(tmp12 - product12) < fabs(tmp12 + product12))
                if (fabs(tmp13) < EPS || fabs(tmp13 - product13) < fabs(tmp13 + product13))
                    if (fabs(tmp23) < EPS || fabs(tmp23 - product23) < fabs(tmp23 + product23)) {
                        rotationVector(0) *= M_PI;
                        rotationVector(1) *= M_PI;
                        rotationVector(2) *= M_PI;
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
                        rotationVector(0) *= M_PI;
                        rotationVector(1) *= M_PI;
                        rotationVector(2) *= M_PI;
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
                        rotationVector(0) *= M_PI;
                        rotationVector(1) *= M_PI;
                        rotationVector(2) *= M_PI;
                        return;
                    }
        }
        assert(false);
    }

    double tmp = angle / 2.0 / sin(angle);
    rotationVector(0) = -(rotationMatrix(1,2) - rotationMatrix(2,1)) * tmp;
    rotationVector(1) =  (rotationMatrix(0,2) - rotationMatrix(2,0)) * tmp;
    rotationVector(2) = -(rotationMatrix(0,1) - rotationMatrix(1,0)) * tmp;
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
}

// Extern template instantiation

template class BeamMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;



}  // namespace Kratos.
