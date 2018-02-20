//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "includes/table.h"
#include "tests/test_element.h"

namespace Kratos
{
    namespace Testing
    {
    //***********************DEFAULT CONSTRUCTOR******************************************//
    //************************************************************************************//

    TestElement::TestElement( 
        IndexType NewId, 
        GeometryType::Pointer pGeometry, 
        const ResidualType TheResidualType
        )
        : Element( NewId, pGeometry )
        , mResidualType( TheResidualType )
    {
        //DO NOT ADD DOFS HERE!!!
    }


    //******************************CONSTRUCTOR*******************************************//
    //************************************************************************************//

    TestElement::TestElement( 
        IndexType NewId, GeometryType::Pointer pGeometry, 
        PropertiesType::Pointer pProperties, 
        const ResidualType TheResidualType
        )
        : Element( NewId, pGeometry, pProperties )
        , mResidualType( TheResidualType )
    {

    }

    //******************************COPY CONSTRUCTOR**************************************//
    //************************************************************************************//

    TestElement::TestElement( TestElement const& rOther)
        :Element(rOther)
        ,mResidualType(rOther.mResidualType)
    {

    }

    //*******************************ASSIGMENT OPERATOR***********************************//
    //************************************************************************************//

    TestElement&  TestElement::operator=(TestElement const& rOther)
    {
        //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

        Element::operator=(rOther);

        return *this;
    }

    //*********************************OPERATIONS*****************************************//
    //************************************************************************************//

    Element::Pointer TestElement::Create( 
        IndexType NewId, 
        NodesArrayType const& rThisNodes, 
        PropertiesType::Pointer pProperties 
        ) const
    {
        //NEEDED TO CREATE AN ELEMENT   
        return Kratos::make_shared<TestElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties, mResidualType );
    }


    //************************************CLONE*******************************************//
    //************************************************************************************//

    Element::Pointer TestElement::Clone( 
        IndexType NewId, 
        NodesArrayType const& rThisNodes 
        ) const
    {
        //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
        //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

        TestElement new_element(NewId, GetGeometry().Create( rThisNodes ), pGetProperties(), mResidualType );

        return Kratos::make_shared<TestElement>(new_element);
    }


    //*******************************DESTRUCTOR*******************************************//
    //************************************************************************************//

    TestElement::~TestElement()
    {
    }

    //************* GETTING METHODS
    //************************************************************************************//
    //************************************************************************************//

    void TestElement::GetDofList( 
        DofsVectorType& rElementalDofList, 
        ProcessInfo& rCurrentProcessInfo
        )
    {
        //NEEDED TO DEFINE THE DOFS OF THE ELEMENT 
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        
        rElementalDofList.resize( 0 );

        rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Y ) );
        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Z ) );    
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestElement::EquationIdVector( 
        EquationIdVectorType& rResult, 
        ProcessInfo& rCurrentProcessInfo
        )
    {
        //NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if ( rResult.size() != dimension )
            rResult.resize( dimension, false );

        rResult[0] = GetGeometry()[0].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[1] = GetGeometry()[0].GetDof( DISPLACEMENT_Y ).EquationId();
        if( dimension == 3)
            rResult[2] = GetGeometry()[0].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    //*********************************DISPLACEMENT***************************************//
    //************************************************************************************//

    void TestElement::GetValuesVector( Vector& rValues, int Step )
    {
        //GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if ( rValues.size() != dimension ) 
            rValues.resize( dimension, false );

        rValues[0] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[1] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
            rValues[2] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }


    //************************************VELOCITY****************************************//
    //************************************************************************************//

    void TestElement::GetFirstDerivativesVector( Vector& rValues, int Step )
    {
        //GIVES THE VECTOR WITH THE TIME DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT VELOCITIES)
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if ( rValues.size() != dimension ) 
            rValues.resize( dimension, false );

        rValues[0] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[1] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
            rValues[2] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Z, Step );
    }

    //*********************************ACCELERATION***************************************//
    //************************************************************************************//

    void TestElement::GetSecondDerivativesVector( Vector& rValues, int Step )
    {
        //GIVES THE VECTOR WITH THE TIME SECOND DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT ACCELERATIONS)
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if ( rValues.size() != dimension ) 
            rValues.resize( dimension, false );

        rValues[0] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[1] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
            rValues[2] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Z, Step );
    }

    //************* COMPUTING  METHODS
    //************************************************************************************//
    //************************************************************************************//

    void TestElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {

        KRATOS_TRY;

        /* Calculate elemental system */

        // Compute RHS (RHS = rRightHandSideVector = Fext - Fint)
        this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

        // Compute LHS
        this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    //***********************************************************************************
    //***********************************************************************************

    void TestElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the RHS
        const unsigned int system_size = dimension;

        if ( rRightHandSideVector.size() != system_size )
            rRightHandSideVector.resize( system_size, false );

        rRightHandSideVector = ZeroVector( system_size ); //resetting RHS

        const array_1d<double, 3 >& current_displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 >& previous_displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);
        const array_1d<double, 3 >& delta_displacement = current_displacement - previous_displacement;
        
        Matrix table_content(2, 34); 
        table_content(0, 0) = 0.0;
        table_content(1, 0) = 0.0;
        table_content(0, 1) = 0.01187446988973706;
        table_content(1, 1) = 0.05863918000028767;
        table_content(0, 2) = 0.02968617472434265;
        table_content(1, 2) = 0.130300687940745;
        table_content(0, 3) = 0.05089058524173028;
        table_content(1, 3) = 0.20586469208541935;
        table_content(0, 4) = 0.0729431721798134;
        table_content(1, 4) = 0.29446429420552855;
        table_content(0, 5) = 0.0916030534351145;
        table_content(1, 5) = 0.36221224757904813;
        table_content(0, 6) = 0.11959287531806619;
        table_content(1, 6) = 0.4338472154489751;
        table_content(0, 7) = 0.1407972858354538;
        table_content(1, 7) = 0.4924620670512766;
        table_content(0, 8) = 0.16539440203562344;
        table_content(1, 8) = 0.5510680719634014;
        table_content(0, 9) = 0.1874469889737066;
        table_content(1, 9) = 0.5992504641747753;
        table_content(0, 10) = 0.21798134011874476;
        table_content(1, 10) = 0.6513220825551012;
        table_content(0, 11) = 0.26293469041560646;
        table_content(1, 11) = 0.7137863502205593;
        table_content(0, 12) = 0.3121289228159458;
        table_content(1, 12) = 0.7592904069809234;
        table_content(0, 13) = 0.351145038167939;
        table_content(1, 13) = 0.7670113558326783;
        table_content(0, 14) = 0.39355385920271424;
        table_content(1, 14) = 0.765596991240671;
        table_content(0, 15) = 0.43172179813401196;
        table_content(1, 15) = 0.7420294086098201;
        table_content(0, 16) = 0.4724342663273962;
        table_content(1, 16) = 0.7067211622781555;
        table_content(0, 17) = 0.49109414758269726;
        table_content(1, 17) = 0.6675590765382461;
        table_content(0, 18) = 0.5190839694656489;
        table_content(1, 18) = 0.6088159479283818;
        table_content(0, 19) = 0.5453774385072095;
        table_content(1, 19) = 0.5578999284523933;
        table_content(0, 20) = 0.5776081424936388;
        table_content(1, 20) = 0.49262683665581863;
        table_content(0, 21) = 0.6098388464800679;
        table_content(1, 21) = 0.4521255831904042;
        table_content(0, 22) = 0.6412213740458016;
        table_content(1, 22) = 0.4181454462215235;
        table_content(0, 23) = 0.6768447837150129;
        table_content(1, 23) = 0.397192060537901;
        table_content(0, 24) = 0.7251908396946566;
        table_content(1, 24) = 0.40749624292126574;
        table_content(0, 25) = 0.7752332485156914;
        table_content(1, 25) = 0.45560564993868147;
        table_content(0, 26) = 0.8176420695504667;
        table_content(1, 26) = 0.5102538668329846;
        table_content(0, 27) = 0.8481764206955048;
        table_content(1, 27) = 0.5571103613541187;
        table_content(0, 28) = 0.888888888888889;
        table_content(1, 28) = 0.6352310589598726;
        table_content(0, 29) = 0.9151823579304497;
        table_content(1, 29) = 0.7081742311396861;
        table_content(0, 30) = 0.9372349448685329;
        table_content(1, 30) = 0.8176343286965619;
        table_content(0, 31) = 0.9491094147582697;
        table_content(1, 31) = 0.8853999754504349;
        table_content(0, 32) = 0.960983884648007;
        table_content(1, 32) = 0.9466467173803182;
        table_content(0, 33) = 1.0;
        table_content(1, 33) = 1.0;
        Table<double, double> this_table(table_content);
        
        switch ( mResidualType )
        {
            case ResidualType::LINEAR:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rRightHandSideVector[j] -= delta_displacement[j] - 1.0;
                break;
            case ResidualType::NON_LINEAR:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rRightHandSideVector[j] -= std::pow(delta_displacement[j], 2) - 1.0;
                break;
            case ResidualType::ARC_LENGTH:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rRightHandSideVector[j] -= this_table.GetValue(current_displacement[j]) - 2.0;
                break;
            default:
                KRATOS_ERROR << "NOT IMPLEMENTED" << std::endl;
        }
    }

    //***********************************************************************************
    //***********************************************************************************

    void TestElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the LHS
        const unsigned int system_size = dimension;

        if ( rLeftHandSideMatrix.size1() != system_size )
            rLeftHandSideMatrix.resize( system_size, system_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

        const array_1d<double, 3 >& delta_displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);
        
        switch ( mResidualType )
        {
            case ResidualType::LINEAR:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rLeftHandSideMatrix(j, j) += 1.0;
                break;
            case ResidualType::NON_LINEAR:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rLeftHandSideMatrix(j, j) += delta_displacement[j] * 2;
                break;
            case ResidualType::ARC_LENGTH:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rLeftHandSideMatrix(j, j) += delta_displacement[j] * 2;
                break;
            default:
                KRATOS_ERROR << "NOT IMPLEMENTED" << std::endl;
        }
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int system_size = dimension;

        if ( rMassMatrix.size1() != system_size )
            rMassMatrix.resize( system_size, system_size, false );

        rMassMatrix = ZeroMatrix( system_size, system_size );

        KRATOS_CATCH( "" );
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestElement::CalculateDampingMatrix( 
        MatrixType& rDampingMatrix, 
        ProcessInfo& rCurrentProcessInfo 
        )
    {
        KRATOS_TRY;

        //0.-Initialize the DampingMatrix:
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the LHS
        const unsigned int system_size = dimension;

        rDampingMatrix = ZeroMatrix( system_size, system_size );

        KRATOS_CATCH( "" );
    }

    //************************************************************************************//
    //************************************************************************************//

    int TestElement::Check( const ProcessInfo& rCurrentProcessInfo )
    {    
        KRATOS_TRY
        
        // Check that all required variables have been registered
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
        
        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for ( std::size_t i = 0; i < this->GetGeometry().size(); ++i ) {
            Node<3>& rnode = this->GetGeometry()[i];
            
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,rnode)
        }
        
        return 0;

        KRATOS_CATCH( "Problem in the Check in the TestElement" )
    }


    //************************************************************************************//
    //************************************************************************************//

    void TestElement::save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void TestElement::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    } // Namespace Testing
} // Namespace Kratos


