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
#include "tests/test_condition.h"

namespace Kratos
{
    namespace Testing
    {
    //***********************DEFAULT CONSTRUCTOR******************************************//
    //************************************************************************************//

    TestCondition::TestCondition( 
        IndexType NewId, 
        GeometryType::Pointer pGeometry, 
        const ResidualType TheResidualType
        )
        : Condition( NewId, pGeometry )
        , mResidualType( TheResidualType )
    {
        //DO NOT ADD DOFS HERE!!!
    }


    //******************************CONSTRUCTOR*******************************************//
    //************************************************************************************//

    TestCondition::TestCondition( 
        IndexType NewId, GeometryType::Pointer pGeometry, 
        PropertiesType::Pointer pProperties, 
        const ResidualType TheResidualType
        )
        : Condition( NewId, pGeometry, pProperties )
        , mResidualType( TheResidualType )
    {

    }

    //******************************COPY CONSTRUCTOR**************************************//
    //************************************************************************************//

    TestCondition::TestCondition( TestCondition const& rOther)
        :Condition(rOther)
        ,mResidualType(rOther.mResidualType)
    {

    }

    //*******************************ASSIGMENT OPERATOR***********************************//
    //************************************************************************************//

    TestCondition&  TestCondition::operator=(TestCondition const& rOther)
    {
        //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

        Condition::operator=(rOther);

        return *this;
    }

    //*********************************OPERATIONS*****************************************//
    //************************************************************************************//

    Condition::Pointer TestCondition::Create( 
        IndexType NewId, 
        NodesArrayType const& rThisNodes, 
        PropertiesType::Pointer pProperties 
        ) const
    {
        //NEEDED TO CREATE AN ELEMENT   
        return Kratos::make_shared<TestCondition>( NewId, GetGeometry().Create( rThisNodes ), pProperties, mResidualType );
    }


    //************************************CLONE*******************************************//
    //************************************************************************************//

    Condition::Pointer TestCondition::Clone( 
        IndexType NewId, 
        NodesArrayType const& rThisNodes 
        ) const
    {
        //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
        //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

        TestCondition new_condition(NewId, GetGeometry().Create( rThisNodes ), pGetProperties(), mResidualType );

        return Kratos::make_shared<TestCondition>(new_condition);
    }


    //*******************************DESTRUCTOR*******************************************//
    //************************************************************************************//

    TestCondition::~TestCondition()
    {
    }

    //************* GETTING METHODS
    //************************************************************************************//
    //************************************************************************************//

    void TestCondition::GetDofList( 
        DofsVectorType& rConditionalDofList, 
        ProcessInfo& rCurrentProcessInfo
        )
    {
        //NEEDED TO DEFINE THE DOFS OF THE ELEMENT 
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        
        rConditionalDofList.resize( 0 );

        rConditionalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_X ) );
        rConditionalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Y ) );
        if( dimension == 3 )
            rConditionalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Z ) );    
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestCondition::EquationIdVector( 
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

    void TestCondition::GetValuesVector( Vector& rValues, int Step )
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

    void TestCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
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

    void TestCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
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

    void TestCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {

        KRATOS_TRY;

        /* Calculate conditional system */

        // Compute RHS (RHS = rRightHandSideVector = Fext - Fint)
        this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

        // Compute LHS
        this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    //***********************************************************************************
    //***********************************************************************************

    void TestCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
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
        
        switch ( mResidualType )
        {
            case ResidualType::LINEAR:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rRightHandSideVector[j] -= 0.0;
                break;
            case ResidualType::NON_LINEAR:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rRightHandSideVector[j] -= 0.0;
                break;
            case ResidualType::ARC_LENGTH:
                
                typedef Table<double,double> DoubleTableType;
                DoubleTableType this_table;
                this_table.PushBack(0.0, 0.0);
                this_table.PushBack(0.01187446988973706, 0.05863918000028767);
                this_table.PushBack(0.02968617472434265, 0.130300687940745);
                this_table.PushBack(0.05089058524173028, 0.20586469208541935);
                this_table.PushBack(0.0729431721798134, 0.29446429420552855);
                this_table.PushBack(0.0916030534351145, 0.36221224757904813);
                this_table.PushBack(0.11959287531806619, 0.4338472154489751);
                this_table.PushBack(0.1407972858354538, 0.4924620670512766);
                this_table.PushBack(0.16539440203562344, 0.5510680719634014);
                this_table.PushBack(0.1874469889737066, 0.5992504641747753);
                this_table.PushBack(0.21798134011874476, 0.6513220825551012);
                this_table.PushBack(0.26293469041560646, 0.7137863502205593);
                this_table.PushBack(0.3121289228159458, 0.7592904069809234);
                this_table.PushBack(0.351145038167939, 0.7670113558326783);
                this_table.PushBack(0.39355385920271424, 0.765596991240671);
                this_table.PushBack(0.43172179813401196, 0.7420294086098201);
                this_table.PushBack(0.4724342663273962, 0.7067211622781555);
                this_table.PushBack(0.49109414758269726, 0.6675590765382461);
                this_table.PushBack(0.5190839694656489, 0.6088159479283818);
                this_table.PushBack(0.5453774385072095, 0.5578999284523933);
                this_table.PushBack(0.5776081424936388, 0.49262683665581863);
                this_table.PushBack(0.6098388464800679, 0.4521255831904042);
                this_table.PushBack(0.6412213740458016, 0.4181454462215235);
                this_table.PushBack(0.6768447837150129, 0.397192060537901);
                this_table.PushBack(0.7251908396946566, 0.40749624292126574);
                this_table.PushBack(0.7752332485156914, 0.45560564993868147);
                this_table.PushBack(0.8176420695504667, 0.5102538668329846);
                this_table.PushBack(0.8481764206955048, 0.5571103613541187);
                this_table.PushBack(0.888888888888889, 0.6352310589598726);
                this_table.PushBack(0.9151823579304497, 0.7081742311396861);
                this_table.PushBack(0.9372349448685329, 0.8176343286965619);
                this_table.PushBack(0.9491094147582697, 0.8853999754504349);
                this_table.PushBack(0.960983884648007, 0.9466467173803182);
                
                for ( unsigned int j = 0; j < dimension; ++j )
                    rRightHandSideVector[j] -= this_table.TableGetNearestValue(current_displacement[j]);
                break;
            default:
                KRATOS_ERROR << "NOT IMPLEMENTED" << std::endl;
        }
    }

    //***********************************************************************************
    //***********************************************************************************

    void TestCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the LHS
        const unsigned int system_size = dimension;

        if ( rLeftHandSideMatrix.size1() != system_size )
            rLeftHandSideMatrix.resize( system_size, system_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestCondition::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
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

    void TestCondition::CalculateDampingMatrix( 
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

    int TestCondition::Check( const ProcessInfo& rCurrentProcessInfo )
    {    
        KRATOS_TRY
        
        // Check that all required variables have been registered
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
        
        // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
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

        KRATOS_CATCH( "Problem in the Check in the TestCondition" )
    }


    //************************************************************************************//
    //************************************************************************************//

    void TestCondition::save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    void TestCondition::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
    }

    } // Namespace Testing
} // Namespace Kratos


