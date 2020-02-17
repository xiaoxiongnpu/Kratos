//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson (thanks Klaus Sautter)
//
//

#if !defined(KRATOS_MPM_EXPLICIT_CENTRAL_DIFFERENCE_SCHEME)
#define KRATOS_MPM_EXPLICIT_CENTRAL_DIFFERENCE_SCHEME

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "containers/array_1d.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_utilities/mpm_boundary_rotation_utility.h"
#include "custom_utilities/mpm_explicit_utilities.h"

namespace Kratos{

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
 * @class MPMExplicitCentralDifferenceScheme
 * @ingroup KratosParticle
 * @brief A MPM explicit scheme
 * @details Scheme options include Forward Euler or Central Difference. 
 * Stress update options include Update Stress First (USF), Update Stress Last (USL) and Modified Update Stress Last (MUSL).
 * 
 * @author Peter
 */
template <class TSparseSpace,
          class TDenseSpace //= DenseSpace<double>
          >
class MPMExplicitCentralDifferenceScheme
    : public MPMExplicitScheme<TSparseSpace, TDenseSpace> {

public:
    /**@name Type Definitions */

    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION(MPMExplicitCentralDifferenceScheme);

    typedef MPMExplicitScheme<TSparseSpace, TDenseSpace>                      BaseType;

    typedef typename BaseType::TDataType                         TDataType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename Element::DofsVectorType                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

    typedef typename BaseType::Pointer                     BaseTypePointer;

    /// The arrays of elements and nodes
    typedef ModelPart::NodesContainerType NodesArrayType;

    /// Definition for the node iterator
    typedef typename ModelPart::NodeIterator NodeIterator;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the numerical limit
    static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{
    void InitializeExplicitScheme (
        ModelPart& rModelPart,
        const SizeType DomainSize = 3
    )override
    {
        KRATOS_TRY

        /// The array of ndoes
        NodesArrayType& r_nodes = rModelPart.Nodes();

        // The first iterator of the array of nodes
        const auto it_node_begin = rModelPart.NodesBegin();

        if (mIsCentralDifference) // TODO can this be cleaned up into one big loop?
        {
            /// Initialise the database of the nodes
            const array_1d<double, 3> zero_array = ZeroVector(3);
            #pragma omp parallel for schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                auto it_node = (it_node_begin + i);
                //it_node->SetValue(NODAL_MASS, 0.0);
                array_1d<double, 3>& r_middle_velocity = it_node->FastGetSolutionStepValue(MIDDLE_VELOCITY);
                r_middle_velocity = ZeroVector(3);
            }
        }


        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            auto it_node = (it_node_begin + i);


            array_1d<double, 3>& r_current_residual = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (IndexType j = 0; j < DomainSize; j++) {
                r_current_residual[j] = 0.0;
            }

            if (mIsCentralDifference)
            {
                array_1d<double, 3>& r_middle_velocity = it_node->FastGetSolutionStepValue(MIDDLE_VELOCITY);
                const array_1d<double, 3>& r_current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
                for (IndexType j = 0; j < DomainSize; j++) {
                    r_middle_velocity[j] = r_current_velocity[j];
                }
            }
        }

        KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************
    void UpdateTranslationalDegreesOfFreedom (
        NodeIterator itCurrentNode,
        const IndexType DisplacementPosition,
        const SizeType DomainSize = 3
    )override
    {

        //PJW
        const double node_X = itCurrentNode->X();
        const double node_Y = itCurrentNode->Y();
        // PJW

//         const double nodal_mass = itCurrentNode->FastGetSolutionStepValue(NODAL_MASS);
         //

//         const double nodal_displacement_damping = itCurrentNode->GetValue(NODAL_DISPLACEMENT_DAMPING);
//         const array_1d<double, 3>& r_current_residual = itCurrentNode->FastGetSolutionStepValue(FORCE_RESIDUAL);

//         array_1d<double, 3>& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY);
//         array_1d<double, 3>& r_current_displacement = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT);
//         array_1d<double, 3>& r_middle_velocity = itCurrentNode->FastGetSolutionStepValue(MIDDLE_VELOCITY);

//         array_1d<double, 3>& r_current_acceleration = itCurrentNode->FastGetSolutionStepValue(ACCELERATION);

//         const array_1d<double, 3>& r_previous_displacement = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
//         const array_1d<double, 3>& r_previous_middle_velocity = itCurrentNode->FastGetSolutionStepValue(MIDDLE_VELOCITY, 1);
//         
         //// Solution of the explicit equation:
//         if (nodal_mass > numerical_limit)
//             // I do this on element lvl
//             //noalias(r_current_acceleration) = (r_current_residual - nodal_displacement_damping * r_current_velocity) / nodal_mass;
//             noalias(r_current_acceleration) = (r_current_residual) / nodal_mass;
//         else
//             noalias(r_current_acceleration) = ZeroVector(3);


//         std::array<bool, 3> fix_displacements = {false, false, false};

//         fix_displacements[0] = (itCurrentNode->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
//         fix_displacements[1] = (itCurrentNode->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
//         if (DomainSize == 3)
//             fix_displacements[2] = (itCurrentNode->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());


//         for (IndexType j = 0; j < DomainSize; j++) {
//             if (fix_displacements[j]) {
//                 r_current_acceleration[j] = 0.0;
//                 r_middle_velocity[j] = 0.0;
//             }

//             r_current_velocity[j] =  r_previous_middle_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * r_current_acceleration[j]; //+ actual_velocity;
//             r_middle_velocity[j] = r_current_velocity[j] + (mTime.Middle - mTime.Previous) * r_current_acceleration[j];
//             r_current_displacement[j] = r_previous_displacement[j] + mTime.Delta * r_middle_velocity[j];

         //	//r_current_velocity[j] = r_middle_velocity[j]; //PJW TESTING

         //	

//         } // for DomainSize



         //PJW integrated momentum form of explicit advance =============================
        std::array<bool, 3> fix_displacements = { false, false, false };
        fix_displacements[0] = (itCurrentNode->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
        fix_displacements[1] = (itCurrentNode->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
        if (DomainSize == 3)
            fix_displacements[2] = (itCurrentNode->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

        array_1d<double, 3>& r_nodal_momenta = itCurrentNode->FastGetSolutionStepValue(NODAL_MOMENTUM);
        array_1d<double, 3>& r_current_residual = itCurrentNode->FastGetSolutionStepValue(FORCE_RESIDUAL);

        //PJW, simple coefficient for central difference stepping
        double alpha = 1.0;
        if (mTime.Previous == 0.0)
        {
            alpha = 0.5;
        }

        // Advance momenta
        for (IndexType j = 0; j < DomainSize; j++) {
            if (fix_displacements[j]) {
                r_nodal_momenta[j] = 0.0;
                r_current_residual[j] = 0.0;
            }
            r_nodal_momenta[j] += alpha * mTime.Delta * r_current_residual[j];

        } // for DomainSize


        // We need to set updated grid velocity here if we are using USL formulation
        if (mStressUpdateOption == 1)
        {
            array_1d<double, 3>& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY);
            r_current_velocity.clear();
            const double nodal_mass = itCurrentNode->FastGetSolutionStepValue(NODAL_MASS);
            if (nodal_mass > numerical_limit)
            {
                for (IndexType j = 0; j < DomainSize; j++)
                {
                    r_current_velocity[j] = r_nodal_momenta[j] / nodal_mass;
                } // for DomainSize
            }
        }
        //PJW integrated momentum form of explicit advance =============================
    }

    
    //***************************************************************************
    //***************************************************************************

    /**
     * initializes time step solution
     * only for reasons if the time step solution is restarted
     * @param r_model_part
     * @param A	LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

            std::cout << "\n\n SCHEME INITIALIZE SS \n\n" << std::endl;

            ProcessInfo CurrentProcessInfo = r_model_part.GetProcessInfo();
        BaseType::InitializeSolutionStep(r_model_part, A, Dx, b);
        // LOOP OVER THE GRID NODES PERFORMED FOR CLEAR ALL NODAL INFORMATION
        #pragma omp parallel for
        for (int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
        {
            auto i = mr_grid_model_part.NodesBegin() + iter;

            // Variables to be cleaned
            double& nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
            double& nodal_density = (i)->FastGetSolutionStepValue(DENSITY);
            array_1d<double, 3 >& nodal_momentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
            array_1d<double, 3 >& nodal_inertia = (i)->FastGetSolutionStepValue(NODAL_INERTIA);
            array_1d<double, 3 >& nodal_force = (i)->FastGetSolutionStepValue(FORCE_RESIDUAL); //PJW
            array_1d<double, 3 >& nodal_displacement = (i)->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 >& nodal_velocity = (i)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 >& nodal_acceleration = (i)->FastGetSolutionStepValue(ACCELERATION); //PJW


            double& nodal_old_pressure = (i)->FastGetSolutionStepValue(PRESSURE, 1);
            double& nodal_pressure = (i)->FastGetSolutionStepValue(PRESSURE);
            if (i->SolutionStepsDataHas(NODAL_MPRESSURE)) {
                double& nodal_mpressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);
                nodal_mpressure = 0.0;
            }

            // Clear
            nodal_mass = 0.0;
            nodal_density = 0.0;
            nodal_momentum.clear();
            nodal_inertia.clear();
            nodal_force.clear(); //PJW

            nodal_displacement.clear();
            nodal_velocity.clear();
            nodal_acceleration.clear();
            nodal_old_pressure = 0.0;
            nodal_pressure = 0.0;

            if (mIsCentralDifference)
            {
                array_1d<double, 3 >& nodal_middle_velocity = (i)->FastGetSolutionStepValue(MIDDLE_VELOCITY); //PJW
                nodal_middle_velocity.clear(); //PJW
            }
        }

        // Extrapolate from Material Point Elements and Conditions
        Scheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);

        KRATOS_CATCH("")
    }


    

    //***************************************************************************
    //***************************************************************************

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    int Check(ModelPart& r_model_part) override
    {
        KRATOS_TRY

            int err = Scheme<TSparseSpace, TDenseSpace>::Check(r_model_part);
        if (err != 0) return err;

        //check that the variables are correctly initialized
        KRATOS_ERROR_IF(DISPLACEMENT.Key() == 0) << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        KRATOS_ERROR_IF(VELOCITY.Key() == 0) << "VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
        KRATOS_ERROR_IF(ACCELERATION.Key() == 0) << "ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;

        //check that variables are correctly allocated
        // TODO add middle velocity
        for (ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin();
            it != r_model_part.NodesEnd(); it++)
        {
            KRATOS_ERROR_IF(it->SolutionStepsDataHas(DISPLACEMENT) == false) << "DISPLACEMENT variable is not allocated for node " << it->Id() << std::endl;
            KRATOS_ERROR_IF(it->SolutionStepsDataHas(VELOCITY) == false) << "VELOCITY variable is not allocated for node " << it->Id() << std::endl;
            KRATOS_ERROR_IF(it->SolutionStepsDataHas(ACCELERATION) == false) << "ACCELERATION variable is not allocated for node " << it->Id() << std::endl;
        }

        //check that dofs exist
        for (ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin();
            it != r_model_part.NodesEnd(); it++)
        {
            KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_X) == false) << "Missing DISPLACEMENT_X dof on node " << it->Id() << std::endl;
            KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_Y) == false) << "Missing DISPLACEMENT_Y dof on node " << it->Id() << std::endl;
            KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_Z) == false) << "Missing DISPLACEMENT_Z dof on node " << it->Id() << std::endl;
        }

        //check for minimum value of the buffer index
        KRATOS_ERROR_IF(r_model_part.GetBufferSize() < 2) << "Insufficient buffer size. Buffer size should be greater than 2. Current size is" << r_model_part.GetBufferSize() << std::endl;

        return 0;
        KRATOS_CATCH("")
    }

    virtual void SchemeCustomInitialization(
        ModelPart& rModelPart,
        const SizeType DomainSize = 3
    ) override
    {
        KRATOS_TRY

            // The array containing the nodes
            NodesArrayType& r_nodes = rModelPart.Nodes();

        // The fisrt node interator
        const auto it_node_begin = rModelPart.NodesBegin();

        // Auxiliar zero array
        const array_1d<double, 3> zero_array = ZeroVector(3);

        // Getting dof position
        const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            // Current step information "N+1" (before step update).
            auto it_node = it_node_begin + i;

            const double nodal_mass = it_node->GetValue(NODAL_MASS);

            array_1d<double, 3>& r_middle_velocity;
            if (mIsCentralDifference)
            {
                r_middle_velocity = it_node->FastGetSolutionStepValue(MIDDLE_VELOCITY);
            }

            const array_1d<double, 3>& r_current_residual = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);

            array_1d<double, 3>& r_current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
            //             array_1d<double,3>& r_current_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3>& r_current_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION);

            // Solution of the explicit equation:
            if (nodal_mass > numerical_limit) {
                r_current_acceleration = r_current_residual / nodal_mass;
            }
            else {
                r_current_acceleration = zero_array;
            }

            std::array<bool, 3> fix_displacements = { false, false, false };

            fix_displacements[0] = (it_node->GetDof(DISPLACEMENT_X, disppos).IsFixed());
            fix_displacements[1] = (it_node->GetDof(DISPLACEMENT_Y, disppos + 1).IsFixed());
            if (DomainSize == 3)
                fix_displacements[2] = (it_node->GetDof(DISPLACEMENT_Z, disppos + 2).IsFixed());

            for (IndexType j = 0; j < DomainSize; j++) {
                if (fix_displacements[j]) {
                    r_current_acceleration[j] = 0.0;

                    if (mIsCentralDifference)
                    {
                        r_middle_velocity[j] = 0.0;
                    }
                    
                }

                if (mIsCentralDifference)
                {
                    r_middle_velocity[j] = 0.0 + (mTime.Middle - mTime.Previous) * r_current_acceleration[j];
                    r_current_velocity[j] = r_middle_velocity[j] + (mTime.Previous - mTime.PreviousMiddle) * r_current_acceleration[j]; //+ actual_velocity;

                }
                // r_current_displacement[j]  = 0.0;

            } // for DomainSize

        }     // for node parallel

        mTime.Previous = mTime.Current;
        mTime.PreviousMiddle = mTime.Middle;
        KRATOS_CATCH("")
    }

   
    /*@} */
    /**@name Operations */
    /*@{ */
    /*@} */
    /**@name Access */
    /*@{ */
    /*@} */
    /**@name Inquiry */
    /*@{ */
    /*@} */
    /**@name Friends */
    /*@{ */

protected:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */

    struct DeltaTimeParameters {
        double PredictionLevel; // 0, 1, 2 // NOTE: Should be a integer?
        double Maximum;         // Maximum delta time
        double Fraction;        // Fraction of the delta time
    };

    /**
     * @brief This struct contains the details of the time variables
     */
    struct TimeVariables {
        double PreviousMiddle; // n-1/2
        double Previous;       // n
        double Middle;         // n+1/2
        double Current;        // n+1

        double Delta;          // Time step
    };

    ///@name Protected static Member Variables
    ///@{

    TimeVariables mTime;            /// This struct contains the details of the time variables

    ModelPart& mr_grid_model_part;

    const int mStressUpdateOption; // 0 = USF, 1 = USL, 2 = MUSL
    const bool mIsCentralDifference;

    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    /*@} */
    /**@name Protected Operations*/
    /*@{ */
    /*@} */
    /**@name Protected  Access */
    /*@{ */
    /*@} */
    /**@name Protected Inquiry */
    /*@{ */
    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */
private:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */
    /*@} */
    /**@name Private Operators*/
    /*@{ */
    /*@} */
    /**@name Private Operations*/
    /*@{ */
    /*@} */
    /**@name Private  Access */
    /*@{ */
    /*@} */
    /**@name Private Inquiry */
    /*@{ */
    /*@} */
    /**@name Unaccessible methods */
    /*@{ */
}; /* Class MPMExplicitScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_MPM_EXPLICIT_SCHEME defined */
