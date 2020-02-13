//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

#if !defined(KRATOS_MPM_EXPLICIT_UTILITIES)
#define KRATOS_MPM_EXPLICIT_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

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
 * @namespace MPMExplicitUtilities
 * @ingroup ParticleMechanicsApplication
 * @brief This namespace includes several utilities necessaries for the computation of the explicit integration
 * @author Peter Wilson
 */
namespace MPMExplicitUtilities
{
    /// The size type definition
    typedef std::size_t SizeType;

    /// The index type definition
    typedef std::size_t IndexType;

    /// The arrays of elements and nodes
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;

    /**
     * @brief This method computes the necessry delta time to avoid numerical instabilities
     * @param rModelPart The model of the problem to solve
     * @param ThisParameters The configuration parameters
     * @return The critical delta time
     */
    double KRATOS_API(PARTICLE_MECHANICS_APPLICATION) CalculateDeltaTime(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /**
     * @brief This method computes the necessry delta time to avoid numerical instabilities (inner method)
     * @param rModelPart The model of the problem to solve
     * @param TimeStepPredictionLevel The prediction level
     * @param MaxDeltaTime The max delta time
     * @param SafetyFactor The safety factor
     * @param MassFactor The factor that multiplies the mass
     * @return The critical delta time
     */
    double KRATOS_API(PARTICLE_MECHANICS_APPLICATION) InnerCalculateDeltaTime(
        ModelPart& rModelPart,
        const double TimeStepPredictionLevel,
        const double MaxDeltaTime,
        const double SafetyFactor,
        const double MassFactor
        );

    // CalculateExplicitKinematics

    // CalcuateExplicitInternalForce

    // FinalizeSolutionStepExplicit

    // UpdateGaussPointExplicit

}; // namespace ExplicitIntegrationUtilities
}  // namespace Kratos
#endif /* KRATOS_MPM_EXPLICIT_UTILITIES defined */
