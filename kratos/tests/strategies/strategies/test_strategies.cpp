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
#include <limits>

// External includes


// Project includes
#include "testing/testing.h"

// Utility includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

namespace Kratos 
{
    namespace Testing 
    {
        /// Tests
       
//         typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//         typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
//         typedef LinearSolver<SparseSpaceType,LocalSpaceType> SolverType;
//         
//         typedef Scheme< SparseSpaceType, LocalSpaceType >  SchemeType;
//         typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, SolverType >  SolvingStrategyType;
//         typedef PointerVectorSet<Dof<double>, SetIdentityFunction<Dof<double>>, std::less<SetIdentityFunction<Dof<double>>::result_type>, std::equal_to<SetIdentityFunction<Dof<double>>::result_type>, Dof<double>* > DofsArrayType;
        
//         static inline DofsArrayType BasicTestStrategyDisplacement(
//             ModelPart& ModelPart,
//             SolvingStrategyType::Pointer pStrategy,
//             std::vector< Dof<double>::Pointer >& DoF
//             )
//         {
//             ModelPart.SetBufferSize(3);
//             
//             ModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
//             
//             auto pnode = ModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
//             
//             pnode->AddDof(DISPLACEMENT_X);
//             pnode->AddDof(DISPLACEMENT_Y);
//             pnode->AddDof(DISPLACEMENT_Z);
//             
//             DoF.reserve(3);
//             DoF.push_back(pnode->pGetDof(DISPLACEMENT_X));
//             DoF.push_back(pnode->pGetDof(DISPLACEMENT_Y));
//             DoF.push_back(pnode->pGetDof(DISPLACEMENT_Z));
//             
//             // Set initial solution
//             const array_1d<double, 3> zero_vector = ZeroVector(3);
//             pnode->FastGetSolutionStepValue(DISPLACEMENT) = zero_vector;
//             pnode->FastGetSolutionStepValue(DISPLACEMENT, 1) = zero_vector;
//             pnode->FastGetSolutionStepValue(DISPLACEMENT, 2) = zero_vector;
//             
//             DofsArrayType Doftemp;
//             Doftemp.reserve(DoF.size());
//             for (auto it= DoF.begin(); it!= DoF.end(); it++)
//                 Doftemp.push_back( it->get() );
//             Doftemp.Sort();
//             
//             CompressedMatrix A = ZeroMatrix(3, 3);
//             Vector Dx = ZeroVector(3);
//             Vector b = ZeroVector(3);
//             
//             pStrategy->Initialize();
//             
//             return Doftemp;
//         }
     
//         /** 
//          * Checks if the Newton Rapshon sstrategy performs correctly the resolution of the system
//          */
//         
//         KRATOS_TEST_CASE_IN_SUITE(DisplacementNRStrategyTest, KratosCoreStrategiesFastSuite) 
//         {
//             constexpr double tolerance = 1e-6;
//             
//             ModelPart model_part("Main");
//             
//             typedef ResidualBasedNewmarkDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedNewmarkDisplacementSchemeType;
//             typename SchemeType::Pointer pscheme = typename SchemeType::Pointer( new ResidualBasedNewmarkDisplacementSchemeType() );
//             
//             ResidualBasedNewtonRaphsonStrategy(
//         ModelPart &model_part,
//         typename TSchemeType::Pointer pScheme,
//         typename TLinearSolver::Pointer pNewLinearSolver,
//         typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria
//             
//             const double delta_time = 1.0e-4;
// 
//             std::vector< Dof<double>::Pointer > DoF;
//             DofsArrayType Doftemp = BasicTestSchemeDisplacement(model_part, pscheme, DoF, delta_time);
//             
//             CompressedMatrix A = ZeroMatrix(3, 3);
//             Vector Dx = ZeroVector(3);
//             Vector b = ZeroVector(3);
//             
//             Node<3>::Pointer pnode = model_part.pGetNode(1);
//             
//             pscheme->Initialize(model_part);
//             
//             const unsigned int number_iterations = 10;
//             for (unsigned int iter = 0; iter < number_iterations; ++iter)
//             {
//                 time += delta_time;
//                 
//                 model_part.CloneTimeStep(time);
//                 
//                 Dx[0] = std::cos(time) - std::cos(time - delta_time);
//                 
//                 pscheme->InitializeSolutionStep(model_part, A, Dx, b);
//                 pscheme->Update(model_part, Doftemp, A, Dx, b);
//                 
//                 const double x = pnode->FastGetSolutionStepValue(DISPLACEMENT_X);
//                 
// //                 // Debug
// //                 std::cout << time << "\t" << x << "\t" << v << "\t" << a << std::endl;
//                 
// //                 KRATOS_CHECK_LESS_EQUAL(std::abs(x - std::cos(time)), tolerance);
// //                 KRATOS_CHECK_LESS_EQUAL(std::abs(v + std::sin(time)), tolerance);
// //                 KRATOS_CHECK_LESS_EQUAL(std::abs(a + std::cos(time)), tolerance);
//             }
//         }
        
    } // namespace Testing
}  // namespace Kratos.

