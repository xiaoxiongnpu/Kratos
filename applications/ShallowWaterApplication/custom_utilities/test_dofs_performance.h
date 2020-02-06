
#ifndef KRATOS_TEST_DOF_SET_H_INCLUDED
#define KRATOS_TEST_DOF_SET_H_INCLUDED

#include "includes/dof.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"

namespace Kratos {

class TestDofsPerformance
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(TestDofsPerformance);

    typedef std::vector<Dof<double>::Pointer>                          DofsVectorType;
    typedef ModelPart::DofsArrayType                                    DofsArrayType;
    typedef Element::EquationIdVectorType                        EquationIdVectorType;
    typedef std::unordered_set < Node<3>::DofType::Pointer, DofPointerHasher> SetType;

    typedef UblasSpace<double, CompressedMatrix, Vector>              SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector>                         LocalSpaceType;
    typedef typename LocalSpaceType::VectorType                       LocalVectorType;
    typedef typename SparseSpaceType::VectorType                     SystemVectorType;
    typedef typename SparseSpaceType::VectorPointerType       SystemVectorPointerType;

    TestDofsPerformance(){}

    ~TestDofsPerformance() = default;

    void TestDofs(ModelPart& rModelPart, int NumTimes)
    {
        double initial_time = OpenMPUtils::GetCurrentTime();

        DofsArrayType dofs;

        SetUpDofs(rModelPart, dofs);

        SetUpSystem(rModelPart, dofs);

        SystemVectorPointerType pDx = SparseSpaceType::CreateEmptyVectorPointer();
        SystemVectorType& rDx = *pDx;

        InitializeVector(dofs, rDx);

        double set_up_dofs_time = OpenMPUtils::GetCurrentTime();

        double assemble_time = 0;
        double updating_time = 0;
        for (int j = 0; j < NumTimes; ++j)
        {
            double loop_initial_time = OpenMPUtils::GetCurrentTime();

            SystemVectorType& rDx = *pDx;
            SparseSpaceType::SetToZero(rDx);

            AssembleDofs(rModelPart, dofs, rDx);

            double assemble_dofs_time = OpenMPUtils::GetCurrentTime();

            UpdateDofs(dofs, rDx);

            double update_dofs_time = OpenMPUtils::GetCurrentTime();

            assemble_time += assemble_dofs_time - loop_initial_time;
            updating_time += update_dofs_time - assemble_dofs_time;
        }

        double final_time = OpenMPUtils::GetCurrentTime();

        int one_time = 1;

        std::cout << std::endl;
        std::cout << "TESTING DOFS" << std::endl;
        std::cout << "Using " << rModelPart.Elements().size() << " elements and iterating " << NumTimes << " times:" << std::endl;
        std::cout << "\tSetting up the system (x" << one_time << ") \t" << set_up_dofs_time - initial_time << std::endl;
        std::cout << "\tAssembling the vector (x" << NumTimes << ") \t" << assemble_time << std::endl;
        std::cout << "\tUpdating the dofs     (x" << NumTimes << ") \t" << updating_time << std::endl;
        std::cout << "\tTotal loop time       (x" << NumTimes << ") \t" << final_time - set_up_dofs_time << std::endl;
    }

    void TestNodes(ModelPart& rModelPart, int NumTimes)
    {
        double initial_time = OpenMPUtils::GetCurrentTime();

        for (int j = 0; j < NumTimes; ++j)
        {
            LoopElements(rModelPart);
            UpdateNodes(rModelPart);
        }

        double final_time = OpenMPUtils::GetCurrentTime();

        std::cout << std::endl;
        std::cout << "TESTING NODES" << std::endl;
        std::cout << "Time for just writing on the nodes using " << rModelPart.Elements().size() << " elements and iterating " << NumTimes << " times:" << std::endl;
        std::cout << "\tWriting on the nodes  (x" << NumTimes << ") \t" << final_time - initial_time << std::endl;
    }

    void SetUpDofs(ModelPart& rModelPart, DofsArrayType& rDofs)
    {
        DofsVectorType tmp_dof_list;
        SetType dof_global_set;
        dof_global_set.reserve(rModelPart.Elements().size()*20);

        #pragma omp parallel firstprivate(tmp_dof_list)
        {
            auto& r_process_info = rModelPart.GetProcessInfo();

            // We create the temporal set and we reserve some space on them
            SetType dofs_tmp_set;
            dofs_tmp_set.reserve(20000);

            // Gets the array of elements from the modeler
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i)
            {
                auto it_elem = rModelPart.ElementsBegin() + i;

                // Gets list of Dof involved on every element
                it_elem->GetDofList(tmp_dof_list, r_process_info);
                dofs_tmp_set.insert(tmp_dof_list.begin(), tmp_dof_list.end());
            }


            // We merge all the sets in one thread
            #pragma omp critical
            {
                dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
            }
        }

        DofsArrayType dof_list;
        dof_list.reserve(dof_global_set.size());
        for (auto it = dof_global_set.begin(); it != dof_global_set.end(); ++it)
        {
            dof_list.push_back( *it );
        }
        dof_list.Sort();
        rDofs = DofsArrayType();
        rDofs = dof_list;
    }

    void SetUpSystem(ModelPart& rModelPart, DofsArrayType& rDofs)
    {
        int n_dofs = static_cast<int>(rDofs.size());

        #pragma omp parallel for firstprivate(n_dofs)
        for (int i = 0; i < n_dofs; i++) {
            typename DofsArrayType::iterator it_dof = rDofs.begin() + i;
            it_dof->SetEquationId(i);
        }
    }

    void InitializeVector(DofsArrayType& rDofs, SystemVectorType& rDx)
    {
        int n_dofs = rDofs.size();
        rDx.resize(n_dofs, false);
    }

    void AssembleDofs(ModelPart& rModelPart, DofsArrayType& rDofs, SystemVectorType& rDx)
    {
        auto& r_process_info = rModelPart.GetProcessInfo();
        
        // Contributions to the system
        LocalVectorType rhs_contribution = LocalVectorType(0);
        
        // Vector containing the localization in the system of the different terms
        EquationIdVectorType equation_id;
        
        const int n_elements = static_cast<int>(rModelPart.NumberOfElements());
        #pragma omp parallel firstprivate(n_elements, rhs_contribution, equation_id)
        {
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < n_elements; i++)
            {
                auto it_elem = rModelPart.ElementsBegin() + i;

                // Calculate the elemental contribution
                it_elem->EquationIdVector(equation_id, r_process_info);
                it_elem->CalculateRightHandSide(rhs_contribution, r_process_info);

                // Assemble the elemental contribution
                AssembleRHS(rDx, rhs_contribution, equation_id);
            }
        }
    }

    void AssembleRHS(SystemVectorType& rb, LocalVectorType& rRHS_Contribution, EquationIdVectorType& rEquationId)
    {
        unsigned int local_size = rRHS_Contribution.size();

        for (unsigned int i_local = 0; i_local < local_size; ++i_local)
        {
            unsigned int i_global = rEquationId[i_local];

            double& b_value = rb[i_global];
            const double& rhs_value = rRHS_Contribution[i_local];

            #pragma omp atomic
            b_value += rhs_value;
        }
    }

    void UpdateDofs(DofsArrayType& rDofs, const SystemVectorType& rDx)
    {
        const int num_dof = static_cast<int>(rDofs.size());

        #pragma omp parallel for
        for(int i = 0; i < num_dof; ++i)
        {
            auto it_dof = rDofs.begin() + i;

            if (it_dof->IsFree())
                it_dof->GetSolutionStepValue() =
                    it_dof->GetSolutionStepValue(1) + SparseSpaceType::GetValue(rDx, it_dof->EquationId());
        }
    }

    void LoopElements(ModelPart& rModelPart)
    {
        const auto& r_process_info = rModelPart.GetProcessInfo();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rModelPart.NumberOfElements()); ++i)
        {
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->AddExplicitContribution(r_process_info);
        }
    }

    void UpdateNodes(ModelPart& rModelPart)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
        {
            auto it_node = rModelPart.NodesBegin() + i;

            auto vn = it_node->FastGetSolutionStepValue(VELOCITY,1);
            auto dv = it_node->FastGetSolutionStepValue(RHS_VELOCITY);
            it_node->FastGetSolutionStepValue(VELOCITY) = vn + dv;

            auto hn = it_node->FastGetSolutionStepValue(HEIGHT,1);
            auto dh = it_node->FastGetSolutionStepValue(RHS_HEIGHT);
            it_node->FastGetSolutionStepValue(HEIGHT) = hn + dh;
        }
    }

}; // class TestDofsPerformance

} // namespace Kratos

#endif