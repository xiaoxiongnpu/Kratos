//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Suneth Warnakulasuriya
//                 Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME)
#define KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME

// Project includes
#include "containers/array_1d.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/openmp_utils.h"
#include "fluid_dynamics_application_variables.h"

// debugging
#include "input_output/vtk_output.h"

#include "custom_strategies/relaxed_dof_updater.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class GenericResidualBasedSimpleSteadyScalarScheme
    : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GenericResidualBasedSimpleSteadyScalarScheme);

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using GeometryType = Element::GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    GenericResidualBasedSimpleSteadyScalarScheme(const double RelaxationFactor)
        : mRelaxationFactor(RelaxationFactor)
    {
        KRATOS_INFO("GenericResidualBasedSimpleSteadyScalarScheme")
            << " Using residual based simple steady scheme with relaxation "
               "factor = "
            << std::scientific << mRelaxationFactor << "\n";
    }

    ~GenericResidualBasedSimpleSteadyScalarScheme() override = default;

    ///@}
    ///@name Operators
    ///@{

    void InitializeSolutionStep(ModelPart& r_model_part,
                                TSystemMatrixType& A,
                                TSystemVectorType& Dx,
                                TSystemVectorType& b) override
    {
        BaseType::InitializeSolutionStep(r_model_part, A, Dx, b);
        mIterationCounter = 0;
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                TSystemMatrixType& rA,
                TSystemVectorType& rDx,
                TSystemVectorType& rb) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet, rDx, mRelaxationFactor);

        KRATOS_CATCH("");
    }

    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    void CalculateSystemContributions(Element::Pointer rCurrentElement,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentElement->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentElement->CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        Matrix SteadyLHS;
        rCurrentElement->CalculateLocalVelocityContribution(
            SteadyLHS, RHS_Contribution, CurrentProcessInfo);
        rCurrentElement->EquationIdVector(EquationId, CurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(LHS_Contribution) += SteadyLHS;

        const double residual_coeff = rCurrentElement->GetValue(ERROR_OVERALL);
        Vector values;
        rCurrentElement->GetFirstDerivativesVector(values, 0);
        AddRelaxation(rCurrentElement->GetGeometry(), LHS_Contribution,
                      RHS_Contribution, residual_coeff, values, CurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
                                                LocalSystemMatrixType& LHS_Contribution,
                                                LocalSystemVectorType& RHS_Contribution,
                                                Condition::EquationIdVectorType& EquationId,
                                                ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentCondition->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentCondition->CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        Matrix SteadyLHS;
        rCurrentCondition->CalculateLocalVelocityContribution(
            SteadyLHS, RHS_Contribution, CurrentProcessInfo);
        rCurrentCondition->EquationIdVector(EquationId, CurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(LHS_Contribution) += SteadyLHS;

        KRATOS_CATCH("");
    }

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
                                    LocalSystemVectorType& rRHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        Matrix LHS_Contribution;
        CalculateSystemContributions(rCurrentElement, LHS_Contribution, rRHS_Contribution,
                                     rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
                                              LocalSystemVectorType& rRHS_Contribution,
                                              Element::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        Matrix LHS_Contribution;
        Condition_CalculateSystemContributions(rCurrentCondition, LHS_Contribution,
                                               rRHS_Contribution, rEquationId,
                                               rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected Operators
    ///@{
    void AddRelaxation(const GeometryType& rGeometry,
                       LocalSystemMatrixType& LHS_Contribution,
                       LocalSystemVectorType& RHS_Contribution,
                       const double ResidualScalar,
                       const Vector& rValues,
                       ProcessInfo& CurrentProcessInfo)
    {
        if (LHS_Contribution.size1() == 0)
            return;

        const unsigned int NumNodes = rGeometry.PointsNumber();
        Matrix Mass;
        this->CalculateLumpedMassMatrix(rGeometry, Mass);

        for (unsigned int iNode = 0; iNode < NumNodes; iNode++)
        {
            const array_1d<double, 3>& r_velocity =
                rGeometry[iNode].FastGetSolutionStepValue(VELOCITY);
            const double length =
                rGeometry[iNode].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
            const double velocity_magnitude = norm_2(r_velocity);

            const double local_limit_dt =
                (velocity_magnitude > 0.0) ? length / velocity_magnitude : 1.0;

            Mass(iNode, iNode) *= ResidualScalar / local_limit_dt;
        }

        noalias(LHS_Contribution) += Mass;
    }

    void CalculateLumpedMassMatrix(const GeometryType& rGeometry,
                                   LocalSystemMatrixType& rLumpedMass) const
    {
        const unsigned int number_of_nodes = rGeometry.PointsNumber();

        if (rLumpedMass.size1() != number_of_nodes)
        {
            rLumpedMass.resize(number_of_nodes, number_of_nodes, false);
        }

        const double size_fraction = rGeometry.DomainSize() / number_of_nodes;
        noalias(rLumpedMass) = IdentityMatrix(number_of_nodes, number_of_nodes);
        noalias(rLumpedMass) = rLumpedMass * size_fraction;
    }
    ///@}

private:
    ///@name Member Variables
    ///@{

    // TSystemVectorType mPreviousB;

    double mPreviousRelaxationFactor;

    unsigned int mIterationCounter = 0;

    VtkOutput* mVtkOutput;

    using DofUpdaterType = RelaxedDofUpdater<TSparseSpace>;
    using DofUpdaterPointerType = typename DofUpdaterType::UniquePointer;

    DofUpdaterPointerType mpDofUpdater = Kratos::make_unique<DofUpdaterType>();

    double mRelaxationFactor;

    ///@}
};

///@}

} // namespace Kratos

#endif /* KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME defined */
