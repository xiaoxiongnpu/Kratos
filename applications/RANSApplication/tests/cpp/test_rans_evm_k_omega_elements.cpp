//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

// System includes
#include <string>
#include <vector>
#include <iostream>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_components.h"
#include "includes/variables.h"
#include "testing/testing.h"

// Application includes
#include "custom_elements/evm_k_omega/rans_evm_omega_element.h"
#include "custom_elements/evm_k_omega/rans_evm_k_omega_k_element.h"
#include "custom_elements/stabilized_convection_diffusion_reaction_utilities.h"
#include "rans_application_variables.h"

namespace Kratos
{
namespace Testing
{
namespace
{
void CreateEvmUnitTestModelPart(const std::string& rElementName,
                                const Variable<double>& rDofVariable,
                                ModelPart& rModelPart)
{
    const auto& r_proto = KratosComponents<Element>::Get(rElementName);
    auto node_ids = std::vector<ModelPart::IndexType>{};
    Matrix coords;
    r_proto.GetGeometry().PointsLocalCoordinates(coords);
    if (coords.size2() == 2)
    {
        for (std::size_t i = 0; i < coords.size1(); ++i)
            rModelPart.CreateNewNode(i + 1, coords(i, 0), coords(i, 1), 0.0);
    }
    else
    {
        for (std::size_t i = 0; i < coords.size1(); ++i)
            rModelPart.CreateNewNode(i + 1, coords(i, 0), coords(i, 1), coords(i, 2));
    }
    for (auto& r_node : rModelPart.Nodes())
    {
        r_node.AddDof(rDofVariable).SetEquationId(r_node.Id());
        node_ids.push_back(r_node.Id());
    }
    auto p_prop = rModelPart.CreateNewProperties(1);
    auto p_element = rModelPart.CreateNewElement(rElementName, 1, node_ids, p_prop);
    rModelPart.SetBufferSize(1);
}

template <typename EvmElement>
double EffectiveKinematicViscosity(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    auto& r_geom = rElement.GetGeometry();
    const Vector N = row(r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1), 0);
    typename Element::GeometryType::ShapeFunctionsGradientsType DN_DX;
    r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::GI_GAUSS_1);
    typename EvmElement::BaseType::ConvectionDiffusionReactionDataType data;
    auto& rRANSElement = dynamic_cast<const typename EvmElement::BaseType&>(rElement);
    rRANSElement.CalculateElementData(data, N, DN_DX[0], rProcessInfo);
    return rRANSElement.CalculateEffectiveKinematicViscosity(data, N, DN_DX[0], rProcessInfo);
}

template <typename EvmElement>
typename EvmElement::BaseType::ConvectionDiffusionReactionDataType ConvectionReactionDiffusionData(
    const Element& rElement, const ProcessInfo& rProcessInfo)
{
    auto& r_geom = rElement.GetGeometry();
    const Vector N = row(r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1), 0);
    typename Element::GeometryType::ShapeFunctionsGradientsType DN_DX;
    r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::GI_GAUSS_1);
    typename EvmElement::BaseType::ConvectionDiffusionReactionDataType data;
    auto& rRANSElement = dynamic_cast<const typename EvmElement::BaseType&>(rElement);
    rRANSElement.CalculateElementData(data, N, DN_DX[0], rProcessInfo);
    return data;
}

template <typename EvmElement>
double ReactionTerm(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const auto data = ConvectionReactionDiffusionData<EvmElement>(rElement, rProcessInfo);
    auto& rRANSElement = dynamic_cast<const typename EvmElement::BaseType&>(rElement);
    auto& r_geom = rElement.GetGeometry();
    const Vector N = row(r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1), 0);
    typename Element::GeometryType::ShapeFunctionsGradientsType DN_DX;
    r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::GI_GAUSS_1);
    return rRANSElement.CalculateReactionTerm(data, N, DN_DX[0], rProcessInfo);
}

template <typename EvmElement>
double DynamicReactionTerm(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const auto data = ConvectionReactionDiffusionData<EvmElement>(rElement, rProcessInfo);
    auto& rRANSElement = dynamic_cast<const typename EvmElement::BaseType&>(rElement);
    auto& r_geom = rElement.GetGeometry();
    const Vector N = row(r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1), 0);
    typename Element::GeometryType::ShapeFunctionsGradientsType DN_DX;
    r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::GI_GAUSS_1);
    const double reaction =
        rRANSElement.CalculateReactionTerm(data, N, DN_DX[0], rProcessInfo);
    const double delta_time = rProcessInfo[DELTA_TIME];
    const double dynamic_tau = rProcessInfo[DYNAMIC_TAU];
    const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
    const double bossak_gamma =
        TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();

    return reaction + dynamic_tau * (1. - bossak_alpha) / (bossak_gamma * delta_time);
}

template <typename EvmElement>
double SourceTerm(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    auto data = ConvectionReactionDiffusionData<EvmElement>(rElement, rProcessInfo);
    auto& rRANSElement = dynamic_cast<const typename EvmElement::BaseType&>(rElement);
    auto& r_geom = rElement.GetGeometry();
    const Vector N = row(r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1), 0);
    typename Element::GeometryType::ShapeFunctionsGradientsType DN_DX;
    r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::GI_GAUSS_1);
    return rRANSElement.CalculateSourceTerm(data, N, DN_DX[0], rProcessInfo);
}

template <typename EvmElement>
double VelocityNorm(const Element& rElement)
{
    auto& r_geom = rElement.GetGeometry();
    const Vector N = row(r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1), 0);
    auto& rRANSElement = dynamic_cast<const typename EvmElement::BaseType&>(rElement);
    const array_1d<double, 3> velocity = rRANSElement.EvaluateInPoint(VELOCITY, N);
    return norm_2(velocity);
}

template <typename EvmElement>
double DivVel(const Element& rElement)
{
    auto& r_geom = rElement.GetGeometry();
    typename Element::GeometryType::ShapeFunctionsGradientsType DN_DX;
    r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::GI_GAUSS_1);
    auto& rRANSElement = dynamic_cast<const typename EvmElement::BaseType&>(rElement);
    const double div = rRANSElement.GetDivergenceOperator(VELOCITY, DN_DX[0]);
    return div;
}

template <typename EvmElement>
void CalculateStreamlineAndCrossWindDiffusionParameters(double& rChi,
                                                        double& rSD,
                                                        double& rCD,
                                                        const Element& rElement,
                                                        const ProcessInfo& rProcessInfo)
{
    // This code was extracted from StabilizedConvectionDiffusionReaction::CalculateDampingMatrix.
    // We needed to add the tests before changing code, but this should be cleaned during refactoring.
    auto& r_geom = rElement.GetGeometry();
    const Vector N = row(r_geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1), 0);
    typename Element::GeometryType::ShapeFunctionsGradientsType DN_DX;
    r_geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::GI_GAUSS_1);
    typename EvmElement::BaseType::ConvectionDiffusionReactionDataType data;
    auto& rRANSElement = dynamic_cast<const typename EvmElement::BaseType&>(rElement);
    rRANSElement.CalculateElementData(data, N, DN_DX[0], rProcessInfo);
    const double effective_kinematic_viscosity =
        rRANSElement.CalculateEffectiveKinematicViscosity(data, N, DN_DX[0], rProcessInfo);
    const array_1d<double, 3> velocity = rRANSElement.EvaluateInPoint(VELOCITY, N);
    const double velocity_magnitude = norm_2(velocity);
    const Matrix r_parameter_derivatives_g =
        rRANSElement.GetGeometryParameterDerivatives()[0];
    Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                       r_parameter_derivatives_g.size2());
    noalias(contravariant_metric_tensor) =
        prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);
    const double reaction =
        rRANSElement.CalculateReactionTerm(data, N, DN_DX[0], rProcessInfo);
    const double delta_time = rProcessInfo[DELTA_TIME];
    const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
    const double bossak_gamma =
        TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
    const double dynamic_tau = rProcessInfo[DYNAMIC_TAU];
    double tau, element_length;
    StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
        tau, element_length, velocity, contravariant_metric_tensor, reaction,
        effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time, dynamic_tau);
    StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
        rChi, rSD, rCD, velocity_magnitude, tau, effective_kinematic_viscosity, reaction,
        bossak_alpha, bossak_gamma, delta_time, element_length, dynamic_tau);
}

void EvmKElement2D3N_SetUp(ModelPart& rModelPart)
{
    // rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1); // relaxed turb kin energy rate
    // rModelPart.AddNodalSolutionStepVariable(RANS_Y_PLUS);
    // rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    CreateEvmUnitTestModelPart("RansEvmKOmegaK2D3N", TURBULENT_KINETIC_ENERGY, rModelPart);
}

void EvmOmegaElement2D3N_SetUp(ModelPart& rModelPart)
{
    // rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    // rModelPart.AddNodalSolutionStepVariable(NORMAL);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2); // relaxed turb kin energy diss. rate 2
    // rModelPart.AddNodalSolutionStepVariable(RANS_Y_PLUS);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    CreateEvmUnitTestModelPart("RansEvmKOmegaOmega2D3N",
                               TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, rModelPart);
}


template <typename EvmElement>
void CheckEvmElementTestData(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    // Ensure several contributions have comparable magnitudes:
    const double effective_kinematic_viscosity =
        EffectiveKinematicViscosity<EvmElement>(rElement, rProcessInfo);
    KRATOS_CHECK(effective_kinematic_viscosity > 0.1 && effective_kinematic_viscosity < 100.0);
    const double reaction = ReactionTerm<EvmElement>(rElement, rProcessInfo);
    KRATOS_CHECK(reaction > 0.1 && reaction < 100.0);
    const double dynamic_reaction = DynamicReactionTerm<EvmElement>(rElement, rProcessInfo);
    KRATOS_CHECK(dynamic_reaction > 0.1 && dynamic_reaction < 100.0);
    const double source = SourceTerm<EvmElement>(rElement, rProcessInfo);
    KRATOS_CHECK(source > 0.1 && source < 100.0);
    const double vel_norm = VelocityNorm<EvmElement>(rElement);
    KRATOS_CHECK(vel_norm > 0.1 && vel_norm < 100.0);
    const double div_vel = DivVel<EvmElement>(rElement);
    KRATOS_CHECK(div_vel > 0.1 && div_vel < 10.0);
    double chi{}, streamline_diffusion{}, crosswind_diffusion{};
    CalculateStreamlineAndCrossWindDiffusionParameters<EvmElement>(
        chi, streamline_diffusion, crosswind_diffusion, rElement, rProcessInfo);
    KRATOS_CHECK(chi > 0.1 && chi < 100.0);
    KRATOS_CHECK(streamline_diffusion > 0.1 && streamline_diffusion < 100.0);
    KRATOS_CHECK(crosswind_diffusion > 0.1 && crosswind_diffusion < 100.0);
}

void EvmKElement2D3N_AssignTestData(ModelPart& rModelPart)
{
    rModelPart.GetProcessInfo()[BOSSAK_ALPHA] = -0.3;
    rModelPart.GetProcessInfo()[DELTA_TIME] = 1.1;
    rModelPart.GetProcessInfo()[DYNAMIC_TAU] = 0.1;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU] = 0.09;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_SIGMA_K] = 0.5;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_GAMMA] = 0.52;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_BETA_ZERO_STAR] = 0.09;
    auto& node1 = rModelPart.GetNode(1);
    node1.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 0.24;
    node1.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) = 0.21;
    node1.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 25.90;
    node1.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 1.29;
    node1.FastGetSolutionStepValue(VELOCITY_X) = -1.85;
    node1.FastGetSolutionStepValue(VELOCITY_Y) = -2.15;
    auto& node2 = rModelPart.GetNode(2);
    node2.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 0.28;
    node2.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) = 0.14;
    node2.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 5.00;
    node2.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 1.04;
    node2.FastGetSolutionStepValue(VELOCITY_X) = -0.78;
    node2.FastGetSolutionStepValue(VELOCITY_Y) = -0.06;
    auto& node3 = rModelPart.GetNode(3);
    node3.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 0.28;
    node3.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) = 0.63;
    node3.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 5.70;
    node3.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 1.11;
    node3.FastGetSolutionStepValue(VELOCITY_X) = 0.01;
    node3.FastGetSolutionStepValue(VELOCITY_Y) = -0.65;
    CheckEvmElementTestData<RansEvmKOmegaKElement<2, 3>>(
        rModelPart.Elements().front(), rModelPart.GetProcessInfo());
}

void EvmOmegaElement2D3N_AssignTestData(ModelPart& rModelPart)
{

    rModelPart.GetProcessInfo()[BOSSAK_ALPHA] = -0.3;
    rModelPart.GetProcessInfo()[DELTA_TIME] = 1.1;
    rModelPart.GetProcessInfo()[DYNAMIC_TAU] = 0.1;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_C1] = 1.44;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_C2] = 1.92;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU] = 0.09;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_SIGMA_OMEGA] = 0.3;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_GAMMA] = 0.52;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_BETA_ZERO] = 0.072;
    rModelPart.GetProcessInfo()[TURBULENCE_RANS_BETA_ZERO_STAR] = 0.9;

    auto& node1 = rModelPart.GetNode(1);
    node1.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 0.28;
    node1.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) = 0.21;
    node1.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE) = 1;
    node1.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 8;
    node1.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 1.29;
    node1.FastGetSolutionStepValue(VELOCITY_X) = -1.85;
    node1.FastGetSolutionStepValue(VELOCITY_Y) = -2.15;

    auto& node2 = rModelPart.GetNode(2);
    node2.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 0.24;
    node2.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) = 0.14;
    node2.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE) = 1;
    node2.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 2.5;
    node2.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 1.04;
    node2.FastGetSolutionStepValue(VELOCITY_X) = -0.78;
    node2.FastGetSolutionStepValue(VELOCITY_Y) = -0.06;

    auto& node3 = rModelPart.GetNode(3);
    node3.FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 0.24;
    node3.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) = 0.63;
    node3.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE) = 1;
    node3.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 2.85;
    node3.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 0.90;
    node3.FastGetSolutionStepValue(VELOCITY_X) = 0.01;
    node3.FastGetSolutionStepValue(VELOCITY_Y) = -0.65;

    CheckEvmElementTestData<RansEvmOmegaElement<2, 3>>(
        rModelPart.Elements().front(), rModelPart.GetProcessInfo());

}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaKElement2D3N_EquationIdVector, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmKElement2D3N_SetUp(model_part);
    // Test:
    auto& r_element = model_part.Elements().front();
    auto eqn_ids = std::vector<std::size_t>{};
    r_element.EquationIdVector(eqn_ids, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(eqn_ids.size(), r_element.GetGeometry().PointsNumber());
    for (std::size_t i = 0; i < eqn_ids.size(); ++i)
    {
        KRATOS_CHECK_EQUAL(eqn_ids[i], i + 1);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaKElement2D3N_GetDofList, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmKElement2D3N_SetUp(model_part);
    // Test:
    auto& r_element = model_part.Elements().front();
    auto dofs = Element::DofsVectorType{};
    r_element.GetDofList(dofs, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(dofs.size(), r_element.GetGeometry().PointsNumber());
    for (std::size_t i = 0; i < dofs.size(); ++i)
    {
        KRATOS_CHECK_EQUAL(dofs[i]->GetVariable(), TURBULENT_KINETIC_ENERGY);
        KRATOS_CHECK_EQUAL(dofs[i]->EquationId(), i + 1);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaKElement2D3N_CalculateLocalSystem, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmKElement2D3N_SetUp(model_part);
    EvmKElement2D3N_AssignTestData(model_part);
    auto& r_element = model_part.Elements().front();
    // Test:
    Matrix LHS;
    Vector RHS;
    r_element.CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(RHS.size(), 3);
    KRATOS_CHECK_NEAR(RHS(0), 10.89437737411318, 1e-12);
    KRATOS_CHECK_NEAR(RHS(1), 2.546968564684157, 1e-12);
    KRATOS_CHECK_NEAR(RHS(2), 2.388045431925697, 1e-12);
    KRATOS_CHECK_EQUAL(LHS.size1(), 3);
    KRATOS_CHECK_EQUAL(LHS.size2(), 3);
    for (std::size_t i = 0; i < LHS.size1(); ++i)
        for (std::size_t j = 0; j < LHS.size1(); ++j)
            KRATOS_CHECK_EQUAL(LHS(i, j), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaKElement2D3N_CalculateRightHandSide, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmKElement2D3N_SetUp(model_part);
    EvmKElement2D3N_AssignTestData(model_part);
    auto& r_element = model_part.Elements().front();
    // Test:
    Vector RHS;
    r_element.CalculateRightHandSide(RHS, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(RHS.size(), 3);
    KRATOS_CHECK_NEAR(RHS(0), 10.89437737411318, 1e-12);
    KRATOS_CHECK_NEAR(RHS(1), 2.546968564684157, 1e-12);
    KRATOS_CHECK_NEAR(RHS(2), 2.388045431925697, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaKElement2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmKElement2D3N_SetUp(model_part);
    EvmKElement2D3N_AssignTestData(model_part);
    auto& r_element = model_part.Elements().front();
    // Test:
    Vector residual = ZeroVector(3);
    Matrix D;
    r_element.CalculateLocalVelocityContribution(D, residual, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(D.size1(), 3);
    KRATOS_CHECK_EQUAL(D.size2(), 3);

    KRATOS_CHECK_NEAR(D(0, 2), -1.187313442742871, 1e-12);
    KRATOS_CHECK_NEAR(D(1, 2), -0.1603968284080569, 1e-12);
    KRATOS_CHECK_NEAR(D(2, 0), -0.7235634427428708, 1e-12);
    KRATOS_CHECK_NEAR(D(2, 1), -0.1474801617413903, 1e-12);
    KRATOS_CHECK_EQUAL(residual.size(), 3);

    KRATOS_CHECK_NEAR(residual(0), -7.487817601538771e+01, 1e-12);
    KRATOS_CHECK_NEAR(residual(1), 1.385629725869158e+01, 1e-12);
    KRATOS_CHECK_NEAR(residual(2), 1.349896599149478e+01, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaKElement2D3N_CalculateMassMatrix, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmKElement2D3N_SetUp(model_part);
    EvmKElement2D3N_AssignTestData(model_part);
    auto& r_element = model_part.Elements().front();
    // Test:
    Matrix M;
    r_element.CalculateMassMatrix(M, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(M.size1(), 3);
    KRATOS_CHECK_EQUAL(M.size2(), 3);
    KRATOS_CHECK_NEAR(M(0, 2), 9.842322258148753e-02, 1e-12);
    KRATOS_CHECK_NEAR(M(1, 2), -1.043933243068875e-02, 1e-12);
    KRATOS_CHECK_NEAR(M(2, 0), -3.197047555027811e-02, 1e-12);
    KRATOS_CHECK_NEAR(M(2, 1), -1.428989480300361e-02, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmKOmegaKElement2D3N_CalculateDampingMatrix, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmKElement2D3N_SetUp(model_part);
    EvmKElement2D3N_AssignTestData(model_part);
    auto& r_element = model_part.Elements().front();
    // Test:
    Matrix D;
    r_element.CalculateDampingMatrix(D, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(D.size1(), 3);
    KRATOS_CHECK_EQUAL(D.size2(), 3);
    KRATOS_CHECK_NEAR(D(0, 2), -1.187313442742871, 1e-12);
    KRATOS_CHECK_NEAR(D(1, 2), -0.1603968284080569, 1e-12);
    KRATOS_CHECK_NEAR(D(2, 0), -0.7235634427428708, 1e-12);
    KRATOS_CHECK_NEAR(D(2, 1), -0.1474801617413903, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmOmega2D3N_EquationIdVector, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmOmegaElement2D3N_SetUp(model_part);
    // Test:
    auto& r_element = model_part.Elements().front();
    auto eqn_ids = std::vector<std::size_t>{};
    r_element.EquationIdVector(eqn_ids, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(eqn_ids.size(), r_element.GetGeometry().PointsNumber());
    for (std::size_t i = 0; i < eqn_ids.size(); ++i)
    {
        KRATOS_CHECK_EQUAL(eqn_ids[i], i + 1);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmOmega2D3N_GetDofList, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmOmegaElement2D3N_SetUp(model_part);
    // Test:
    auto& r_element = model_part.Elements().front();
    auto dofs = Element::DofsVectorType{};
    r_element.GetDofList(dofs, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(dofs.size(), r_element.GetGeometry().PointsNumber());
    for (std::size_t i = 0; i < dofs.size(); ++i)
    {
        KRATOS_CHECK_EQUAL(dofs[i]->GetVariable(), TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        KRATOS_CHECK_EQUAL(dofs[i]->EquationId(), i + 1);
    }
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmOmega2D3N_CalculateLocalSystem, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmOmegaElement2D3N_SetUp(model_part);
    EvmOmegaElement2D3N_AssignTestData(model_part);
    auto& r_element = model_part.Elements().front();
    // Test:
    Matrix LHS;
    Vector RHS;
    r_element.CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(RHS.size(), 3);
    // std::cout<<std::scientific<<std::setprecision(15)<<RHS(0)<<std::endl;
    // std::cout<<std::scientific<<std::setprecision(15)<<RHS(1)<<std::endl;
    // std::cout<<std::scientific<<std::setprecision(15)<<RHS(2)<<std::endl;
    KRATOS_CHECK_NEAR(RHS(0), 2.920434481882852e+00, 1e-12);
    KRATOS_CHECK_NEAR(RHS(1), 4.890436151478367e-01, 1e-12);
    KRATOS_CHECK_NEAR(RHS(2), 4.116607565847451e-01, 1e-12);
    KRATOS_CHECK_EQUAL(LHS.size1(), 3);
    KRATOS_CHECK_EQUAL(LHS.size2(), 3);
    for (std::size_t i = 0; i < LHS.size1(); ++i)
        for (std::size_t j = 0; j < LHS.size1(); ++j)
            KRATOS_CHECK_EQUAL(LHS(i, j), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmOmega2D3N_CalculateRightHandSide, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmOmegaElement2D3N_SetUp(model_part);
    EvmOmegaElement2D3N_AssignTestData(model_part);
    auto& r_element = model_part.Elements().front();
    // Test:
    Vector RHS;
    r_element.CalculateRightHandSide(RHS, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(RHS.size(), 3);
    // std::cout<<std::scientific<<std::setprecision(15)<<RHS(0)<<std::endl;
    // std::cout<<std::scientific<<std::setprecision(15)<<RHS(1)<<std::endl;
    // std::cout<<std::scientific<<std::setprecision(15)<<RHS(2)<<std::endl;
    KRATOS_CHECK_NEAR(RHS(0), 2.920434481882852e+00, 1e-12);
    KRATOS_CHECK_NEAR(RHS(1), 4.890436151478367e-01, 1e-12);
    KRATOS_CHECK_NEAR(RHS(2), 4.116607565847451e-01, 1e-12);
}
//ask suneth about the effect of changing kinetic energy as it was making source term negative
KRATOS_TEST_CASE_IN_SUITE(RansEvmOmega2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmOmegaElement2D3N_SetUp(model_part);
    EvmOmegaElement2D3N_AssignTestData(model_part);
    auto& r_element = model_part.Elements().front();
    // Test:
    Vector residual = ZeroVector(3);
    Matrix D;
    r_element.CalculateLocalVelocityContribution(D, residual, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(D.size1(), 3);
    KRATOS_CHECK_EQUAL(D.size2(), 3);   
// std::cout<<std::scientific<<std::setprecision(15)<< D(0, 2)<<std::endl;
    
//     std::cout<<std::scientific<<std::setprecision(15)<<D(1, 2)<<std::endl;
    
//     std::cout<<std::scientific<<std::setprecision(15)<<D(2, 0)<<std::endl;

//     std::cout<<std::scientific<<std::setprecision(15)<<D(2, 1)<<std::endl;
    KRATOS_CHECK_NEAR(D(0, 2), -7.946736751223777e-01, 1e-12);
    KRATOS_CHECK_NEAR(D(1, 2), 2.147821365647794e-02, 1e-12);
    KRATOS_CHECK_NEAR(D(2, 0), -3.309236751223776e-01, 1e-12);
    KRATOS_CHECK_NEAR(D(2, 1), 3.439488032314460e-02, 1e-12);
    KRATOS_CHECK_EQUAL(residual.size(), 3);
    // std::cout<<std::scientific<<std::setprecision(15)<<residual(0)<<std::endl;
    // std::cout<<std::scientific<<std::setprecision(15)<<residual(1)<<std::endl;
    // std::cout<<std::scientific<<std::setprecision(15)<<residual(2)<<std::endl;

    KRATOS_CHECK_NEAR(residual(0), -5.660890252604245e-01, 1e-12);
    KRATOS_CHECK_NEAR(residual(1), -5.977480398818545e-02, 1e-12);
    KRATOS_CHECK_NEAR(residual(2), -4.612108612915561e-02, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmOmega2D3N_CalculateMassMatrix, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmOmegaElement2D3N_SetUp(model_part);
    EvmOmegaElement2D3N_AssignTestData(model_part);
    auto& r_element = model_part.Elements().front();
    // Test:
    Matrix M;
    r_element.CalculateMassMatrix(M, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(M.size1(), 3);
    KRATOS_CHECK_EQUAL(M.size2(), 3);
    // std::cout<<std::scientific<<std::setprecision(15)<< M(0, 2)<<std::endl;
    
    // std::cout<<std::scientific<<std::setprecision(15)<<M(1, 2)<<std::endl;
    
    // std::cout<<std::scientific<<std::setprecision(15)<<M(2, 0)<<std::endl;

    // std::cout<<std::scientific<<std::setprecision(15)<<M(2, 1)<<std::endl;
    KRATOS_CHECK_NEAR(M(0, 2), 1.261323485870145e-01, 1e-12);
    KRATOS_CHECK_NEAR(M(1, 2), -2.815106894849080e-02, 1e-12);
    KRATOS_CHECK_NEAR(M(2, 0), -5.912919589821915e-02, 1e-12);
    KRATOS_CHECK_NEAR(M(2, 1), -3.303130474714071e-02, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansEvmOmega2D3N_CalculateDampingMatrix, KratosRansFastSuite_1)
{
    // Setup:
    Model model;
    auto& model_part = model.CreateModelPart("test");
    EvmOmegaElement2D3N_SetUp(model_part);
    EvmOmegaElement2D3N_AssignTestData(model_part);
    auto& r_element = model_part.Elements().front();
    // Test:
    Matrix D;
    r_element.CalculateDampingMatrix(D, model_part.GetProcessInfo());
    KRATOS_CHECK_EQUAL(D.size1(), 3);
    KRATOS_CHECK_EQUAL(D.size2(), 3);
    // std::cout<<std::scientific<<std::setprecision(15)<< D(0, 2)<<std::endl;
    
    // std::cout<<std::scientific<<std::setprecision(15)<<D(1, 2)<<std::endl;
    
    // std::cout<<std::scientific<<std::setprecision(15)<<D(2, 0)<<std::endl;

    // std::cout<<std::scientific<<std::setprecision(15)<<D(2, 1)<<std::endl;
    KRATOS_CHECK_NEAR(D(0, 2), -7.946736751223777e-01, 1e-12);
    KRATOS_CHECK_NEAR(D(1, 2), 2.147821365647794e-02, 1e-12);
    KRATOS_CHECK_NEAR(D(2, 0), -3.309236751223776e-01, 1e-12);
    KRATOS_CHECK_NEAR(D(2, 1), 3.439488032314460e-02, 1e-12);
}

} // namespace Testing
} // namespace Kratos.
