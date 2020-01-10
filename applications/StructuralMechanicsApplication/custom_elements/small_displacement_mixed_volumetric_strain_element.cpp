// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/checks.h"

// Application includes
#include "custom_elements/small_displacement_mixed_volumetric_strain_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

Element::Pointer SmallDisplacementMixedVolumetricStrainElement::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainElement>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedVolumetricStrainElement::Clone(
    IndexType NewId,
    NodesArrayType const &rThisNodes) const
{
    KRATOS_TRY

    SmallDisplacementMixedVolumetricStrainElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedVolumetricStrainElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::EquationIdVector(
    EquationIdVectorType &rResult,
    ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const IndexType n_nodes = r_geometry.PointsNumber();
    const IndexType dim = r_geometry.WorkingSpaceDimension();
    const IndexType dof_size = n_nodes * (dim + 1);

    if (rResult.size() != dof_size)
    {
        rResult.resize(dof_size);
    }

    const IndexType disp_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType eps_vol_pos = r_geometry[0].GetDofPosition(VOLUMETRIC_STRAIN);

    IndexType aux_index = 0;
    if (dim == 2)
    {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node)
        {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN, eps_vol_pos).EquationId();
        }
    }
    else
    {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node)
        {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Z, disp_pos + 2).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN, eps_vol_pos).EquationId();
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::GetDofList(
    DofsVectorType &rElementalDofList,
    ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const IndexType n_nodes = r_geometry.PointsNumber();
    const IndexType dim = r_geometry.WorkingSpaceDimension();
    const IndexType dof_size = n_nodes * (dim + 1);

    if (rElementalDofList.size() != dof_size)
    {
        rElementalDofList.resize(dof_size);
    }

    if (dim == 2)
    {
        for (IndexType i = 0; i < n_nodes; ++i)
        {
            rElementalDofList[i * (dim + 1)] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * (dim + 1) + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * (dim + 1) + 2] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
        }
    }
    else if (dim == 3)
    {
        for (IndexType i = 0; i < n_nodes; ++i)
        {
            rElementalDofList[i * (dim + 1)] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * (dim + 1) + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * (dim + 1) + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[i * (dim + 1) + 3] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::Initialize()
{
    KRATOS_TRY

    // Integration method initialization
    mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
    const auto &r_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    // Constitutive Law Vector initialisation
    if (mConstitutiveLawVector.size() != r_integration_points.size())
    {
        mConstitutiveLawVector.resize(r_integration_points.size());
    }

    // Initialize material
    InitializeMaterial();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto &r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node)
    {
        const auto &r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d)
        {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Set te constitutive law values
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss)
    {
        // Recompute the kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Set the constitutive variables
        SetConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->InitializeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto &r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node)
    {
        const auto &r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d)
        {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Set te constitutive law values
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss)
    {
        // Recompute the kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Set the constitutive variables
        SetConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->FinalizeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size)
    {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size)
    {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    
    bool compute_lhs = true;
    CalculateLocalSystemAuxiliary(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo, compute_lhs);
}





void SmallDisplacementMixedVolumetricStrainElement::CalculateLocalSystemAuxiliary(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo,
    bool compute_lhs
    )
{
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const SizeType n_nodes = GetGeometry().PointsNumber();
    const SizeType block_size = dim + 1;

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node)
    {
        const auto &r_disp = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d)
        {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = GetGeometry()[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Set auxiliary arrays
    Vector voigt_identity = ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d)
    {
        voigt_identity[d] = 1.0;
    }

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    Matrix Bi(strain_size, dim);
    Matrix Bj(strain_size, dim);
    Vector Gi(dim);
    Vector Gj(dim);
    

    Matrix Luu(dim, dim);
    Vector Lue(dim);
    Vector Leu(dim);
    double Lee = 0.0;

    const SizeType n_gauss = GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
    const auto &r_integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
    {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, GetGeometry().IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        Vector C_m = prod(constitutive_variables.D, voigt_identity);

        const auto body_3 = GetBodyForce(GetGeometry().IntegrationPoints(GetIntegrationMethod()), i_gauss);
        Vector body_force(dim);
        for(IndexType l = 0; l < dim; ++l) 
            body_force[l] = body_3[l];

        // Calculate the linearised bulk and shear moduli
        const double bulk_modulus = CalculateLinearisedBulkModulus(constitutive_variables);
        const double shear_modulus = CalculateLinearisedShearModulus(constitutive_variables);

        //tau_2 stab constant
        const double tau_2 = 0.0 * std::min(0.01, 4.0 * shear_modulus / bulk_modulus);

        // Calculate tau_1 stabilization constant
        Matrix aux = ZeroMatrix(dim, dim);
        Matrix Ddev = prod(kinematic_variables.DevStrainOp, constitutive_variables.D);
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node)
        {
            for (IndexType k = 0; k < strain_size; ++k)
                for (IndexType l = 0; l < dim; ++l)
                    Bi(k, l) = kinematic_variables.B(k, i_node * dim + l);

            aux += prod(trans(Bi), Matrix(prod(constitutive_variables.D, Bi)));
            //aux += prod(trans(Bi), Matrix(prod(Ddev, Bi)));
        }
        double det;
        Matrix tau_1(dim, dim);
        MathUtils<double>::InvertMatrix(aux, tau_1, det);

        //TODO: remove!
        double h = sqrt(GetGeometry().Area());
        tau_1 = 0.0* h*h/(shear_modulus + bulk_modulus) * identity_matrix<double>(dim,dim);
        

        Vector tau1_body_force = prod(tau_1, body_force);

        for (IndexType i = 0; i < n_nodes; ++i)
        {
            IndexType iblock = i * block_size;

            noalias(Gi) = row(kinematic_variables.DN_DX, i);
            const double Ni = kinematic_variables.N[i];

            for (IndexType k = 0; k < strain_size; ++k)
                for (IndexType l = 0; l < dim; ++l)
                    Bi(k, l) = kinematic_variables.B(k, i * dim + l);

            double eps_v_gauss = 0.0;
            double div_u_gauss = 0.0;

            for (IndexType j = 0; j < n_nodes; ++j)
            {
                eps_v_gauss += kinematic_variables.N[j] * kinematic_variables.VolumetricNodalStrains[j];
                for (IndexType d = 0; d < dim; ++d)
                    div_u_gauss += kinematic_variables.DN_DX(j, d) * kinematic_variables.Displacements(j * dim + d);
            }

            double eps_minus_div = eps_v_gauss - div_u_gauss;

            Vector Ri_disp = Ni*body_force;
            noalias(Ri_disp) -= prod(trans(Bi), constitutive_variables.StressVector);
            noalias(Ri_disp) += (eps_minus_div / dim) * prod(trans(Bi), C_m);
            noalias(Ri_disp) -= ( (1-tau_2) * bulk_modulus * eps_minus_div) * Gi;

            double Ri_mass = (1-tau_2) * bulk_modulus * eps_minus_div * Ni;
            Ri_mass -= bulk_modulus * inner_prod(Gi, tau1_body_force);

            for (IndexType j = 0; j < n_nodes; ++j)
            { 
                IndexType jblock = j * block_size;

                noalias(Gj) = row(kinematic_variables.DN_DX, j);
                const double Nj = kinematic_variables.N[j];

                Vector tau_1_Gj = prod(tau_1, Gj);
                Ri_mass += bulk_modulus * bulk_modulus * kinematic_variables.VolumetricNodalStrains[j] *inner_prod(Gi,  tau_1_Gj); //todo OPTIMIZE

                if(compute_lhs)
                {
                    for (IndexType k = 0; k < strain_size; ++k)
                        for (IndexType l = 0; l < dim; ++l)
                            Bj(k, l) = kinematic_variables.B(k, j * dim + l);

                    noalias(Luu) = prod(trans(Bi), Matrix( prod(constitutive_variables.D, Bj ) ) );
                    noalias(Luu) -= (1-tau_2) * bulk_modulus * outer_prod(Gi, Gj);
                    
                    noalias(Lue) = (1-tau_2) * bulk_modulus * Nj * Gi;
                    noalias(Leu) = (1-tau_2) * bulk_modulus * Ni * Gj;

                    Lee =  -(1-tau_2) * bulk_modulus * Ni * Nj;
                    Lee -= bulk_modulus * bulk_modulus * inner_prod(Gi, tau_1_Gj); //todo OPTIMIZE

                    //sum LHS into the correct positions of the matrix and multiply by the integration weight
                    rLeftHandSideMatrix(iblock + dim  , jblock + dim ) += w_gauss * Lee; //eps-eps block

                    for (IndexType k = 0; k < dim; ++k) 
                    {
                        rLeftHandSideMatrix(iblock + k    , jblock + dim ) += w_gauss * Lue(k); // u - eps block
                        rLeftHandSideMatrix(iblock + dim  , jblock + k   ) += w_gauss * Leu(k); // eps - u block
                    }
                    
                    for (IndexType k = 0; k < dim; ++k) 
                        for (IndexType l = 0; l < dim; ++l)
                            rLeftHandSideMatrix(iblock + k    , jblock + l   ) += w_gauss * Luu(k, l); //uu block

                }
            }

            //sum RHS into the correct positions of the matrix and multiply by the integration weight
            for (IndexType k = 0; k < dim; ++k)
                rRightHandSideVector[iblock+k   ] += w_gauss * Ri_disp[k];

            rRightHandSideVector[iblock+dim ] += w_gauss * Ri_mass;
        }
    }

    if(Id() == 1)
    {
        KRATOS_WATCH(rRightHandSideVector)
        KRATOS_WATCH(rLeftHandSideMatrix)
    }


}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix,
    ProcessInfo &rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size)
    {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node)
    {
        const auto &r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d)
        {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Set auxiliary arrays
    Vector voigt_identity = ZeroVector(strain_size);
    for (IndexType d = 0; d < dim; ++d)
    {
        voigt_identity[d] = 1.0;
    }

    // Calculate the LHS contributions
    rLeftHandSideMatrix.clear();

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto &r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
    {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_geometry.IntegrationPoints(this->GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Calculate the linearised bulk modulus
        const double bulk_modulus = CalculateLinearisedBulkModulus(constitutive_variables);

        // Calculate tau_1 stabilization constant
        Matrix B_i(strain_size, dim);
        Matrix aux = ZeroMatrix(dim, dim);
        Matrix Ddev = prod(kinematic_variables.DevStrainOp, constitutive_variables.D);
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node)
        {
            for (IndexType k = 0; k < strain_size; ++k)
            {
                for (IndexType l = 0; l < dim; ++l)
                {
                    B_i(k, l) = kinematic_variables.B(k, i_node * dim + l);
                }
            }
            aux += prod(trans(B_i), Matrix(prod(Ddev, B_i)));
        }
        double det;
        Matrix tau_1_mat_w_k(dim, dim);
        MathUtils<double>::InvertMatrix(aux, tau_1_mat_w_k, det);
        tau_1_mat_w_k *= (bulk_modulus * w_gauss);

        // Calculate tau_2 stabilization constant
        const double shear_modulus = CalculateLinearisedShearModulus(constitutive_variables);
        const double tau_2 = std::min(1.0e-2, 4.0 * shear_modulus / bulk_modulus);
        
        // Calculate and add the LHS contributions
        const double w_1_tau_2_k = w_gauss * (1 - tau_2) * bulk_modulus;

        const Matrix transB_C = prod(trans(kinematic_variables.B), cons_law_values.GetConstitutiveMatrix());
        const Matrix transB_C_B = prod(transB_C, kinematic_variables.B);

        const Vector C_m_voigt = prod(cons_law_values.GetConstitutiveMatrix(), voigt_identity);
        const Matrix C_m = MathUtils<double>::StressVectorToTensor(C_m_voigt);

        for (IndexType i = 0; i < n_nodes; ++i)
        {
            const Vector G_I = row(kinematic_variables.DN_DX, i);
            const Vector tau_1_G_I_C_m = prod(G_I, Matrix(prod(tau_1_mat_w_k, C_m)));
            for (IndexType j = 0; j < n_nodes; ++j)
            {
                const Vector G_J = row(kinematic_variables.DN_DX, j);
                const Matrix G_I_transG_J = outer_prod(G_I, G_J);
                // Add mass conservation volumetric strain contribution
                double &r_LHS_mass_mass = rLeftHandSideMatrix(i * block_size + dim, j * block_size + dim);
                r_LHS_mass_mass -= w_1_tau_2_k * kinematic_variables.N[i] * kinematic_variables.N[j];
                r_LHS_mass_mass -= inner_prod(tau_1_G_I_C_m, G_J);
                for (IndexType d = 0; d < dim; ++d)
                {
                    // Add mass conservation displacement divergence contribution
                    rLeftHandSideMatrix(i * block_size + dim, j * block_size + d) += w_1_tau_2_k * kinematic_variables.N[i] * G_J(d);
                    // Add momentum internal force contribution
                    for (IndexType d2 = 0; d2 < dim; ++d2)
                    {
                        double &r_LHS_mom_mom = rLeftHandSideMatrix(i * block_size + d, j * block_size + d2);
                        r_LHS_mom_mom += w_gauss * transB_C_B(i * dim + d, j * dim + d2);
                        r_LHS_mom_mom -= w_1_tau_2_k * G_I_transG_J(d, d2);
                    }
                    // Add momentum volumetric strain contribution
                    rLeftHandSideMatrix(i * block_size + d, j * block_size + dim) += w_1_tau_2_k * G_I(d) * kinematic_variables.N(j);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateRightHandSide(
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size)
    {
        rRightHandSideVector.resize(matrix_size, false);
    }
    rRightHandSideVector.clear();

    Matrix dummy(0,0);
    bool compute_lhs = false;
    CalculateLocalSystemAuxiliary(dummy,rRightHandSideVector,rCurrentProcessInfo, compute_lhs);


 }

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::InitializeMaterial()
{
    KRATOS_TRY

    const auto &r_properties = GetProperties();
    if (r_properties[CONSTITUTIVE_LAW] != nullptr)
    {
        const auto &r_geometry = GetGeometry();
        const auto &r_N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        IndexType aux = 0;
        for (auto &it_gauss_pt : mConstitutiveLawVector)
        {
            it_gauss_pt = (r_properties[CONSTITUTIVE_LAW])->Clone();
            (it_gauss_pt)->InitializeMaterial(r_properties, r_geometry, row(r_N_values, aux));
            aux++;
        }
    }
    else
    {
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacementMixedVolumetricStrainElement::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::SetConstitutiveVariables(
    KinematicVariables &rThisKinematicVariables,
    ConstitutiveVariables &rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters &rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType &IntegrationPoints)
{
    // Here we essentially set the input parameters
    rValues.SetStrainVector(rThisKinematicVariables.EquivalentStrain); // equivalent total strain
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N);        // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF);             // assuming that det(F) is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F);        // assuming that F is computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D);      //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateConstitutiveVariables(
    KinematicVariables &rThisKinematicVariables,
    ConstitutiveVariables &rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters &rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType &IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure)
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> SmallDisplacementMixedVolumetricStrainElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber) const
{
    return StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateKinematicVariables(
    KinematicVariables &rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod &rIntegrationMethod) const
{
    const auto &r_geometry = GetGeometry();
    const auto &r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);

    // Shape functions
    rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    // Calculate the inverse Jacobian
    GeometryUtils::JacobianOnInitialConfiguration(
        r_geometry,
        r_integration_points[PointNumber],
        rThisKinematicVariables.J0);
    MathUtils<double>::InvertMatrix(
        rThisKinematicVariables.J0,
        rThisKinematicVariables.InvJ0,
        rThisKinematicVariables.detJ0);
    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0)
        << "Element ID: " << this->Id() << " is inverted. det(J0) = " << rThisKinematicVariables.detJ0 << std::endl;

    // Calculate the shape functions gradients
    GeometryUtils::ShapeFunctionsGradients(
        r_geometry.ShapeFunctionsLocalGradients(rIntegrationMethod)[PointNumber],
        rThisKinematicVariables.InvJ0,
        rThisKinematicVariables.DN_DX);

    // Compute B
    CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.DN_DX);

    // Calculate the equivalent total strain
    CalculateEquivalentStrain(rThisKinematicVariables);

    // Compute equivalent F
    ComputeEquivalentF(rThisKinematicVariables.F, rThisKinematicVariables.EquivalentStrain);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateB(
    Matrix &rB,
    const Matrix &rDN_DX) const
{
    KRATOS_TRY;

    StructuralMechanicsElementUtilities::CalculateB(*this, rDN_DX, rB);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateEquivalentStrain(KinematicVariables &rThisKinematicVariables) const
{
    const auto &r_geom = GetGeometry();
    const SizeType n_nodes = r_geom.PointsNumber();
    const SizeType dim = r_geom.WorkingSpaceDimension();

    // Add the deviatoric contribution to the equivalent strain
    const Vector total_strain = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);
    rThisKinematicVariables.EquivalentStrain = prod(rThisKinematicVariables.DevStrainOp, total_strain);

    // Interpolate and add the nodal volumetric strain contribution
    double gauss_vol_strain = 0.0;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node)
    {
        gauss_vol_strain += rThisKinematicVariables.N[i_node] * rThisKinematicVariables.VolumetricNodalStrains[i_node];
    }
    for (IndexType d = 0; d < dim; ++d)
    {
        rThisKinematicVariables.EquivalentStrain[d] += (1.0 / dim) * gauss_vol_strain;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::ComputeEquivalentF(
    Matrix &rF,
    const Vector &rStrainTensor) const
{
    StructuralMechanicsElementUtilities::ComputeEquivalentF(*this, rStrainTensor, rF);
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateElementSize(const KinematicVariables &rThisKinematicVariables) const
{
    const auto &r_geometry = GetGeometry();
    switch (r_geometry.GetGeometryType())
    {
    case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
    {
        const BoundedMatrix<double, 3, 2> aux_grads(rThisKinematicVariables.DN_DX);
        return ElementSizeCalculator<2, 3>::GradientsElementSize(aux_grads);
    }
    case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
    {
        const BoundedMatrix<double, 4, 3> aux_grads(rThisKinematicVariables.DN_DX);
        return ElementSizeCalculator<3, 4>::GradientsElementSize(aux_grads);
    }
    case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4:
    {
        return ElementSizeCalculator<2, 4>::AverageElementSize(r_geometry);
    }
    case GeometryData::KratosGeometryType::Kratos_Hexahedra3D8:
    {
        return ElementSizeCalculator<3, 8>::AverageElementSize(r_geometry);
    }
    default:
    {
        KRATOS_ERROR << "Asking for a non-implemented geometry element size calculation." << std::endl;
    }
    }
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateLinearisedBulkModulus(const ConstitutiveVariables &rThisConstitutiveVariables) const
{
    const auto &r_geom = GetGeometry();
    const SizeType dim = r_geom.WorkingSpaceDimension();

    double bulk_modulus = 0.0;
    for (IndexType i = 0; i < dim; ++i)
    {
        for (IndexType j = 0; j < dim; ++j)
        {
            bulk_modulus += rThisConstitutiveVariables.D(i, j);
        }
    }

    return bulk_modulus / std::pow(dim, 2);
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedVolumetricStrainElement::CalculateLinearisedShearModulus(const ConstitutiveVariables &rThisConstitutiveVariables) const
{
    const auto &r_geom = GetGeometry();
    const SizeType dim = r_geom.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto &rD = rThisConstitutiveVariables.D;

    double tmp=0;
    double shear_modulus = (3*(rD(0,0)+rD(1,1) - rD(0,1)) + 9 * rD(2,2) ) /30.0; //modified from https://arxiv.org/pdf/cond-mat/0509705.pdf
    // double shear_modulus = std::numeric_limits<double>::max();
    
    // for (unsigned int i = dim; i < strain_size; ++i)
    // {
    //     shear_modulus = std::min(shear_modulus, rD(i, i));
    // }

    return shear_modulus;
}

/***********************************************************************************/
/***********************************************************************************/

int SmallDisplacementMixedVolumetricStrainElement::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    int check = SmallDisplacementMixedVolumetricStrainElement::BaseType::Check(rCurrentProcessInfo);

    // Base check
    
    check = StructuralMechanicsElementUtilities::SolidElementCheck(*this, rCurrentProcessInfo, mConstitutiveLawVector);

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(VOLUMETRIC_STRAIN)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const auto &r_geometry = this->GetGeometry();
    for (IndexType i = 0; i < r_geometry.size(); i++)
    {
        const NodeType &r_node = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUMETRIC_STRAIN, r_node)

        KRATOS_CHECK_DOF_IN_NODE(VOLUMETRIC_STRAIN, r_node)
    }

    return check;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const auto &r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss)
    {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has(rVariable))
    {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    }
    else
    {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::CalculateOnIntegrationPoints(
    const Variable<Vector> &rVariable,
    std::vector<Vector> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const auto &r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss)
    {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has(rVariable))
    {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    }
    else if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR)
    {
        // Create and initialize element variables:
        const SizeType n_nodes = r_geometry.PointsNumber();
        const SizeType dim = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node)
        {
            const auto &r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            for (IndexType d = 0; d < dim; ++d)
            {
                kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
            }
            kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
        }

        // Create the constitutive variables and values containers
        ConstitutiveVariables constitutive_variables(strain_size);
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
        auto &r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss)
        {
            // Calculate kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Call the constitutive law to update material variables
            if (rVariable == CAUCHY_STRESS_VECTOR)
            {
                // Compute material reponse
                CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
            }
            else
            {
                // Compute material reponse
                CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_PK2);
            }

            // Check sizes and save the output stress
            if (rOutput[i_gauss].size() != strain_size)
            {
                rOutput[i_gauss].resize(strain_size, false);
            }
            rOutput[i_gauss] = constitutive_variables.StressVector;
        }
    }
    else
    {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::GetValueOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::GetValueOnIntegrationPoints(
    const Variable<Vector> &rVariable,
    std::vector<Vector> &rValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::save(Serializer &rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainElement::BaseType);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod", IntMethod);
    rSerializer.save("mConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedVolumetricStrainElement::load(Serializer &rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementMixedVolumetricStrainElement::BaseType);
    int IntMethod;
    rSerializer.load("IntegrationMethod", IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

} // Namespace Kratos
