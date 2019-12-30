//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez, based on Iñigo Lopez and Riccardo Rossi work
//
#include "embedded_incompressible_potential_flow_element.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    const EmbeddedIncompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);
    const int kutta = r_this.GetValue(KUTTA);

    BoundedVector<double,NumNodes> distances;
    for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
        distances[i_node] = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    }
    const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<Dim,NumNodes>(distances);

    if (is_embedded && wake == 0)
        CalculateEmbeddedLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    else
        BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateEmbeddedLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);
    rLeftHandSideMatrix.clear();

    array_1d<double, NumNodes> potential;
    Vector distances(NumNodes);
    const EmbeddedIncompressiblePotentialFlowElement& r_this = *this;
    const int kutta = r_this.GetValue(KUTTA);
    for(unsigned int i_node = 0; i_node<NumNodes; i_node++) {
        if ((kutta == 1) && (this->GetGeometry()[i_node].GetValue(TRAILING_EDGE)))
            distances(i_node) = this->GetGeometry()[i_node].GetValue(TEMPERATURE);
        else
            distances(i_node) = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    }

    // Vector nodal_distances(NumNodes);
    // for(unsigned int i_node = 0; i_node<NumNodes; i_node++) {
    //     nodal_distances(i_node) = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    // }
    bool is_te = false;
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        if (this->GetGeometry()[i].GetValue(TRAILING_EDGE))
        {
            is_te = true;
        }
    }

    potential = PotentialFlowUtilities::GetPotentialOnNormalElement<Dim, NumNodes>(*this);

    ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(distances);
    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;
    pModifiedShFunc -> ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::GI_GAUSS_1);

    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    // Computing Normal
    std::vector<Vector> cut_normal;
    pModifiedShFunc -> ComputePositiveSideInterfaceAreaNormals(cut_normal,GeometryData::GI_GAUSS_1);
    double norm_normal = sqrt(inner_prod(cut_normal[0],cut_normal[0]));
    auto unit_normal = cut_normal[0]/norm_normal;
    this->SetValue(VELOCITY_LOWER,unit_normal);

    BoundedMatrix<double, 2, 1 > n_kutta;
    n_kutta(0,0)=cut_normal[0][0]/norm_normal;
    n_kutta(1,0)=cut_normal[0][1]/norm_normal;
    BoundedMatrix<double, 2, 1 > n_angle;

    double angle_in_deg = rCurrentProcessInfo[ROTATION_ANGLE];
    n_angle(0,0)=sin(angle_in_deg*Globals::Pi/180);
    n_angle(1,0)=cos(angle_in_deg*Globals::Pi/180);
    double projection = -n_kutta(0,0)*n_angle(1,0)+n_kutta(1,0)*n_angle(0,0);

    bool is_neighbour = false;
    for (unsigned int i = 0; i < NumNodes; ++i){
        const GlobalPointersVector<Element>& r_node_elem_candidates = this -> GetGeometry()[i].GetValue(NEIGHBOUR_ELEMENTS);
        for (std::size_t j = 0; j < r_node_elem_candidates.size(); j++) {
            auto r_elem = r_node_elem_candidates(j);
            if (r_elem->Is(STRUCTURE))
                is_neighbour = true;
        }
    }
    bool is_projection = false;
    if (std::abs(projection)>0.5 && is_te)
    // if ((std::abs(projection)>0.5 && is_te) || (is_te))
    // if ((std::abs(projection)>0.5 && is_te) || (is_neighbour && kutta == 1))
    // if ((std::abs(projection)>0.5 && is_te) || (is_neighbour && kutta == 0))
    // if ((std::abs(projection)>0.5 && is_te) || (is_neighbour))
    // if ((std::abs(projection)>0.5 && is_te) || (is_neighbour) || (kutta==1 && is_te))
    // if ((std::abs(projection)>0.5 && is_te) || (kutta==1 && is_te))
    {
        KRATOS_WATCH(this->Id())
        KRATOS_WATCH(projection)
        n_kutta = n_angle;
        is_projection = true;
    }

    // }else{
    // }

    BoundedMatrix<double,NumNodes,Dim> DN_DX;
    BoundedMatrix<double, NumNodes, NumNodes> lhs_kutta = ZeroMatrix(NumNodes, NumNodes);
    for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
        DN_DX=positive_side_sh_func_gradients(i_gauss);
        Matrix test=prod(DN_DX,n_kutta);
        noalias(rLeftHandSideMatrix) += free_stream_density*prod(DN_DX,trans(DN_DX))*positive_side_weights(i_gauss);;
        noalias(lhs_kutta) += free_stream_density * positive_side_weights(i_gauss) * prod(test,trans(test));
    }


    auto penalty = rCurrentProcessInfo[INITIAL_PENALTY];
    // if (kutta==1){
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        if (this->GetGeometry()[i].GetValue(TRAILING_EDGE))
        // if (this->GetGeometry()[i].GetValue(TRAILING_EDGE))
        {
            for (unsigned int j = 0; j < NumNodes; ++j)
            {
                // rLeftHandSideMatrix(i, j) = lhs_kutta(i, j);
                rLeftHandSideMatrix(i, j) += penalty*lhs_kutta(i, j);
            }
        }

        if (is_projection && (this->GetGeometry()[i].GetValue(WING_TIP)   ))
        {
            for (unsigned int j = 0; j < NumNodes; ++j)
            {
                // rLeftHandSideMatrix(i, j) = lhs_kutta(i, j);
                rLeftHandSideMatrix(i, j) += penalty*lhs_kutta(i, j);
            }
        }
    // }
    }
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, potential);
}

template <>
ModifiedShapeFunctions::Pointer EmbeddedIncompressiblePotentialFlowElement<2,3>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}

template <>
ModifiedShapeFunctions::Pointer EmbeddedIncompressiblePotentialFlowElement<3,4>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Generic geometry check
    int out = BaseType::Check(rCurrentProcessInfo);
    if (out != 0)
    {
        return out;
    }

    for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(GEOMETRY_DISTANCE,this->GetGeometry()[i]);
    }

    return out;

    KRATOS_CATCH("");
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int Dim, int NumNodes>
std::string EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "EmbeddedIncompressiblePotentialFlowElement #" << this->Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EmbeddedIncompressiblePotentialFlowElement #" << this->Id();
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    this->pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class EmbeddedIncompressiblePotentialFlowElement<2, 3>;
template class EmbeddedIncompressiblePotentialFlowElement<3, 4>;


} // namespace Kratos
