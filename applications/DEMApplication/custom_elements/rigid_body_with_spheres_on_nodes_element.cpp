// Created by: Miguel Angel Celigueta, maceli@cimne.upc.edu

// System includes

// Project includes
#include "rigid_body_with_spheres_on_nodes_element.h"

namespace Kratos {

    RigidBodyWithSpheresOnNodesElement3D::RigidBodyWithSpheresOnNodesElement3D()
    : RigidBodyElement3D() {
    }

    RigidBodyWithSpheresOnNodesElement3D::RigidBodyWithSpheresOnNodesElement3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : RigidBodyElement3D(NewId, pGeometry) {
    }

    RigidBodyWithSpheresOnNodesElement3D::RigidBodyWithSpheresOnNodesElement3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : RigidBodyElement3D(NewId, pGeometry, pProperties) {
    }

    RigidBodyWithSpheresOnNodesElement3D::RigidBodyWithSpheresOnNodesElement3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : RigidBodyElement3D(NewId, ThisNodes) {
    }

    Element::Pointer RigidBodyWithSpheresOnNodesElement3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
        return Element::Pointer(new RigidBodyWithSpheresOnNodesElement3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    // Destructor
    RigidBodyWithSpheresOnNodesElement3D::~RigidBodyWithSpheresOnNodesElement3D() {
    }

    void RigidBodyWithSpheresOnNodesElement3D::CustomInitialize(ModelPart& rigid_body_element_sub_model_part, ParticleCreatorDestructor::Pointer p_creator_destructor, const bool do_seed_nodes_with_spheres, ModelPart& spheres_model_part) {
        RigidBodyElement3D::CustomInitialize(rigid_body_element_sub_model_part, p_creator_destructor, do_seed_nodes_with_spheres, spheres_model_part);
        if(do_seed_nodes_with_spheres) {
            const double radius = rigid_body_element_sub_model_part[RADIUS];
            Properties::Pointer props = rigid_body_element_sub_model_part.pGetProperties(1);
            for (size_t k = 0; k < mListOfNodes.size(); k++) {
                NodesArrayType::iterator node_k = mListOfNodes.begin() + k;
                const array_1d<double, 3>& coords = node_k->Coordinates();
                p_creator_destructor->CreateSphericParticle(spheres_model_part, coords, props, radius, "PolyhedronSkinParticle3D");
            }
        }
    }

    void RigidBodyWithSpheresOnNodesElement3D::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void RigidBodyWithSpheresOnNodesElement3D::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

} // namespace Kratos