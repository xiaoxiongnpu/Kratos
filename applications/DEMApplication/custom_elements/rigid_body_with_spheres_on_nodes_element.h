// Created by: Miguel Angel Celigueta, maceli@cimne.upc.edu

#if !defined KRATOS_RIGID_BODY_WITH_SPHERES_ON_NODES_ELEMENT_H_INCLUDED
#define KRATOS_RIGID_BODY_WITH_SPHERES_ON_NODES_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "rigid_body_element.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) RigidBodyWithSpheresOnNodesElement3D : public RigidBodyElement3D {

    public:
        /// Pointer definition of RigidBodyWithSpheresOnNodesElement3D
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(RigidBodyWithSpheresOnNodesElement3D);

        RigidBodyWithSpheresOnNodesElement3D();
        RigidBodyWithSpheresOnNodesElement3D(IndexType NewId, GeometryType::Pointer pGeometry);
        RigidBodyWithSpheresOnNodesElement3D(IndexType NewId, NodesArrayType const& ThisNodes);
        RigidBodyWithSpheresOnNodesElement3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor
        virtual ~RigidBodyWithSpheresOnNodesElement3D();

        void CustomInitialize(ModelPart& rigid_body_element_sub_model_part, ParticleCreatorDestructor::Pointer p_creator_destructor, const bool do_seed_nodes_with_spheres, ModelPart& spheres_model_part) override;

        virtual std::string Info() const override
        {
	    std::stringstream buffer;
	    buffer << "Rigid Body With Spheres On Nodes Element #" << Id();
	    return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const override
        {
	    rOStream << "Rigid Body With Spheres On Nodes Element #" << Id();
        }

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const override
        {
	    //mpGeometry->PrintData(rOStream);
        }

    protected:


    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override;

        virtual void load(Serializer& rSerializer) override;

    }; // Class RigidBodyWithSpheresOnNodesElement3D

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, RigidBodyWithSpheresOnNodesElement3D& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const RigidBodyWithSpheresOnNodesElement3D& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos

#endif // KRATOS_RIGID_BODY_WITH_SPHERES_ON_NODES_ELEMENT_H_INCLUDED  defined
