//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter
//
//

#ifndef CABLE_NET_MPC_PROCESS_H
#define CABLE_NET_MPC_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "spatial_containers/spatial_containers.h"


// Application includes
#include "custom_processes/apply_multi_point_constraints_process.h"

namespace Kratos
{

class CableNetMpcProcess : public ApplyMultipointConstraintsProcess
{
  public:

    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef std::vector<NodeTypePointer>::iterator NodeIterator;
    typedef std::vector<double> DoubleVector;
    typedef DoubleVector::iterator DoubleVectorIterator;
    typedef std::size_t SizeType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;
    

    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(CableNetMpcProcess);

    /// Constructor.
    CableNetMpcProcess(ModelPart &model_part,
                                      Parameters rParameters) : ApplyMultipointConstraintsProcess(model_part,rParameters)
    {}


    /**
     * @brief This function finds neighbor nodes for each slave node and couples the nodes
     *        this is the main function of this process
     */
    void CoupleModelParts()
    {
        ModelPart &master_model_part    = mr_model_part.GetSubModelPart(m_parameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part     = mr_model_part.GetSubModelPart(m_parameters["slave_sub_model_part_name"].GetString());
        const double neighbor_search_radius      = m_parameters["neighbor_search_radius"].GetDouble();
        const int bucket_size           =  m_parameters["bucket_size"].GetInt();
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

        NodeVector master_node_list(r_nodes_master.size());
        this->CreateListOfNodesOfMasterSubModelPart(master_node_list);

        KDTree::Pointer search_tree =
         Kratos::shared_ptr<KDTree>(new KDTree(master_node_list.begin(), master_node_list.end(), bucket_size));


        const int max_number_of_neighbors = 2;
        for(NodeType& node_i : r_nodes_slave)
        {

            //1.) find nodal neighbors
            NodeVector neighbor_nodes( max_number_of_neighbors );
            DoubleVector resulting_squared_distances( max_number_of_neighbors );

            SizeType number_of_neighbors = search_tree->SearchInRadius( node_i,
                                                                    neighbor_search_radius,
                                                                    neighbor_nodes.begin(),
                                                                    resulting_squared_distances.begin(),
                                                                    max_number_of_neighbors );

            
            DoubleVector list_of_weights( number_of_neighbors, 0.0 );

            this->CalculateNodalWeights(resulting_squared_distances,list_of_weights,number_of_neighbors);
            this->CoupleSlaveToNeighborMasterNodes(node_i,neighbor_nodes,list_of_weights,number_of_neighbors);

            this->SetmIsInitialized(true);

            //DoubleVector list_of_weights2( number_of_neighbors, 0.0 );
            //test new function to calculate weight
            //this->ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights2);
            

            //if(m_parameters["debug_info"].GetBool()) KRATOS_WATCH(list_of_weights);

            //std::cout << "slave: " << node_i.Id() << " has " << number_of_neighbors << " masters " << std::endl;
            //std::cout << "###################################################" << std::endl;
        }
    }


    /**
     * @brief This function couples nodes by calling the parent class functions
     */
    void CoupleSlaveToNeighborMasterNodes(const NodeType& rCurrentSlaveNode,
     const NodeVector& rNeighborNodes, const DoubleVector& rNodalNeighborWeights,
     const SizeType& rNumberOfNeighbors)
    {
        ModelPart &master_model_part    = mr_model_part.GetSubModelPart(m_parameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part     = mr_model_part.GetSubModelPart(m_parameters["slave_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();


        for(SizeType dof_iterator=0;dof_iterator<m_parameters["variable_names"].size();++dof_iterator)
        {
            VariableComponentType current_dof = KratosComponents<VariableComponentType>::Get(m_parameters["variable_names"][dof_iterator].GetString());
            //KRATOS_WATCH(current_dof);

            for(SizeType master_iterator =0;master_iterator<rNumberOfNeighbors;++master_iterator)
            {
                ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodesAndVariableComponents(
                    r_nodes_master[rNeighborNodes[master_iterator]->Id()],current_dof,
                    r_nodes_slave[rCurrentSlaveNode.Id()],current_dof,rNodalNeighborWeights[master_iterator],0);

                if(m_parameters["debug_info"].GetBool()){
                    std::cout << rNeighborNodes[master_iterator]->Id() << "-----" << rCurrentSlaveNode.Id() << "-----" << rNodalNeighborWeights[master_iterator] << std::endl;
                }
                
            } // each master node
        }  // each dof

    }


    /**
     * @brief This function creates a NodeVector of the master nodes to be used for the Kd tree
     */
    void CreateListOfNodesOfMasterSubModelPart(NodeVector& MasterNodeList)
    {
        ModelPart &master_model_part = mr_model_part.GetSubModelPart(m_parameters["master_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes = master_model_part.Nodes();

        MasterNodeList.resize(r_nodes.size());
        auto i_begin = master_model_part.NodesBegin();

        for (SizeType i(0);i<r_nodes.size();++i)
        {
            NodeTypePointer pnode = *((i_begin+i).base());
            MasterNodeList[i] = pnode;
        }
    }


    /**
     * @brief This function re-calculates the weights used in mpc
     */
    void CalculateNodalWeights(const DoubleVector& rResultingSquaredDistances, DoubleVector& rNodalNeighborWeights, const SizeType& rNumberOfNeighbors)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();
        double total_nodal_distance = 0.00;

        if((rNumberOfNeighbors==1)&&(std::abs(rResultingSquaredDistances[0])<numerical_limit))
         {rNodalNeighborWeights[0] = 1.00;}
        else
        {   
            for (SizeType i=0;i<rNumberOfNeighbors;++i) total_nodal_distance+=std::sqrt(rResultingSquaredDistances[i]);   
            for (SizeType i=0;i<rNumberOfNeighbors;++i)
            {
                rNodalNeighborWeights[i] = std::sqrt(rResultingSquaredDistances[rNumberOfNeighbors-(i+1)])/total_nodal_distance;
            }
        }
    }



    void ExecuteInitializeSolutionStep() override
    {
        if (this->GetmIsInitialized()) 
            {if (m_parameters["reform_every_step"].GetBool()) 
                {this->CoupleModelParts();}
            }
        else this->CoupleModelParts();
    }




    /////////////////////////////////////
    /////////----> test functions

    void ComputeWeightForAllNeighbors(  ModelPart::NodeType& design_node,
                                        NodeVector& neighbor_nodes,
                                        SizeType number_of_neighbors,
                                        DoubleVector& list_of_weights)
    {


        double total_length(0.0);
        DoubleVector temp_vector(number_of_neighbors,0.00);
        KRATOS_WATCH(number_of_neighbors);
        KRATOS_WATCH(temp_vector);
        for(SizeType neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            double current_length(this->CalculateCurrentLength(design_node,neighbor_node));

            temp_vector[neighbor_itr] = current_length;
            total_length += current_length;
        }

        for(SizeType i=0;i<number_of_neighbors;++i) temp_vector[i] /= total_length;
        for(SizeType i=0;i<number_of_neighbors;++i) list_of_weights[i] = temp_vector[number_of_neighbors-(i+1)];
    }

    double CalculateCurrentLength(ModelPart::NodeType& rNodeI,ModelPart::NodeType& rNodeJ) {
        const double du =
            rNodeJ.FastGetSolutionStepValue(DISPLACEMENT_X) -
            rNodeI.FastGetSolutionStepValue(DISPLACEMENT_X);
        const double dv =
            rNodeJ.FastGetSolutionStepValue(DISPLACEMENT_Y) -
            rNodeI.FastGetSolutionStepValue(DISPLACEMENT_Y);
        const double dw =
            rNodeJ.FastGetSolutionStepValue(DISPLACEMENT_Z) -
            rNodeI.FastGetSolutionStepValue(DISPLACEMENT_Z);
        const double dx = rNodeJ.X0() - rNodeI.X0();
        const double dy = rNodeJ.Y0() - rNodeI.Y0();
        const double dz = rNodeJ.Z0() - rNodeI.Z0();
        const double l = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy) +
                                    (dw + dz) * (dw + dz));
        return l;
    }
    //<------- ////////////////
    /////////////////////////////////////



    void SetmIsInitialized(const bool& check) {this->mIsInitialized = check;}
    bool GetmIsInitialized() const {return this->mIsInitialized;} 



  protected:

  private:

    bool mIsInitialized = false;

}; // Class 

}; // namespace 

#endif 
