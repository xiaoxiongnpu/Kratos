//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Marc Nunez
//


#include "define_embedded_wake_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"


namespace Kratos
{
// Constructor for DefineEmbeddedWakeProcess Process
DefineEmbeddedWakeProcess::DefineEmbeddedWakeProcess(ModelPart& rModelPart,
                    ModelPart& rWakeModelPart
                ):
    Process(),
    mrModelPart(rModelPart),
    mrWakeModelPart(rWakeModelPart)
{}

void DefineEmbeddedWakeProcess::Execute()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]>2) << "DOMAIN_SIZE is greater than 2. DefineEmbeddedWakeProcess is only implemented for 2D cases!" << std::endl;

    ComputeDistanceToWake();
    MarkWakeElements();
    ComputeTrailingEdgeNode();
    MarkKuttaWakeElements();

    KRATOS_CATCH("");
}


void DefineEmbeddedWakeProcess::ComputeDistanceToWake(){

    CalculateDiscontinuousDistanceToSkinProcess<2> distance_calculator(mrModelPart, mrWakeModelPart);
    distance_calculator.Execute();

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = mrModelPart.NodesBegin() + i;
        it_node->SetValue(WAKE_DISTANCE, 0.0);
    }

    // for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
    //     auto it_node = mrModelPart.NodesBegin() + i;
    //     std::size_t current_id = it_node->Id();
    //     double& distance = it_node->FastGetSolutionStepValue(GEOMETRY_DISTANCE);
    //     if (distance<0.0){
    //         const GlobalPointersVector<Element>& r_node_elem_candidates = it_node -> GetValue(NEIGHBOUR_ELEMENTS);
    //         bool all_positive = true;
    //         for (std::size_t j = 0; j < r_node_elem_candidates.size(); j++) {
    //             auto r_geometry = r_node_elem_candidates(j)->GetGeometry();
    //             for (std::size_t j_node = 0; j_node<r_geometry.size(); j_node++ ){
    //                 double neighbour_distance = r_geometry[j_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE);
    //                 if (neighbour_distance<0.0 && current_id != r_geometry[j_node].Id()){
    //                     all_positive = false;
    //                     break;
    //                 }
    //             }
    //             // if (!all_positive)
    //         }
    //         if (all_positive) {
    //             KRATOS_WATCH("HEY")
    //             distance = -1*distance;
    //         }
    //     }
    // }
    // std::vector<std::vector<std::size_t>> list_of_lists;
    // std::vector<std::size_t> list_of_ids;
    // for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
    //     auto it_node = mrModelPart.NodesBegin() + i;
    //     list_of_ids.push_back(it_node->Id());
    // }

    // std::size_t initial_id = list_of_ids[0];
    // for (auto it=list_of_ids.begin();
    //                           it!=list_of_ids.end();
    //                           /*it++*/) <----------- I commented it.
    // {

    //     if(list_of_ids[it])
    //         it = allPlayers.erase(it);
    //     else
    //         ++it;
    // }

    //     double& distance = it_node->FastGetSolutionStepValue(GEOMETRY_DISTANCE);
    //     if (distance<0.0){
    //         const GlobalPointersVector<Element>& r_node_elem_candidates = it_node -> GetValue(NEIGHBOUR_ELEMENTS);
    //         bool all_positive = true;
    //         for (std::size_t j = 0; j < r_node_elem_candidates.size(); j++) {
    //             auto r_geometry = r_node_elem_candidates(j)->GetGeometry();
    //             for (std::size_t j_node = 0; j_node<r_geometry.size(); j_node++ ){
    //                 double neighbour_distance = r_geometry[j_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE);
    //                 if (neighbour_distance<0.0 && current_id != r_geometry[j_node].Id()){
    //                     all_positive = false;
    //                     break;
    //                 }
    //             }
    //             // if (!all_positive)
    //         }
    //         if (all_positive) {
    //             KRATOS_WATCH("HEY")
    //             distance = -1*distance;
    //         }
    //     }
    // }

}

void DefineEmbeddedWakeProcess::MarkWakeElements(){

    ModelPart& deactivated_model_part = mrModelPart.CreateSubModelPart("deactivated_model_part");
    std::vector<std::size_t> deactivated_elements_id_list;
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin() + i;

        BoundedVector<double, 3> nodal_distances_to_wake = it_elem->GetValue(ELEMENTAL_DISTANCES);
        it_elem->SetValue(WAKE_ELEMENTAL_DISTANCES, nodal_distances_to_wake);

        // Selecting the cut (wake) elements
        const bool is_wake_element = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(nodal_distances_to_wake);

        BoundedVector<double,3> geometry_distances;
        Vector distances(3);
        for(unsigned int i_node = 0; i_node<3; i_node++){
            geometry_distances[i_node] = it_elem->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
            distances(i_node) = geometry_distances[i_node];
        }
        const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(geometry_distances);
        if (is_embedded){
            ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(it_elem->pGetGeometry(), Vector(geometry_distances));
            // Computing Normal
            std::vector<Vector> cut_normal;
            pModifiedShFunc -> ComputePositiveSideInterfaceAreaNormals(cut_normal,GeometryData::GI_GAUSS_1);
            double norm_normal = sqrt(inner_prod(cut_normal[0],cut_normal[0]));
            auto unit_normal = cut_normal[0]/norm_normal;
            it_elem->SetValue(VELOCITY_LOWER,unit_normal);
        }
        // Mark wake element and save their nodal distances to the wake
        if (is_wake_element) {
            if (is_embedded){
                #pragma omp critical
                {
                    deactivated_elements_id_list.push_back(it_elem->Id());
                }
                it_elem->Set(SOLID, true);
                KRATOS_WATCH(it_elem->Id())
                KRATOS_WATCH(geometry_distances)
                // it_elem->Set(ACTIVE, false);
            }
            else{
                it_elem->SetValue(WAKE, true);
            }
            auto r_geometry = it_elem->GetGeometry();
            for (unsigned int i = 0; i < it_elem->GetGeometry().size(); i++) {
                r_geometry[i].SetLock();
                r_geometry[i].SetValue(WAKE_DISTANCE, nodal_distances_to_wake(i));
                r_geometry[i].UnSetLock();
            }
        }
    }
    deactivated_model_part.AddElements(deactivated_elements_id_list);
}

ModifiedShapeFunctions::Pointer DefineEmbeddedWakeProcess::pGetModifiedShapeFunctions(const GeomPointerType pGeometry, const Vector& rDistances) const {
        return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

void DefineEmbeddedWakeProcess::ComputeTrailingEdgeNode(){

    double max_distance = 0.0;
    ModelPart& deactivated_model_part = mrModelPart.GetSubModelPart("deactivated_model_part");
    Element::Pointer p_max_elem;

    auto wake_origin = mrModelPart.GetProcessInfo()[WAKE_ORIGIN];

    // Find furthest deactivated element to the wake origin
    for (int i = 0; i < static_cast<int>(deactivated_model_part.Elements().size()); i++) {
        ModelPart::ElementIterator it_elem = deactivated_model_part.ElementsBegin() + i;

        BoundedVector<double,2> distance_vector;
        distance_vector[0] = wake_origin[0] - it_elem->GetGeometry().Center().X();
        distance_vector[1] = wake_origin[1] - it_elem->GetGeometry().Center().Y();
        double norm = norm_2(distance_vector);
        if(norm>max_distance){
            max_distance = norm;
            p_max_elem = mrModelPart.pGetElement(it_elem->Id());
        }
    }
    p_max_elem->SetValue(KUTTA, true);
    KRATOS_WATCH(p_max_elem->Id())
///////////////////////////////////////////////////////////////////////////////////////////////////

    // Mark nodes of the furthest deactivated element and store its neighbour elements
    for (unsigned int i_node= 0; i_node < p_max_elem->GetGeometry().size(); i_node++) {
        p_max_elem->GetGeometry()[i_node].Set(MARKER,true);
        const GlobalPointersVector<Element>& r_node_elem_candidates = p_max_elem -> GetGeometry()[i_node].GetValue(NEIGHBOUR_ELEMENTS);
        for (std::size_t j = 0; j < r_node_elem_candidates.size(); j++) {
            mKuttaWakeElementCandidates.push_back(r_node_elem_candidates(j));
        }
    }

    mrModelPart.RemoveSubModelPart("deactivated_model_part");
}

void DefineEmbeddedWakeProcess::MarkKuttaWakeElements(){

    // Find elements that touch the furthest deactivated element and that are part of the wake.
    std::vector<std::size_t> trailing_edge_node_list;
    for (std::size_t i = 0; i < mKuttaWakeElementCandidates.size(); i++)
    {
        auto& r_geometry = mKuttaWakeElementCandidates[i].GetGeometry();
        if (mKuttaWakeElementCandidates[i].GetValue(WAKE) && mKuttaWakeElementCandidates[i].Is(ACTIVE)) {
            int counter = 0;
            for (std::size_t i_node= 0; i_node < r_geometry.size(); i_node++) {
                if(r_geometry[i_node].Is(MARKER)){
                    trailing_edge_node_list.push_back(r_geometry[i_node].Id());
                    counter++;
                }
                if (counter>1)
                    mKuttaWakeElementCandidates[i].Set(STRUCTURE);
            }
        }

    }


    if (mrModelPart.HasSubModelPart("trailing_edge_sub_model_part")){
        mrModelPart.RemoveSubModelPart("trailing_edge_sub_model_part");
    }
    mrModelPart.CreateSubModelPart("trailing_edge_sub_model_part");

    std::sort(trailing_edge_node_list.begin(),
              trailing_edge_node_list.end());
    mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").AddNodes(trailing_edge_node_list);

    std::size_t max_number_of_structure_elements = 0;
    auto p_max_node = mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").NodesBegin();

    bool is_found = false;

    for (int i = 0; i < static_cast<int>(mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").Nodes().size()); i++) {
        auto it_node = mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").NodesBegin() + i;

        // const GlobalPointersVector<Element>& r_node_elem_candidates = it_node -> GetValue(NEIGHBOUR_ELEMENTS);
        // std::size_t counter = 0;
        // for (std::size_t j = 0; j < r_node_elem_candidates.size(); j++) {
        //     if (r_node_elem_candidates(j)->Is(STRUCTURE)){
        //         counter++;
        //     }
        // }
        // if (counter > max_number_of_structure_elements) {
        //     max_number_of_structure_elements = counter;
        //     p_max_node = it_node;
        // }
        if (it_node->GetValue(WAKE_DISTANCE) > 0.0) {
            it_node->SetValue(TRAILING_EDGE, true);
            is_found = true;
        }
    }
    if (!is_found)
        p_max_node->SetValue(TRAILING_EDGE, true);

}
}// Namespace Kratos




// ////////////////////////////////////////////////////////////////////////////////




//     // Mark nodes of the furthest deactivated element and store its neighbour elements
//     for (unsigned int i_node= 0; i_node < p_max_elem->GetGeometry().size(); i_node++) {
//         p_max_elem->GetGeometry()[i_node].SetValue(TRAILING_EDGE,true);
//         const GlobalPointersVector<Element>& r_node_elem_candidates = p_max_elem -> GetGeometry()[i_node].GetValue(NEIGHBOUR_ELEMENTS);
//         for (std::size_t j = 0; j < r_node_elem_candidates.size(); j++) {
//             mKuttaWakeElementCandidates.push_back(r_node_elem_candidates(j));
//         }
//     }

//     mrModelPart.RemoveSubModelPart("deactivated_model_part");
// }

// void DefineEmbeddedWakeProcess::MarkKuttaWakeElements(){

//     // Find elements that touch the furthest deactivated element and that are part of the wake.
//     std::vector<std::size_t> trailing_edge_node_list;
//     for (std::size_t i = 0; i < mKuttaWakeElementCandidates.size(); i++)
//     {
//         auto& r_geometry = mKuttaWakeElementCandidates[i].GetGeometry();
//         if (mKuttaWakeElementCandidates[i].GetValue(WAKE) && mKuttaWakeElementCandidates[i].Is(ACTIVE)) {
//             for (std::size_t i_node= 0; i_node < r_geometry.size(); i_node++) {
//                 if(r_geometry[i_node].GetValue(TRAILING_EDGE)){
//                     trailing_edge_node_list.push_back(r_geometry[i_node].Id());
//                     mKuttaWakeElementCandidates[i].Set(STRUCTURE);
//                 }
//             }
//         }

//     }

//     if (mrModelPart.HasSubModelPart("trailing_edge_sub_model_part")){
//         mrModelPart.RemoveSubModelPart("trailing_edge_sub_model_part");
//     }
//     mrModelPart.CreateSubModelPart("trailing_edge_sub_model_part");

//     std::sort(trailing_edge_node_list.begin(),
//               trailing_edge_node_list.end());
//     mrModelPart.GetSubModelPart("trailing_edge_sub_model_part").AddNodes(trailing_edge_node_list);

// }
// }// Namespace Kratos
