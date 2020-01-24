//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Massimiliano Zecchetto
//  Collaborators:
//
//-------------------------------------------------------------
//

#if !defined( UPDATE_THERMAL_MODEL_PART_H )
#define  UPDATE_THERMAL_MODEL_PART_H



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
//#include "utilities/openmp_utils.h" //TODO: check if it is possible to do in parallel
#include "processes/process.h"


namespace Kratos
{
///@}
///@name Kratos Classes
///@{


/// A tool to generate a copy of a ModelPart, sharing the same nodes as the original.
class UpdateThermalModelPartProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    //typedef Element BaseType;
    //typedef Element::ElementType ElementType;
    //typedef ModelPart::ElementType ElementType;
    //typedef Table<double,double> TableType;
    typedef Element BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(UpdateThermalModelPartProcess);

    ///@}
    ///@name Life Cycle
    ///@{
//const Element& two_dimension_reference_element, const Element& three_dimension_reference_element,
    /// Default Constructor
    UpdateThermalModelPartProcess(ModelPart& origin_model_part, ModelPart& destination_model_part,
                                Parameters rParameters
                                ) : Process(Flags()), rOriginModelPart(origin_model_part), rDestinationModelPart(destination_model_part)
                                //TwoDReferenceElement(two_dimension_reference_element), ThreeDReferenceElement(three_dimension_reference_element)//mrModelPart(model_part)
    {
        KRATOS_TRY
        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "domain_size": 2,
                "three": "unknown_name",
                "two": "unknown_name"
            }  )" );

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        DomainSize = rParameters["domain_size"].GetInt();
        TwoDElement = rParameters["two"].GetString();
        ThreeDElement = rParameters["three"].GetString();
        //unsigned int TableId = rParameters["table"].GetInt();
        //mpTable = model_part.pGetTable(TableId);

        //const Element& TwoDReferenceElement = KratosComponents<Element>::Get(TwoDElement);
        //const Element& ThreeDReferenceElement = KratosComponents<Element>::Get(ThreeDElement);
        KRATOS_CATCH("");
    }

    /// Copy constructor.
    //UpdateThermalModelPartProcess(UpdateThermalModelPartProcess const& rOther) = delete; TRY TO DECOMMENT THIS!!!!! IdP COMMENTED

    /// Destructor.
    //~UpdateThermalModelPartProcess() override = default; THIS IS THE ORIGINAL FROM CONNECTIVITY PRESERVER MODELER
    ~UpdateThermalModelPartProcess() override {}

    ///@}
    ///@name Operators
    ///@{

	//void operator()()
	//{
	//	ExecuteInitializeSolutionStep();
	//}
    /// Assignment operator.
    //UpdateThermalModelPartProcess & operator=(UpdateThermalModelPartProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Generate a copy of rOriginModelPart in rDestinationModelPart, using the given element and condtion types.
    /** This function fills rDestinationModelPart using data obtained from rOriginModelPart. The elements
     *  and conditions of rDestinationModelPart part use the same connectivity (and id) as in
     *  rOriginModelPart but their type is determined by rReferenceElement and rReferenceBoundaryCondition.
     *  Note that both ModelParts will share the same nodes, as well as ProcessInfo and tables.
     *  SubModelParts and, in MPI, communicator data will be replicated in DestinationModelPart.
     *  @param rOriginModelPart The source ModelPart.
     *  @param rDestinationModelPart The ModelPart to be filled by this function.
     *  @param rReferenceElement The Element type for rDestinationModelPart.
     *  @param rReferenceBoundaryCondition The Condition type for rDestinationModelPart.
     */

    void Execute() override {}
    void ExecuteInitialize() override {}
    void ExecuteInitializeSolutionStep() override
        //ModelPart& rOriginModelPart,
        //ModelPart& rDestinationModelPart,
        //Element const& rReferenceElement,
        //Condition const& rReferenceBoundaryCondition) override
    {
        KRATOS_TRY;
        //KRATOS_WATCH("BEFORE RESET DESTINATION MODEL PART");
        //KRATOS_WATCH(rOriginModelPart);
        //KRATOS_WATCH(rDestinationModelPart);
        this->ResetDestinationModelPart(rOriginModelPart, rDestinationModelPart);
        //KRATOS_WATCH("AFTER RESET DESTINATION MODEL PART");
        //KRATOS_WATCH(rOriginModelPart);
        //KRATOS_WATCH(rDestinationModelPart);
        this->CopyNodes(rOriginModelPart, rDestinationModelPart);
        //KRATOS_WATCH("AFTER COPY NODES");
        //KRATOS_WATCH(rOriginModelPart);
        //KRATOS_WATCH(rDestinationModelPart);
        this->DuplicateElements();//(rOriginModelPart, rDestinationModelPart);
        //KRATOS_WATCH("AFTER COPY ELEMENTS");
        //KRATOS_WATCH(rOriginModelPart);
        //KRATOS_WATCH(rDestinationModelPart);

        //this->DuplicateConditions(rOriginModelPart, rDestinationModelPart);

        //this->DuplicateSubModelParts();//(rOriginModelPart, rDestinationModelPart);
        //KRATOS_WATCH("AFTER TRANSFER ELEMENTS");
        //KRATOS_WATCH(rOriginModelPart);
        //KRATOS_WATCH(rDestinationModelPart);

        KRATOS_CATCH("");
}

    ///@}
protected:
    ModelPart& rOriginModelPart;
    ModelPart& rDestinationModelPart;
    //Condition const& rReferenceBoundaryCondition;
    std::string TwoDElement;
    //ElementType::Pointer mpTable;
    //Element& TwoDReferenceElement;
    //Element TwoDReferenceElement; this does not work without being a reference &
    std::string ThreeDElement;
    //Element ThreeDReferenceElement;
    //TableType::Pointer mpTable;
    //Element const& TwoDReferenceElement;
    //Element const& ThreeDReferenceElement;
    unsigned int DomainSize;

private:
    ///@name Private Operations
    ///@{
    //UpdateThermalModelPartProcess &operator=(UpdateThermalModelPartProcess const &rOther);
    void ResetDestinationModelPart(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart) const
    {
        rOriginModelPart.RemoveNodesFromAllLevels(TO_ERASE);

        VariableUtils().SetFlag(TO_ERASE, true, rDestinationModelPart.Nodes());
        VariableUtils().SetFlag(TO_ERASE, true, rDestinationModelPart.Elements());
        //VariableUtils().SetFlag(TO_ERASE, true, rDestinationModelPart.Conditions());

        rDestinationModelPart.RemoveNodesFromAllLevels(TO_ERASE);
        rDestinationModelPart.RemoveElementsFromAllLevels(TO_ERASE);
        //rDestinationModelPart.RemoveConditionsFromAllLevels(TO_ERASE);

        VariableUtils().SetFlag(TO_ERASE, false, rOriginModelPart.Nodes());
    }

    //void UpdateThermalModelPartProcess::CopyCommonData(
    void CopyNodes(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart) const
    {
        // Assign the nodes to the new model part
        rDestinationModelPart.AddNodes(rOriginModelPart.NodesBegin(), rOriginModelPart.NodesEnd());
        for (auto i_part = rOriginModelPart.SubModelPartsBegin(); i_part != rOriginModelPart.SubModelPartsEnd(); ++i_part) {
            ModelPart& destination_part = rDestinationModelPart.GetSubModelPart(i_part->Name());
            destination_part.AddNodes(i_part->NodesBegin(), i_part->NodesEnd());
        }
    }

    //void DuplicateElements(
    //    ModelPart& rOriginModelPart,
    //    ModelPart& rDestinationModelPart,
    //    const Element& rReferenceElement
    //) const;

    //void UpdateThermalModelPartProcess::DuplicateElements(
    void DuplicateElements() const//(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart) const
        //const Element& rReferenceElement) const
    {
        const Element& TwoDReferenceElement = KratosComponents<Element>::Get(TwoDElement);
        const Element& ThreeDReferenceElement = KratosComponents<Element>::Get(ThreeDElement);
        // Generate the elements
        for (ModelPart::SubModelPartIterator i_mp = rOriginModelPart.SubModelPartsBegin(); i_mp != rOriginModelPart.SubModelPartsEnd(); i_mp++)
        {
            //KRATOS_WATCH(i_mp->Name());
            //#pragma omp parallel
                //{
                    // skipping sub model parts without elements
                    if (i_mp->NumberOfElements()) {
                    // transfering the elements
                    std::vector<ModelPart::IndexType> ids;
                    ids.reserve(i_mp->Elements().size());
                    ModelPart& destination_part = rDestinationModelPart.GetSubModelPart(i_mp->Name());
                    //ModelPart::ElementIterator ElemBegin;
                    //ModelPart::ElementIterator ElemEnd;
                    //OpenMPUtils::PartitionedIterators(rOriginModelPart.Elements(), ElemBegin, ElemEnd);
                    // selecting only the fluid and solid sub model parts
                    //ModelPart::ElementsContainerType temp_elements; original from connectivity modeler
                    ModelPart::ElementsContainerType temp_elements;
                    //temp_elements.reserve(rOriginModelPart.NumberOfElements());
                    temp_elements.reserve(i_mp->NumberOfElements());
                    if ((i_mp->Is(SOLID) && i_mp->IsNot(ACTIVE)) || (i_mp->Is(FLUID) && i_mp->IsNot(ACTIVE)) || (i_mp->Is(BOUNDARY) && i_mp->Is(RIGID)))
                    {
                        //KRATOS_WATCH(i_mp->Name());
                        if (DomainSize == 3)
                        {
                            // Generate the elements
                            for (auto i_elem = i_mp->ElementsBegin(); i_elem != i_mp->ElementsEnd(); ++i_elem) {
                                Properties::Pointer properties = i_elem->pGetProperties();

                                // Reuse the geometry of the old element (to save memory)
                                Element::Pointer p_element = ThreeDReferenceElement.Create(i_elem->Id(), i_elem->pGetGeometry(), properties);

                                temp_elements.push_back(p_element);
                            }
                            // la riga seguente non e' corretta in quanto gli elementi vanno aggiunti al destination submodelpart
                            //i_mp->AddElements(temp_elements.begin(), temp_elements.end()); //commenting because I can do that with the funtion already provided by the connectivity preserve modeler
                            rDestinationModelPart.AddElements(temp_elements.begin(), temp_elements.end());
                        } else {
                            // Generate the elements
                            for (auto i_elem = i_mp->ElementsBegin(); i_elem != i_mp->ElementsEnd(); ++i_elem) {
                                Properties::Pointer properties = i_elem->pGetProperties();

                                // Reuse the geometry of the old element (to save memory)
                                Element::Pointer p_element = TwoDReferenceElement.Create(i_elem->Id(), i_elem->pGetGeometry(), properties);

                                temp_elements.push_back(p_element);
                            }
                            // la riga seguente non e' corretta in quanto gli elementi vanno aggiunti al destination submodelpart
                            //i_mp->AddElements(temp_elements.begin(), temp_elements.end()); //commenting because I can do that with the funtion already provided by the connectivity preserve modeler
                            //KRATOS_WATCH("rDestinationModelPart");
                            //KRATOS_WATCH(rDestinationModelPart);
                            rDestinationModelPart.AddElements(temp_elements.begin(), temp_elements.end());
                        }
                        //adding elements by index
                        for(auto it=i_mp->ElementsBegin(); it!=i_mp->ElementsEnd(); ++it)
                            ids.push_back(it->Id());
                        destination_part.AddElements(ids, 0);
                    } //else {
                        // Skipping the Computing Model Part of the pfem fluid model part
                        // How to skip the computing model part of the Thermal Model part???
                        //if (!((i_mp->Is(ACTIVE) && i_mp->Is(SOLID)) || (i_mp->Is(ACTIVE) && i_mp->Is(FLUID))))
                        //{
                        //    if (DomainSize == 3)
                        //    {
                        //        //KRATOS_WATCH(i_mp->Name());
                        //        for (auto i_elem = i_mp->ElementsBegin(); i_elem != i_mp->ElementsEnd(); ++i_elem) {
                        //        Properties::Pointer properties = i_elem->pGetProperties();

                        //        // Reuse the geometry of the old element (to save memory)
                        //        Element::Pointer p_element = TwoDReferenceElement.Create(i_elem->Id(), i_elem->pGetGeometry(), properties);

                        //        temp_elements.push_back(p_element);
                        //    }
                        //    // la riga seguente non e' corretta in quanto gli elementi vanno aggiunti al destination submodelpart
                        //    //i_mp->AddElements(temp_elements.begin(), temp_elements.end()); //commenting because I can do that with the funtion already provided by the connectivity preserve modeler
                        //    rDestinationModelPart.AddElements(temp_elements.begin(), temp_elements.end());
                        //    for(auto it=i_mp->ElementsBegin(); it!=i_mp->ElementsEnd(); ++it)
                        //        ids.push_back(it->Id());
                        //    destination_part.AddElements(ids, 0);
                        //    }
                        //}

                    //}
                    }
                //}

        }
    }

    //void DuplicateConditions(
    //    ModelPart& rOriginModelPart,
    //    ModelPart& rDestinationModelPart,
    //    const Condition& rReferenceBoundaryCondition
    //) const;

    //void UpdateThermalModelPartProcess::DuplicateConditions(
    //void DuplicateConditions(
    //    ModelPart& rOriginModelPart,
    //    ModelPart& rDestinationModelPart,
    //    const Condition& rReferenceBoundaryCondition) const
    //{
    //    // Generate the conditions
    //    ModelPart::ConditionsContainerType temp_conditions;
    //    temp_conditions.reserve(rOriginModelPart.NumberOfConditions());
    //    for (auto i_cond = rOriginModelPart.ConditionsBegin(); i_cond != rOriginModelPart.ConditionsEnd(); ++i_cond) {
    //        Properties::Pointer properties = i_cond->pGetProperties();

    //        // Reuse the geometry of the old element (to save memory)
    //        Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(i_cond->Id(), i_cond->pGetGeometry(), properties);

    //        temp_conditions.push_back(p_condition);
    //    }

    //    rDestinationModelPart.AddConditions(temp_conditions.begin(), temp_conditions.end());
    //}

    //void UpdateThermalModelPartProcess::DuplicateSubModelParts(
    //void DuplicateSubModelParts() const//(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart) const
    //{
        //for (auto i_part = rOriginModelPart.SubModelPartsBegin(); i_part != rOriginModelPart.SubModelPartsEnd(); ++i_part) {
            //if(!rDestinationModelPart.HasSubModelPart(i_part->Name())) {
            //    rDestinationModelPart.CreateSubModelPart(i_part->Name()); I can get rid of this because I already cloned the pfem model part
            //}
            //ModelPart& destination_part = rDestinationModelPart.GetSubModelPart(i_part->Name());
            //destination_part.AddNodes(i_part->NodesBegin(), i_part->NodesEnd());
            //if ((i_part->Is(SOLID) && i_part->IsNot(ACTIVE)) || (i_part->Is(FLUID) && i_part->IsNot(ACTIVE)) || (i_part->Is(BOUNDARY) && i_part->Is(RIGID)))
            //{

            //ModelPart& destination_part = rDestinationModelPart.GetSubModelPart(i_part->Name());

            //destination_part.AddNodes(i_part->NodesBegin(), i_part->NodesEnd());

            //std::vector<ModelPart::IndexType> ids;
            //ids.reserve(i_part->Elements().size());

            // Execute only if we created elements in the destination
            //if (rDestinationModelPart.NumberOfElements() > 0)
            //{
            //    //adding by index
            //    for(auto it=i_part->ElementsBegin(); it!=i_part->ElementsEnd(); ++it)
            //        ids.push_back(it->Id());
            //    destination_part.AddElements(ids, 0); //adding by index
            //}
            //} else {
            //    //SKIPPING PFEM COMPUTING MODEL PART
            //    if (!((i_part->Is(ACTIVE) && i_part->Is(SOLID)) || (i_part->Is(ACTIVE) && i_part->Is(FLUID))))
            //    {
            //        if (DomainSize == 3)
            //        {
            //            ModelPart& destination_part = rDestinationModelPart.GetSubModelPart(i_part->Name());

            //            destination_part.AddNodes(i_part->NodesBegin(), i_part->NodesEnd());

            //            std::vector<ModelPart::IndexType> ids;
            //            ids.reserve(i_part->Elements().size());

            //            // Execute only if we created elements in the destination
            //            if (rDestinationModelPart.NumberOfElements() > 0)
            //            {
            //                //adding by index
            //                for(auto it=i_part->ElementsBegin(); it!=i_part->ElementsEnd(); ++it)
            //                    ids.push_back(it->Id());
            //                destination_part.AddElements(ids, 0); //adding by index
            //            }
            //        }
            //    }
            //}
            // Execute only if we created conditions in the destination
            //if (rDestinationModelPart.NumberOfConditions() > 0)
            //{
            //    ids.clear();
            //    for(auto it=i_part->ConditionsBegin(); it!=i_part->ConditionsEnd(); ++it)
            //        ids.push_back(it->Id());
            //    destination_part.AddConditions(ids, 0);
            //}

            // Recursively call this function to duplicate any child SubModelParts
            //this->DuplicateSubModelParts(*i_part, destination_part);
        //}/
    //}

    ///@}
};//Class UpdateThermalModelPartProcess

///@}
/// input stream function
//inline std::istream &operator>>(std::istream &rIStream,
//								UpdateThermalModelPartProcess &rThis);

///// output stream function
//inline std::ostream &operator<<(std::ostream &rOStream,
//								const UpdateThermalModelPartProcess &rThis)
//{
//	rThis.PrintInfo(rOStream);
//	rOStream << std::endl;
//	rThis.PrintData(rOStream);
//
//	return rOStream;
//}
///@}


} // namespace Kratos.

#endif //UPDATE_THERMAL_MODEL_PART_INCLUDED  defined