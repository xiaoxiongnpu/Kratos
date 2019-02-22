// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "adjoint_lift_response_function_coordinates_jump.h"
// #include "node.h"
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{
    AdjointLiftJumpCoordinatesResponseFunction::AdjointLiftJumpCoordinatesResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
     : AdjointPotentialResponseFunction(rModelPart, ResponseSettings)
    {
        // This response function currently only works in 2D!
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const int domain_size = r_current_process_info[DOMAIN_SIZE];
        KRATOS_ERROR_IF(domain_size != 2) << "Invalid DOMAIN_SIZE: " << domain_size << std::endl;

        // Get pointer to element that contains the traced node
        this->GetNeighboringElementPointer();
    }

    AdjointLiftJumpCoordinatesResponseFunction::~AdjointLiftJumpCoordinatesResponseFunction(){}

    void AdjointLiftJumpCoordinatesResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);
        rResponseGradient.clear();
        if( rAdjointElement.Id() == mpNeighboringElement->Id() )
        {
            const array_1d<double, 3> v_inf = rProcessInfo.GetValue(VELOCITY_INFINITY);
            double v_norm = norm_2(v_inf);
            double derivative= 2.0/v_norm;
            unsigned int NumNodes = rAdjointElement.GetGeometry().size();
            for(IndexType i = 0; i < NumNodes; ++i)
            {
                if(rAdjointElement.GetGeometry()[i].GetValue(TRAILING_EDGE))
                {   
                    rResponseGradient[i] = derivative;
                    rResponseGradient[i+NumNodes] = -derivative;
                }
            }
        }
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculateGradient(const Condition& rAdjointCondition,
                                   const Matrix& rResidualGradient,
                                   Vector& rResponseGradient,
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rResponseGradient.size() != rResidualGradient.size1())
            rResponseGradient.resize(rResidualGradient.size1(), false);

        rResponseGradient.clear();
        KRATOS_CATCH("");
    }
  

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("")
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);

        KRATOS_CATCH("")
    }

    void AdjointLiftJumpCoordinatesResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;


        if (rSensitivityGradient.size() != 0)
            rSensitivityGradient.resize(0, false);
     
        KRATOS_CATCH("");
    }

    void AdjointLiftJumpCoordinatesResponseFunction::GetNeighboringElementPointer()
    {
        KRATOS_TRY;

        for (auto elem_it = mrModelPart.Elements().ptr_begin(); elem_it != mrModelPart.Elements().ptr_end(); ++elem_it)
        {   
            if ((*elem_it)->Is(STRUCTURE)){
                mpNeighboringElement = (*elem_it);
                return;
            }
        }
        KRATOS_ERROR << "No neighboring element is available for the traced node." << std::endl;

        KRATOS_CATCH("");
    }
} // namespace Kratos.


