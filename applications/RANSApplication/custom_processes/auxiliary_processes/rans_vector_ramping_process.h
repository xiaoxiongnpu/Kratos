//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_VECTOR_RAMPING_PROCESS_H_INCLUDED)
#define KRATOS_RANS_VECTOR_RAMPING_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "utilities/python_function_callback_utility.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @brief Checks lower and upper bounds of a scalar
 *
 * This process checks lower and upper bounds of a variable in nodes in a given modelpart
 *
 */

class KRATOS_API(RANS_APPLICATION) RansVectorRampingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    /// Pointer definition of RansVectorRampingProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansVectorRampingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansVectorRampingProcess(Model& rModel, Parameters rParameters);

    /// Destructor.
    ~RansVectorRampingProcess() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    PythonGenericFunctionUtility::Pointer mpFunction;

    Model& mrModel;
    Parameters mrParameters;
    std::string mModelPartName;

    std::string mVariableName;
    Vector mDirection;

    double mModulus = 0.0;
    int mRampingIterations;
    int mCurrentRampingIteration;
    int mEchoLevel;
    bool mIsFunction = false;
    bool mIsConstrained;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RansVectorRampingProcess& operator=(RansVectorRampingProcess const& rOther);

    /// Copy constructor.
    RansVectorRampingProcess(RansVectorRampingProcess const& rOther);

    ///@}

}; // Class RansVectorRampingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const RansVectorRampingProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_VECTOR_RAMPING_PROCESS_H_INCLUDED defined
