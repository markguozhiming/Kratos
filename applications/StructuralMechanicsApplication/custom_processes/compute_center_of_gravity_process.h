// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Philipp Bucher
//                   Salman Yousaf
//

#if !defined(KRATOS_COMPUTE_CENTER_OF_GRAVITY_PROCESS)
#define KRATOS_COMPUTE_CENTER_OF_GRAVITY_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{
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

/**
 * @class ComputeCenterOfGravityProcess
 *
 * @ingroup StructuralMechanicsApplication
 *
 * @brief This method computes the tonal mass of a structure
 * @details It takes into account the beam, shells and solid elements
 *
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ComputeCenterOfGravityProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeCenterOfGravityProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeCenterOfGravityProcess);

    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Point                                           PointType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef Geometry<PointType>                     GeometryPointType;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    typedef ModelPart::ElementsContainerType        ElementsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeCenterOfGravityProcess(
        ModelPart& rThisModelPart
        ):mrThisModelPart(rThisModelPart)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~ComputeCenterOfGravityProcess() override
    = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

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
    std::string Info() const override
    {
        return "ComputeCenterOfGravityProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeCenterOfGravityProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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

    ModelPart& mrThisModelPart;              // The main model part

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
    ComputeCenterOfGravityProcess& operator=(ComputeCenterOfGravityProcess const& rOther) = delete;

    /// Copy constructor.
    //ComputeCenterOfGravityProcess(ComputeCenterOfGravityProcess const& rOther);


    ///@}

}; // Class ComputeCenterOfGravityProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ComputeCenterOfGravityProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ComputeCenterOfGravityProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
#endif /* KRATOS_COMPUTE_CENTER_OF_GRAVITY_PROCESS defined */
