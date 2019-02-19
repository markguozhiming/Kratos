//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					     Kratos default license: kratos/license.txt
//
//  Main author: Suneth Warnakulasuriya
//

#if !defined(KRATOS_TURBULENCE_EDDY_VISCOSITY_PROCESS_H_INCLUDED)
#define KRATOS_TURBULENCE_EDDY_VISCOSITY_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "../custom_strategies/residual_based_bossak_velocity_scheme.h"
#include "../rans_constitutive_laws_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "processes/variational_distance_calculation_process.h"

namespace Kratos
{
///@addtogroup RANSConstitutiveLawsApplication
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

/// Auxiliary process to set Boussinesq buoyancy forces in variable temperature flows.
/** This process modifies the BODY_FORCE variable according to the Boussinesq hypothesis
    so that the fluid element can take natural convection into account.

    This process makes use of the following data:
    - TEMPERATURE from the nodal solution step data: current temperature for the node (mandatory).
    - AMBIENT_TEMPERATURE from ProcessInfo: The reference temperature for the simulation (mandatory).
    - gravity from the Parameters passed in the constructor: an array that defines the gravity vector (mandatory).
    - thermal_expansion_coefficient from the Parameters: a double defining the thermal expansion coefficient for the fluid (optional).

    With this, the process calculates the Boussinesq force and assings it to the BODY_FORCE solution step variable of each node.
    The force is set to (1 + thermal_expansion_coefficient*(temperature - ambient_temperature) ) * g

    If the thermal expansion coefficient is not provided, it is assumed to be (1/ambient_temperature).
    This is the usual value for perfect gases (if the temperature is given in Kelvin).
 */
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class TurbulenceEddyViscosityModelProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodesContainerType NodesArrayType;

    /// Pointer definition of TurbulenceEddyViscosityModelProcess
    KRATOS_CLASS_POINTER_DEFINITION(TurbulenceEddyViscosityModelProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    TurbulenceEddyViscosityModelProcess(ModelPart& rModelPart,
                                        Parameters& rParameters,
                                        typename TLinearSolver::Pointer pLinearSolver);

    /// Destructor.
    ~TurbulenceEddyViscosityModelProcess() override
    {
        // delete mpDistanceCalculator;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

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
    ModelPart& mrModelPart;
    Parameters& mrParameters;
    typename TLinearSolver::Pointer mpLinearSolver;

    bool mIsMeshMoving;
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void AddSolutionStepVariables()
    {
        KRATOS_TRY

        mrModelPart.GetNodalSolutionStepVariablesList().push_back(DISTANCE);
        mrModelPart.GetNodalSolutionStepVariablesList().push_back(FLAG_VARIABLE);
        mrModelPart.GetNodalSolutionStepVariablesList().push_back(KINEMATIC_VISCOSITY);
        mrModelPart.GetNodalSolutionStepVariablesList().push_back(TURBULENT_VISCOSITY);

        KRATOS_INFO("TurbulenceModel")<<"Added eddy viscosity turbulence model solution step variables.\n";

        KRATOS_CATCH("");
    }

    virtual void AddDofs()
    {
        KRATOS_INFO("TurbulenceModel")<<"Added eddy viscosity turbulence model dofs.\n";
    }

    virtual void InitializeTurbulenceModelPart()
    {
        KRATOS_THROW_ERROR(std::runtime_error, "Calling base class InitializeTurbulenceModelPart.", "");
    }

    virtual void UpdateFluidViscosity()
    {
        KRATOS_TRY

        NodesArrayType& nodes = mrModelPart.Nodes();

// Modifying viscosity of the nodes with the calculated turbulent viscosity
#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            auto it_node = nodes.begin() + i;
            const double kinematic_viscosity =
                it_node->FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
            const double turbulent_viscosity =
                it_node->FastGetSolutionStepValue(TURBULENT_VISCOSITY);

            double& effective_viscosity = it_node->FastGetSolutionStepValue(VISCOSITY);
            effective_viscosity = kinematic_viscosity + turbulent_viscosity;
        }

        KRATOS_CATCH("");
    }


    virtual void InitializeConditionFlags(const Flags& rFlag)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "Calling base class InitializeConditionFlags.", "");
    }

    void GenerateModelPart(ModelPart& rOriginModelPart,
                           ModelPart& rDestinationModelPart,
                           const Element& rReferenceElement,
                           const Condition& rReferenceCondition);

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

    VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>* mpDistanceCalculator;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateWallDistances();

    void InitializeNodeFlags(const Parameters& rParameters, const Flags& rFlag)
    {
        KRATOS_TRY

        for (std::size_t i = 0; i < rParameters.size(); ++i)
        {
            std::string model_part_name = rParameters.GetArrayItem(i).GetString();
            KRATOS_ERROR_IF(!mrModelPart.HasSubModelPart(model_part_name))
                << "TurbulenceEddyViscosityModelProcess: Wall condition "
                << model_part_name << " not found." << std::endl;
            ModelPart& current_model_part = mrModelPart.GetSubModelPart(model_part_name);

            NodesArrayType& nodes_array = current_model_part.Nodes();

#pragma omp parallel for
            for (int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
            {
                auto it_node = nodes_array.begin() + i;
                it_node->Set(rFlag, true);
            }
        }

        KRATOS_CATCH("");
    }

    // void AssignBoundaryConditions();

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
    TurbulenceEddyViscosityModelProcess& operator=(
        TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    TurbulenceEddyViscosityModelProcess(
        TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    ///@}

}; // Class TurbulenceEddyViscosityModelProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TURBULENCE_EDDY_VISCOSITY_PROCESS_H_INCLUDED  defined