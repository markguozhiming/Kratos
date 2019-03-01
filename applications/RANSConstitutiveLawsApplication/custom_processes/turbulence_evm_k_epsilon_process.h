//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_TURBULENCE_EVM_K_EPSILON_PROCESS_H_INCLUDED)
#define KRATOS_TURBULENCE_EVM_K_EPSILON_PROCESS_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// Project includes
#include "custom_strategies/general_convergence_criteria.h"
#include "includes/cfd_variables.h"
#include "includes/communicator.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residual_based_bossak_velocity_scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/brent_iteration.h"
#include "utilities/derivatives_extension.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_conditions/evm_epsilon_wall_condition.h"
#include "custom_elements/evm_k_epsilon/evm_epsilon_element.h"
#include "custom_elements/evm_k_epsilon/evm_k_element.h"
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.h"
#include "custom_processes/turbulence_eddy_viscosity_model_process.h"
#include "rans_constitutive_laws_application_variables.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "utilities/color_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class ResidualBasedBossakTurbulentKineticEnergyScheme
    : public ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>
{
    class NodalDerivativesExtension : public DerivativesExtension
    {
        Node<3>* mpNode;

    public:
        explicit NodalDerivativesExtension(Node<3>* pNode) : mpNode(pNode)
        {
        }

        void GetFirstDerivativesVector(std::vector<IndirectScalar<double>>& rVector,
                                       std::size_t Step,
                                       ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            rVector[0] = MakeIndirectScalar(*this->mpNode, TURBULENT_KINETIC_ENERGY, Step);
        }

        void GetSecondDerivativesVector(std::vector<IndirectScalar<double>>& rVector,
                                        std::size_t Step,
                                        ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            rVector[0] =
                MakeIndirectScalar(*this->mpNode, TURBULENT_KINETIC_ENERGY_RATE, Step);
        }

        void GetFirstDerivativesDofsVector(std::vector<Dof<double>::Pointer>& rVector,
                                           ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            rVector[0] = this->mpNode->pGetDof(TURBULENT_KINETIC_ENERGY);
        }
    };

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBossakTurbulentKineticEnergyScheme);

    typedef ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace> BaseType;

    /// Constructor.

    ResidualBasedBossakTurbulentKineticEnergyScheme(const double AlphaBossak)
        : ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>(AlphaBossak, true, false)
    {
    }

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        const int number_of_nodes = rModelPart.NumberOfNodes();

#pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++)
        {
            Node<3>& r_current_node = *(rModelPart.NodesBegin() + i);
            r_current_node.SetValue(
                DERIVATIVES_EXTENSION,
                Kratos::make_shared<NodalDerivativesExtension>(&r_current_node));
        }

        KRATOS_CATCH("");
    }

    ///@}
};

template <class TSparseSpace, class TDenseSpace>
class ResidualBasedBossakTurbulentEnergyDissipationRateScheme
    : public ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>
{
    class NodalDerivativesExtension : public DerivativesExtension
    {
        Node<3>* mpNode;

    public:
        explicit NodalDerivativesExtension(Node<3>* pNode) : mpNode(pNode)
        {
        }

        void GetFirstDerivativesVector(std::vector<IndirectScalar<double>>& rVector,
                                       std::size_t Step,
                                       ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            rVector[0] = MakeIndirectScalar(
                *this->mpNode, TURBULENT_ENERGY_DISSIPATION_RATE, Step);
        }

        void GetSecondDerivativesVector(std::vector<IndirectScalar<double>>& rVector,
                                        std::size_t Step,
                                        ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            rVector[0] = MakeIndirectScalar(
                *this->mpNode, TURBULENT_ENERGY_DISSIPATION_RATE_2, Step);
        }

        void GetFirstDerivativesDofsVector(std::vector<Dof<double>::Pointer>& rVector,
                                           ProcessInfo& rCurrentProcessInfo) override
        {
            rVector.resize(1);
            rVector[0] = this->mpNode->pGetDof(TURBULENT_ENERGY_DISSIPATION_RATE);
        }
    };

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBossakTurbulentEnergyDissipationRateScheme);

    typedef ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace> BaseType;

    /// Constructor.

    ResidualBasedBossakTurbulentEnergyDissipationRateScheme(const double AlphaBossak)
        : ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>(AlphaBossak, true, false)
    {
    }

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        const int number_of_nodes = rModelPart.NumberOfNodes();

#pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++)
        {
            Node<3>& r_current_node = *(rModelPart.NodesBegin() + i);
            r_current_node.SetValue(
                DERIVATIVES_EXTENSION,
                Kratos::make_shared<NodalDerivativesExtension>(&r_current_node));
        }

        KRATOS_CATCH("");
    }

    ///@}
};

template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class TurbulenceEvmKEpsilonProcess
    : public TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TurbulenceEvmKEpsilonProcess
    KRATOS_CLASS_POINTER_DEFINITION(TurbulenceEvmKEpsilonProcess);

    typedef TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> StrategyType;
    typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;
    typedef typename SchemeType::Pointer SchemePointerType;
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BuilderAndSolverType;
    typedef typename BuilderAndSolverType::Pointer BuilderAndSolverPointerType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> ConvergenceCriteriaType;
    typedef typename ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    TurbulenceEvmKEpsilonProcess(ModelPart& rModelPart,
                                 Parameters& rParameters,
                                 typename TLinearSolver::Pointer pDistanceLinearSolver,
                                 typename TLinearSolver::Pointer pKLinearSolver,
                                 typename TLinearSolver::Pointer pEpsilonLinearSolver)
        : TurbulenceEddyViscosityModelProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>(
              rModelPart, rParameters, pDistanceLinearSolver),
          mpKLinearSolver(pKLinearSolver),
          mpEpsilonLinearSolver(pEpsilonLinearSolver)
    {
        KRATOS_TRY

        KRATOS_INFO("TurbulenceModel")
            << "Initializing k-epsilon turbulence model" << std::endl;

        Parameters default_parameters = Parameters(R"(
        {
            "convergence_tolerances":
            {
                "k_relative_tolerance": 1e-3,
                "k_absolute_tolerance": 1e-5,
                "epsilon_relative_tolerance": 1e-3,
                "epsilon_absolute_tolerance": 1e-5,
                "turbulent_viscosity_relative_tolerance": 1e-3,
                "turbulent_viscosity_absolute_tolerance": 1e-5,
                "k_max_iterations": 10,
                "epsilon_max_iterations": 10,
                "maximum_coupling_iterations": 10,
                "echo_level": 2
            },
            "echo_level"     : 0,
            "scheme_settings": {},
            "constants":
            {
                "wall_smoothness_beta"    : 5.2,
                "von_karman"              : 0.41,
                "c_mu"                    : 0.09,
                "c1"                      : 1.44,
                "c2"                      : 1.92,
                "sigma_k"                 : 1.0,
                "sigma_epsilon"           : 1.3
            },
            "flow_parameters":
            {
                "ramp_up_time"                : 0.03,
                "turbulent_mixing_length"     : 0.01,
                "turbulent_viscosity_fraction": 0.05,
                "free_stream_velocity"        : 1.0,
                "free_stream_k"               : 1.0,
                "free_stream_epsilon"         : 1.0,
                "velocity_ratio_constant"     : 0.003
            }
        })");

        this->mrParameters["model_properties"].ValidateAndAssignDefaults(default_parameters);

        const Parameters& model_properties =
            this->mrParameters["model_properties"]["constants"];

        rModelPart.GetProcessInfo()[WALL_SMOOTHNESS_BETA] =
            model_properties["wall_smoothness_beta"].GetDouble();
        rModelPart.GetProcessInfo()[WALL_VON_KARMAN] =
            model_properties["von_karman"].GetDouble();
        rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU] =
            model_properties["c_mu"].GetDouble();
        rModelPart.GetProcessInfo()[TURBULENCE_RANS_C1] =
            model_properties["c1"].GetDouble();
        rModelPart.GetProcessInfo()[TURBULENCE_RANS_C2] =
            model_properties["c2"].GetDouble();
        rModelPart.GetProcessInfo()[TURBULENT_KINETIC_ENERGY_SIGMA] =
            model_properties["sigma_k"].GetDouble();
        rModelPart.GetProcessInfo()[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA] =
            model_properties["sigma_epsilon"].GetDouble();

        const Parameters& flow_parameters =
            this->mrParameters["model_properties"]["flow_parameters"];
        rModelPart.GetProcessInfo()[TURBULENT_MIXING_LENGTH] =
            flow_parameters["turbulent_mixing_length"].GetDouble();
        rModelPart.GetProcessInfo()[TURBULENT_VISCOSITY_FRACTION] =
            flow_parameters["turbulent_viscosity_fraction"].GetDouble();

        mFreestreamVelocity = flow_parameters["free_stream_velocity"].GetDouble();
        mFreestreamK = flow_parameters["free_stream_k"].GetDouble();
        mFreestreamEpsilon = flow_parameters["free_stream_epsilon"].GetDouble();
        mVelocityRatio = flow_parameters["velocity_ratio_constant"].GetDouble();
        mRampUpTime = flow_parameters["ramp_up_time"].GetDouble();

        KRATOS_ERROR_IF(
            this->mrModelPart.HasSubModelPart("TurbulenceModelPartRANSEVMK"))
            << "TurbulenceEddyViscosityModelProcess: "
               "TurbulenceModelPartRANSEVMK is "
               "already found."
            << std::endl;
        this->mrModelPart.CreateSubModelPart("TurbulenceModelPartRANSEVMK");
        KRATOS_ERROR_IF(this->mrModelPart.HasSubModelPart(
            "TurbulenceModelPartRANSEVMEpsilon"))
            << "TurbulenceEddyViscosityModelProcess: "
               "TurbulenceModelPartRANSEVMEpsilon is "
               "already found."
            << std::endl;
        this->mrModelPart.CreateSubModelPart(
            "TurbulenceModelPartRANSEVMEpsilon");

        mpTurbulenceKModelPart =
            &(this->mrModelPart.GetSubModelPart("TurbulenceModelPartRANSEVMK"));
        mpTurbulenceEpsilonModelPart = &(this->mrModelPart.GetSubModelPart(
            "TurbulenceModelPartRANSEVMEpsilon"));

        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~TurbulenceEvmKEpsilonProcess()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Execute() override
    {
        ProcessInfo& r_current_process_info = this->mrModelPart.GetProcessInfo();

        const int process_step = r_current_process_info[RANS_MODELLING_PROCESS_STEP];

        if (process_step == 1)
            this->AddSolutionStepVariables();
        else if (process_step == 2)
            this->AddDofs();
        else if (process_step == 3)
        {
            const double current_time = r_current_process_info[TIME];

            if (current_time < mRampUpTime)
                return;
            if (this->mEchoLevel > 0)
                KRATOS_INFO("TurbulenceModel") << "Solving RANS equations...\n";

            this->SolveStep();
            this->UpdateFluidViscosity();
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,
                               "Provided step in RANS_MODELLING_PROCESS_STEP "
                               "is not supported.",
                               "");
        }
    }

    /// this function is designed for being called at the beginning of the
    /// computations right after reading the model and the groups
    virtual void ExecuteInitialize() override
    {
        BaseType::ExecuteInitialize();

        this->GenerateSolutionStrategies();
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        BaseType::ExecuteInitializeSolutionStep();

        if (IsBoundaryConditionsAssigned)
            return;

        ProcessInfo& r_current_process_info = this->mrModelPart.GetProcessInfo();
        const double current_time = r_current_process_info[TIME];
        if (current_time < mRampUpTime)
            return;

        this->AssignBoundaryConditions();

        IsBoundaryConditionsAssigned = true;

        KRATOS_CATCH("");
    }

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
    virtual std::string Info() const override
    {
        return std::string("TurbulenceEvmKEpsilonProcess");
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void InitializeTurbulenceModelPart() override
    {
        // KRATOS_TRY;

        KRATOS_INFO("TurbulenceModel")
            << "Initializing turbulence model part\n";

        std::stringstream element_postfix;
        element_postfix << TDim << "D" << TDim + 1 << "N";

        std::stringstream condition_postfix;
        condition_postfix << TDim << "D";

        const Element& k_element =
            KratosComponents<Element>::Get("RANSEVMK" + element_postfix.str());
        const Element& epsilon_element =
            KratosComponents<Element>::Get("RANSEVMEPSILON" + element_postfix.str());
        const Condition& k_cond =
            KratosComponents<Condition>::Get("Condition" + condition_postfix.str());
        // const Condition& epsilon_cond = KratosComponents<Condition>::Get(
        //     "EpsilonWallCondition" + condition_postfix.str());
        const Condition& epsilon_cond =
            KratosComponents<Condition>::Get("Condition" + condition_postfix.str());

        this->GenerateModelPart(this->mrModelPart, *mpTurbulenceKModelPart, k_element, k_cond);
        this->GenerateModelPart(this->mrModelPart, *mpTurbulenceEpsilonModelPart,
                                epsilon_element, epsilon_cond);

        KRATOS_INFO("TurbulenceModel") << *mpTurbulenceKModelPart;
        KRATOS_INFO("TurbulenceModel") << *mpTurbulenceEpsilonModelPart;

        // KRATOS_CATCH("");
    }

    void InitializeConditionFlags(const Flags& rFlag) override
    {
        KRATOS_TRY;

        InitializeConditionFlagsForModelPart(mpTurbulenceKModelPart, rFlag);
        InitializeConditionFlagsForModelPart(mpTurbulenceEpsilonModelPart, rFlag);

        KRATOS_CATCH("");
    }

    void AddSolutionStepVariables() override
    {
        KRATOS_TRY

        this->mrModelPart.GetNodalSolutionStepVariablesList().push_back(TURBULENT_KINETIC_ENERGY);
        this->mrModelPart.GetNodalSolutionStepVariablesList().push_back(
            TURBULENT_ENERGY_DISSIPATION_RATE);
        this->mrModelPart.GetNodalSolutionStepVariablesList().push_back(TURBULENT_KINETIC_ENERGY_RATE);
        this->mrModelPart.GetNodalSolutionStepVariablesList().push_back(
            TURBULENT_ENERGY_DISSIPATION_RATE_2);
        this->mrModelPart.GetNodalSolutionStepVariablesList().push_back(FRICTION_VELOCITY);

        BaseType::AddSolutionStepVariables();

        KRATOS_CATCH("");
    }

    void AddDofs() override
    {
        VariableUtils().AddDof<Variable<double>>(TURBULENT_KINETIC_ENERGY, this->mrModelPart);
        VariableUtils().AddDof<Variable<double>>(
            TURBULENT_ENERGY_DISSIPATION_RATE, this->mrModelPart);
        VariableUtils().AddDof<Variable<double>>(TURBULENT_KINETIC_ENERGY_RATE,
                                                 this->mrModelPart);
        VariableUtils().AddDof<Variable<double>>(
            TURBULENT_ENERGY_DISSIPATION_RATE_2, this->mrModelPart);

        BaseType::AddDofs();
    }
    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::map<std::size_t, double> mRatioTolerance;
    std::map<std::size_t, double> mAbsTolerance;

    typename TLinearSolver::Pointer mpKLinearSolver;

    typename StrategyType::Pointer mpKStrategy;

    typename TLinearSolver::Pointer mpEpsilonLinearSolver;

    typename StrategyType::Pointer mpEpsilonStrategy;

    double mFreestreamK;

    double mFreestreamEpsilon;

    double mFreestreamVelocity;

    double mVelocityRatio;

    int mMaximumCouplingIterations;

    double mTurbulentViscosityRelativeTolerance;

    double mTurbulentViscosityAbsoluteTolerance;

    double mRampUpTime;

    bool IsBoundaryConditionsAssigned = false;

    ConvergenceCriteriaPointerType mpConvergenceCriteria;

    BuilderAndSolverPointerType mpKBuilderAndSolver;

    BuilderAndSolverPointerType mpEpsilonBuilderAndSolver;

    ModelPart* mpTurbulenceKModelPart;

    ModelPart* mpTurbulenceEpsilonModelPart;

    ///@}
    ///@name Private Operators
    ///@{

    void InitializeConditionFlagsForModelPart(ModelPart* pModelPart, const Flags& rFlag)
    {
        KRATOS_TRY

        auto& conditions_array = pModelPart->Conditions();

        for (auto& it_cond : conditions_array)
        {
            bool is_flag_true = true;

            auto& r_geometry = it_cond.GetGeometry();

            for (auto& it_node : r_geometry.Points())
                if (!(it_node.Is(rFlag)))
                    is_flag_true = false;
            it_cond.Set(rFlag, is_flag_true);
        }

        KRATOS_CATCH("");
    }

    void GenerateSolutionStrategies()
    {
        KRATOS_TRY

        KRATOS_INFO("TurbulenceModel")
            << "Generating turbulence modelling strategies.\n";

        const Parameters& r_convergence_parameters =
            this->mrParameters["model_properties"]["convergence_tolerances"];

        double k_relative_tolerance =
            r_convergence_parameters["k_relative_tolerance"].GetDouble();
        double k_absolute_tolerance =
            r_convergence_parameters["k_absolute_tolerance"].GetDouble();
        int k_max_iterations = r_convergence_parameters["k_max_iterations"].GetInt();

        double epsilon_relative_tolerance =
            r_convergence_parameters["epsilon_relative_tolerance"].GetDouble();
        double epsilon_absolute_tolerance =
            r_convergence_parameters["epsilon_absolute_tolerance"].GetDouble();
        int epsilon_max_iterations =
            r_convergence_parameters["epsilon_max_iterations"].GetInt();

        this->mTurbulentViscosityRelativeTolerance =
            r_convergence_parameters["turbulent_viscosity_relative_tolerance"].GetDouble();
        this->mTurbulentViscosityAbsoluteTolerance =
            r_convergence_parameters["turbulent_viscosity_absolute_tolerance"].GetDouble();
        this->mMaximumCouplingIterations =
            r_convergence_parameters["maximum_coupling_iterations"].GetInt();

        bool CalculateReactions = false;
        bool ReformDofSet = false;

        const double alpha_bossak =
            this->mrParameters["model_properties"]["scheme_settings"]
                              ["alpha_bossak"]
                                  .GetDouble();

        // K solution strategy
        BuilderAndSolverPointerType pKBuilderAndSolver = BuilderAndSolverPointerType(
            new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(
                mpKLinearSolver));
        mpKBuilderAndSolver = pKBuilderAndSolver;

        SchemePointerType pKScheme = SchemePointerType(
            new ResidualBasedBossakTurbulentKineticEnergyScheme<TSparseSpace, TDenseSpace>(
                alpha_bossak));

        ConvergenceCriteriaPointerType pKConvergenceCriteria = ConvergenceCriteriaPointerType(
            new GeneralConvergenceCriteria<TSparseSpace, TDenseSpace>(
                k_relative_tolerance, k_absolute_tolerance));

        pKConvergenceCriteria->SetEchoLevel(
            r_convergence_parameters["echo_level"].GetInt());

        mpKStrategy = typename StrategyType::Pointer(
            new ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
                *this->mpTurbulenceKModelPart, pKScheme, mpKLinearSolver,
                pKConvergenceCriteria, pKBuilderAndSolver, k_max_iterations,
                CalculateReactions, ReformDofSet, this->mIsMeshMoving));

        // Epsilon solution strategy
        BuilderAndSolverPointerType pEpsilonBuilderAndSolver = BuilderAndSolverPointerType(
            new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(
                mpEpsilonLinearSolver));

        mpEpsilonBuilderAndSolver = pEpsilonBuilderAndSolver;

        SchemePointerType pEpsilonScheme = SchemePointerType(
            new ResidualBasedBossakTurbulentEnergyDissipationRateScheme<TSparseSpace, TDenseSpace>(
                alpha_bossak));

        ConvergenceCriteriaPointerType pEpsilonConvergenceCriteria =
            ConvergenceCriteriaPointerType(new GeneralConvergenceCriteria<TSparseSpace, TDenseSpace>(
                epsilon_relative_tolerance, epsilon_absolute_tolerance));

        pEpsilonConvergenceCriteria->SetEchoLevel(
            r_convergence_parameters["echo_level"].GetInt());

        mpEpsilonStrategy = typename StrategyType::Pointer(
            new ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
                // new LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
                *this->mpTurbulenceEpsilonModelPart, pEpsilonScheme, mpEpsilonLinearSolver,
                pEpsilonConvergenceCriteria, pEpsilonBuilderAndSolver, epsilon_max_iterations,
                CalculateReactions, ReformDofSet, this->mIsMeshMoving));

        // mpConvergenceCriteria = ConvergenceCriteriaPointerType(
        //     new TurbulenceViscosityConvergenceCriteria<TSparseSpace, TDenseSpace>(
        //         turbulent_viscosity_relative_tolerance, turbulent_viscosity_absolute_tolerance));
        // mpConvergenceCriteria->Initialize(*mpTurbulenceModelPart);

        mpKBuilderAndSolver->SetEchoLevel(
            this->mrParameters["model_properties"]["echo_level"].GetInt());
        mpEpsilonBuilderAndSolver->SetEchoLevel(
            this->mrParameters["model_properties"]["echo_level"].GetInt());
        mpKStrategy->SetEchoLevel(
            this->mrParameters["model_properties"]["echo_level"].GetInt());
        mpEpsilonStrategy->SetEchoLevel(
            this->mrParameters["model_properties"]["echo_level"].GetInt());

        KRATOS_CATCH("");
    }

    void AssignBoundaryConditions()
    {
        const ProcessInfo& r_current_process_info = this->mrModelPart.GetProcessInfo();

        const double mixing_length = r_current_process_info[TURBULENT_MIXING_LENGTH];
        const double C_mu = r_current_process_info[TURBULENCE_RANS_C_MU];

        // Set Boundary Conditions
        const int NumNodes = this->mrModelPart.NumberOfNodes();
#pragma omp parallel for
        for (int i = 0; i < NumNodes; i++)
        {
            ModelPart::NodeIterator iNode = this->mrModelPart.NodesBegin() + i;
            if (iNode->Is(STRUCTURE))
            {
                this->SetUpWallNode(*iNode);
                this->InitializeValues(*iNode, 0.0, 0.0);
            }
            else if (iNode->Is(INLET))
            {
                array_1d<double, 3>& velocity =
                    iNode->FastGetSolutionStepValue(VELOCITY, 0);
                const double velocity_mag = norm_2(velocity);
                const double k = velocity_mag * mVelocityRatio;
                const double epsilon = C_mu * std::pow(k, 1.5) / mixing_length;

                this->SetUpFreestreamNode(*iNode);
                this->InitializeValues(*iNode, k, epsilon);
            }
            else // internal node, set initial value
            {
                array_1d<double, 3>& velocity =
                    iNode->FastGetSolutionStepValue(VELOCITY, 0);
                const double velocity_mag = norm_2(velocity);
                const double k = std::pow(velocity_mag / mixing_length, 2);
                const double epsilon = C_mu * std::pow(k, 1.5) / mixing_length;

                this->InitializeValues(*iNode, k, epsilon);
            }
        }
    }

    void SetUpWallNode(Node<3>& rNode)
    {
        rNode.Fix(TURBULENT_KINETIC_ENERGY);
        rNode.Fix(TURBULENT_ENERGY_DISSIPATION_RATE);
    }

    void SetUpFreestreamNode(Node<3>& rNode)
    {
        rNode.Fix(TURBULENT_KINETIC_ENERGY);
        rNode.Fix(TURBULENT_ENERGY_DISSIPATION_RATE);
    }

    void InitializeValues(Node<3>& rNode, double K, double Epsilon)
    {
        rNode.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = K;
        rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = Epsilon;
    }

    void AssignFrictionVelocities()
    {
        const ProcessInfo& r_current_process_info = this->mrModelPart.GetProcessInfo();

        const double von_karman = r_current_process_info[WALL_VON_KARMAN];
        const double beta = r_current_process_info[WALL_SMOOTHNESS_BETA];
        const double c_mu = r_current_process_info[TURBULENCE_RANS_C_MU];

        auto& nodes_array = this->mrModelPart.Nodes();
#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
        {
            auto it_node = nodes_array.begin() + i;
            it_node->FastGetSolutionStepValue(FRICTION_VELOCITY) = 0.0;
        }

        auto& element_array = this->mrModelPart.Elements();
        // #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(element_array.size()); ++i)
        {
            auto it_element = element_array.begin() + i;
            auto& r_geometry = it_element->GetGeometry();

            int wall_points = 0;
            int non_structure_index = 0;
            for (long unsigned int i = 0; i < r_geometry.PointsNumber(); ++i)
            {
                if (r_geometry[i].Is(STRUCTURE))
                    wall_points++;
                else
                    non_structure_index = static_cast<int>(i);
            }

            if (wall_points == TDim)
            {
                double& friction_velocity =
                    r_geometry[non_structure_index].FastGetSolutionStepValue(FRICTION_VELOCITY);
                // #pragma omp atomic
                friction_velocity = -1.0;
                // it_element->SetValue(STRUCTURE, true);
            }
        }

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
        {
            auto it_node = nodes_array.begin() + i;
            double& u_tau = it_node->FastGetSolutionStepValue(FRICTION_VELOCITY);
            if (!(it_node->Is(STRUCTURE)) && (u_tau == -1.0))
            {
                const array_1d<double, 3> velocity =
                    it_node->FastGetSolutionStepValue(VELOCITY);
                const double velocity_magnitude = norm_2(velocity);
                const double kinematic_viscosity =
                    it_node->FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
                const double wall_distance = it_node->FastGetSolutionStepValue(DISTANCE);
                u_tau = EvmKepsilonModelUtilities::CalculateUTau(
                    velocity_magnitude, wall_distance, kinematic_viscosity, beta, von_karman);
            }
            else if (!(it_node->Is(STRUCTURE)) && (u_tau == 0.0))
            {
                const double tke =
                    it_node->FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
                u_tau = std::pow(c_mu, 0.25) * std::sqrt(tke);
            }
        }
    }

    //     void AssignFrictionVelocities()
    //     {
    //         const ProcessInfo& r_current_process_info = this->mrModelPart.GetProcessInfo();

    //         const double von_karman = r_current_process_info[WALL_VON_KARMAN];
    //         const double beta = r_current_process_info[WALL_SMOOTHNESS_BETA];

    //         auto& nodes_array = this->mrModelPart.Nodes();
    // #pragma omp parallel for
    //         for (int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
    //         {
    //             auto it_node = nodes_array.begin() + i;
    //             it_node->FastGetSolutionStepValue(FRICTION_VELOCITY) = 0.0;
    //         }

    //         auto& element_array = this->mrModelPart.Elements();
    //         // #pragma omp parallel for
    //         for (int i = 0; i < static_cast<int>(element_array.size()); ++i)
    //         {
    //             auto it_element = element_array.begin() + i;
    //             auto& r_geometry = it_element->GetGeometry();

    //             int wall_points = 0;
    //             int non_structure_index = 0;
    //             for (long unsigned int i = 0; i < r_geometry.PointsNumber(); ++i)
    //             {
    //                 if (r_geometry[i].Is(STRUCTURE))
    //                     wall_points++;
    //                 else
    //                     non_structure_index = static_cast<int>(i);
    //             }

    //             if (wall_points == TDim)
    //             {
    //                 double& friction_velocity =
    //                     r_geometry[non_structure_index].FastGetSolutionStepValue(FRICTION_VELOCITY);
    //                 // #pragma omp atomic
    //                 friction_velocity = -1.0;
    //                 // it_element->SetValue(STRUCTURE, true);
    //             }
    //         }

    // // Calculate friction velocity for near wall nodes
    // #pragma omp parallel for
    //         for (int i = 0; i < static_cast<int>(nodes_array.size()); ++i)
    //         {
    //             auto it_node = nodes_array.begin() + i;
    //             double& u_tau = it_node->FastGetSolutionStepValue(FRICTION_VELOCITY);
    //             if (!(it_node->Is(STRUCTURE)) && (u_tau == -1.0))
    //             {
    //                 const array_1d<double, 3> velocity =
    //                     it_node->FastGetSolutionStepValue(VELOCITY);
    //                 const double velocity_magnitude = norm_2(velocity);
    //                 const double kinematic_viscosity =
    //                     it_node->FastGetSolutionStepValue(KINEMATIC_VISCOSITY);
    //                 const double wall_distance = it_node->FastGetSolutionStepValue(DISTANCE);
    //                 u_tau = EvmKepsilonModelUtilities::CalculateUTau(
    //                     velocity_magnitude, wall_distance, kinematic_viscosity, beta, von_karman);
    //             }
    //         }

    //         // #pragma omp parallel for
    //         bool filled_domain = true;
    //         while (filled_domain)
    //         {
    //             filled_domain = false;
    //             for (int i = 0; i < static_cast<int>(element_array.size()); ++i)
    //             {
    //                 auto it_element = element_array.begin() + i;
    //                 auto& r_geometry = it_element->GetGeometry();

    //                 double element_tau = -1.0;
    //                 for (long unsigned int i = 0; i < r_geometry.PointsNumber(); ++i)
    //                 {
    //                     double u_tau = r_geometry[i].FastGetSolutionStepValue(FRICTION_VELOCITY);
    //                     if (!(r_geometry[i].Is(STRUCTURE)) && (u_tau != 0.0))
    //                     {
    //                         element_tau = u_tau;
    //                         break;
    //                     }
    //                 }

    //                 if (element_tau != -1.0)
    //                 {
    //                     for (long unsigned int i = 0; i < r_geometry.PointsNumber(); ++i)
    //                     {
    //                         double& u_tau = r_geometry[i].FastGetSolutionStepValue(FRICTION_VELOCITY);
    //                         if (!(r_geometry[i].Is(STRUCTURE)) && (u_tau == 0.0))
    //                         {
    //                             u_tau = element_tau;
    //                             filled_domain = true;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }

    void SolveStep()
    {
        const int NumberOfNodes = this->mrModelPart.NumberOfNodes();

        this->AssignFrictionVelocities();

        this->UpdateTurbulentViscosity();

        mpKStrategy->Initialize();
        mpKStrategy->InitializeSolutionStep();
        mpKStrategy->Predict();

        mpEpsilonStrategy->Initialize();
        mpEpsilonStrategy->InitializeSolutionStep();
        mpEpsilonStrategy->Predict();

        Vector old_turbulent_viscosity = ZeroVector(NumberOfNodes);
        Vector new_turbulent_viscosity = ZeroVector(NumberOfNodes);

        bool is_converged = false;
        int step = 0;
        while (!is_converged && (step < mMaximumCouplingIterations))
        {
#pragma omp parallel for
            for (int i = 0; i < NumberOfNodes; ++i)
            {
                old_turbulent_viscosity[i] =
                    (this->mrModelPart.NodesBegin() + i)->FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            }

            if (this->mEchoLevel > 0)
                KRATOS_INFO("TurbulenceModel") << "Solving for K...\n";
            mpKStrategy->SolveSolutionStep();

            // this->CorrectTurbulentKineticEnergy();

            if (this->mEchoLevel > 0)
                KRATOS_INFO("TurbulenceModel") << "Solving for Epsilon...\n";
            mpEpsilonStrategy->SolveSolutionStep();

            this->UpdateTurbulentViscosity();

#pragma omp parallel for
            for (int i = 0; i < NumberOfNodes; ++i)
            {
                new_turbulent_viscosity[i] =
                    (this->mrModelPart.NodesBegin() + i)->FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            }

            double increase_norm(0.0);

#pragma omp parallel for reduction(+ : increase_norm)
            for (int i = 0; i < NumberOfNodes; ++i)
            {
                increase_norm += new_turbulent_viscosity[i] - old_turbulent_viscosity[i];
            }

            const double ratio = std::abs(increase_norm / NumberOfNodes);
            increase_norm = std::abs(increase_norm);

            if (this->mEchoLevel > 0)
                std::cout << "[" << step << "] CONVERGENCE CHECK: TURBULENT_VISCOSITY: ratio = "
                          << std::scientific << ratio << "; exp.ratio = " << std::scientific
                          << mTurbulentViscosityRelativeTolerance
                          << ": abs = " << std::scientific << increase_norm
                          << "; exp.abs = " << std::scientific
                          << mTurbulentViscosityAbsoluteTolerance << std::endl;

            is_converged = (increase_norm < mTurbulentViscosityAbsoluteTolerance) ||
                           (ratio < mTurbulentViscosityRelativeTolerance);

            if (this->mEchoLevel > 0 && is_converged)
                std::cout << "[" << step << "] "
                          << "CONVERGENCE CHECK: TURBULENT_VISCOSITY: "
                             "*** CONVERGENCE IS ACHIEVED ***\n";

            step++;
        }

        if (!is_converged)
        {
            std::cout
                << "|----------------------------------------------------|"
                << std::endl;
            std::cout << "|    " << BOLDFONT(FRED("ATTENTION: Max coupling iterations exceeded"))
                      << "     |" << std::endl;
            std::cout
                << "|----------------------------------------------------|"
                << std::endl;
        }

        mpKStrategy->FinalizeSolutionStep();
        mpEpsilonStrategy->FinalizeSolutionStep();
    }

    void CorrectTurbulentKineticEnergy()
    {
        const int NumNodes = this->mrModelPart.NumberOfNodes();

#pragma omp parallel for
        for (int i = 0; i < NumNodes; i++)
        {
            ModelPart::NodeIterator iNode = this->mrModelPart.NodesBegin() + i;
            double& tke = iNode->FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            tke = std::max<double>(0.0, tke);
        }
    }

    void UpdateTurbulentViscosity()
    {
        const int NumNodes = this->mrModelPart.NumberOfNodes();

        const ProcessInfo& r_current_process_info = this->mrModelPart.GetProcessInfo();

        const double C_mu = r_current_process_info[TURBULENCE_RANS_C_MU];
        const double nu_fraction = r_current_process_info[TURBULENT_VISCOSITY_FRACTION];
        const double mixing_length = r_current_process_info[TURBULENT_MIXING_LENGTH];

#pragma omp parallel for
        for (int i = 0; i < NumNodes; i++)
        {
            ModelPart::NodeIterator iNode = this->mrModelPart.NodesBegin() + i;
            const double tke = iNode->FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            const double epsilon =
                iNode->FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
            const double wall_distance = iNode->FastGetSolutionStepValue(DISTANCE);
            const double u_tau = iNode->FastGetSolutionStepValue(FRICTION_VELOCITY);
            const double nu = iNode->FastGetSolutionStepValue(KINEMATIC_VISCOSITY);

            const double y_plus =
                EvmKepsilonModelUtilities::CalculateYplus(u_tau, wall_distance, nu);
            const double f_mu = EvmKepsilonModelUtilities::CalculateFmu(y_plus);
            const double nu_min = nu_fraction * nu;

            const double nu_t = EvmKepsilonModelUtilities::CalculateTurbulentViscosity(
                C_mu, f_mu, tke, epsilon, mixing_length, nu_min);

            iNode->FastGetSolutionStepValue(TURBULENT_VISCOSITY) = nu_t;
        }
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TurbulenceEvmKEpsilonProcess& operator=(TurbulenceEvmKEpsilonProcess const& rOther);

    /// Copy constructor.
    TurbulenceEvmKEpsilonProcess(TurbulenceEvmKEpsilonProcess const& rOther);

    ///@}

}; // Class TurbulenceEvmKEpsilonProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator>>(
    std::istream& rIStream,
    TurbulenceEvmKEpsilonProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// output stream function
template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TurbulenceEvmKEpsilonProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_TURBULENCE_EVM_K_EPSILON_PROCESS_H_INCLUDED  defined