from __future__ import print_function, absolute_import, division

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_coupled_solver import CoSimulationBaseCouplingSolver

# Other imports
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import couplingsolverprint, red, green, cyan, bold
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(model, cosim_solver_settings, solver_name):
    return GaussSeidelStrongCouplingSolver(model, cosim_solver_settings, solver_name)

class GaussSeidelStrongCouplingSolver(CoSimulationBaseCouplingSolver):
    def __init__(self, model, cosim_solver_settings, solver_name):
        if not cosim_solver_settings['coupling_sequence'].size() == 2:
            raise Exception("Exactly two solvers have to be specified for the " + self.__class__.__name__ + "!")

        super(GaussSeidelStrongCouplingSolver, self).__init__(model, cosim_solver_settings, solver_name)

        self.convergence_accelerators_list = cs_tools.CreateConvergenceAccelerators(
            self.cosim_solver_settings["convergence_accelerators"],
            self.participating_solvers,
            self.echo_level)

        self.convergence_criteria_list = cs_tools.CreateConvergenceCriteria(
            self.cosim_solver_settings["convergence_criteria"],
            self.participating_solvers,
            self.echo_level)

        self.num_coupling_iterations = self.cosim_solver_settings["num_coupling_iterations"].GetInt()

    def Initialize(self):
        super(GaussSeidelStrongCouplingSolver, self).Initialize()
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.Initialize()
        for conv_crit in self.convergence_criteria_list:
            conv_crit.Initialize()

    def Finalize(self):
        super(GaussSeidelStrongCouplingSolver, self).Finalize()
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.Finalize()
        for conv_crit in self.convergence_criteria_list:
            conv_crit.Finalize()

    def InitializeSolutionStep(self):
        super(GaussSeidelStrongCouplingSolver, self).InitializeSolutionStep()
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.InitializeSolutionStep()
        for conv_crit in self.convergence_criteria_list:
            conv_crit.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(GaussSeidelStrongCouplingSolver, self).FinalizeSolutionStep()
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.FinalizeSolutionStep()
        for conv_crit in self.convergence_criteria_list:
            conv_crit.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        for k in range(self.num_coupling_iterations):
            if self.echo_level > 0:
                couplingsolverprint(self._Name(),
                                    cyan("Coupling iteration:"), bold(str(k+1)+" / " + str(self.num_coupling_iterations)))

            for conv_acc in self.convergence_accelerators_list:
                conv_acc.InitializeNonLinearIteration()
            for conv_crit in self.convergence_criteria_list:
                conv_crit.InitializeNonLinearIteration()

            for solver_name, solver in self.participating_solvers.items():
                self._SynchronizeInputData(solver_name)
                solver.SolveSolutionStep()
                self._SynchronizeOutputData(solver_name)

            for conv_acc in self.convergence_accelerators_list:
                conv_acc.FinalizeNonLinearIteration()
            for conv_crit in self.convergence_criteria_list:
                conv_crit.FinalizeNonLinearIteration()

            is_converged = True
            for conv_crit in self.convergence_criteria_list:
                is_converged = is_converged and conv_crit.IsConverged()
            if is_converged:
                if self.echo_level > 0:
                    couplingsolverprint(self._Name(), green("### CONVERGENCE WAS ACHIEVED ###"))
                break
            else:
                for conv_acc in self.convergence_accelerators_list:
                    conv_acc.ComputeUpdate()

            if k+1 >= self.num_coupling_iterations and self.echo_level > 0:
                couplingsolverprint(self._Name(), red("XXX CONVERGENCE WAS NOT ACHIEVED XXX"))

    def PrintInfo(self):
        super(GaussSeidelStrongCouplingSolver, self).PrintInfo()

        couplingsolverprint(self._Name(), "Uses the following objects:")
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.PrintInfo()
        for conv_crit in self.convergence_criteria_list:
            conv_crit.PrintInfo()

    def Check(self):
        super(GaussSeidelStrongCouplingSolver, self).Check()
        for conv_acc in self.convergence_accelerators_list:
            conv_acc.Check()
        for conv_crit in self.convergence_criteria_list:
            conv_crit.Check()

    def _Name(self):
        return self.__class__.__name__