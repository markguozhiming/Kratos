from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import co_simulation_tools as tools
try :
    import KratosMultiphysics
except ModuleNotFoundError:
    print(tools.bcolors.FAIL + 'KRATOS is not available ! Please ensure that Kratos is available for usage !'+ tools.bcolors.ENDC)
    exit()

try :
    from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
except ModuleNotFoundError:
    print(tools.bcolors.FAIL + 'KRATOS_FLUID_DYNAMICS is not available ! Please ensure that it is available for usage !'+ tools.bcolors.ENDC)
    exit()
# Importing the base class
from . import kratos_base_field_solver

def Create(name, cosim_solver_settings):
    return KratosCoSimulationFluidSolver(name, cosim_solver_settings)

class KratosCoSimulationFluidSolver(kratos_base_field_solver.KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        return FluidDynamicsAnalysis(self.model, self.project_parameters)

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__