from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
import swimming_DEM_solver
from swimming_DEM_solver import Say
BaseSolver = swimming_DEM_solver.SwimmingDEMSolver

class InterpolationTestSolver(BaseSolver):
    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        super(InterpolationTestSolver, self).__init__(model,
                                                      project_parameters,
                                                      field_utility,
                                                      fluid_solver,
                                                      dem_solver,
                                                      variables_manager)

    def ReturnExactVelocity(self, t, x, y, z):
        interpolate_process_data = self.project_parameters['processes']['check_interpolated_fluid_velocity'][0]
        interpolate_process_parameters = interpolate_process_data['Parameters']
        field_def = [entry.GetString() for entry in interpolate_process_parameters['value']]
        field = [eval(field_def[0]),
                 eval(field_def[1]),
                 eval(field_def[2])]

        return Vector(field)

    def CannotIgnoreFluidNow(self):
        return self.calculating_fluid_in_current_step

    def SolveFluidSolutionStep(self):

        for node in self.fluid_solver.main_model_part.Nodes:
            velocity = self.ReturnExactVelocity(self.next_time_to_solve_fluid, node.X, node.Y, node.Z)
            node.SetSolutionStepValue(VELOCITY, velocity)

    def SolveDEM(self):
        super(InterpolationTestSolver, self).SolveDEM()

    def SolveDEMSolutionStep(self):
        pass

