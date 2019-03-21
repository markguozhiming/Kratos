import KratosMultiphysics
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle
from KratosMultiphysics.ParticleMechanicsApplication.apply_mpm_particle_dirichlet_condition_process import ApplyMPMParticleDirichletConditionProcess

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMCouplingInterfaceConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMCouplingInterfaceConditionProcess(ApplyMPMParticleDirichletConditionProcess):
    def __init__(self, Model, settings ):

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "particles_per_condition"   : 0,
                "imposition_type"           : "penalty",
                "penalty_factor"            : 0,
                "constrained"               : "fixed",
                "option"                    : ""
            }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)
        self.model = Model
        self.model_part_name = settings["model_part_name"].GetString()

        # Initiate base class - Dirichlet condition
        super(ApplyMPMCouplingInterfaceConditionProcess, self).__init__(Model, settings)

        # Set INTERFACE flag active
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, self.model_part.Conditions)


    def ExecuteBeforeSolutionLoop(self):
        # Get updated model_part
        if (self.model_part_name.startswith('Background_Grid.')):
            self.model_part_name = self.model_part_name.replace('Background_Grid.','')
        mpm_material_model_part_name = "MPM_Material." + self.model_part_name
        self.model_part = self.model[mpm_material_model_part_name]

        ### Translate conditions with INTERFACE flag into a new model part "MPM_Coupling_Interface" responsible for coupling with structure
        # Create coupling model part
        if not self.model.HasModelPart("MPM_Coupling_Interface"):
            self.model.CreateModelPart("MPM_Coupling_Interface")
        self.coupling_model_part = self.model.GetModelPart("MPM_Coupling_Interface").CreateSubModelPart(self.model_part_name)

        # Add variables to the coupling model part
        self._add_coupling_variables_to_model_part(self.coupling_model_part)

        # Create nodes and fill coupling model part
        for mpc in self.model_part.Conditions:
            if (mpc.Is(KratosMultiphysics.INTERFACE)):
                node_id         = mpc.Id
                node_coordinate = mpc.GetValue(KratosParticle.MPC_COORD)
                coupling_node   = self.coupling_model_part.CreateNewNode(node_id, node_coordinate[0], node_coordinate[1], node_coordinate[2])

                ## Set Displacement and Normal
                normal = mpc.GetValue(KratosParticle.MPC_NORMAL)
                coupling_node.SetValue(KratosMultiphysics.NORMAL,normal)


    def ExecuteInitializeSolutionStep(self):
        ### Send displacement from coupling_mp to mp
        for coupling_node in self.coupling_model_part.Nodes:
            coupling_id  = coupling_node.Id

            ## IMPOSED DISPLACEMENT
            total_displacement = coupling_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0)
            old_displacement = self.model_part.GetCondition(coupling_id).GetValue(KratosParticle.MPC_DISPLACEMENT)
            incremental_displacement = total_displacement - old_displacement
            self.model_part.GetCondition(coupling_id).SetValue(KratosParticle.MPC_IMPOSED_DISPLACEMENT,incremental_displacement)

            ## ADD VELOCITY
            current_velocity = coupling_node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            self.model_part.GetCondition(coupling_id).SetValue(KratosParticle.MPC_VELOCITY, current_velocity)

            ## ADD NORMAL
            normal = coupling_node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            self.model_part.GetCondition(coupling_id).SetValue(KratosParticle.MPC_NORMAL, normal)


    def ExecuteFinalizeSolutionStep(self):
        ### Get contact force from mp to coupling_mp
        for mpc in self.model_part.Conditions:
            if (mpc.Is(KratosMultiphysics.INTERFACE)):
                coupling_id   = mpc.Id
                contact_force = mpc.GetValue(KratosParticle.MPC_CONTACT_FORCE)
                self.coupling_model_part.GetNode(coupling_id).SetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE,contact_force)


    # Local functions
    def _add_coupling_variables_to_model_part(self, coupling_model_part):
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
        coupling_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
