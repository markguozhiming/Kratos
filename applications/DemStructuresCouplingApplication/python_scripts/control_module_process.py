import KratosMultiphysics
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ControlModuleProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ControlModuleProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]
        variable_name = settings["variable_name"].GetString()

        print('holaaaaaaaaaaaaaaaaaaaaaaaaaaaa')
        print(model_part)

        self.components_process_list = []

        if settings["fixed"][0].GetBool() == True:
            x_params = KratosMultiphysics.Parameters("{}")
            x_params.AddValue("model_part_name",settings["model_part_name"])
            x_params.AddEmptyValue("variable_name").SetString(variable_name+"_X")
            x_params.AddEmptyValue("reaction_variable_name").SetString(variable_name+"_X")
            x_params.AddValue("target_stress_table",settings["target_stress_table"][0])
            x_params.AddValue("initial_velocity",settings["initial_velocity"][0])
            x_params.AddValue("compression_length",settings["compression_length"])
            x_params.AddValue("young_modulus",settings["young_modulus"])
            x_params.AddValue("start_time",settings["start_time"])
            self.components_process_list.append(DemFem.ControlModuleProcess(model_part, x_params))

        if settings["fixed"][1].GetBool() == True:
            y_params = KratosMultiphysics.Parameters("{}")
            y_params.AddValue("model_part_name",settings["model_part_name"])
            y_params.AddEmptyValue("variable_name").SetString(variable_name+"_Y")
            y_params.AddEmptyValue("reaction_variable_name").SetString(variable_name+"_Y")
            y_params.AddValue("target_stress_table",settings["target_stress_table"][1])
            y_params.AddValue("initial_velocity",settings["initial_velocity"][1])
            y_params.AddValue("compression_length",settings["compression_length"])
            y_params.AddValue("young_modulus",settings["young_modulus"])
            y_params.AddValue("start_time",settings["start_time"])
            self.components_process_list.append(DemFem.ControlModuleProcess(model_part, y_params))

        if settings["fixed"][2].GetBool() == True:
            z_params = KratosMultiphysics.Parameters("{}")
            z_params.AddValue("model_part_name",settings["model_part_name"])
            z_params.AddEmptyValue("variable_name").SetString(variable_name+"_Z")
            z_params.AddEmptyValue("reaction_variable_name").SetString(variable_name+"_Z")
            z_params.AddValue("target_stress_table",settings["target_stress_table"][2])
            z_params.AddValue("initial_velocity",settings["initial_velocity"][2])
            z_params.AddValue("compression_length",settings["compression_length"])
            z_params.AddValue("young_modulus",settings["young_modulus"])
            z_params.AddValue("start_time",settings["start_time"])
            self.components_process_list.append(DemFem.ControlModuleProcess(model_part, z_params))

    def ExecuteInitialize(self):

        for component in self.components_process_list:
            component.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        for component in self.components_process_list:
            component.ExecuteFinalizeSolutionStep()