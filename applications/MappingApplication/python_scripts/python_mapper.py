import KratosMultiphysics as KM

# Base Class for all python-based mappers
# same interface as the C++ version (custom_mappers/mapper.h)
# see "custom_mappers/mapper.h" for more details

class PythonMapper(object):
    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        self.model_part_origin = model_part_origin
        self.model_part_destination = model_part_destination
        self.mapper_settings = mapper_settings

    def Map(variable_origin, variable_destination, mapper_flags=KM.Flags()):
        raise NotImplementedError

    def InverseMap(variable_origin, variable_destination, mapper_flags=KM.Flags()):
        raise NotImplementedError

    def UpdateInterface(mapper_flags=KM.Flags(), search_radius=0.0):
        raise NotImplementedError