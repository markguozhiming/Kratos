import KratosMultiphysics as KM

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
if not CheckIfApplicationsAvailable("FSIApplication"):
    raise Exception("The FSIApplication needs to be compiled to use this mapper")
from KratosMultiphysics.FSIApplication import AdvancedNMPointsMapper

from KratosMultiphysics.MappingApplication.python_mapper import PythonMapper

def CreateMapper(model_part_origin, model_part_destination, mapper_settings):
    return ApproximateMortarMapper(model_part_origin, model_part_destination, mapper_settings)

class ApproximateMortarMapper(PythonMapper):
    # Wrapper for the AdvancedNMPointsMapper from the FSIApplication
    def __init(self, model_part_origin, model_part_destination, mapper_settings):
        super(ApproximateMortarMapper, self).__init__(model_part_origin, model_part_destination, mapper_settings)

        # TODO save additional settings needed by this mapper