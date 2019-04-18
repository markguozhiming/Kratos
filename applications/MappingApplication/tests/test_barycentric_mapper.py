from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import basic_mapper_tests
import blade_mapping_test

class BarycentricBasicTestsLine(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "barycentric",
            "interpolation_type" : "line",
            "interface_submodel_part_origin": "line_tri",
            "interface_submodel_part_destination": "line_quad",
            "echo_level" : 0
        }""")
        super(BarycentricBasicTestsLine, cls).setUpMapper(mapper_params)

class BarycentricBasicTestsLineSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "barycentric",
            "interpolation_type" : "line",
            "interface_submodel_part_origin": "line_quad",
            "interface_submodel_part_destination": "line_tri",
            "echo_level" : 0
        }""")
        super(BarycentricBasicTestsLineSwitchedSides, cls).setUpMapper(mapper_params, switch_sides=True)

class BarycentricBasicTestsSurface(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "barycentric",
            "interpolation_type" : "triangle",
            "interface_submodel_part_origin": "surface_tri",
            "interface_submodel_part_destination": "surface_quad",
            "echo_level" : 0
        }""")
        super(BarycentricBasicTestsSurface, cls).setUpMapper(mapper_params)

class BarycentricBasicTestsSurfaceSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "barycentric",
            "interpolation_type" : "triangle",
            "interface_submodel_part_origin": "surface_quad",
            "interface_submodel_part_destination": "surface_tri",
            "echo_level" : 0
        }""")
        super(BarycentricBasicTestsSurfaceSwitchedSides, cls).setUpMapper(mapper_params, switch_sides=True)

class BarycentricBasicTestsVolume(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "barycentric",
            "interpolation_type" : "tetrahedra",
            "echo_level" : 0
        }""")
        super(BarycentricBasicTestsVolume, cls).setUpMapper(mapper_params)

class BarycentricBasicTestsVolumeSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "barycentric",
            "interpolation_type" : "tetrahedra",
            "echo_level" : 0
        }""")
        super(BarycentricBasicTestsVolumeSwitchedSides, cls).setUpMapper(mapper_params, switch_sides=True)

class BarycentricBladeMapping(blade_mapping_test.BladeMappingTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "barycentric",
            "interpolation_type" : "triangle",
            "echo_level" : 0
        }""")
        super(BarycentricBladeMapping, cls).setUpMapper(mapper_params)
        cls.print_output = False

if __name__ == '__main__':
    import sys
    if "--using-mpi" in sys.argv:
        from KratosMultiphysics import mpi # initialize MPI
        sys.argv.remove("--using-mpi") # has to be removed bcs Unittest cannot parse it
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()
