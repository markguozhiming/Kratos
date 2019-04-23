from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7
"""
This is a map of the name of the available IO types to be specified in
JSON file and their python module (file) name.
New IOs should be registered here with an additional entry.
eg : "name_in_JSON" : "python module(file) name"
"""
available_ios = {
    "kratos" : "kratos_interfaces.kratos_io",
    "su2"    : "su2_interface.su2_io",
    "dummy"  : "dummy_interface.dummy_co_simulation_io"
}

def CreateIO(model, io_type, settings):
    """
    This function creates and returns the IO used for CoSimulation
    New IOs have to be registered by adding them to "available_ios"
    """
    if io_type in available_ios:
        module_full_name = 'KratosMultiphysics.CoSimulationApplication.solver_interfaces.'+available_ios[io_type]
        io_module = __import__(module_full_name,fromlist=[available_ios[io_type]])
        return io_module.Create(model, settings)
    else:
        err_msg  = 'The requested IO "' + io_name + '" is not available!\n'
        err_msg += 'The available IOs are :\n'
        for avail_io in available_ios:
            err_msg += "\t" + avail_io + "\n"
        raise NameError(err_msg)
