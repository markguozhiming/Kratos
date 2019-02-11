//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "includes/define_python.h"
#include "includes/kratos_parameters.h"

#include "spaces/ublas_space.h"

// strategies
#include "custom_strategies/residual_based_bossak_velocity_scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"

// linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
namespace Python
{
void AddCustomStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    // typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    // typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> BaseSolvingStrategyType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    //********************************************************************
    //********************************************************************
    typedef ResidualBasedBossakVelocityScheme<SparseSpaceType, LocalSpaceType> ResidualBasedBossakVelocityScheme;
    py::class_<ResidualBasedBossakVelocityScheme, typename ResidualBasedBossakVelocityScheme::Pointer, BaseSchemeType>(
        m, "ResidualBasedBossakVelocityScheme")
        .def(py::init<Kratos::Parameters>());
}

} // namespace Python.
} // Namespace Kratos
