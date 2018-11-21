//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// System includes

#if defined(KRATOS_PYTHON)

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "DEM_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(KratosDEMApplication,m)
{
    class_<KratosDEMApplication, KratosDEMApplication::Pointer, KratosApplication>(m, "KratosDEMApplication")
        .def(init<>());
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomConstitutiveLawsToPython(m);

    //Constitutive law
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER )

    //scheme
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_TRANSLATIONAL_INTEGRATION_SCHEME_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_ROTATIONAL_INTEGRATION_SCHEME_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER )

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PROBABILITY_DISTRIBUTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXCENTRICITY_PROBABILITY_DISTRIBUTION )

    // OPTIONS AND FLAGS
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TOP)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BOTTOM)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FORCE_INTEGRATION_GROUP)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TABLE_NUMBER_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TABLE_NUMBER_ANGULAR_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TABLE_NUMBER_FORCE)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TABLE_NUMBER_MOMENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TABLE_NUMBER) // JIG: To erase (1 January 2019)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BOUNDING_BOX_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ROTATION_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CRITICAL_TIME_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VIRTUAL_MASS_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SEARCH_CONTROL)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_TIME_TO_PRINT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COORDINATION_NUMBER)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAX_AMPLIFICATION_RATIO_OF_THE_SEARCH_RADIUS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CLEAN_INDENT_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TRIHEDRON_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ROLLING_FRICTION_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POISSON_EFFECT_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEIGH_INITIALIZED)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TRIAXIAL_TEST_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPUTE_STRESS_TENSOR_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FIX_VELOCITIES_FLAG)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTAINS_CLUSTERS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANDOM_ORIENTATION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_RESOLUTION_METHOD)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPUTE_FEM_RESULTS_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BREAKABLE_CLUSTER)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTINUUM_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FLOATING_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_ENGINE_POWER)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_MAX_ENGINE_FORCE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_THRESHOLD_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_ENGINE_PERFORMANCE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_DRAG_CONSTANT_X)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_DRAG_CONSTANT_Y)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_DRAG_CONSTANT_Z)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CLUSTER_INFORMATION) //Dangerous. Requires adding ClusterInformation to python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CLUSTER_FILE_NAME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INJECTOR_ELEMENT_TYPE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_GHOST)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_VELOCITY_X_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_VELOCITY_Y_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_VELOCITY_Z_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_ANGULAR_VELOCITY_X_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_ANGULAR_VELOCITY_Y_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INITIAL_ANGULAR_VELOCITY_Z_VALUE)

    // *************** Continuum only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DELTA_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CASE_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SKIN_SPHERE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_COHESION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_TENSION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COHESIVE_GROUP)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PROPERTIES_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_MESH_OPTION)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FAILURE_CRITERION_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONCRETE_TEST_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IF_BOUNDARY_ELEMENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_PRECONSOLIDATION_PRESSURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_M_CAMCLAY_SLOPE)
    // *************** Continuum only END ***************

    // MATERIAL PARAMETERS

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_MASS_COEFF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_MOMENT_OF_INERTIA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ROLLING_FRICTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ROLLING_FRICTION_WITH_WALLS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HISTORICAL_MIN_K)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_INERTIA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_DENSITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FRICTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_STATIC_FRICTION_COEF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_DYNAMIC_FRICTION_COEF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COEFFICIENT_OF_RESTITUTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_ROTATION_DAMP_RATIO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DAMPING_GAMMA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, K_NORMAL)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, K_TANGENTIAL)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_RADIUS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAX_STRESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, GAMMA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXCENTRICITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXCENTRICITY_STANDARD_DEVIATION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FABRIC_COEFFICIENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POISSON_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INTERNAL_COHESION)

    // *************** Continuum only BEGIN *************
    // *************** Dempack Constitutive Law only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SLOPE_FRACTION_N1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SLOPE_FRACTION_N2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SLOPE_FRACTION_N3)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SLOPE_LIMIT_COEFF_C1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SLOPE_LIMIT_COEFF_C2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SLOPE_LIMIT_COEFF_C3)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, YOUNG_MODULUS_PLASTIC)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PLASTIC_YIELD_STRESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DAMAGE_FACTOR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHEAR_ENERGY_COEF)
    // *************** Dempack Constitutive Law only END ***************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DONZE_G1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DONZE_G2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DONZE_G3)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DONZE_MAX_DEF)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_FAILURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_ORIENTATION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_SIGMA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_TAU)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FAILURE_CRITERION_STATE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, UNIDIMENSIONAL_DAMAGE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_SIGMA_MIN)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TENSION_LIMIT_INCREASE_SLOPE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_TAU_ZERO)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_INTERNAL_FRICC)
    // *************** Continuum only END *************

    // GEOMETRIC PARAMETERS

    // *************** Continuum only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_CONTACT_AREA_HIGH)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_CONTACT_AREA_LOW)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MEAN_CONTACT_AREA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, REPRESENTATIVE_VOLUME)
    // *************** Continuum only END ***************

    // INLET PARAMETERS

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INLET_START_TIME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INLET_STOP_TIME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INLET_NUMBER_OF_PARTICLES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, STANDARD_DEVIATION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAX_RAND_DEVIATION_ANGLE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IMPOSED_MASS_FLOW_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MASS_FLOW)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LINEAR_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, INLET_INITIAL_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, INLET_INITIAL_PARTICLES_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INLET_MAX_PARTICLES_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DENSE_INLET)

    // KINEMATICS

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, PARTICLE_ROTATION_ANGLE)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EULER_ANGLES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DOMAIN_IS_PERIODIC)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DOMAIN_MIN_CORNER)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DOMAIN_MAX_CORNER)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ORIENTATION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ORIENTATION_REAL) // JIG: SHOULD BE REMOVED IN THE FUTURE
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ORIENTATION_IMAG) // JIG: SHOULD BE REMOVED IN THE FUTURE
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, PRINCIPAL_MOMENTS_OF_INERTIA)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DELTA_DISPLACEMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DELTA_ROTA_DISPLACEMENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITY_START_TIME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITY_STOP_TIME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ANGULAR_VELOCITY_START_TIME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ANGULAR_VELOCITY_STOP_TIME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_BODY_MOTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FREE_BODY_MOTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_BODY_MASS)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RIGID_BODY_CENTER_OF_MASS)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RIGID_BODY_INERTIAS)
    // ****************** Quaternion Integration BEGIN ******************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUX_ORIENTATION)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LOCAL_AUX_ANGULAR_VELOCITY)
    // ******************* Quaternion Integration END *******************

    // FORCE AND MOMENTUM

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, PARTICLE_MOMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MAX_ROTA_MOMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ROLLING_RESISTANCE_MOMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ELASTIC_FORCES)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONTACT_FORCES)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RIGID_ELEMENT_FORCE)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, TANGENTIAL_ELASTIC_FORCES)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, FORCE_REACTION)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MOMENT_REACTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_PRESSURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_NODAL_AREA)

    // ENERGY

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_ELASTIC_ENERGY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_ROTATIONAL_KINEMATIC_ENERGY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_GRAVITATIONAL_ENERGY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_INELASTIC_VISCODAMPING_ENERGY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_INELASTIC_FRICTIONAL_ENERGY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPUTE_ENERGY_OPTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, GLOBAL_DAMPING)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NORMAL_IMPACT_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TANGENTIAL_IMPACT_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FACE_NORMAL_IMPACT_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FACE_TANGENTIAL_IMPACT_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LINEAR_IMPULSE)

    // *************** Continuum only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, INITIAL_ROTA_MOMENT)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_BLOCK_CONTACT_FORCE)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LOCAL_CONTACT_FORCE)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_CONTACT_FORCES )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEIGHBOUR_SIZE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEIGHBOUR_RATIO)

    // CONCRETE TEST

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FIXED_VEL_TOP)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FIXED_VEL_BOT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AREA_VERTICAL_TAPA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AREA_VERTICAL_CENTRE)

    // TENSION

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DEM_STRESS_TENSOR)

    // APPLIED LOADS

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BLAST_RADIUS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BLAST_CURVE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BLAST_PRESSURE_MAX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BLAST_TIME_PRESSURE_MAX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BLAST_SHAPE_FACTOR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BLAST_TIME_DELAY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BLAST_BOREHOLE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BLAST_NPOINTS)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BLAST_COORDINATES_1)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BLAST_COORDINATES_2)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BLAST_COORDINATES_3)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BLAST_COORDINATES_4)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BLAST_COORDINATES_5)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BLAST_COORDINATES_6)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BLAST_COORDINATES_7)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BLAST_COORDINATES_8)
    // *************** Continuum only END *************

    // DUMMIES INT AND DOUBLE VARIABLES
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DUMMY_SWITCH)

    //EXPORTS
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXPORT_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXPORT_PARTICLE_FAILURE_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRINT_EXPORT_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRINT_STRESS_TENSOR_OPTION)

    //For DEM_FEM Element
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, LOCAL_DAMP_RATIO)

    //SLS DEM-FEM
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLE_FRICTION) //deprecated since April 11th, 2018
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_FRICTION) //deprecated since April 11th, 2018
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHEAR_STRESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NON_DIMENSIONAL_VOLUME_WEAR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IMPACT_WEAR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SEVERITY_OF_WEAR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BRINELL_HARDNESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPUTE_WEAR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IMPACT_WEAR_SEVERITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WALL_COHESION)

    //DEM_CLUSTERS
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CLUSTER_VOLUME)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LOCAL_ANGULAR_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CHARACTERISTIC_LENGTH)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SPRAYED_MATERIAL)

    //BOUNDING BOX
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BOUNDING_BOX_START_TIME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BOUNDING_BOX_STOP_TIME)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_FACE_ROTA_SPEED)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_FACE_AXIAL_SPEED)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_FACE_PROP_ID)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RIGID_FACE_ROTA_ORIGIN_COORD)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RIGID_FACE_ROTA_AXIAL_DIR)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, RIGID_FACE_ROTA_GLOBAL_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_FACE_BEGIN_TIME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_FACE_END_TIME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RIGID_FACE_FLAG)

    //OPTIMIZATION

    // *************** Thermal only BEGIN *************
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HEATFLUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, THERMAL_CONDUCTIVITY)
    // *************** Thermal only END ***************

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
