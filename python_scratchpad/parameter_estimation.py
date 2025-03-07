# Load required standard modules
import numpy as np
from matplotlib import pyplot as plt

# Load required tudatpy modules
from tudatpy import constants
from tudatpy.interface import spice
from tudatpy import numerical_simulation
from tudatpy.numerical_simulation import environment
from tudatpy.numerical_simulation import environment_setup
from tudatpy.numerical_simulation import propagation_setup
from tudatpy.numerical_simulation import estimation, estimation_setup
from tudatpy.numerical_simulation.estimation_setup import observation
from tudatpy.astro.time_conversion import DateTime
from tudatpy.astro import element_conversion

# Load spice kernels
spice.load_standard_kernels()

# Set simulation start and end epochs
simulation_start_epoch = DateTime(2000, 1, 1).epoch()
simulation_end_epoch = DateTime(2000, 1, 4).epoch()

bodies_to_create = ["Sun", "Mars", "Jupiter"]

# Create default body settings for bodies_to_create
global_frame_origin = "Mars"
global_frame_orientation = "J2000"
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create, global_frame_origin, global_frame_orientation)

# Create empty body settings for the satellite
body_settings.add_empty_settings("MRO")

body_settings.get("MRO").constant_mass = 2000

# Create aerodynamic coefficient interface settings
reference_area_drag = 20  # Average projection area of a 3U CubeSat
drag_coefficient = 1.5
aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
    reference_area_drag, [drag_coefficient, 0.0, 0.0]
)

# Add the aerodynamic interface to the body settings
body_settings.get("MRO").aerodynamic_coefficient_settings = aero_coefficient_settings

# Create radiation pressure settings
reference_area_radiation = 20  # Average projection area of a 3U CubeSat
radiation_pressure_coefficient = 1.2
occulting_bodies_dict = dict()
occulting_bodies_dict["Sun"] = ["Mars"]
vehicle_target_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
    reference_area_radiation, radiation_pressure_coefficient, occulting_bodies_dict )


# Add the radiation pressure interface to the body settings
body_settings.get("MRO").radiation_pressure_target_settings = vehicle_target_settings

# Create system of bodies
bodies = environment_setup.create_system_of_bodies(body_settings)

# Define bodies that are propagated
bodies_to_propagate = ["MRO"]

# Define central bodies of propagation
central_bodies = ["Mars"]

# Define the accelerations acting on Delfi-C3
accelerations_settings_MRO = dict(
    Sun=[
        propagation_setup.acceleration.radiation_pressure(),
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Mars=[
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Jupiter=[
        propagation_setup.acceleration.point_mass_gravity()
    ]
)

# Create global accelerations dictionary
acceleration_settings = {"MRO": accelerations_settings_MRO}

# Create acceleration models
acceleration_models = propagation_setup.create_acceleration_models(
    bodies,
    acceleration_settings,
    bodies_to_propagate,
    central_bodies)