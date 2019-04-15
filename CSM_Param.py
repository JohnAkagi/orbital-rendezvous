# -*- coding: utf-8 -*-
"""
Parameters related to the spacecraft
"""

# SPS - Service Propultions System
# Main engine on Apollo Service Module

SPS_mass = 18410 #kg propellants (https://en.wikipedia.org/wiki/Apollo_command_and_service_module)
SPS_thrust = 91000 #N (https://en.wikipedia.org/wiki/Apollo_command_and_service_module)

# Mass flow is calculated using the given thrust and specific impulse found at https://en.wikipedia.org/wiki/Apollo_command_and_service_module
# See https://en.wikipedia.org/wiki/Specific_impulse
SPS_mass_flow = 29.5 #kg/s 


# RCS - Reaction Control System
# Used for changing orientation of spacecraft
RCS_mass = 155 #kg total propellent, (https://en.wikipedia.org/wiki/Apollo_command_and_service_module)
RCS_mass_flow = .15 #kg/s, calculated using specific impulse found at http://www.astronautix.com/a/apollosm.html
RCS_thrust = 440 #N, (https://en.wikipedia.org/wiki/Apollo_command_and_service_module)


# CSM - Command and Service Module
# Physical dimensions of CSM

CSM_mass = 11900 #kg, dry mass (https://en.wikipedia.org/wiki/Apollo_command_and_service_module)
CSM_d = 3.9 #m, diameter of CSM (https://en.wikipedia.org/wiki/Apollo_command_and_service_module)
CSM_r = CSM_d/2 #m, radius of CSM

#CSM moment of inertia (https://www.hq.nasa.gov/alsj/SNA-8-D-027III-Rev2-CsmLmSpacecraftOperationalDataBook-Volume3-MassProperties.pdf)
# Since simulation is 2D, we use one of the two higher moments of inertia to spin the craft to change direction of thruster
CSM_I = 108465 # m^2kg

#Gravitional value for G*Earth mass
mu = 3.986e14