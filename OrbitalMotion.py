# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:53:17 2019

@author: johna
"""

import numpy as np
import CSM_Param as PARAM
from scipy.optimize import minimize
from scipy.integrate import odeint


def motion(state, t, mass, J, thrust):
    """
    This function propogates forward the spacecraft in conjunction with ODEINT
    
    Parameters:
        state (array of floats): x position, y position, x velocity, y velocity, heading, angular velocity, SPS propellant mass, RCS propellant mass
        t (float): time
        mass (floats): spacecraft dry mass
        J (float): moment of inertia
        thrust (array of floats): SPS thrust (range of 0 to 1), RCS thrust (range -1 to 1)
        bias (list of floats): Decimal describing how far off of nominal the thrusts are
    Positions in meters, velocities in meters/second
    Positions are relative to the center of the Earth
    Angular velocity in radians/second
    Heading is in radians, y axis is zero, positive towards x axis
    Mass is in kg
    SPS is main thrust in direction of heading
    RCS is rotation thrust
    """
    
    #Initialize derivative vector
    state_dot = np.zeros(8)
    
    #Seperate attitude measurements
    pos = state[0:2]
    pos_x = state[0]
    pos_y = state[1]
    vel_x = state[2]
    vel_y = state[3]
    heading = state[4]
    ang_vel = state[5]
    SPS_mass = state[6]
    RCS_mass = state[7]
    
    
    SPS_command = thrust[0]
    RCS_command = thrust[1]
    
    
    total_mass = SPS_mass + RCS_mass + mass
    
    #Find distance between Earth and Spacecraft
    d = np.linalg.norm(pos)
    
    #Calculate magnitude of main thruster
    SPS_thrust_force = SPS_command*PARAM.SPS_thrust
    #Calculate magnitude of RCS thruster
    #The RCS thruster fire in pairs, hence the factor of 2 added
    RCS_thrust_force = RCS_command*PARAM.RCS_thrust*2
    
    #Calculate acceleration due to SPS
    x_accel_SPS = SPS_thrust_force*np.sin(heading)/total_mass
    y_accel_SPS = SPS_thrust_force*np.cos(heading)/total_mass
    
    #Calculate acceleration due to RCS
    #Torque is FxR
    torque = RCS_thrust_force * PARAM.CSM_r
    ang_accel = torque/PARAM.CSM_I
    
    #Calculate mass flow of thrusters
    SPS_flow = -PARAM.SPS_mass_flow*SPS_command
    RCS_flow = -PARAM.RCS_mass_flow*np.abs(RCS_command)
    
    #If there is no propellant mass, there is no thrust, acceleration, etc.
    if SPS_mass <= 0:
        SPS_thrust_force = 0
        SPS_mass = 0
        SPS_flow = 0
        x_accel_SPS = 0
        y_accel_SPS = 0
    if RCS_mass <= 0:
        RCS_thrust_force = 0
        RCS_flow = 0
        RCS_mass = 0
        ang_accel = 0
    
    #Calculate acceleration due to gravity
    r_dd = -PARAM.mu/(d**3)
    x_accel_g = pos_x*r_dd
    y_accel_g = pos_y*r_dd
    
    
    #Update derivatives
    state_dot[0] = vel_x
    state_dot[1] = vel_y
    state_dot[2] = x_accel_g + x_accel_SPS
    state_dot[3] = y_accel_g + y_accel_SPS
    state_dot[4] = ang_vel
    state_dot[5] = ang_accel
    state_dot[6] = SPS_flow
    state_dot[7] = RCS_flow

    
    return(state_dot)
    
    
#def motion_fsolve(bias, state, t, mass, J, thrust):
#    """
#    This function propogates forward the spacecraft in conjunction with ODEINT
#    
#    Parameters:
#        state (array of floats): x position, y position, x velocity, y velocity, heading, angular velocity, SPS propellant mass, RCS propellant mass
#        t (float): time
#        mass (floats): spacecraft dry mass
#        J (float): moment of inertia
#        thrust (array of floats): SPS thrust (range of 0 to 1), RCS thrust (range -1 to 1)
#        bias (list of floats): Decimal describing how far off of nominal the thrusts are
#    Positions in meters, velocities in meters/second
#    Positions are relative to the center of the Earth
#    Angular velocity in radians/second
#    Heading is in radians, y axis is zero, positive towards x axis
#    Mass is in kg
#    SPS is main thrust in direction of heading
#    RCS is rotation thrust
#    """
#    
#    #Initialize derivative vector
#    state_dot = np.zeros(8)
#    
#    #Seperate attitude measurements
#    pos = state[0:2]
#    pos_x = state[0]
#    pos_y = state[1]
#    vel_x = state[2]
#    vel_y = state[3]
#    heading = state[4]
#    ang_vel = state[5]
#    SPS_mass = state[6]
#    RCS_mass = state[7]
#    
#    
#    SPS_command = thrust[0]*(bias[0])
#    RCS_command = thrust[1]*(bias[1])
#    
#    
#    total_mass = SPS_mass + RCS_mass + mass
#    
#    #Find distance between Earth and Spacecraft
#    d = np.linalg.norm(pos)
#    
#    #Calculate magnitude of main thruster
#    SPS_thrust_force = SPS_command*PARAM.SPS_thrust
#    #Calculate magnitude of RCS thruster
#    #The RCS thruster fire in pairs, hence the factor of 2 added
#    RCS_thrust_force = RCS_command*PARAM.RCS_thrust*2
#    
#    #Calculate acceleration due to SPS
#    x_accel_SPS = SPS_thrust_force*np.sin(heading)/total_mass
#    y_accel_SPS = SPS_thrust_force*np.cos(heading)/total_mass
#    
#    #Calculate acceleration due to RCS
#    #Torque is FxR
#    torque = RCS_thrust_force * PARAM.CSM_r
#    ang_accel = torque/PARAM.CSM_I
#    
#    #Calculate mass flow of thrusters
#    SPS_flow = -PARAM.SPS_mass_flow*SPS_command
#    RCS_flow = -PARAM.RCS_mass_flow*np.abs(RCS_command)
#    
#    #If there is no propellant mass, there is no thrust, acceleration, etc.
#    if SPS_mass <= 0:
#        SPS_thrust_force = 0
#        SPS_mass = 0
#        SPS_flow = 0
#        x_accel_SPS = 0
#        y_accel_SPS = 0
#    if RCS_mass <= 0:
#        RCS_thrust_force = 0
#        RCS_flow = 0
#        RCS_mass = 0
#        ang_accel = 0
#    
#    #Calculate acceleration due to gravity
#    r_dd = -PARAM.mu/(d**3)
#    x_accel_g = pos_x*r_dd
#    y_accel_g = pos_y*r_dd
#    
#    #Update derivatives
#    state_dot[0] = vel_x
#    state_dot[1] = vel_y
#    state_dot[2] = x_accel_g + x_accel_SPS
#    state_dot[3] = y_accel_g + y_accel_SPS
#    state_dot[4] = ang_vel
#    state_dot[5] = ang_accel
#    state_dot[6] = SPS_flow
#    state_dot[7] = RCS_flow
#
#    
#    return(state_dot)
    
def estimator_main(bias,starting_point,commands,match_data,spacing,init_bias):

    sim_data = np.empty((len(commands),len(match_data[0])))
    sol = starting_point
    for time in np.arange(len(commands)):
        sol_temp = odeint(motion, sol, [0,spacing], args=(PARAM.CSM_mass,PARAM.CSM_I,[commands[time,0]*bias[0], commands[time,1]*bias[1]]))
        sol = sol_temp[-1,:]
        sim_data[time,:] = [sol[0], sol[1], sol[5]]
    error = np.sum(((match_data - sim_data)**2)@np.array([[1],[1],[1e3]])) + (1000*np.sum((bias - init_bias)**2))

    return(error)
    
def estimator(match_data, starting_point, commands, bias0, spacing):
    
    bias = minimize(estimator_main, bias0, args=(starting_point,commands,match_data,spacing,bias0),bounds=[(bias0[0]-.02,bias0[0]+.02),(bias0[1]-.02,bias0[1]+.02)])

    return(bias.x)
    
def min_change(bias0, match_data, starting_point, commands, spacing):
    
    bias = minimize(estimator_main, bias0, args=(starting_point,commands,match_data,spacing),bounds=[(.5,1.5),(.5,1.5)])
    return(1000*np.sum((bias.x - bias0)**2))
    
def opt_estimator(match_data, starting_point, commands, bias0, spacing):
    
    bias = minimize(min_change, bias0, args=(match_data, starting_point,commands,spacing),bounds=[(.5,1.5),(.5,1.5)])
    return(bias.x)
    
    
    
    
    