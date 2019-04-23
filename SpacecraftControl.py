# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 11:43:03 2019

@author: johna
"""

from gekko_model import gekko_model
import numpy as np
from OrbitalMotion import motion, estimator, opt_estimator
import matplotlib.pyplot as plt
import matplotlib as mpl
import CSM_Param as PARAM
from scipy.integrate import odeint

#Test the model
if __name__ == "__main__":
    #Initial state
    r = [0, 8410000, 6900, 0, 0, 0, PARAM.SPS_mass, PARAM.RCS_mass]
    #target = [8400000, 0, 0, -6900, 0, 0, PARAM.SPS_mass, PARAM.RCS_mass]
    target = [0, 8400000, 6910, 0, 0, 0, PARAM.SPS_mass, PARAM.RCS_mass]
    
    spacing = 1.0
    time_simulation = np.arange(0,3600*48,1)
    time_simulation = np.arange(0,600,spacing)
    time_n = len(time_simulation)
    full_time = time_simulation
    
    thrust = [0, 0] #SPS and RCS
    
    r_all = np.empty((time_n, len(r)))
    target_all = np.empty((time_n, len(target)))
    thrust_all = np.empty((time_n, len(thrust)))
    estimate_all = np.empty((time_n, 3))
    mhe_all = np.empty((time_n, 3))
    bias_estimate_all = np.empty((time_n, 2))
    noise_all = np.empty((time_n, 2))
    estimate_2 = np.empty((time_n,3))
    model_all = np.empty((time_n,len(r)))
    
    r_all[0,:] = r
    target_all[0,:] = target
    
    param = [PARAM.CSM_mass, PARAM.SPS_thrust, PARAM.RCS_thrust, PARAM.SPS_mass_flow, PARAM.RCS_mass_flow,
             PARAM.CSM_I, PARAM.CSM_r]
    physics = [PARAM.mu]
    
    model = gekko_model(r, param, physics, target)
    #mpc_time = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 30, 60, 90, 180, 300, 600, 1200, 3600, 7200]
    mpc_time = np.arange(0,30,spacing)
    
    sim = gekko_model(r, param, physics, target)
    sim_time = [0., 1.]
    sim.sim_setup(sim_time)

    model.mpc_setup(mpc_time)
    

    

    
    for time in np.arange(len(time_simulation)):
        
        model.mpc_update(r, target)
        thrust = model.mpc_solve()
        r_all[time, :] = r
        target_all[time,:] = target
        thrust_all[time, :] = thrust
        
        sol = odeint(motion, r, [0,spacing], args=(PARAM.CSM_mass,PARAM.CSM_I,[thrust[0], thrust[1]]))
        r = sol[-1,:]
        
        sol = odeint(motion, target, [0,spacing], args=(PARAM.CSM_mass,PARAM.CSM_I,[0, 0]))
        target = sol[-1,:]

        #model.mpc_update(r, target)
        if np.mod(time, 10) == 0:
            print("{:.3f}%".format(time*100/len(time_simulation)))
            

    
    fig, ax = plt.subplots()
    #Plot Orbit
    ax.plot(r_all[:,0],r_all[:,1],'r-',label="Chaser")
    ax.plot(target_all[:,0], target_all[:,1],'b-',label="Target",alpha=.5)
    #Plot Earth
    earth = mpl.patches.Circle((0,0),6371000)
    ax.add_artist(earth)
    ax.legend()
    #ax.scatter((0),[0],c='b')
    #ax.set_aspect('equal', 'box')
    
    fig, ax = plt.subplots(nrows=2, ncols=2)
    
    distance = np.sqrt((r_all[:,0] - target_all[:,0])**2 + (r_all[:,1] - target_all[:,1])**2)
    vel_diff = np.sqrt(r_all[:,2]**2 + r_all[:,3]**2) - np.sqrt(target_all[:,2]**2 + target_all[:,3]**2)
    radius_r = np.sqrt(r_all[:,0]**2 + r_all[:,1]**2)
    radius_t = np.sqrt(target_all[:,0]**2 + target_all[:,1]**2)
    phase_r = np.rad2deg(np.arctan2(r_all[:,0], r_all[:,1]))
    phase_t = np.rad2deg(np.arctan2(target_all[:,0], target_all[:,1]))
    
    ax[0,0].plot(full_time,distance)
    ax[0,1].plot(full_time,vel_diff)
    ax[1,0].plot(full_time,phase_r-phase_t,label="Chaser")
    #ax[1,0].plot(full_time,phase_t,label="Target",alpha=.5)
    
    ax[1,1].plot(full_time,radius_r,label="Chaser")
    ax[1,1].plot(full_time,radius_t,label="Target")
    
    fig, ax = plt.subplots(2)
    ax[0].plot(full_time, thrust_all[:,0])
    ax[0].set_ylabel("Thrust Fraction")
    ax[0].set_title("Main Thruster")
    
    ax[1].plot(full_time, thrust_all[:,1])
    ax[1].set_ylabel("Thrust Fraction")
    ax[1].set_xlabel("Time (s)")
    ax[1].set_title("Control Thruster")
    
    fig.tight_layout()
    
    fig, ax = plt.subplots(2)
    ax[0].plot(full_time,distance)
    ax[0].set_ylabel("Separation (m)")
    ax[0].set_xlabel("Time (s)")
    ax[0].set_title("Separation")
    
    ax[1].plot(full_time, vel_diff)
    ax[1].set_ylabel("Velocity (m/s)")
    ax[1].set_xlabel("Time (s)")
    ax[1].set_title("Relative Velocity")
    
    fig.tight_layout()
    
    fig, ax = plt.subplots()
    ax.plot(full_time,radius_r,'r',label="Chaser",alpha=.75)
    ax.plot(full_time,radius_t,'b',label="Target",alpha=.75)
    ax.set_ylabel('Radius (m)')
    ax.set_xlabel('Time (s)')
    ax.legend()
    fig.tight_layout()
    
    fig, ax = plt.subplots(1)
    ax.plot(full_time,phase_r-phase_t)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Angle (degree)")
    ax.set_title("Angular Separation")
    
    fig, ax = plt.subplots(1)
    ax.plot(full_time,r_all[:,5])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Angular Velocity")
    ax.set_title("Angular Velocity")
    
    
#    fig, ax = plt.subplots(nrows=1, ncols=2)
#    
#    ax[0].plot(full_time, r_all[:,6], 'k', label='Actual',linewidth=2.0, alpha=.5)
#    ax[0].plot(full_time, model_all[:,6],'r:', label='Estimate',linewidth=2.0)
#    ax[0].set_ylabel("Fuel (kg)")
#    ax[0].set_xlabel("Time (s)")
#    ax[0].set_title("Main Thrusters")
#    ax[0].legend()
#    
#    
#    ax[1].plot(full_time, r_all[:,7], 'k', label='Acutal',linewidth=2.0, alpha=.5)
#    ax[1].plot(full_time, model_all[:,7],'r:',label='Estimate',linewidth=2.0)
#    ax[1].set_ylabel("Fuel (kg)")
#    ax[1].set_xlabel("Time (s)")
#    ax[1].set_title("Control Thrusters")
#    ax[1].legend()
#    
#    fig.tight_layout()
#    
#    fig2, ax2 = plt.subplots(1)
#    ax2.plot(full_time, bias_estimate_all[:,1],'r', label="RCS Estimate")
#    ax2.plot(full_time, bias_estimate_all[:,0],'b', label="SPS Estimate")
#    ax2.plot(full_time, noise_all[:,1],'r:', label="RCS Input")
#    ax2.plot(full_time, noise_all[:,0],'b:', label="SPS Input")
#    ax2.plot([0, full_time[-1]],[bias[1], bias[1]], 'r-', label="RCS Nominal", alpha=.5)
#    ax2.plot([0, full_time[-1]],[bias[0], bias[0]], 'b-', label="SPS Nominal", alpha=.5)
#    ax2.legend()
#    ax2.set_ylabel('Input Noise')
#    ax2.set_xlabel('Time (s)')
#    ax2.set_ylim([.7,1.3])