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
    r = [100, 8400100, 6900, 0, 0, 0, PARAM.SPS_mass, PARAM.RCS_mass]
    #target = [8400000, 0, 0, -6900, 0, 0, PARAM.SPS_mass, PARAM.RCS_mass]
    target = [0, 8400000, 6910, 0, 0, 0, PARAM.SPS_mass, PARAM.RCS_mass]
    
    spacing = 1.0
    time_simulation = np.arange(0,3600*48,1)
    time_simulation = np.arange(0,160,spacing)
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
    
    del_SPS = np.arange(0,1.1,.2)
    del_RCS = np.arange(-1,1.1,.2)
    
    sens_radius = np.zeros((len(del_SPS),len(del_RCS)))
    sens_phase = np.zeros((len(del_SPS),len(del_RCS)))
    sens_angvel = np.zeros((len(del_SPS), len(del_RCS)))
    
    sol = odeint(motion, r, [0,spacing], args=(PARAM.CSM_mass,PARAM.CSM_I,[.2, 0]))

    for i in np.arange(len(del_SPS)):
        thrust[0] = del_SPS[i]
        for j in np.arange(len(del_RCS)):
            thrust[1] = del_RCS[j]
            sol = odeint(motion, r, [0,spacing], args=(PARAM.CSM_mass,PARAM.CSM_I,[thrust[0], thrust[1]]))
            result = sol[-1,:]
            radius = np.sqrt(result[0]**2 + result[1]**2)
            phase = np.arctan2(result[0],result[1])
            sens_radius[i,j] = radius
            sens_phase[i,j] = np.rad2deg(phase)
            sens_angvel[i,j] = np.rad2deg(result[5])
        
    
    sol = odeint(motion, r, [0,spacing], args=(PARAM.CSM_mass,PARAM.CSM_I,[0,0]))
    result = sol[-1,:]
    radius = np.sqrt(result[0]**2 + result[1]**2)
    phase = np.arctan2(result[0],result[1])
    
    sens_radius = sens_radius - radius
    sens_phase = sens_phase - np.rad2deg(phase)
    sens_angvel = sens_angvel - np.rad2deg(result[5])


            

    
#    fig, ax = plt.subplots()
#    #Plot Orbit
#    ax.plot(r_all[:,0],r_all[:,1],'r-',label="Chaser")
#    ax.plot(target_all[:,0], target_all[:,1],'b-',label="Target",alpha=.5)
#    #Plot Earth
#    earth = mpl.patches.Circle((0,0),6371000)
#    ax.add_artist(earth)
#    ax.legend()
#    #ax.scatter((0),[0],c='b')
#    #ax.set_aspect('equal', 'box')
#    
#    fig, ax = plt.subplots(nrows=2, ncols=2)
#    
#    distance = np.sqrt((r_all[:,0] - target_all[:,0])**2 + (r_all[:,1] - target_all[:,1])**2)
#    vel_diff = np.sqrt(r_all[:,2]**2 + r_all[:,3]**2) - np.sqrt(target_all[:,2]**2 + target_all[:,3]**2)
#    radius_r = np.sqrt(r_all[:,0]**2 + r_all[:,1]**2)
#    radius_t = np.sqrt(target_all[:,0]**2 + target_all[:,1]**2)
#    phase_r = np.rad2deg(np.arctan2(r_all[:,0], r_all[:,1]))
#    phase_t = np.rad2deg(np.arctan2(target_all[:,0], target_all[:,1]))
#    
#    ax[0,0].plot(full_time,distance)
#    ax[0,1].plot(full_time,vel_diff)
#    ax[1,0].plot(full_time,phase_r,label="Chaser")
#    ax[1,0].plot(full_time,phase_t,label="Target")
#    
#    ax[1,1].plot(full_time,radius_r,label="Chaser")
#    ax[1,1].plot(full_time,radius_t,label="Target")
#    
#    fig, ax = plt.subplots(2)
#    ax[0].plot(full_time, thrust_all[:,0])
#    
#    ax[1].plot(full_time, thrust_all[:,1])
    
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