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
    r = [0, 8400010, 6900, 0, 0, 0, PARAM.SPS_mass, PARAM.RCS_mass]
    target = [8400000, 0, 0, -6900, 0, 0, PARAM.SPS_mass, PARAM.RCS_mass]
    target = [0, 8400000, 6910, 0, 0, 0, PARAM.SPS_mass, PARAM.RCS_mass]
    
    bias = [1.05, .95]
    
    spacing = 1.0
    time_simulation = np.arange(0,3600*48,1)
    time_simulation = np.arange(0,160,spacing)
#    free_fall_time = np.arange(0,0,spacing)
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
    
    #model = gekko_model(r, param, physics, target)
    #mpc_time = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 30, 60, 90, 180, 300, 600, 1200, 3600, 7200]
    #mpc_time = np.arange(0,30,spacing)
    
    #mhe = gekko_model(r, param, physics, target)
    #mhe_time = np.arange(0,60,spacing)
    #mhe.mhe_setup(mhe_time)
    
    mhe_time = np.arange(0,30,spacing)
    
    sim = gekko_model(r, param, physics, target)
    sim_time = [0., 1.]
    sim.sim_setup(sim_time)
    
    sim_per = gekko_model(r, param, physics, target)
    sim_time = [0., 1.]
    sim_per.sim_setup(sim_time)
    
    commands = np.zeros((len(mhe_time),2))
    start_point = np.zeros((len(mhe_time),8))
    match_data = np.zeros((len(mhe_time),3))
    
    #mpc_time = np.append(mpc_time,[100, 200, 300])
    #model.mpc_setup(mpc_time)
    #model.sim_setup(full_time)
    
    #model.sim_solve(np.zeros(len(full_time)), np.zeros(len(full_time)))
    
    state_carryover = [r[0], r[1], r[5]]
    
    bias_guess = [1., 1.]
    
    model_r = r
    
    for time in np.arange(len(time_simulation)):
        
        
        #thrust = model.mpc_solve()
        #mhe.mhe_on()
        if time > 0:
            thrust = [1., .25]
        if time > 30:
            thrust = [.25, 1.]
        if time > 60:
            thrust = [.25, -1.]
        if time > 90:
            thrust = [.5, 0.25]
        if time > 120:
            thrust = [.75, -.6]
#        if time > 50:
#            thrust = [.75, .75]
        r_all[time, :] = r
        model_all[time, :] = model_r
        target_all[time,:] = target
        thrust_all[time, :] = thrust
        
        #thrust = [1.,0.]
        

        
        
        
        noise = np.random.normal(scale=.05,size=(2))
        
        noise_total = bias + noise
        noise_all[time, :] = noise_total[0]
        
        #estimate_all[time,:] = sim.sim_solve(thrust[0]*noise_total[0], thrust[1]*noise_total[1])
        #state_carryover = estimate_all[time,:]
        
        #estimate_2[time,:] = sim_per.sim_solve(thrust[0], thrust[1])
        
        #if time > mhe_time[-1]:
            #bias_model, state_model = mhe.mhe_solve(state_carryover[0], state_carryover[1], state_carryover[2], thrust[0], thrust[1])
            #mhe_all[time,:] = state_model
            #bias_estimate_all[time,:] = bias_model
            
            
        #elif time < mhe_time[-1]:
            
            #measurements.append([sim.x.MODEL, sim.y.MODEL, sim.x_vel.MODEL ,sim.y_vel.MODEL, sim.heading.MODEL, sim.ang_vel.MODEL, sim.SPS_mass.VALUE, sim.RCS_mass.VALUE])
            #commands.append(thrust)
#            
#            mhe.mhe_add(state_carryover[0], state_carryover[1], state_carryover[2], thrust[0], thrust[1])
#            mhe_all[time,:] = [state_carryover[0], state_carryover[1], state_carryover[2]]
#            bias_estimate_all[time,:] = bias_model
            
#            mhe.x.VALUE = sim.x.MODEL
#            mhe.y.VALUE = sim.y.MODEL
#            mhe.y_vel.VALUE = sim.y_vel.MODEL
#            mhe.x_vel.VALUE = sim.x_vel.MODEL
#            mhe.heading.VALUE = sim.heading.MODEL
#            mhe.ang_vel.VALUE = sim.ang_vel.MODEL
#            mhe.SPS_mass.VALUE = sim.SPS_mass.VALUE
#            mhe.RCS_mass.VALUE = sim.RCS_mass.VALUE
        #else:
            #mhe = gekko_model(measurements[-1], param, physics, target)
            #mhe.mhe_setup(mhe_time)
            #for value_list in measurements:
        
        start_point = np.roll(start_point,-1,axis=0)
        start_point[-1,:] = model_r
            
        
        sol = odeint(motion, r, [0,spacing], args=(PARAM.CSM_mass,PARAM.CSM_I,[thrust[0]*noise_total[0], thrust[1]*noise_total[1]]))
        r = sol[-1,:]
        
        match_data = np.roll(match_data,-1,axis=0)
        match_data[-1,:] = [r[0], r[1], r[5]]
        
        commands = np.roll(commands,-1,axis=0)
        commands[-1,:] = thrust
        
        model = odeint(motion, model_r, [0,spacing], args=(PARAM.CSM_mass,PARAM.CSM_I,[thrust[0]*bias_guess[0], thrust[1]*bias_guess[1]]))
        model = odeint(motion, r_all[time,:], [0,spacing], args=(PARAM.CSM_mass,PARAM.CSM_I,[thrust[0]*bias_guess[0], thrust[1]*bias_guess[1]]))
        model_r = model[-1,:]
        
        
        if (time >= mhe_time[-1]):# and (np.mod(time,5) == 0):
            bias_guess = estimator(match_data, start_point[0], commands, bias_guess, spacing)
            model_temp = start_point[0]
            for i in np.arange(len(start_point)):
                updated_model = odeint(motion, model_temp, [0,spacing], args=(PARAM.CSM_mass,PARAM.CSM_I,[commands[i,0]*bias_guess[0], commands[i,1]*bias_guess[1]]))
                model_temp = updated_model[-1,:]
                if i < len(start_point)-1:
                    start_point[i+1,:] = model_temp
            model_r = model_temp
            print(bias_guess)
            
        bias_estimate_all[time,:] = bias_guess
        noise_all[time,:] = noise_total
#        
        #sol = odeint(motion, target, [0,spacing], args=(PARAM.CSM_mass,0,[0,0], bias))
        #target = sol[-1,:]
        
        
        
        #model.mpc_update(r, target)
        if np.mod(time, 10) == 0:
            print("{:.3f}%".format(time*100/len(time_simulation)))
            

    
    fig, ax = plt.subplots()
    #Plot Orbit
    ax.plot(r_all[:,0],r_all[:,1],'ro',label="Chaser")
    ax.plot(model_all[:,0], model_all[:,1],'bx',label="Model")
    #Plot Earth
    earth = mpl.patches.Circle((0,0),6371000)
    ax.add_artist(earth)
    ax.legend()
    ax.scatter((0),[0],c='b')
    #ax.set_aspect('equal', 'box')
    
    fig, ax = plt.subplots(nrows=3, ncols=3)
    
    ax[0,0].plot(full_time, r_all[:,0],'r',label="Truth", alpha=.5)
    ax[0,0].plot(full_time, model_all[:,0],'b:',label="Estimate")
    ax[0,0].set_title('X Position')
    
    ax[0,1].plot(full_time, r_all[:,1],'r',label="Truth", alpha=.5)
    ax[0,1].plot(full_time, model_all[:,1],'b:',label="Estimate")
    ax[0,1].set_title('Y Position')
    
    ax[0,2].plot(full_time, r_all[:,2],'r',label="Truth", alpha=.5)
    ax[0,2].plot(full_time, model_all[:,2],'b:',label="Estimate")
    ax[0,2].set_title('X Velocity')
    
    ax[1,0].plot(full_time, r_all[:,3],'r',label="Truth", alpha=.5)
    ax[1,0].plot(full_time, model_all[:,3],'b:',label="Estimate")
    ax[1,0].set_title('Y Veolcity')
    
    ax[1,1].plot(full_time, r_all[:,4],'r',label="Truth", alpha=.5)
    ax[1,1].plot(full_time, model_all[:,4],'b:',label="Estimate")
    ax[1,1].set_title('Heading')
    
    ax[1,2].plot(full_time, r_all[:,5],'r',label="Truth", alpha=.5)
    ax[1,2].plot(full_time, model_all[:,5],'b:',label="Estimate")
    ax[1,2].set_title('Angular Velocity')
    
    ax[2,0].plot(full_time, r_all[:,6],'r',label="Truth", alpha=.5)
    ax[2,0].plot(full_time, model_all[:,6],'b:',label="Estimate")
    ax[2,0].set_title('Main Fuel')
    
    ax[2,1].plot(full_time, r_all[:,7],'r',label="Truth", alpha=.5)
    ax[2,1].plot(full_time, model_all[:,7],'b:',label="Estimate")
    ax[2,1].set_title('Control Fuel')
    
    ax[2,2].plot([np.NaN],[np.NaN],'r',label='Truth', alpha=.5)
    ax[2,2].plot([np.NaN],[np.NaN],'b:',label='Estimate')
    ax[2,2].get_yaxis().set_visible(False)
    ax[2,2].get_xaxis().set_visible(False)
    ax[2,2].legend()
    
    fig.tight_layout()
    
    fig, ax = plt.subplots(nrows=1, ncols=2)
    
    ax[0].plot(full_time, r_all[:,6], 'k', label='Actual',linewidth=2.0, alpha=.5)
    ax[0].plot(full_time, model_all[:,6],'r:', label='Estimate',linewidth=2.0)
    ax[0].set_ylabel("Fuel (kg)")
    ax[0].set_xlabel("Time (s)")
    ax[0].set_title("Main Thrusters")
    ax[0].legend()
    
    
    ax[1].plot(full_time, r_all[:,7], 'k', label='Acutal',linewidth=2.0, alpha=.5)
    ax[1].plot(full_time, model_all[:,7],'r:',label='Estimate',linewidth=2.0)
    ax[1].set_ylabel("Fuel (kg)")
    ax[1].set_xlabel("Time (s)")
    ax[1].set_title("Control Thrusters")
    ax[1].legend()
    
    fig.tight_layout()
    
    fig2, ax2 = plt.subplots(1)
    ax2.plot(full_time, bias_estimate_all[:,1],'r', label="RCS Estimate")
    ax2.plot(full_time, bias_estimate_all[:,0],'b', label="SPS Estimate")
    ax2.plot(full_time, noise_all[:,1],'r:', label="RCS Input")
    ax2.plot(full_time, noise_all[:,0],'b:', label="SPS Input")
    ax2.plot([0, full_time[-1]],[bias[1], bias[1]], 'r-', label="RCS Nominal", alpha=.5)
    ax2.plot([0, full_time[-1]],[bias[0], bias[0]], 'b-', label="SPS Nominal", alpha=.5)
    ax2.legend()
    ax2.set_ylabel('Input Noise')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylim([.7,1.3])