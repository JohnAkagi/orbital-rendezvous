# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 15:16:31 2019

@author: johna
"""

from gekko import GEKKO
import numpy as np

class gekko_model():
    
    def __init__(self, x0, param0, physics0, x_target):
        
        x_pos = x0[0]
        y_pos = x0[1]
        x_vel = x0[2]
        y_vel = x0[3]
        heading = x0[4]
        ang_vel = x0[5]
        SPS_mass = x0[6]
        RCS_mass = x0[7]
        
        dry_mass = param0[0]
        SPS_thrust = param0[1]
        RCS_thrust = param0[2]
        SPS_mass_flow = param0[3]
        RCS_mass_flow = param0[4]
        J = param0[5]
        r = param0[6]
        
        x_pos_t = x_target[0]
        y_pos_t = x_target[1]
        x_vel_t = x_target[2]
        y_vel_t = x_target[3]
        
        mu = physics0[0]
        
        self.m = GEKKO(remote=False)
        
        self.x_t = self.m.Var(value=x_pos_t)
        self.y_t = self.m.Var(value=y_pos_t)
        self.x_vel_t = self.m.Var(value=x_vel_t)
        self.y_vel_t = self.m.Var(value=y_vel_t)
        
        self.dry_mass = self.m.Param(value=dry_mass)
        self.mu = self.m.Param(value=mu)
        self.SPS_flow = self.m.FV(value=SPS_mass_flow, lb=0, ub=SPS_mass)
        self.RCS_flow = self.m.FV(value=RCS_mass_flow, lb=0, ub=RCS_mass)
        self.r = self.m.Param(value=r)
        
        self.RCS_bias = self.m.FV(value=1.0,ub=1.5,lb=.5)
        self.SPS_bias = self.m.FV(value=1.0,ub=1.5,lb=.5)
        
        self.RCS_com = self.m.MV(value=0,lb=-1,ub=1, integer=True)
        self.SPS_com = self.m.MV(value=0,lb=0,ub=1, integer=True)
        
        self.x = self.m.CV(value=x_pos, name='x')
        self.y = self.m.CV(value=y_pos, name='y')
        self.x_vel = self.m.CV(value=x_vel, name='x_vel')
        self.y_vel = self.m.CV(value=y_vel, name='y_vel')
        self.heading = self.m.CV(value=heading, name='heading')
        self.ang_vel = self.m.CV(value=ang_vel, name='ang_vel')
        self.SPS_mass = self.m.Var(value=SPS_mass, lb=0, ub=SPS_mass, name='SPS_mass')
        self.RCS_mass = self.m.Var(value=RCS_mass, lb=0, ub=RCS_mass, name='RCS_mass')
        
        self.J = self.m.FV(value=J)
        self.SPS_thrust = self.m.FV(value=SPS_thrust)
        self.RCS_thrust = self.m.FV(value=RCS_thrust)
        
        self.SPS_thrust_f = self.m.Intermediate(self.SPS_com*self.SPS_thrust*self.SPS_bias)
        self.RCS_thrust_f = self.m.Intermediate(self.RCS_com*self.RCS_thrust*self.RCS_bias)
        
        self.d_t = self.m.Intermediate(self.m.sqrt(self.x_t**2 + self.y_t**2))
        self.r_dd_t = self.m.Intermediate(-self.mu/self.d_t**3)
        
        self.mass_t = self.m.Intermediate(self.dry_mass + self.RCS_mass + self.SPS_mass)
        
        self.d = self.m.Intermediate(self.m.sqrt(self.x**2 + self.y**2))
        self.r_dd = self.m.Intermediate(-self.mu/(self.d**3))
        
        self.x_accel_g = self.m.Intermediate(self.x*self.r_dd)
        self.y_accel_g = self.m.Intermediate(self.y*self.r_dd)
        
        self.x_accel_SPS = self.m.Intermediate(self.SPS_thrust_f*self.m.sin(self.heading)/self.mass_t)
        self.y_accel_SPS = self.m.Intermediate(self.SPS_thrust_f*self.m.cos(self.heading)/self.mass_t)
        
        self.phase = self.m.Intermediate(self.m.atan(self.x/self.y))
        self.phase_t = self.m.Intermediate(self.m.atan(self.x_t/self.y_t))
        
        self.m.Equation(self.x.dt() == self.x_vel)
        self.m.Equation(self.y.dt() == self.y_vel)
        self.m.Equation(self.x_vel.dt() == self.x_accel_g + self.x_accel_SPS)
        self.m.Equation(self.y_vel.dt() == self.y_accel_g + self.y_accel_SPS)
        self.m.Equation(self.heading.dt() == self.ang_vel)
        self.m.Equation(self.ang_vel.dt() == self.RCS_thrust_f * self.r / self.J)
        self.m.Equation(self.SPS_mass.dt() == -self.SPS_flow * self.SPS_com * self.SPS_bias)
        self.m.Equation(self.RCS_mass.dt() == -self.RCS_flow * self.m.abs(self.RCS_com) * self.RCS_bias)
        
        self.m.Equation(self.x_t.dt() == self.x_vel_t)
        self.m.Equation(self.y_t.dt() == self.y_vel_t)
        self.m.Equation(self.x_vel_t.dt() == self.x_t * self.r_dd_t)
        self.m.Equation(self.y_vel_t.dt() == self.y_t * self.r_dd_t)
        
        
        
    def mpc_setup(self, time):
        self.m.time = time
        
        p = np.zeros(len(time))
        p[-2:] = 1
        p[:] = 1
        
        self.final = self.m.Param(value=p)
        
        self.RCS_com.STATUS = 1
        self.SPS_com.STATUS = 1
        
        
        self.m.options.MAX_ITER = 1000
        
        self.m.Equation(self.d >= 6471000) #Distance must be larger than the radius of the Earth + 100 km
        self.m.Equation(self.RCS_mass >= 0) #Must have positive mass in tanks
        self.m.Equation(self.SPS_mass >= 0)
        self.m.Equation(self.m.abs(self.ang_vel) <= 2*3.14)
        
        self.m.Obj((self.final*(self.d - self.d_t)**2))
        self.m.Obj(1e12*self.final*(self.phase - self.phase_t)**2)
        self.m.Obj((self.m.sqrt(self.x_vel**2 + self.y_vel**2) - self.m.sqrt(self.x_vel_t**2 + self.y_vel_t**2))**2)
        
        self.m.options.IMODE = 6
        
        
    def mpc_update(self, x0, x_target):
        x_pos = x0[0]
        y_pos = x0[1]
        x_vel = x0[2]
        y_vel = x0[3]
        heading = x0[4]
        ang_vel = x0[5]
        SPS_mass = x0[6]
        RCS_mass = x0[7]
        
        x_pos_t = x_target[0]
        y_pos_t = x_target[1]
        x_vel_t = x_target[2]
        y_vel_t = x_target[3]
        
        self.x.VALUE = x_pos
        self.y.VALUE = y_pos
        self.x_vel.VALUE = x_vel
        self.y_vel.VALUE = y_vel
        self.heading.VALUE = heading
        self.ang_vel.VALUE = ang_vel
        self.SPS_mass.VALUE = SPS_mass
        self.RCS_mass.VALUE = RCS_mass
        
        self.x_t.VALUE = x_pos_t
        self.y_t.VALUE = y_pos_t
        self.x_vel_t.VALUE = x_vel_t
        self.y_vel_t.VALUE = y_vel_t
        
    def mpc_solve(self):
        try:
            self.m.solve(disp=False)
            print("Finished")
            RCS_com = self.RCS_com.NEWVAL
            SPS_com = self.SPS_com.NEWVAL
        except:
            print("Did not finish")
            RCS_com = 0
            SPS_com = 0
                
        return([SPS_com, RCS_com])
        
    def mhe_setup(self, time):
        self.m.time = time
        self.m.options.EV_TYPE = 1
        self.m.options.IMODE = 5
        self.m.options.ICD_CALC = 1
        
    def mhe_on(self):
        self.x.FSTATUS = 1 #Receive measurements
        self.y.FSTATUS = 1
        self.ang_vel.FSTATUS = 1
        self.RCS_com.FSTATUS = 1
        self.SPS_com.FSTATUS = 1
        self.SPS_bias.FSTATUS = 0
        self.RCS_bias.FSTATUS = 0
        
        
        self.x.STATUS = 1 
        self.y.STATUS = 1
        self.ang_vel.STATUS = 1
        self.RCS_com.STATUS = 0
        self.SPS_com.STATUS = 0
        self.SPS_bias.STATUS = 1
        self.RCS_bias.STATUS = 1
        
        self.SPS_bias.DMAX = .05
        self.RCS_bias.DMAX = .05
        
    def mhe_add(self, x, y, ang_vel, SPS_com, RCS_com):
        self.x.MEAS = x
        self.y.MEAS = y
        self.ang_vel.MEAS = ang_vel
        self.SPS_com.MEAS = SPS_com
        self.RCS_com.MEAS = RCS_com
        
        
    def mhe_solve(self, x, y, ang_vel, SPS_com, RCS_com):
        
        self.x.MEAS = x
        self.y.MEAS = y
        self.ang_vel.MEAS = ang_vel
        self.SPS_com.MEAS = SPS_com
        self.RCS_com.MEAS = RCS_com
        
        self.m.solve(disp=False)
        
        return([self.SPS_bias.NEWVAL, self.RCS_bias.NEWVAL], [self.x.MODEL, self.y.MODEL, self.ang_vel.MODEL])
        
    def sim_setup(self,time):
        self.m.time = time
        
        self.RCS_com.FSTATUS = 1
        self.SPS_com.FSTATUS = 1
        
        self.m.options.SOLVER = 3
        self.m.options.IMODE = 4
        
    def sim_solve(self, SPS_com, RCS_com):
        self.RCS_com.MEAS = RCS_com
        self.SPS_com.MEAS = SPS_com
        self.SPS_com.STATUS = 0
        self.RCS_com.STATUS = 0
        self.m.solve(disp=False)
        
        return([self.x.MODEL, self.y.MODEL, self.ang_vel.MODEL])