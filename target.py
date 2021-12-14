# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 11:34:55 2021

@author: abhil
"""

from math import sin, cos, radians
from scipy.optimize import minimize, differential_evolution
from drag import solve_ode_scipy
import numpy as np
import matplotlib.pyplot as plt

def objective(initial, t):
    
    timestep = 0.1
    
    h = 0
    v0 = initial[0]
    theta = initial[1]
    tx = t[0]
    ty = t[1]
    theta_r = radians(theta)    
    
    v0_bnds = (100, 800)
    theta_bnds = (10, 80)
    
    bnds = [v0_bnds, theta_bnds]

    # sol = differential_evolution(solve_ode_scipy, bounds = bnds, args = [v0, theta, h, timestep])
    
    # y0 = np.array([v0, theta])
    
    f0 = solve_ode_scipy(v0, theta, h, timestep)
    # f1 = drag_ode(t, [v0*cos(theta_r), v0*sin(theta_r), 0, 0])
                      
    sol = minimize(f0, initial, bounds = bnds)
    
    return sol
    
c0 = objective(initial = (100, 80), t = (1000,50))
