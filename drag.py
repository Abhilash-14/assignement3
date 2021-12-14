# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 13:49:20 2021

@author: Abhilash Sharma
"""


#hello


from math import sqrt, radians, pi, sin, cos, tan 
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def drag_ode(t,y = []):
    
    gx, gy = 0, 9.81
    
    vx, vy, rx, ry = y[0], y[1], y[2], y[3]
    
    drxdt, drydt = vx, vy
    
    modV = sqrt((vx**2) + (vy**2))
    
    M = 5.5
    A = pi*((0.1/2)**2)
    cof = 0.479
    rho = 1.225
        
    dvxdt = ((-0.5*rho*cof*A*modV)*vx-M*gx)/M
    dvydt = ((-0.5*rho*cof*A*modV)*vy-M*gy)/M

    sol = [dvxdt, dvydt, drxdt, drydt]
    
    return sol

def solve_ode_euler(v0, theta, h, timestep, int_time = 100):
            
    theta_r = radians(theta)
    
    vx = v0*cos(theta_r)
    vy = v0*sin(theta_r)
    rx = 0
    ry = h
    t = 0
        
    i = 0
    n = float(int_time)/timestep 
    
    vx_L = [vx]
    vy_L = [vy]
    rx_L = [rx]
    ry_L = [ry]
    t_L = [t]
    
    while ry_L[-1] >= 0 and i < n:
        
        if i > 0:
            y_dot = drag_ode(t, [vx_L[i - 1], vy_L[i - 1], rx_L[i - 1], ry_L[i - 1]])
            
            new_vx = vx_L[i - 1] + timestep*y_dot[0]
            new_vy = vy_L[i - 1] + timestep*y_dot[1]
            new_rx = rx_L[i - 1] + timestep*y_dot[2]
            new_ry = ry_L[i - 1] + timestep*y_dot[3]
                
            vx_L.append(new_vx)
            vy_L.append(new_vy)
            rx_L.append(new_rx)
            ry_L.append(new_ry)
            t_L.append(t)
            
        t += timestep
        i += 1
        
    vx_arr = np.array(vx_L)
    vy_arr = np.array(vy_L)
    rx_arr = np.array(rx_L)
    ry_arr = np.array(ry_L)
    t_arr = np.array(t_L)
    
    a = np.column_stack((vx_arr, vy_arr))
    b = np.column_stack((rx_arr, ry_arr))
    c = np.column_stack((a,b))
    sol = np.column_stack((c,t_arr))
        
    return sol   

def solve_ode_scipy(v0, theta, h, timestep, int_time = 100):
    
    theta_r = radians(theta)
    
    vx = v0*cos(theta_r)
    vy = v0*sin(theta_r)
    rx = 0
    ry = h
    
    n = int(int_time)/timestep
    t = np.linspace(0, 100, int(n))
    
    def event(t, y):
        return y[3]
    
    event.terminal = True
    sol = solve_ivp(fun = drag_ode, t_span =[t[0], t[-1]], y0 = [vx, vy, rx, ry],\
                    t_eval = t, events = [event])   
    arr = sol.y    
   
    return arr
    
def no_drag(v0, theta, h):
    
    g = 9.81
    
    theta_r = radians(theta)
    
    d1 = v0*cos(theta_r)/g
    d2 = v0*sin(theta_r)
    d3 = sqrt(((v0*sin(theta_r))**2) + 2*g*h)
    
    d = (d1)*(d2 + d3)
        
    t = d/(v0*cos(theta_r))
    
    x = np.linspace(0, int(d), 10000)
    y = h + x*tan(theta_r) - (x**2)*(g/(2*((v0*cos(theta_r))**2)))       
   
    return x, y, t

if __name__ == '__main__':
    print("The initial conditions for the projectile will be needed to predict its motion.")
    print("It should be noted that the timestep is set by default to 0.1s.")
    
    while True:
        try:
            v0 = input("Enter an initial speed [m/s]: ")
            
            if float(v0) > 0 and v0.isdigit() == True:        
                print("Valid initial speed detected! \n")
                break
        
        except ValueError:
            print("Invalid initial speed. Try again!")
    
    print("The angle of attack is between 0 and 90 degrees.")
            
    while True:
        try:
            theta = input("Enter an initial angle of attack [deg]: ")
            
            if (float(theta) > 0 and float(theta) <= 90) and theta.isdigit() == True:
                print("Valid angle of attack detected!")
                break
        except ValueError:
            print("Invalid angle of attack detected. Try again!")
            
    while True:
        try:
            h = input("Enter an initial height [m]: ")
            
            if float(h) >= 0 and h.isdigit() == True:
                print("Valid height detected")
                break
        
        except ValueError:
            print("Invalid height deteced. Try again!")
    
    print()
    print("The default integration time has been set to 100s. This can be changed below.")
    print("However, if you wish to keep it at 100s, enter 'n' to keep it at the default setting.")
    
    while True:
        try:
            int_time = input("Enter a maximum integration time: ")
            
            if int_time == 'n':
                print("The integration time will be set to 100s.")
                int_time = 100
                break
            
            elif float(int_time) > 0 and int_time.isdigit() == True:
                print("Valid integration time detected!")
                break
                
        except ValueError:
            print("Invalid integration time detected. Try again!")
        
    timestep = 0.1
    sol1 = no_drag(float(v0), float(theta), float(h))
    sol2 = solve_ode_euler(float(v0), float(theta), float(h), timestep)
    
    plt.plot(sol1[0], sol1[1], 'blue', label = 'Without drag')
    plt.plot(sol2[:,2], sol2[:,3], 'orange', label = 'With drag')
    
    plt.xlabel("Downrange Distance [m]")
    plt.ylabel("Height [m]")
    plt.title("Projectile Motion")
    plt.legend()
    plt.xlim(0,)
    plt.ylim(0,)
     
    plt.show()
    
    L = [float(sol2[-1:,4]), float(sol2[-1:,2]), max(sol2[:,3])]
    
    print()
    print("Basic Statistics")
    print("Travel time [s]:", round(L[0], 4), end = "s \n")
    print("Distance travelled [m]:", round(L[1], 4), end = "m \n")
    print("Maximum height [m]:", round(L[2], 4), end = "m") 
