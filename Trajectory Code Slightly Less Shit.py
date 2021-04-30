import numpy as np
import matplotlib.pyplot as plt
import ambiance

"""Initial Conditions"""
init_alt = 400E3 #[m]
init_v = 8000    #[m/s]
init_fpa = 1     #[degrees]

"""Constants"""
G = 6.67430E-11  #[Nm^2/kg^2]
Me = 5.97219E24  #[kg]
Re = 6371E3      #[m]
y0 = 6620        #[m] - Scale Height - AERO4800 Textbook
beta = 1/y0
surface_rho = 1.225 #[kg/m^3]

"""Vehicle Properties"""
mass = 10        #[kg]
noser = 1.0      #[m]
Cd = 1
Cl = 0
S = np.pi*noser**2 #[m^2] - Reference Area
BC = mass/(S*Cd)

"""Loop Properties"""
dt = 0.5 #[s]
time = 86400 #[s]
steps = time/dt + 1


def g_acc(alt): 
    #This function calculates the acceleration due to gravity at the
    #specified altitude
    g = (G*Me)/((Re + alt)**2)
    return g

def density(alt):
    #calculates density at a certain altitude
    rho = surface_rho*np.exp(-beta*alt)
    return rho

def drag(v,rho):
    #calculates drag based upon flow conditions
    Fd = 0.5*Cd*(v**2)*rho*S
    return Fd

def wackycalcs():
    vx = init_v*np.cos(init_fpa)
    vy = init_v*np.sin(init_fpa)
    vx_prev = init_v*np.cos(init_fpa)
    vy_prev = init_v*np.sin(init_fpa)
    fpa = init_fpa
    alt = init_alt
    for t in np.linspace(0,time,steps):

        rho = density(alt)
        g = g_acc(alt)
        Fx = -drag(vx,rho)
        Fy = mass*g - drag(vy,rho)
        vx = vx_prev + ((Fx/mass) * t)
        vy = vy_prev + ((Fy/mass) * t)
        fpa = np.arctan(vy/vx)
        vx_prev = vx
        vy_prev = vy
        


