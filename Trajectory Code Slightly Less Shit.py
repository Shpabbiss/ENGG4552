import numpy as np
import matplotlib.pyplot as plt

"""Initial Conditions"""
init_alt = 400E3 #[m]
init_v = 8000    #[m/s]
init_fpa = 1     #[degrees]

"""Constants"""
G = 6.67430E-11  #[Nm^2/kg^2]
Me = 5.97219E24  #[kg]
Re = 6371E3      #[m]
y0 = 6620
beta = 1/y0

"""Vehicle Parameters"""
mass = 10        #[kg]
noser = 1.0      #[m]
Cd = 1
Cl = 0
S = np.pi*noser**2 #[m^2] - Reference Area
BC = mass/(S*Cd)

g = (G*Me)/(Re**2)

