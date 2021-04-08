"""Initial Trajectory Bullshit"""
import numpy as np
import matplotlib.pyplot as plt

init_altitude = 400*10**3 #[m]
init_velocity = 8000    #[m/s]
surface_density = 1.225 #[kg/m^3]
gas_constant = 287.05   #[J/kg K]
abs_temp = 0 + 273.15   #[K] - Not 100% sure what this should be??
g = 9.81
beta = 1/6620       #Beta^-1 - retrieved from AERO4800 textbook
beta_op2 = 0.000139 #Beta^-1 - retrieved from Lit Review (Portfolio)
US_atmo_vals = [1.225,0.9093,0.7364,0.59,0.4135,0.1948,0.04008,0.001027,\
                0.00001846]  #US Standard Atmospheric Densities
US_atmo_heights = [0,3,5,7,10,15,25,50,80] #Corresponding altitudes


def scaleheight():   #Need to check this formula - conflicting information
    b = - (g/gas_constant)/(abs_temp)
    return(b)

def atmo_density(h):  #Calculates density based on altitude
    rho = surface_density*np.exp(scaleheight()*h)
    return rho

def atmo_density_plot():#Plots density vs altitude + compares with US Standard
    x_val = np.linspace(0,100*10**3,1000)
    y_val = []
    for x in x_val:
        y_val.append(atmo_density(x))
    plt.plot((x_val/1000),y_val)
    plt.plot(US_atmo_heights,US_atmo_vals)
    plt.title("Atmospheric Density vs Altitude")
    plt.xlabel("Altitude (km)")
    plt.ylabel("Density (kg/m^3)")
    
atmo_density_plot()
