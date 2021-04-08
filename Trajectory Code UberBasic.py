"""Initial Trajectory Bullshit"""

import numpy as np
import matplotlib.pyplot as plt


"""Initial Conditions"""
init_altitude = 400*10**3 #[m]
init_v = 8000    #[m/s]
gamma = np.radians(12)               #[degrees]


"""Constants"""
surface_density = 1.225 #[kg/m^3]
gas_constant = 287.05   #[J/kg K]
abs_temp = 0 + 273.15   #[K] - Not 100% sure what this should be??
g = 9.81
y0 = 6620           #[m] Scale height
beta = 1/6620       #Beta^-1 - retrieved from AERO4800 textbook
#beta_op2 = 0.000139 #Beta^-1 - retrieved from Lit Review (Portfolio)


"""US Atmospheric Standards"""
US_atmo_vals = [1.225,0.9093,0.7364,0.59,0.4135,0.1948,0.04008,0.001027,\
                0.00001846]  #US Standard Atmospheric Densities
US_atmo_heights = [0,3,5,7,10,15,25,50,80] #Corresponding altitudes


"""Vehicle Properties"""
mass = 10      #[kg]
Cd = 1          #Coefficient of Drag
Cl = 0          #Coefficient of Lift
S = 2.5           #[m^2] Reference Area - Based on CAD model
BC = mass/(S*Cd)

"""Functions/Calculations"""

#def scaleheight():   #Need to check this formula - conflicting information
 #   b = - (g/gas_constant)/(abs_temp)
  #  return(b)

def atmo_density(h):  #Calculates density based on altitude
    rho = surface_density*np.exp(beta*h)
    return rho

def atmo_density_plot():#Plots density vs altitude + compares with US Standard
    x_vals = np.linspace(0,100*10**3,1000)
    y_vals = []
    for x in x_vals:
        y_vals.append(atmo_density(x))
    plt.plot((x_vals/1000),y_vals,label = "Calculated")
    plt.plot(US_atmo_heights,US_atmo_vals,label = "US Standard")
    plt.title("Atmospheric Density vs Altitude")
    plt.legend()
    plt.xlabel("Altitude (km)")
    plt.ylabel("Density (kg/m^3)")
    plt.show()

def velocity(init_alt): #calculates the velocity at a given altitude
    alt_vals = np.linspace(0,init_alt,init_alt+1)
    v_vals = np.zeros(init_alt+1)
    i = 0
    for h in alt_vals:
        v_vals[i] = init_v*np.exp(-(((surface_density*y0)/(2*BC*np.sin(gamma)))*\
                             np.exp(-h/y0)))
        i += 1
    return v_vals

def deceleration(init_alt): #cqlculates acceleration based on altitude
    alt_vals = np.linspace(0,init_alt,init_alt+1)
    a_vals = np.zeros(init_alt+1)
    i = 0
    for h in alt_vals:
        a_vals[i] = ((surface_density*(init_v**2))/(2*BC*g)*np.exp(-h/y0))*np.exp(-\
        (((surface_density*y0)/(BC*np.sin(gamma)))*np.exp(-h/y0)))
        i += 1
    return a_vals

def vs_altitude_plot(v_vals,a_vals,init_alt): 
    #plots velocity and acceleration vs altitude
    x_vals = np.linspace(0,init_alt,init_alt+1)
    fig,ax = plt.subplots()
    ax.plot((x_vals/1000),v_vals,color="blue")
    ax.set_xlabel("Altitude (km)")
    ax.set_ylabel("Velocity (m/s)",color = "blue")
    ax2 = ax.twinx()
    ax2.plot((x_vals/1000),a_vals,color="red")
    ax2.set_ylabel("Deceleration (m/s^2)",color = "red")
    plt.title("Spacecraft Velocity vs Altitude")
    plt.show()
    
        
            


"""Running the Program"""
#atmo_density_plot()
v_vals = velocity(150000)
a_vals = deceleration(150000)
vs_altitude_plot(v_vals,a_vals,150000)




    
