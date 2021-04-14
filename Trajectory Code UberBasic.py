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
r = 1.5
mass = 10      #[kg]
Cd = 1          #Coefficient of Drag
Cl = 0          #Coefficient of Lift
Cfavg = 0.02         #Body Averaged Skin Friction Coefficient (PLACEHOLDER 
# UNTIL CFD- value taken from 4800 powerpoint)         
S = np.pi*r**2           #[m^2] Reference Area - Based on CAD model
Sw = S          #[m^2] Whole Exposed Surface - Based on CAD model
BC = mass/(S*Cd) #Ballistic Coefficient 

"""Functions/Calculations"""

#def scaleheight():   #Need to check this formula - conflicting information
 #   b = - (g/gas_constant)/(abs_temp)
  #  return(b)

def atmo_density(init_alt):  #Calculates density based on altitude
    alt_vals = np.linspace(0,init_alt,init_alt+1)
    rho_vals = np.zeros(init_alt+1)
    i = 0
    for h in alt_vals:
        rho_vals[i] = surface_density*np.exp(-beta*h)
        i += 1
    return rho_vals

def atmo_density_plot(init_alt,rho):
    #Plots density vs altitude + compares with US Standard
    alt_vals = np.linspace(0,init_alt,init_alt+1) 
    plt.plot((alt_vals/1000),rho,label = "Calculated")
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
        v_vals[i] = init_v*np.exp(-(((surface_density*y0)/\
                                     (2*BC*np.sin(gamma)))*np.exp(-h/y0)))
        i += 1
    return v_vals

def deceleration(init_alt): #cqlculates acceleration based on altitude
    alt_vals = np.linspace(0,init_alt,init_alt+1)
    a_vals = np.zeros(init_alt+1)
    i = 0
    for h in alt_vals:
        a_vals[i] = ((surface_density*(init_v**2))/(2*BC*g)*np.exp(-h/y0))*\
            np.exp(-(((surface_density*y0)/(BC*np.sin(gamma)))*np.exp(-h/y0)))
        i += 1
    return a_vals

def heat_flux(init_alt,v_vals,rho_vals):
    #Calculating heat flux at each altitude
    alt_vals = np.linspace(0,init_alt,init_alt+1)
    q_vals1 = np.zeros(init_alt+1)
    q_vals2 = np.zeros(init_alt+1)
    i = 0
    for h in alt_vals:
        #q_vals 1 is the formula found in lit review, 2 is from AERO4800
        q_vals1[i] = 1.83*(10**-4)*(v_vals[i]**3)*np.sqrt(rho_vals[i]/r)
        q_vals2[i] = np.sqrt(surface_density)*(init_v**3)*np.exp(-h/(2*y0))*\
        np.exp(((-3*surface_density*y0)/(2*BC*np.sin(gamma)))*np.exp(-h/y0))
        i+=1
    return q_vals1, q_vals2

def heat_flux_plot(q1,q2,init_alt):
    #This Works its just the same line on the graph - wack
    x_vals = np.linspace(0,init_alt,init_alt+1)
    fig,ax = plt.subplots()
    ax.plot(x_vals/1000,q1,color='red',label='Lit Review Formula')
    ax.set_xlabel("Altitude (km)")
    ax.set_ylabel("Heat Flux (J/s) - Lit Review")
    ax2 = ax.twinx()
    ax2.plot(x_vals/1000,q2,color='blue',label='AERO4800 Formula')
    ax2.set_ylabel("Heat Flux (J/s) - AERO4800")
    plt.title("Heat Flux vs Altitude")
    #plt.legend()
    plt.show()
    
    #plt.plot(x_vals/1000,q1,color='red',label='Lit Review Formula')
    #plt.plot(x_vals/1000,q2,color='blue',label='AERO4800 Formula')
    #plt.legend()
    #plt.title("Heat Flux vs Altitude")
    #plt.xlabel("Altitude (km)")
    #plt.ylabel("Heat Flux q (J/s)")
    #plt.show()        

def vs_altitude_plot(v_vals,a_vals,init_alt): 
    #plots velocity and acceleration vs altitude
    x_vals = np.linspace(0,init_alt,init_alt+1)
    fig,ax = plt.subplots()
    ax.plot((x_vals/1000),v_vals,color="blue")
    ax.set_xlabel("Altitude (km)")
    ax.set_ylabel("Velocity (m/s)",color = "blue")
    ax2 = ax.twinx()
    ax2.plot((x_vals/1000),a_vals,color="red")
    ax2.set_ylabel("Deceleration (g's)",color = "red")
    plt.title("Spacecraft Velocity vs Altitude")
    plt.show()
    
def total_heat():
    #Calculates total heat load over trajectory Q
    Qtot = (mass/4)*((Cfavg*Sw)/(Cd*S))*(init_v)**2
    return Qtot        


"""Running the Program"""
Graphing_Altitude = 150000
rho_vals = atmo_density(Graphing_Altitude)
atmo_density_plot(Graphing_Altitude,rho_vals)
v_vals = velocity(Graphing_Altitude)
a_vals = deceleration(Graphing_Altitude)
q_vals1,q_vals2 = heat_flux(Graphing_Altitude,v_vals,rho_vals)
vs_altitude_plot(v_vals,a_vals,Graphing_Altitude)
heat_flux_plot(q_vals1,q_vals2,Graphing_Altitude)
Qtot = total_heat()
print("Total heat load absorbed =",Qtot*10**-6,"MJ/kg")


    
