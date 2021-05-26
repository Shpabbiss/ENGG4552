#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 17:55:49 2021

@author: matt
"""


import scipy.integrate as sp
import numpy as np
import matplotlib.pyplot as plt

"""Initial Conditions"""
init_alt = 300E3            #[m]
init_disp = 0               #[m]
init_v = 8000               #[m/s]
init_fpa = np.radians(1)    #[radians]

"""Constants"""
G = 6.67430E-11             #[Nm^2/kg^2]
Me = 5.97219E24             #[kg]
Re = 6371E3                 #[m]
y0 = 6620                   #[m] - Scale Height - AERO4800 Textbook
beta = 1/y0
surface_rho = 1.225         #[kg/m^3]
pica_rho = 280              #[kg/m^3]  (ASSUMED PICA DENSITY - WILL NEED CHANGING)

"""Vehicle Properties"""
mass = 10                   #[kg] Total vehicle mass            (ASSUMED - WILL NEED CHANGING)
noser = 0.3                 #[m] Nose radius                    (ASSUMED - WILL NEED CHANGING)
Cd = 1.4                      #Coefficient of drag                (ASSUMED - WILL NEED CHANGING)
Cl = 0                      #Coefficient of lift                (ASSUMED - WILL NEED CHANGING)
radius = 0.3
S = 0.05106375          #[m^2] - Reference Area
BC = mass/(S*Cd)            #Ballistic coefficient
em = 0.9                    #Emissivity                         (ASSUMED - WILL NEED CHANGING)
sig = 5.67*(10**-8)         #[W / (m^2 k^4)] Stefan-Boltzmann constant

            #Cp's:      Carbon = 710 [J/kgK]            Aluminium = 900 [J/kgK]   
            
cpc = 710  # 900            #[J/kgK] Specific heat of carbon    (ASSUMED - WILL NEED CHANGING)
thicc = 0.01                #[m] Heat shield thickness          (ASSUMED - WILL NEED CHANGING) (Potentially do this as a function later on as likely will be thicker at nose)



"""Loop Properties"""
dt = 0.05                   #[s]
time = 2000                 #[s]
steps = time/dt + 1
T_initial = 300             #K


###############################################################################
###############################################################################


"""Area Calculations"""

#Area Slicing Parameters

ang = 90 # Angle to in degrees
scale = 2 # Steps per degree
step = ((ang * scale) + 1)
angles = np.linspace(0, ang, step)

# Area iterative generator - Calculates sliced area based on angle and radius
def area_slice(r, t0, t1):
    A = 2 * np.pi * r * r * (1 - np.cos(np.radians(t1)) - (1 - np.cos(np.radians(t0))))
    return(A)

#Area Calculation Sanity Check

ar_it = 0
i = 0
areas = []

for i in range(step - 1):
    ar_it += area_slice(noser, angles[i], angles[i+1])
    areas.append(area_slice(noser, angles[i], angles[i+1]))
    
ar_ac = 2 * np.pi * noser * noser * (1 - np.cos(np.radians(ang)))

print("\nSanity Check: Analytical Area = Iterative Area (Exactly)")
print("Analytical Area =", round(ar_ac, 4), "m^2")
print("Iterative Area =", round(ar_it, 4), "m^2")
print()

###############################################################################
###############################################################################


def g_acc(alt): 
    
    """Calculates acceleration due to gravity at a certain altitude.
    params
    input
    alt [float] Altitude in m
    output
    acceleration in m/s^2
    """
    g = (G*Me)/((Re + alt)**2)
    return g

def density(alt):
    
    """Calculates density at a certain altitude.
    params
    input
    alt [float] Altitude in m
    output
    density in kg/m^3
    """
    rho = surface_rho*np.exp(-beta*alt)
    return rho

def drag(v,rho):
    
    """Calculates drag force at certain flow conditions.
    params
    input
    v [float] velocity in m/s
    rho [float] density in kg/m^3
    output
    drag force in N
    """
    
    Fd = 0.5*Cd*(v**2)*rho*S
    return Fd

def lift(v,rho):
    
    """Calculates lift force at certain flow conditions.
    params
    input
    v [float] velocity in m/s
    rho [float] density in kg/m^3
    output
    lfit force in N
    """
    
    Fl = 0.5*Cl*(v**2)*rho*S
    return Fl

def TrajectorySolver(gamma):
    
    """Calculates Trajectory of Re-Entry Vehicle based on vehicle paramaters.
    params
    input
    gamma - initial flight path angle in degrees
    variables at top of code
    output
    altitude array in m
    velocity array in m/s
    time at which parachute deployed in s
    """
    
    init_fpa = np.radians(gamma)
    t_vals = np.linspace(0,time,int(steps))
    
    #sets up required values for loop with initial parameters
    vx = vx_prev = init_v*np.cos(init_fpa)
    vy = vy_prev = init_v*np.sin(init_fpa)
    fpa = init_fpa
    alt = alt_prev = init_alt
    disp = disp_prev = init_disp
    
    #initialises arrays
    v_vals = np.zeros(int(steps))
    alt_vals = np.zeros(int(steps))
    disp_vals = np.zeros(int(steps))
    i = 0
    
    #for loop calculates trajectory
    for i,t in enumerate(t_vals):

        v = np.sqrt((vx**2)+(vy**2)) #Gets total velocity vector
        alt_vals[i] = alt            #Adds altitude val to array
        disp_vals[i] = disp
        v_vals[i] = v                #Adds velocity val to array
        rho = density(alt)           #Calculates density at this step
        g = g_acc(alt)               #Calculates grav acc at this step
        liftval = lift(v,rho)
        dragval = drag(v,rho)        #Calculates Drag at this step
        Fx = liftval*np.sin(fpa) - dragval*np.cos(fpa)#X Direction Force Calcs
        Fy = mass*g - liftval*np.cos(fpa) \
            - dragval*np.sin(fpa)   #Y Direction Force Calcs
        
        vx = vx_prev + ((Fx/mass) * dt) #Calcs v in x direction
        if vx <= 0:     
            vx = 0
        vy = vy_prev + ((Fy/mass) * dt) #Calcs v in y direction
        
        if vx == 0:         #Calcs new flight path angle based on vx and vy
            fpa = np.radians(90)
        else:
            fpa = np.arctan(vy/vx)
        
        alt = alt_prev - (vy*dt + 0.5*(Fy/mass)*(dt**2)) #Calcs new altitude
        disp = disp_prev + (vx*dt + 0.5*(Fx/mass)*(dt**2)) #Calcs new x disp
        vx_prev = vx   #Sets vx val for next loop
        vy_prev = vy   #Sets vy val for next loop
        alt_prev = alt #Sets alt val for next loop
        disp_prev = disp
        tlim = t
        if alt <= 2000: #Breaks loop at parachute deployment
            print("Parachute deployed at",round(alt,2), "m and Velocity = ",\
                  round(v,2),"m/s after", t, "seconds of flight time.")
            break 
    return alt_vals,disp_vals,v_vals,tlim

    
def a_val(v_vals):
    
    """Produces acceleration array based on velocity array.
    params
    input
    v_vals array [floats] velocity in m/s
    output
    Acceleration array in g's'
    """
    
    a = np.zeros(int(steps))
    for i in range(int(steps)):
        if i == 0:
            a[i] = 0
        else:
            a[i] = v_vals[i-1] - v_vals[i]
    a_vals = a/(9.81*dt)
    return a_vals

def array_cleaner(gamma):
    
    """Cleans extra zero values from TrajectorySolver arrays.
    params
    input
    gamma - initial flight path angle in degrees
    variables at top of code
    output
    TrajectorySolver and Acc arrays pruned of excess end values.
    """
    
    t_vals = np.linspace(0,time,int(steps))
    alt_vals,disp_vals,v_vals,tlim = TrajectorySolver(gamma)   
    a_vals = a_val(v_vals)
    i = 0
    check = 0
    

    for i in range(int(steps)):
        
        #If the parachute hasn't deployed, do nothing and go to next timestep
        if i <= (tlim/dt): 
            check += 1    
            continue
        
        #If parachute has deployed, delete unnecessary values
        else:
            alt_vals = np.delete(alt_vals,check)
            disp_vals = np.delete(disp_vals,check)
            v_vals = np.delete(v_vals,check)
            a_vals = np.delete(a_vals,check)
            t_vals = np.delete(t_vals,check)
        
    rho_vals = density(alt_vals) #Creates density array 
        
    return alt_vals,disp_vals,rho_vals,v_vals,a_vals,t_vals,gamma

"""Running the Code"""
alt_vals,disp_vals,rho_vals,v_vals,a_vals,t_vals,gamma = array_cleaner(0)
alt_vals1,disp_vals1,rho_vals1,v_vals1,a_vals1,t_vals1,gamma1=array_cleaner(5)
alt_vals2,disp_vals2,rho_vals2,v_vals2,a_vals2,t_vals2,gamma2=array_cleaner(10)
alt_vals3,disp_vals3,rho_vals3,v_vals3,a_vals3,t_vals3,gamma3=array_cleaner(15)



###############################################################################
###############################################################################
###############################################################################
###############################################################################


"""Heat Flux Setup and iterative solution"""

#Setup
#Velocity_List = velocity(G_alt)
#Rho_List = atmo_density(G_alt)
  
Heat_flux_master_list = [0, 0, 0, 0]
Heat_load_master_list = [0, 0, 0, 0]
q_conv_master_list = [0, 0, 0, 0]
q_rerad_master_list = [0, 0, 0, 0]
q_net_master_list = [0, 0, 0, 0]
Ts_master_list = [0, 0, 0, 0]
Heat_Load_in_master_list = [0, 0, 0, 0]
velocity_master_list = [v_vals, v_vals1, v_vals2, v_vals3]
rho_master_list = [rho_vals, rho_vals1, rho_vals2, rho_vals3]
t_vals_master_list = [t_vals, t_vals1, t_vals2, t_vals3]

# [columns, rows]

# Temperature array
Ts = np.ones((len(v_vals), step))
Ts[:,0] = 300





#Itertive Numerical Process
"""The for-loop cycles first through 'p' lists of
velocity and density values, each list being a 
different flight path angle (FPA). For each FPA, 
the loop then cycles through time steps, and for 
each time step the area slicing method is conducted
and summed."""
for p in range(len(velocity_master_list)):
    Heat_flux_list = []
    Heat_load_list = []
    q_conv_list = []
    q_rerad_list = []
    q_net_list = []
    Ts = np.ones((len(v_vals), step))
    Ts[:,0] = 300
    for i in range(len(velocity_master_list[p])):
        Heat_flux_slices = 0
        Heat_load_slices = 0
        q_conv_stag = 7.455*(10**(-9))*(rho_master_list[p][i]**0.4705)*(velocity_master_list[p][i]**(3.089))*(noser**(-0.52))*10000
        q_conv_prev = q_conv_stag.copy()
        q_conv_list.append(q_conv_stag/(10**6))
        for j in range(step - 1):
            area = areas[j]
            q_conv_ang = q_conv_stag*(np.cos(np.radians(angles[j+1])))
            q_conv_slice = (q_conv_ang + q_conv_prev)/2
            if i == 0:
                q_rerad = sig*em*((300)**4)
            if i != 0:
                Temp = (Ts[i-1, j] + Ts[i-1,j+1])/2
                q_rerad = sig*em*(Temp**4)
            if j == 0:    
                q_rerad_list.append(-q_rerad*(10**(-6)))
                q_net_list.append((q_conv_slice-q_rerad)*10**(-6))
            Heat_flux_conv = q_conv_slice * area
            Heat_Load_in_master_list[p] += Heat_flux_conv * dt
            Heat_flux_rerad = q_rerad * area
            Heat_flux_slices += area * (q_conv_slice - q_rerad)
            Heat_load_slices += area * (q_conv_slice - q_rerad) * dt
            q_conv_prev = q_conv_ang.copy()
            Ts[i,j] = Ts[i-1,j] + (((Heat_flux_conv - Heat_flux_rerad) * dt) / (area * pica_rho * thicc * cpc))
        Heat_flux_list.append(Heat_flux_slices/(10**6))
        Heat_load_list.append(Heat_load_slices)
    Heat_flux_master_list[p] = Heat_flux_list
    Heat_load_master_list[p] = Heat_load_list
    q_conv_master_list[p] = q_conv_list
    q_rerad_master_list[p] = q_rerad_list
    q_net_master_list[p] = q_net_list
    Ts_master_list[p] = Ts[:,0]




#Stagnation Point Convective Heat Flux - Time
plt.figure()
plt.plot(t_vals,q_conv_master_list[0], label="FPA = 0")
plt.plot(t_vals1,q_conv_master_list[1], label="FPA = 5")
plt.plot(t_vals2,q_conv_master_list[2], label="FPA = 10")
plt.plot(t_vals3,q_conv_master_list[3], label="FPA = 15")
plt.title("Stagnation Point Convective Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.legend()
plt.show()

#Stagnation Point Convective Heat Flux - Altitude
plt.figure()
plt.plot(q_conv_master_list[0], alt_vals,label="FPA = 0")
plt.plot(q_conv_master_list[1], alt_vals1,label="FPA = 5")
plt.plot(q_conv_master_list[2], alt_vals2,label="FPA = 10")
plt.plot(q_conv_master_list[3], alt_vals3,label="FPA = 15")
plt.title("Stagnation Point Convective Heat Flux Plotted as a Function of Altitude")
plt.xlabel('Heat Flux (MW/m^2)')
plt.ylabel('Altitude (m)')
plt.legend()
plt.show()




#Stagnation Reradiative Heat Flux - Time
plt.figure()
plt.plot(t_vals,q_rerad_master_list[0], label="FPA = 0")
plt.plot(t_vals1,q_rerad_master_list[1], label="FPA = 5")
plt.plot(t_vals2,q_rerad_master_list[2], label="FPA = 10")
plt.plot(t_vals3,q_rerad_master_list[3], label="FPA = 15")
plt.title("Stagnation Point Reradiative Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.legend()
plt.show()

#Stagnation Reradiative Heat Flux - Altitude
plt.figure()
plt.plot(q_rerad_master_list[0], alt_vals,label="FPA = 0")
plt.plot(q_rerad_master_list[1], alt_vals1,label="FPA = 5")
plt.plot(q_rerad_master_list[2], alt_vals2,label="FPA = 10")
plt.plot(q_rerad_master_list[3], alt_vals3,label="FPA = 15")
plt.title("Stagnation Point Reradiative Heat Flux Plotted as a Function of Altitude")
plt.xlabel('Heat Flux (MW/m^2)')
plt.ylabel('Altitude (m)')
plt.legend()
plt.show()




#Stagnation Point Total Plotting - Time 
plt.figure()
plt.plot(t_vals,q_net_master_list[0], label="FPA = 0")
plt.plot(t_vals1,q_net_master_list[1], label="FPA = 5")
plt.plot(t_vals2,q_net_master_list[2], label="FPA = 10")
plt.plot(t_vals3,q_net_master_list[3], label="FPA = 15")
plt.title("Stagnation Point Net Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.legend()
plt.show()

#Stagnation Point Total Plotting - Altitude 
plt.figure()
plt.plot(q_net_master_list[0], alt_vals,label="FPA = 0")
plt.plot(q_net_master_list[1], alt_vals1,label="FPA = 5")
plt.plot(q_net_master_list[2], alt_vals2,label="FPA = 10")
plt.plot(q_net_master_list[3], alt_vals3,label="FPA = 15")
plt.title("Stagnation Point Net Heat Flux Plotted as a Function of Altitude")
plt.xlabel('Heat Flux (MW/m^2)')
plt.ylabel('Altitude (m)')
plt.legend()
plt.show()




#Entire Surface Plotting - Time 
plt.figure()
plt.plot(t_vals,Heat_flux_master_list[0], label="FPA = 0")
plt.plot(t_vals1,Heat_flux_master_list[1], label="FPA = 5")
plt.plot(t_vals2,Heat_flux_master_list[2], label="FPA = 10")
plt.plot(t_vals3,Heat_flux_master_list[3], label="FPA = 15")
plt.title("Total Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW)')
plt.legend()
plt.show()

#Entire Surface Plotting - Altitude  
plt.figure()
plt.plot(Heat_flux_master_list[0], alt_vals,label="FPA = 0")
plt.plot(Heat_flux_master_list[1], alt_vals1,label="FPA = 5")
plt.plot(Heat_flux_master_list[2], alt_vals2,label="FPA = 10")
plt.plot(Heat_flux_master_list[3], alt_vals3,label="FPA = 15")
plt.title("Total Heat Flux Plotted as a Function of Altitude ")
plt.xlabel('Heat Flux (MW)')
plt.ylabel('Altitude  (m)')
plt.legend()
plt.show()



FPA_List = [0, 5, 10, 15]

for i in range(len(FPA_List)): 
    print("\n\n\nFPA =", FPA_List[i], ":")
    print("\nThe maximum heat flux over the enitre surface =", round(max(Heat_flux_master_list[i]), 4), "MW")
    print("\nThe maximum convective heat flux at the stagnation point =", round(max(q_conv_master_list[i]),4), "MW/m^2")
    print("\nThe maximum reradiative heat flux at the stagnation point =", round(min(q_rerad_master_list[i]),4), "MW/m^2")
    print("\nThe maximum net heat flux at the stagnation point =", round(max(q_net_master_list[i]), 4), "MW/m^2")
    print("\nThe maximum temperature =", round(np.max(Ts_master_list[i]), 4), "K")
    print("\nTotal heat load (SIMPS) =", round(sp.simps(t_vals_master_list[i] , Heat_flux_master_list[i]), 4), "MJ")
    # print("\nThe maximum dynamic pressure =", round(np.max(dp), 4), "kPa/m^2")
    print("\nThe total heat load =", round(Heat_Load_in_master_list[i], 4), "MJ")
    print("\nThe total heat load =", round(dt*sum(q_conv_master_list[i]), 4), "MJ/m^2")

