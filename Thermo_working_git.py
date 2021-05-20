#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 17:55:49 2021

@author: matt
"""



import numpy as np
import matplotlib.pyplot as plt

"""Initial Conditions"""
init_alt = 200E3            #[m]
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
Cd = 1                      #Coefficient of drag                (ASSUMED - WILL NEED CHANGING)
Cl = 0                      #Coefficient of lift                (ASSUMED - WILL NEED CHANGING)
radius = 0.3
S = np.pi*radius**2          #[m^2] - Reference Area
BC = mass/(S*Cd)            #Ballistic coefficient
em = 0.9                    #Emissivity                         (ASSUMED - WILL NEED CHANGING)
sig = 5.67*(10**-8)         #[W / (m^2 k^4)] Stefan-Boltzmann constant
cpc = 710                   #[J/kgK] Specific heat of carbon    (ASSUMED - WILL NEED CHANGING)
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
    A1 = 2 * np.pi * r * r * (1 - np.cos(np.radians(t1)))
    A0 = 2 * np.pi * r * r * (1 - np.cos(np.radians(t0)))
    return (A1 - A0)

#Area Calculation Sanity Check

ar_it = 0
i = 0

for i in range(step - 1):
    ar_it += area_slice(noser, angles[i], angles[i+1])
    
ar_ac = 2 * np.pi * noser * noser * (1 - np.cos(np.radians(ang)))

print("\nSanity Check: Analytical Area = Iterative Area (Exactly)")
print("Analytical Area =", round(ar_ac, 4), "m^2")
print("Iterative Area =", round(ar_it, 4), "m^2")

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

def TrajectorySolver():
    
    """Calculates Trajectory of Re-Entry Vehicle based on vehicle paramaters.
    params
    input
    variables at top of code
    output
    altitude array in m
    velocity array in m/s
    time at which parachute deployed in s
    """
    
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
            print("\nParachute deployed at",round(alt,2), "m and Velocity = ",\
                  round(v,2),"m/s \nafter", t, "seconds of flight time.")
            break 
    return alt_vals,disp_vals,v_vals,tlim

def plotter(t_vals,alt_vals,disp_vals,rho_vals,v_vals,a_vals):
    
    """Plots Graphs
    params
    input
    alt_vals array [floats] Altitude in m
    v_vals array [floats] velocity in m/s
    tlim [float] Parachute deployment time in s
    a_vals array [floats] acceleration in g's'
    
    output
    Altitude, Velocity, Acceleration vs Time
    """
    
    #Plot Altitude vs Time
    plt.plot(t_vals, (alt_vals/1000))
    plt.title("Altitude vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Altitude (km)")
    #plt.xlim(0,tlim)
    plt.show()
    
    #Plot Altitude vs Displacement over Ground
    plt.plot(disp_vals/1000,alt_vals/1000)
    plt.title("Altitude vs Displacement over Ground")
    plt.xlabel("Ground Displacement (km)")
    plt.ylabel("Altitude (km)")
    plt.show()
    
    #Plot Density vs Time
    plt.plot(t_vals,rho_vals)
    plt.title("Air Density vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Air Density (kg/m^3)")
    plt.show()
    
    #Plot Velocity
    plt.plot(t_vals, (v_vals/1000))
    plt.title("Velocity vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (km/s)")
    #plt.xlim(0,tlim)
    plt.show()
    
    #Plot Decceleration
    plt.plot(t_vals,a_vals)
    plt.title("Deceleration vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Deceleration (g's)")
    #plt.xlim(0,tlim)
    plt.show()
    
    
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

def array_cleaner():
    
    """Cleans extra zero values from TrajectorySolver arrays.
    params
    input
    variables at top of code
    output
    TrajectorySolver and Acc arrays pruned of excess end values.
    """
    
    t_vals = np.linspace(0,time,int(steps))
    alt_vals,disp_vals,v_vals,tlim = TrajectorySolver()   
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
    
    rho_vals = density(alt_vals) #creates density array for cliff
        
    return alt_vals,disp_vals,rho_vals,v_vals,a_vals,t_vals

"""Running the Code"""
print(steps)
alt_vals,disp_vals,rho_vals,v_vals,a_vals,t_vals = array_cleaner()
plotter(t_vals,alt_vals,disp_vals,rho_vals,v_vals,a_vals)


###############################################################################
###############################################################################
###############################################################################
###############################################################################


"""Heat Flux Setup and iterative solution"""

#Setup
#Velocity_List = velocity(G_alt)
#Rho_List = atmo_density(G_alt)
q_conv_list = []
q_rerad_list = []
q_net_list = []  
Heat_flux_list = []
Heat_load_list = []
Drag_list = []
Lift_list = []
dp = []
# [columns, rows]

# Temperature array
Ts = np.ones((len(v_vals), step))
Ts[:,0] = 300


Heat_load_in = 0
#Itertive Numerical Process
for i in range(len(v_vals)):
    Heat_flux_slices = 0
    Heat_load_slices = 0
    q_conv_stag = 7.455*(10**(-9))*(rho_vals[i]**0.4705)*(v_vals[i]**(3.089))*(noser**(-0.52))*10000
    dynam_press = (0.5 * (rho_vals[i]) * ((v_vals[i]) ** 2))
    dp.append(dynam_press/1000)
    q_conv_prev = q_conv_stag.copy()
    q_conv_list.append(q_conv_stag/(10**6))
    dragg = 0.5 * rho_vals[i] * v_vals[i] * v_vals[i] * Cd * S
    liftt  =0.5 * rho_vals[i] * v_vals[i] * v_vals[i] * Cl * S
    Drag_list.append(dragg/1000)
    Lift_list.append(liftt)
    for j in range(step - 1):
        area = area_slice(noser, angles[j], angles[j+1])
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
        Heat_load_in += Heat_flux_conv * dt
        Heat_flux_rerad = q_rerad * area
        Heat_flux_slices += area * (q_conv_slice - q_rerad)
        Heat_load_slices += area * (q_conv_slice - q_rerad) * dt
        q_conv_prev = q_conv_ang.copy()
        Ts[i,j] = Ts[i-1,j] + (((Heat_flux_conv - Heat_flux_rerad) * dt) / (area * pica_rho * thicc * cpc))
    Heat_flux_list.append(Heat_flux_slices/1000)
    Heat_load_list.append(Heat_load_slices)


#F_drag = 0.5*rho*v*v*Cd*A

#Stagnation Point Plotting
plt.figure()
plt.plot(t_vals,q_conv_list)
plt.title("Stagnation Convective Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.show()

#Heat Flux
plt.figure()
plt.plot(t_vals,q_rerad_list)
plt.title("Stagnation Reradiative Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.show()

#Entire Surface Plotting
plt.figure()
plt.plot(t_vals, Heat_flux_list)
plt.title("Total Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (kW)')
plt.show()

#Nose Temperature
plt.figure()
plt.plot(t_vals, Ts[:,0], 'b', label='Nose Temperature')
plt.xlabel('Time (s)')
plt.ylabel('Nose Temperature (K)')
plt.title("Nose Temperature Plotted as a Function of Time")
plt.show()

#Drag
plt.figure()
plt.plot(t_vals, Drag_list, 'b', label='Drag')
plt.xlabel('Time (s)')
plt.ylabel('Drag (kN)')
plt.title("Drag Force on Capsule Plotted as a Function of Time")
plt.show()

#Lift
plt.figure()
plt.plot(t_vals, Lift_list, 'b', label='Lift')
plt.xlabel('Time (s)')
plt.ylabel('Lift (N)')
plt.title("Lift Force on Capsule Plotted as a Function of Time")
plt.show()

#Dynamic Pressure
plt.figure()
plt.plot(t_vals, dp, 'b', label='Dynamic Pressure')
plt.xlabel('Time (s)')
plt.ylabel('Dynamic Pressure (kPa/m^2)')
plt.title("Dynamic Pressure at Nose Plotted as a Function of Time")
plt.show()

#88 Degree Temperature
#plt.figure()
#plt.plot(t_vals, Ts[:,176], 'b', label='88 Degree Temperature')
#plt.xlabel('Time (s)')
#plt.ylabel('88 Degree Temperature (K)')
#plt.title("88 Degree Temperature Plotted as a Function of Time")
#plt.show()

#Maximum Temperature wrt to angle
plt.figure()
plt.plot(np.linspace(0, ang, step), Ts[3273, :], 'b', label='Max temperature')
plt.xlabel('Angle')
plt.ylabel('Nose Temperature (K)')
plt.title("Maximum Temperature Plotted as a Function of Angle")
plt.show()



print("\n\nThe maximum heat flux over the enitre surface =", round(max(Heat_flux_list)*10**-3, 4), "MW")
print("\nThe maximum convective heat flux at the stagnation point =", round(max(q_conv_list),4), "MW/m^2")
print("\nThe maximum reradiative heat flux at the stagnation point =", round(min(q_rerad_list),4), "MW/m^2")
print("\nThe maximum net heat flux at the stagnation point =", max(q_net_list), "MW/m^2")
print("\nThe maximum temperature =", round(np.max(Ts), 4), "K")
print("\nThe maximum dynamic pressure =", round(np.max(dp), 4), "kPa/m^2")
print("\nThe total heat load =", Heat_load_in*(10**-6), "MJ")
print("\nThe total heat load =", dt*sum(q_conv_list), "MJ/m^2")

