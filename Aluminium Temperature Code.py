# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 12:00:27 2021

@author: cliff
"""

"""New Thermal Analysis Code - Thermal_Analysis_Final may be bugged"""

import scipy.integrate as sp
import numpy as np
import matplotlib.pyplot as plt

"""Initial Conditions"""
init_alt = 200E3            #[m]
init_disp = 0               #[m]
init_v = 7800               #[m/s]
init_fpa = np.radians(1)    #[radians]

"""Constants"""
G = 6.67430E-11             #[Nm^2/kg^2]
Me = 5.97219E24             #[kg]
Re = 6371E3                 #[m]
y0 = 6620                   #[m] - Scale Height - AERO4800 Textbook
beta = 1/y0
surface_rho = 1.225         #[kg/m^3]
pica_rho = 280              #[kg/m^3]  (ASSUMED PICA DENSITY - WILL NEED CHANGING)
al_rho = 2710

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
         
cpc = 1925.56  # 900            #[J/kgK] Specific heat of carbon    (ASSUMED - WILL NEED CHANGING)
#cpa = 1050            #[m] Heat shield thickness          (ASSUMED - WILL NEED CHANGING) (Potentially do this as a function later on as likely will be thicker at nose)
#thicc_a = 0.012
thicc = 0.034

"""Loop Properties"""
dt = 0.5                   #[s]
time = 360                 #[s]
steps = time/dt + 1


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
    ar_it += area_slice(radius, angles[i], angles[i+1])
    areas.append(area_slice(radius, angles[i], angles[i+1]))

ar_ac = 2 * np.pi * radius * radius * (1 - np.cos(np.radians(ang)))

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

def g_mod(v):
    if v >= 7800:
        g_mod = 0
    else:
        g_mod = (7800 - v)/7800
    return g_mod
    
        

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
        #print(v)
        #print(alt)
        alt_vals[i] = alt            #Adds altitude val to array
        disp_vals[i] = disp
        v_vals[i] = v                #Adds velocity val to array
        rho = density(alt)           #Calculates density at this step
        ga = g_acc(alt)               #Calculates grav acc at this step
        gm = g_mod(v)
        g = ga * gm
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

def plotter(t_vals,alt_vals,disp_vals,v_vals,a_vals):
    
    """Plots Graphs
    params
    input
    alt_vals array [floats] Altitude in m
    v_vals array [floats] velocity in m/s
    a_vals array [floats] acceleration in g's'
    output
    Altitude, X Displacement, Velocity, Decceleration vs Time
    """
    
    #Plot Altitude vs Time
    plt.plot(t_vals, alt_vals)
    plt.title("Altitude vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Altitude (m)")
    plt.show()
    
    #Plot Altitude vs Displacement over Ground
    plt.plot(disp_vals/1000,alt_vals/1000)
    plt.title("Altitude vs Displacement over Ground")
    plt.xlabel("Ground Displacement (km)")
    plt.ylabel("Altitude (km)")
    plt.show()
    
    #Plot Velocity
    plt.plot(t_vals, v_vals)
    plt.title("Velocity vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (m/s)")
    plt.show()
    
    #Plot Decceleration
    plt.plot(t_vals,a_vals)
    plt.title("Deceleration vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Deceleration (g's)")
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
        
    rho_vals = density(alt_vals) #creates density array for cliff
        
    return alt_vals,disp_vals,rho_vals,v_vals,a_vals,t_vals,gamma

def plot_comparisons(alt_vals,disp_vals,v_vals,a_vals,t_vals,gamma,alt_vals1,\
                     disp_vals1,v_vals1,a_vals1,t_vals1,gamma1,alt_vals2,\
                         disp_vals2,v_vals2,a_vals2,t_vals2,gamma2,alt_vals3,\
                             disp_vals3,v_vals3,a_vals3,t_vals3,gamma3):
    
    """Plots Graphs
    params
    input
    4 alt_vals array [floats] Altitude in m
    4 v_vals array [floats] velocity in m/s
    a_vals array [floats] acceleration in g's'
    
    output
    Altitude, X Displacement, Velocity, Decceleration vs Time
    """
    
    #Plot Altitude vs Time
    plt.plot(t_vals, alt_vals,label="Initial FPA = "+str(gamma))
    #plt.plot(t_vals1, alt_vals1,label="Initial FPA = " + str(gamma1))
    #plt.plot(t_vals2, alt_vals2,label="Initial FPA = " + str(gamma2))
    #plt.plot(t_vals3, alt_vals3,label="Initial FPA = " + str(gamma3))
    plt.title("Altitude vs Time")
    plt.legend()
    plt.xlabel("Time (s)")
    plt.ylabel("Altitude (m)")
    plt.show()
    
    #Plot Altitude vs Displacement over Ground
    plt.plot(disp_vals/1E3, alt_vals/1E3,label="Initial FPA = "+str(gamma))
    #plt.plot(disp_vals1/1E3, alt_vals1/1E3,label="Initial FPA = "+str(gamma1))
    #plt.plot(disp_vals2/1E3, alt_vals2/1E3,label="Initial FPA = "+str(gamma2))
    #plt.plot(disp_vals3/1E3, alt_vals3/1E3,label="Initial FPA = "+str(gamma3))
    plt.title("Altitude vs Displacement over Ground")
    plt.legend()
    plt.xlabel("Ground Displacement (km)")
    plt.ylabel("Altitude (km)")
    plt.show()
    
    #Plot Velocity
    plt.plot(t_vals, v_vals,label="Initial FPA = "+str(gamma))
    #plt.plot(t_vals1, v_vals1,label="Initial FPA = " + str(gamma1))
    #plt.plot(t_vals2, v_vals2,label="Initial FPA = " + str(gamma2))
    #plt.plot(t_vals3, v_vals3,label="Initial FPA = " + str(gamma3))
    plt.title("Velocity vs Time")
    plt.legend()
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity (m/s)")
    plt.show()
    
    #Plot Decceleration
    plt.plot(t_vals, a_vals,label="Initial FPA = "+str(gamma))
    #plt.plot(t_vals1, a_vals1,label="Initial FPA = " + str(gamma1))
    #plt.plot(t_vals2, a_vals2,label="Initial FPA = " + str(gamma2))
    #plt.plot(t_vals3, a_vals3,label="Initial FPA = " + str(gamma3))
    plt.title("Deceleration vs Time")
    plt.legend()
    plt.xlabel("Time (s)")
    plt.ylabel("Deceleration (g's)")
    plt.show()
    
    #Plot Decceleration vs Altitude
    plt.plot(alt_vals/1E3, a_vals,label="Initial FPA = "+str(gamma))
    #plt.plot(alt_vals1/1E3, a_vals1,label="Initial FPA = " + str(gamma1))
    #plt.plot(alt_vals2/1E3, a_vals2,label="Initial FPA = " + str(gamma2))
    #plt.plot(alt_vals3/1E3, a_vals3,label="Initial FPA = " + str(gamma3))
    plt.title("Decceleration vs Altitude")
    plt.legend()
    plt.xlabel("Altitude (km)")
    plt.ylabel("Deceleration (g's)")
    plt.show()


"""Running the Code"""
print(steps)
alt_vals,disp_vals,rho_vals,v_vals,a_vals,t_vals,gamma=array_cleaner(10)
#alt_vals1,disp_vals1,rho_vals1,v_vals1,a_vals1,t_vals1,gamma1=array_cleaner(5)
#alt_vals2,disp_vals2,rho_vals2,v_vals2,a_vals2,t_vals2,gamma2=array_cleaner(10)
#alt_vals3,disp_vals3,rho_vals3,v_vals3,a_vals3,t_vals3,gamma3=array_cleaner(15)
plotter(t_vals,alt_vals,disp_vals,v_vals,a_vals)
#plot_comparisons(alt_vals,disp_vals,v_vals,a_vals,t_vals,gamma,alt_vals1,\
                      #disp_vals1,v_vals1,a_vals1,t_vals1,gamma1,alt_vals2,\
                       #   disp_vals2,v_vals2,a_vals2,t_vals2,gamma2,alt_vals3,\
                        #      disp_vals3,v_vals3,a_vals3,t_vals3,gamma3)





#def g_acc(alt): 
#    
#    """Calculates acceleration due to gravity at a certain altitude.
#    params
#    input
#    alt [float] Altitude in m
#    output
#    acceleration in m/s^2
#    """
#    g = (G*Me)/((Re + alt)**2)
#    return g
#
#def density(alt):
#    
#    """Calculates density at a certain altitude.
#    params
#    input
#    alt [float] Altitude in m
#    output
#    density in kg/m^3
#    """
#    rho = surface_rho*np.exp(-beta*alt)
#    return rho
#
#def drag(v,rho):
#    
#    """Calculates drag force at certain flow conditions.
#    params
#    input
#    v [float] velocity in m/s
#    rho [float] density in kg/m^3
#    output
#    drag force in N
#    """
#    
#    Fd = 0.5*Cd*(v**2)*rho*S
#    return Fd
#
#def lift(v,rho):
#    
#    """Calculates lift force at certain flow conditions.
#    params
#    input
#    v [float] velocity in m/s
#    rho [float] density in kg/m^3
#    output
#    lfit force in N
#    """
#    
#    Fl = 0.5*Cl*(v**2)*rho*S
#    return Fl
#
#def TrajectorySolver(gamma):
#    
#    """Calculates Trajectory of Re-Entry Vehicle based on vehicle paramaters.
#    params
#    input
#    gamma - initial flight path angle in degrees
#    variables at top of code
#    output
#    altitude array in m
#    velocity array in m/s
#    time at which parachute deployed in s
#    """
#    
#    init_fpa = np.radians(gamma)
#    t_vals = np.linspace(0,time,int(steps))
#    
#    #sets up required values for loop with initial parameters
#    vx = vx_prev = init_v*np.cos(init_fpa)
#    vy = vy_prev = init_v*np.sin(init_fpa)
#    fpa = init_fpa
#    alt = alt_prev = init_alt
#    disp = disp_prev = init_disp
#    
#    #initialises arrays
#    v_vals = np.zeros(int(steps))
#    alt_vals = np.zeros(int(steps))
#    disp_vals = np.zeros(int(steps))
#    i = 0
#    
#    #for loop calculates trajectory
#    for i,t in enumerate(t_vals):
#
#        v = np.sqrt((vx**2)+(vy**2)) #Gets total velocity vector
#        alt_vals[i] = alt            #Adds altitude val to array
#        disp_vals[i] = disp
#        v_vals[i] = v                #Adds velocity val to array
#        rho = density(alt)           #Calculates density at this step
#        g = g_acc(alt)               #Calculates grav acc at this step
#        liftval = lift(v,rho)
#        dragval = drag(v,rho)        #Calculates Drag at this step
#        Fx = liftval*np.sin(fpa) - dragval*np.cos(fpa)#X Direction Force Calcs
#        Fy = mass*g - liftval*np.cos(fpa) \
#            - dragval*np.sin(fpa)   #Y Direction Force Calcs
#        
#        vx = vx_prev + ((Fx/mass) * dt) #Calcs v in x direction
#        if vx <= 0:     
#            vx = 0
#        vy = vy_prev + ((Fy/mass) * dt) #Calcs v in y direction
#        
#        if vx == 0:         #Calcs new flight path angle based on vx and vy
#            fpa = np.radians(90)
#        else:
#            fpa = np.arctan(vy/vx)
#        
#        alt = alt_prev - (vy*dt + 0.5*(Fy/mass)*(dt**2)) #Calcs new altitude
#        disp = disp_prev + (vx*dt + 0.5*(Fx/mass)*(dt**2)) #Calcs new x disp
#        vx_prev = vx   #Sets vx val for next loop
#        vy_prev = vy   #Sets vy val for next loop
#        alt_prev = alt #Sets alt val for next loop
#        disp_prev = disp
#        tlim = t
#        if alt <= 2000: #Breaks loop at parachute deployment
#        #    print("Parachute deployed at",round(alt,2), "m and Velocity = ",\
#        #          round(v,2),"m/s after", t, "seconds of flight time.")
#            break 
#    return alt_vals,disp_vals,v_vals,tlim
#
#    
#def a_val(v_vals):
#    
#    """Produces acceleration array based on velocity array.
#    params
#    input
#    v_vals array [floats] velocity in m/s
#    output
#    Acceleration array in g's'
#    """
#    
#    a = np.zeros(int(steps))
#    for i in range(int(steps)):
#        if i == 0:
#            a[i] = 0
#        else:
#            a[i] = v_vals[i-1] - v_vals[i]
#    a_vals = a/(9.81*dt)
#    return a_vals
#
#def array_cleaner(gamma):
#    
#    """Cleans extra zero values from TrajectorySolver arrays.
#    params
#    input
#    gamma - initial flight path angle in degrees
#    variables at top of code
#    output
#    TrajectorySolver and Acc arrays pruned of excess end values.
#    """
#    
#    t_vals = np.linspace(0,time,int(steps))
#    alt_vals,disp_vals,v_vals,tlim = TrajectorySolver(gamma)   
#    a_vals = a_val(v_vals)
#    i = 0
#    check = 0
#    
#
#    for i in range(int(steps)):
#        
#        #If the parachute hasn't deployed, do nothing and go to next timestep
#        if i <= (tlim/dt): 
#            check += 1    
#            continue
#        
#        #If parachute has deployed, delete unnecessary values
#        else:
#            alt_vals = np.delete(alt_vals,check)
#            disp_vals = np.delete(disp_vals,check)
#            v_vals = np.delete(v_vals,check)
#            a_vals = np.delete(a_vals,check)
#            t_vals = np.delete(t_vals,check)
#        
#    rho_vals = density(alt_vals) #Creates density array 
#        
#    return alt_vals,disp_vals,rho_vals,v_vals,a_vals,t_vals,gamma
#
#"""Running the Code"""
#alt_vals,disp_vals,rho_vals,v_vals,a_vals,t_vals,gamma = array_cleaner(10)
## alt_vals1,disp_vals1,rho_vals1,v_vals1,a_vals1,t_vals1,gamma1=array_cleaner(5)
## alt_vals2,disp_vals2,rho_vals2,v_vals2,a_vals2,t_vals2,gamma2=array_cleaner(10)
## alt_vals3,disp_vals3,rho_vals3,v_vals3,a_vals3,t_vals3,gamma3=array_cleaner(15)



###############################################################################
###############################################################################
###############################################################################
###############################################################################


"""Heat Flux Setup and iterative solution"""

#Setup
#Velocity_List = velocity(G_alt)
#Rho_List = atmo_density(G_alt)
  
# Heat_flux_master_list = [0, 0, 0, 0]
# Heat_load_master_list = [0, 0, 0, 0]
# q_conv_master_list = [0, 0, 0, 0]
# q_rerad_master_list = [0, 0, 0, 0]
# q_net_master_list = [0, 0, 0, 0]
# Ts_master_list = [0, 0, 0, 0]
# Heat_Load_in = []
# velocity_master_list = [v_vals, v_vals1, v_vals2, v_vals3]
# rho_master_list = [rho_vals, rho_vals1, rho_vals2, rho_vals3]
# t_vals_master_list = [t_vals, t_vals1, t_vals2, t_vals3]

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

Heat_Load_in = []
Heat_flux_list = []
Heat_load_list = []
q_conv_list = []
q_rerad_list = []
q_net_list = []
Ts = np.ones((len(v_vals), step))
Tc = np.ones((len(v_vals), step))
#Tal = np.ones((len(v_vals), step))
Ts[0,:] = 300
Tc[0,:] = 300
#Tal[:,:] = 300
#radius_in = noser - thicc
#pica_mass = pica_rho*((2*np.pi*(noser**3)/3)-(2*np.pi*(radius_in**3)/3))
#radius_al = noser - thicc - thicc_a
#pica_mass = pica_rho*((2*np.pi*(radius_in**3)/3)-(2*np.pi*(radius_al**3)/3))
k_pica = 0.577 #Thermal conductivity for carbon
k_glue = 0.55
Tc_List = []
Ts_List = []
#Tal_List = []
heatrad = []
heatcond = []
for i in range(len(v_vals)):
    Heat_flux_slices = 0
    Heat_load_slices = 0
    q_conv_stag = 7.455*(10**(-9))*(rho_vals[i]**0.4705)*(v_vals[i]**(3.089))*(noser**(-0.52))*10000
    q_conv_prev = q_conv_stag.copy()
    q_conv_list.append(q_conv_stag/(10**6))
    for j in range(step - 1):
        area = areas[j]
        q_conv_ang = q_conv_stag*(np.cos(np.radians(angles[j+1])))
        q_conv_slice = (q_conv_ang + q_conv_prev)/2
        if i == 0:
            q_rerad = sig*em*((300)**4)
            Ts[i,j] = 300
            q_cond = 0
            Tc[i,j] = 300 
            #Tal[i,j] = 300
            #Tal[i+1,j] = 300
        if i != 0:
            Temp = (Ts[i-1, j] + Ts[i-1,j+1])/2
            q_rerad = sig*em*((Temp)**4)
            q_cond = (((Ts[i-1,j]+Ts[i-1,j+1])/2) - ((Tc[i-1,j]+Tc[i-1,j+1])/2))*((k_pica * area)/thicc)
            Tc[i,j] = Tc[i-1,j] + (((q_cond) * dt) / (area * pica_rho * thicc * cpc))
            #if j == 0:
                #print("Pica Conduction =", round(q_cond, 4))
        """if i > 1:
            q_cond2 = (((Tc[i-1,j]+Tc[i-1,j+1])/2) - ((Tal[i-1,j]+Tal[i-1,j+1])/2)) * ((k_glue * area)/thicc_a)
            Tal[i,j] = Tal[i-1,j] + (((q_cond2) * dt) / (area * al_rho * thicc_a * cpa))
            #Tal[i,j] = Tal[i-1,j] + ((q_cond2*thicc_a*dt)/(0.55 * area))"""
        """if i > 1:
            q_cond2 = (((Tc[i-1,j]+Tc[i-1,j+1])/2) - ((Tal[i-1,j]+Tal[i-1,j+1])/2))*((0.55 * area)/thicc_a)
            Tal[i,j] = Tal[i-1,j] + (((q_cond2) * dt) / (area * al_rho * thicc_a * cpa))
            #Tal[i,j] = Tal[i-1,j] + ((q_cond2*thicc_a*dt)/(0.55 * area))
        """    
            #if j == 0:
                #print("Aluminium Conduction =", round(q_cond2, 4))
        Heat_flux_conv = q_conv_slice * area
        Heat_flux_rerad = q_rerad * area
        Heat_flux_cond = q_cond
        
        if i != 0:
            Ts[i,j] = Ts[i-1,j] + (((Heat_flux_conv - Heat_flux_rerad - Heat_flux_cond) * dt) / (area * pica_rho * thicc * cpc))
        
        Heat_Load_in += Heat_flux_conv * dt
        Heat_flux_slices += area * (q_conv_slice - q_rerad - q_cond)
        Heat_load_slices += area * (q_conv_slice - q_rerad - q_cond) * dt
        
        q_conv_prev = q_conv_ang.copy()
        if j == 0:    
            q_rerad_list.append(-q_rerad*(10**(-6)))
            q_net_list.append((q_conv_slice-q_rerad)*10**(-6))
            Tc_List.append(Tc[i,j])
            Ts_List.append(Ts[i,j])
            #Tal_List.append(Tal[i,j])
    Heat_flux_list.append(Heat_flux_slices/(10**6))
    Heat_load_list.append(Heat_load_slices)
    heatcond.append(Heat_flux_cond)

plt.figure()
plt.plot(t_vals,Ts_List, "r", label="Ts")
#plt.plot(t_vals,Tal_List, "b", label="Tal")
plt.plot(t_vals,Tc_List, "k", label="Tc")
plt.title("Nose Temperatures Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()

Talll = max(Tc_List)
Hottt = max(Ts_List)

print("\n\nFor a heat shield thickness of {h:0.5f} m".format(h = thicc))
print("\nThe maximum temperature of the aluminium = {g:0.3f} K".format(g = Talll))
print("\nThe maximum surface temperature of the heat shield = {c:0.3f} K\n\n".format(c = Hottt))



    
"""
#Stagnation Point Convective Heat Flux - Time
plt.figure()
plt.plot(t_vals,q_conv_list)
plt.title("Stagnation Point Convective Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.legend()
plt.show()

#Stagnation Reradiative Heat Flux - Time
plt.figure()
plt.plot(t_vals,q_rerad_list)
plt.title("Stagnation Point Reradiative Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.legend()
plt.show()

#Entire Surface Plotting - Time 
plt.figure()
plt.plot(t_vals,Heat_flux_list)
plt.title("Total Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW)')
plt.legend()
plt.show()
"""
"""
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
"""
"""
FPA_List = [0, 5, 10, 15]

# for i in range(len(FPA_List)): 
#     print("\n\n\nFPA =", FPA_List[i], ":")
#     print("\nThe maximum heat flux over the enitre surface =", round(max(Heat_flux_master_list[i]), 4), "MW")
#     print("\nThe maximum convective heat flux at the stagnation point =", round(max(q_conv_master_list[i]),4), "MW/m^2")
#     print("\nThe maximum reradiative heat flux at the stagnation point =", round(min(q_rerad_master_list[i]),4), "MW/m^2")
#     print("\nThe maximum net heat flux at the stagnation point =", round(max(q_net_master_list[i]), 4), "MW/m^2")
#     print("\nThe maximum temperature =", round(np.max(Ts_master_list[i]), 4), "K")
#     print("\nTotal heat load (SIMPS) =", round(sp.simps(t_vals_master_list[i] , Heat_flux_master_list[i]), 4), "MJ")
#     # print("\nThe maximum dynamic pressure =", round(np.max(dp), 4), "kPa/m^2")
#     print("\nThe total heat load =", round(Heat_Load_in_master_list[i], 4), "MJ")
#     print("\nThe total heat load =", round(dt*sum(q_conv_master_list[i]), 4), "MJ/m^2")
"""