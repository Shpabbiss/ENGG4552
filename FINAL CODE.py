
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 12:00:27 2021
Last updated on Tue Oct 19 16:00:23 2021

ENGG4552 Heat Shield Trajectory Code

Authors: Matthew Rickerby, Daniel Burgess, Clifford Asmussen, Nathan Holyoak,
         Daniel Kearney, Shanil Panchal, Tharang Addepalli

"""

# Import required modules
import numpy as np
import matplotlib.pyplot as plt

"""Initial Conditions"""
init_alt = 200E3            #[m]            Initial altitude
init_disp = 0               #[m]            Initial displacement
init_v = 7800               #[m/s]          Initial velocity

"""Constants"""
G = 6.67430E-11             #[Nm^2/kg^2]    Gravitational constant
Me = 5.97219E24             #[kg]           Mass of Earth
Re = 6371E3                 #[m]            Radius of Earth
y0 = 6620                   #[m]            Scale Height - AERO4800 Textbook
beta = 1/y0                 #               Scaling factor               
surface_rho = 1.225         #[kg/m^3]       Surface air density
pica_rho = 1155.54          #[kg/m^3]       Density of PICA
char_rho = 280              #[kg/m^3]       Density of Char (from research)
ins_rho = 280               #[kg/m^3]       Carbon felt density
al_rho = 2710               #[kg/m^3]       Aluminium density
fo_rho = 96                 #[kg/m^3]       Internal foam density
mass_cube = 1.33            #[kg]           CubeSat mass
mass_struc = 2              #[kg]           Aluminium structure mass (both sheets included)
mass_bolts = 1              #[kg]           Assumed bolt mass
mass_other = 5              #[kg]           Assumed mass of parachute and other components 

"""Vehicle Properties"""
radius = 0.14586            #[m]            Radius of capsule
ang = 87                    #[deg]          Angle prototype wraps to
scale = 1                   #               Steps per degree
step = ((ang * scale) + 1)  #               Steps in angles
angles = np.linspace(0, ang, step) #        Angle steps generation

noser = 0.07927             #[m]            Nose radius
Cd = 1.4                    #               Coefficient of drag (from cart 3d)
Cl = 0                      #               Coefficient of lift (ballistic entry)
S = 0.05106375              #[m^2]          Reference area

# Radiative properties
em = 0.9                    #               Emissivity (approximated as carbon emissivity)
sig = 5.67*(10**-8)         #[W/(m^2 k^4)]  Stefan-Boltzmann constant

# Specific Heat
cpc = 1925.56               #[J/kgK]        PICA (approximated as specific heat of carbon)
cpch = 4605.48              #[J/kgK]        Charr (from research)
cpin = 730                  #[J/kgK]        Carbon felt
cpfo = 1050                 #[J/kgK]        Insulative foam
cpal = 1050                  #[J/kgK]        Aluminium

# Thermal conductivity
k_pica = 0.577              #[W/mK]         PICA
k_char = 0.577              #[W/mK]         Charr (from research)
k_ins = 0.13                #[W/mK]         Carbon felt
k_foa = 0.05                #[W/mK]         Insulative foam


"""Thicknesses"""
pthc = 0.010                #[m]            PICA
sthc = 0.026                #[m]            Carbon felt
fthc = 0.0135               #[m]            Insulative foam
althc = 0.016               #[m]            Aluminium

#Open lists and append thickness around front heat shield
thic_p = []
thic_s = []
thic_f = fthc # Insulative foam
thic_al = althc
for i in range(step):
    thicccc = pthc # Constant PICA thickness
    thic_p.append(thicccc)
    
    # Calculate extra thickness in nose
    if i < 30.5: # In nose radius
        thicss = (-6.9921 * (10**-6) * i**2) + (1.9423*i*10**-6) + (0.032418 - sthc)
    else: # In flat part
        thicss = 0.00
    thic_s.append(thicss + sthc)
    
    
"""Loop Properties"""
dt = 0.05                   #[s]            Timestep size
time = 2000                 #[s]            Time of simulation
steps = time/dt + 1         #               Number of timesteps


###############################################################################
###############################################################################


"""Area Calculations"""
# Area iterative generator - Calculates sliced spherical area based on angle and radius
# Used for nose
def area_slice(r, t0, t1):
    A = 2 * np.pi * r * r * (1 - np.cos(np.radians(t1)) - (1 - np.cos(np.radians(t0))))
    return(A)

# Area iterative generator - Calculates sliced conical area based on angle and radius
# Used for flat front sides of shield
def areacone_slice(t1, t2):
    r1 = radius * (1 - np.cos(np.radians(t1)))
    r2 = radius * (1 - np.cos(np.radians(t2)))
    h1 = (r1 / np.tan(np.radians(t1))) - 0.0032
    h2 = (r2 / np.tan(np.radians(t2))) - 0.0032
    A1 = np.pi * r1  * np.sqrt((h1**2) + (r1**2))
    A2 = np.pi * r2  * np.sqrt((h2**2) + (r2**2))
    return(A2-A1)

# Start iterations at 0
ar_it = 0
vol_it = 0
vol_s = 0
mass_it = 0
i = 0
areas = []
conestart = 30.502          #[mm]           Heigh where curved nose transitiones to flat sides

# Iterate through angles appending area, volume and masses
for i in range(step - 1):
    if i <= conestart:
        x = area_slice(radius, angles[i], angles[i+1])
        y = x * thic_p[i]
        ar_it += x
        vol_it += y
        mass_it += y * pica_rho
        vol_s += x * thic_s[i]
        areas.append(x)
    if i > conestart:
        x = areacone_slice(angles[i], angles[i+1])
        y = x * thic_p[i]
        ar_it += x
        vol_it += y
        mass_it += y * pica_rho
        vol_s += x * thic_s[i]
        areas.append(x)

# Extract masses and areas
mass_in = vol_s * ins_rho
ar_ac_cir = 2 * np.pi * noser * noser * (1 - np.cos(np.radians(30.52)))
cone1 = np.pi * radius * np.sqrt((0.0859181**2) + (radius**2))
cone2 = np.pi * 0.04023 * np.sqrt((0.04023**2) + (0.0208181**2))
ar_ac_con = cone1 - cone2
ar_ac = ar_ac_cir + ar_ac_con

###############################################################################
###############################################################################

# Sum masses for total vehicle mass
mass = 4.32 #mass_it + mass_cube + mass_struc + mass_bolts + mass_in + mass_other
BC = mass/(S*Cd)            #               Ballistic coefficient

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
        #    print("Parachute deployed at",round(alt,2), "m and Velocity = ",\
        #          round(v,2),"m/s after", t, "seconds of flight time.")
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
alt_vals,disp_vals,rho_vals,v_vals,a_vals,t_vals,gamma = array_cleaner(10)

# Itertive Numerical Process
"""The for-loop cycles first through 'p' lists of
velocity and density values, each list being a 
different flight path angle (FPA). For each FPA, 
the loop then cycles through time steps, and for 
each time step the area slicing method is conducted
and summed."""

# Open arrays for heat loads and heat fluxes
Heat_Load_in = []
Heat_flux_list = []
Heat_load_list = []
q_conv_list = []
q_rerad_list = []
q_net_list = []

# Open matrices from temperature of each material
# Appends temperature for all angles for each timestep
Ts = np.ones((len(v_vals), step))
Tc = np.ones((len(v_vals), step))
Ti = np.ones((len(v_vals), step))
Ta = np.ones((len(v_vals), step))
Tal = np.ones((len(v_vals), step))

# Open matrix for ablation
# Appends ablation for all angles for each timestep
ab_mount = np.zeros((len(v_vals), step)) #[m]  
              
# Set intial temperatures              
Ts[:,:] = 300               #[K]    PICA
Tc[:,:] = 300               #[K]    Cubesat
Ti[:,:] = 300               #[K]    Carbon felt
Ta[:,:] = 300               #[K]    Foam
Tal[:,:] = 300              #[K]    Aluminium

# Open lists to append temperatures at different points
# Nose
Ti0_List = []
Tc0_List = []
Ts0_List = []
Ta0_List = []
Tal0_List = []

# 21 degrees
Ti21_List = []
Tc21_List = []
Ts21_List = []
Ta21_List = []
Tal21_List = []

# 42 degrees
Ti42_List = []
Tc42_List = []
Ts42_List = []
Ta42_List = []
Tal42_List = []

# 64 degrees
Ti64_List = []
Tc64_List = []
Ts64_List = []
Ta64_List = []
Tal64_List = []

# 86 degrees (Edge)
Ti86_List = []
Tc86_List = []
Ts86_List = []
Ta86_List = []
Tal86_List = []

# Open arrays for recorded values
# heatrad = []
# heatcond = []
Q_abl = []              #[MW/m^2]           Ablated heat flux
Q_abtot = []            #[MW]               Ablated heat flux
mass_abl = []           #[kg]               Ablated mass
line_ab = []            #[m]                Linear ablation
ab_rate = []            #[m/s]              Linear ablation rate at any angle
ab_rate1 = []           #[m/s]              Linear ablation rate at nose
ab_tot = np.zeros(step) #[m]                Total ablation
ab_time = []            #[s]                Ablation time

# Start max heating time
heattime = 0            #[s]
maxtime = 3 * dt        #[s]
mgtime = 3 * dt         #[s]
mTtime = 3 * dt         #[s]

# Iterate through each timestep
for i in range(len(v_vals)): # i is timestep
    # Set flux and load for this slice to 0
    Heat_flux_slices = 0
    Heat_load_slices = 0
    
    # Calculate stagnation point convection heat flux (W/m^2)
    q_conv_stag = 7.455*(10**(-9))*(rho_vals[i]**0.4705)*(v_vals[i]**(3.089))*(noser**(-0.52))*10000
    
    # If stagnation point greater than 0.5 MW/m^2, add dt to max heating time
    if q_conv_stag >= 0.5 * (10**6):
        heattime += dt

    # Copy and append stagnation convective heat flux
    q_conv_prev = q_conv_stag.copy()
    q_conv_list.append(q_conv_stag/(10**6)) # MW/m^2
            
        
    
    # For each timestep, iterate through each angle and calcualte fluxes and temperatures
    for j in range(step - 1): # j is angle
        # Set area for this angle
        area = areas[j]
        
        # Set heat flux at each angle using trend from stardust data
        if j < 15.75:
            q_conv_ang = q_conv_stag * (((1.01*(10**-5))*(j**4)) - (0.000293 * (j**3)) + (0.0022571 * (j**2)) - (0.019333 * j) + 1)
        if j >= 15.75:
            q_conv_ang = q_conv_stag * ((-0.001289 * j) + 0.75435)

        # Set convective heat flux using average over this area
        q_conv_slice = (q_conv_ang + q_conv_prev)/2
        
        # If first timestep, set all loads to starting
        if i == 0:
            q_rerad = sig*em*((300)**4)
            q_cond = 0
            q_cond_c = 0
            Q_abla = 0
            thic_char = 0
            abrate = 0
            if j == 0:
                Q_abl.append(Q_abla)
                mass_abl.append(0)

        # If not first timestep, calculate all loads and fluxes
        if i != 0:
            # Starting temperature is set as the last temperature for this slice
            Temp = Ts[i-1, j]
            
            # If temperature is less than ablation temperature, no ablation
            if Temp < (1100 + 273.15):
                Q_abla = 0                      #[MW/m^2] No ablated heat flux
                if j == 0:
                    Q_abl.append(Q_abla)        # Append heat flux per area
                mass_abl.append(0)              #[kg]     No mass ablated
                ab_tot[j] += 0                  #[m]      No ablation
                abrate = 0                      #[m/s]    No ablation rate
                
                # If temperature is more than ablation temperature, find ablation
            if Temp >= (1100 + 273.15):
                abrate = ( 0.00000001 * Temp - 0.000013732) #[m/s]   Ablation rate
                ab_rate.append(abrate)                      #        Append ablation rate
                abm = abrate * dt                           #[m]     Ablation amount
                ab_mount[i,j] = abm                         #[m]     Append ablation amount
                ab_tot[j] += abm                            #[m]     Add to this slices ablation amount
                Q_abla = pica_rho * abrate * 116299895.55   #[W/m^2] Ablative heat flux per unit area (from research)
                
                # Append values at nose
                if j == 0:
                    ab_time.append(dt)                      #[s]      Ablation time
                    Q_abl.append(Q_abla * 10**-6)           #[MW/m^2] Ablation rate
                    line_ab.append(abrate*dt)               #[m]      Linear ablation

                # Append total ablation amount
                Q_abtot.append(area*Q_abla*dt)                    #[J]  Total heat ablated
                mass_abl.append(dt * area* Q_abla / 116299895.55) #[kg] Mass ablated

            # Charr thickness
            thic_char = ab_tot[j]       

            # Calculate other heat fluxes
            # Reradiated heat flux
            q_rerad = sig*em*((Temp)**4) #[W/m^2]           
            
            # Conduction from top of PICA to top of carbon felt
            q_cond = ((Ts[i-1,j]) - (Ti[i-1,j]))*(((k_pica * (1 - thic_char)) + (k_char * thic_char)) * area)/(thic_p[j])
            
            # Conduction from top of carbon felt to top of insulative foam
            q_cond_ins = ((Ti[i-1,j]) - (Ta[i-1,j]))*((k_ins * area)/thic_s[j])

            # Conduction from top of insulative foam to cubesat
            q_cond_fo = ((Ta[i-1,j]) - (Tc[i-1,j]))*((k_foa * area)/thic_f)
            
            # Calculate internal temperatures at end of timestep from resulting net heat flux
            # Carbon felt
            Ti[i,j] = Ti[i-1,j] + (((q_cond - q_cond_ins) * dt) / (area * pica_rho * thic_p[j] * cpc))
            
            # Insulative foam 
            Ta[i,j] = Ta[i-1,j]  + (((q_cond_ins) * dt) / (area * ins_rho * thic_s[j] * cpin))
            
            # Aluminium structure
            Tal[i,j] = Tal[i-1,j]  + (((q_cond_ins) * dt) / (area * ins_rho * thic_al * cpal))
            
            # CubeSat
            Tc[i,j] = Tc[i-1,j]  + (((q_cond_fo) * dt) / (area * fo_rho * thic_f * cpfo))
            
            # To account for overreradiation at beginning set all temperatures to always above 300K
            if Ti[i,j] < 300:
                Ti[i,j] = 300
            if Ta[i,j] < 300:
                Ta[i,j] = 300
            if Tc[i,j] < 300:
                Tc[i,j] = 300
            if Tal[i,j] < 300:
                Tal[i,j] = 300
                
        # Calculate heat fluxes for area slice
        Heat_flux_conv = q_conv_slice * area    # W
        Heat_flux_rerad = q_rerad * area        # W
        Heat_flux_cond = q_cond                 # W
        Heat_flux_abla = Q_abla * area          # W
        
        # If not first timestep, find resulting surface temperature for next timestep from net heat flux
        if i != 0:
            Ts[i,j] = Ts[i-1,j] + (((Heat_flux_conv - Heat_flux_abla - Heat_flux_rerad - Heat_flux_cond) * dt) / (area * pica_rho * thic_p[j] * cpc))
            if Ts[i,j] < 300:
                Ts[i,j] = 300
        
        # Add loads and fluxes at surface
        Heat_Load_in += Heat_flux_conv * dt                                 #[J] Convective heat load 
        Heat_flux_slices += area * (q_conv_slice - q_rerad - Q_abla)        #[W] Net heat flux
        Heat_load_slices += area * (q_conv_slice - q_rerad - Q_abla) * dt   #[J] Net heat load
        
        # Set previous convection for calculating average flux over area at next angle
        q_conv_prev = q_conv_ang.copy()
        
        # Append temperatures at different angles for plotting
        # Nose (Stagnation point)
        if j == 0:    
            q_rerad_list.append(-q_rerad*(10**(-6)))                    #[MW/m^2] Reradiated heat flux
            q_net_list.append((q_conv_slice-q_rerad-Q_abla)*10**(-6))   #[MW/m^2] Net heat flux
            ab_rate1.append(abrate)                                     #[m/s]    Ablation rate
            Ti0_List.append(Ti[i,j])
            Tc0_List.append(Tc[i,j])
            Ts0_List.append(Ts[i,j])
            Ta0_List.append(Ta[i,j])
            Tal0_List.append(Tal[i,j])
        if j == 21:    
            Ti21_List.append(Ti[i,j])
            Tc21_List.append(Tc[i,j])
            Ts21_List.append(Ts[i,j])
            Ta21_List.append(Ta[i,j])
            Tal21_List.append(Tal[i,j])
        if j == 42:    
            Ti42_List.append(Ti[i,j])
            Tc42_List.append(Tc[i,j])
            Ts42_List.append(Ts[i,j])
            Ta42_List.append(Ta[i,j])
            Tal42_List.append(Tal[i,j])
        if j == 64:    
            Ti64_List.append(Ti[i,j])
            Tc64_List.append(Tc[i,j])
            Ts64_List.append(Ts[i,j])
            Ta64_List.append(Ta[i,j])
            Tal64_List.append(Tal[i,j])
        if j == 86:    
            Ti86_List.append(Ti[i,j])
            Tc86_List.append(Tc[i,j])
            Ts86_List.append(Ts[i,j])
            Ta86_List.append(Ta[i,j])
            Tal86_List.append(Tal[i,j])
    
    if i > 3:        
        if q_conv_list[-1] > q_conv_list[-2]:
            maxtime += dt
        if a_vals[i] > a_vals[i-1]:
            # print("AAAA")
            mgtime += dt
        if Ts[i,0] > Ts[i-1,0]:
            # print("BBBB")
            mTtime += dt
    # Append heat flux and load for each timesep
    Heat_flux_list.append(Heat_flux_slices/(10**6)) #[MW/m^2]
    Heat_load_list.append(Heat_load_slices)         #[J]


# Open lists to plot temperatures and values of interest
surftemp_plot = []      #[K]    Surface temperatures
sttemp_plot = []        #[K]    Structure temperatures
cubetemp_plot = []      #[K]    Cubesat temperatures
intemp_plot = []        #[K]    Carbon felt temperatures
altemp_plot = []        #[K]    Foam temperatures
ablation_plot = []      #[mm]   Ablation amount
angle_plot = []         #[deg]  Angles

# For each angle append: the max temperature, ablation amount, and angle
for n in range(step-1):
    surftemp_plot.append(max(Ts[:,n]))
    cubetemp_plot.append(max(Tc[:,n]))
    intemp_plot.append(max(Ti[:,n]))
    altemp_plot.append(max(Ta[:,n]))
    sttemp_plot.append(max(Tal[:,n]))
    ablation_plot.append(sum(ab_mount[:,n]) * 10**3)
    angle_plot.append(n)

# Set edge values in accordance with Stardust data
for n in range(1):
    angle_plot.append(86)
    surftemp_plot.append(300)
    cubetemp_plot.append(300)
    intemp_plot.append(300)
    altemp_plot.append(300)
    sttemp_plot.append(300)
    ablation_plot.append(0)


# Extract relevant data for results
altop = np.max(Ta) - 273.15       #[C]      Maximum insulative foam temperature
Talll = np.max(Ti) - 273.15       #[C]      Maximum carbon felt temperature
stTalll = np.max(Tal) - 273.15    #[C]      Maximum structure temperature
Hottt = np.max(Ts) - 273.15       #[C]      Maximum surface temperature
Tcube = np.max(Tc) - 273.15       #[C]      Maximum CubeSat temperature
maxcon = max(q_conv_list)         #[MW/m^2] Maximum convective heat flux
QMAX = max(q_net_list)            #[MW/m^2] Maximum net heat flux
HLOAD = (sum(Heat_flux_list)*dt)  #[MJ]     Total Heat Load 
ab_len = sum(line_ab) * 10**3     #[mm/s]   Length of ablation 
ab_tim = sum(ab_time)             #[s]      Time ablating 
ab_rates = Q_abl                  #[MW/m^2] Ablation heat fluxes at nose 


# Plot PICA thickness as a function of angle
plt.figure(453)
plt.plot(angle_plot, thic_p)
plt.title("PICA Thickness Plotted as a Function of Angle")
plt.xlabel('Angle (degree)')
plt.ylabel('PICA Thickness (mm)')
plt.show()

# Plot carbon felt thickness as a function of angle
plt.figure(457)
plt.plot(angle_plot, thic_s)
plt.title("Carbon Felt Thickness Plotted as a Function of Angle")
plt.xlabel('Angle (degree)')
plt.ylabel('PICA Thickness (mm)')
plt.show()
    
# Plot linear ablation amount as a function of angle
plt.figure(124242)
plt.plot(angle_plot, ablation_plot)
plt.title("Linear Ablation Amount Plotted as a Function of Angle")
plt.xlabel('Angle (degree)')
plt.ylabel('Linear Ablation (mm)')
plt.show()

# Plot maximum surface temperature as a function of angle
plt.figure(124243)
plt.plot(angle_plot, surftemp_plot)
plt.title("Maximum Surface Temperature Plotted as a Function of Angle")
plt.xlabel('Angle (degree)')
plt.ylabel('Temperature (K)')
plt.show()

# Plot maximum carbon felt temperature as a function of angle
plt.figure(124249)
plt.plot(angle_plot, intemp_plot)
plt.title("Maximum Carbon Felt Temperature Plotted as a Function of Angle")
plt.xlabel('Angle (degree)')
plt.ylabel('Temperature (K)')
plt.show()

# Plot maximum insulative foam temperature as a function of angle
plt.figure(124244)
plt.plot(angle_plot, altemp_plot)
plt.title("Maximum Insulative Foam Temperature Plotted as a Function of Angle")
plt.xlabel('Angle (degree)')
plt.ylabel('Temperature (K)')
plt.show()

# Plot maximum CubeSat temperature as a function of angle
plt.figure(1243)
plt.plot(angle_plot, cubetemp_plot)
plt.title("Maximum CubeSat Temperature Plotted as a Function of Angle")
plt.xlabel('Angle (degree)')
plt.ylabel('Temperature (K)')
plt.show()

# Plot nose temperatures as a function of time
plt.figure(0)
plt.plot(t_vals,Ts0_List, "r", label="PICA")
plt.plot(t_vals,Tal0_List, "k", label="Structure")
plt.plot(t_vals,Ti0_List, "b", label="Carbon Felt")
plt.plot(t_vals,Ta0_List, "c-.", label="Foam")
plt.plot(t_vals,Tc0_List, "g", label="CubeSat")
plt.title("Nose Temperatures Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()

# Plot temperatures at and anlge 21 degrees from nose as a function of time
plt.figure(21)
plt.plot(t_vals,Ts21_List, "r", label="PICA")
plt.plot(t_vals,Tal21_List, "k", label="Structure")
plt.plot(t_vals,Ta21_List, "c-.", label="Foam")
plt.plot(t_vals,Ti21_List, "b", label="Carbon Felt")
plt.plot(t_vals,Tc21_List, "g", label="CubeSat")
plt.title("21 Degrees Temperatures Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()

# Plot temperatures at and anlge 42 degrees from nose as a function of time
plt.figure(42)
plt.plot(t_vals,Ts42_List, "r", label="PICA")
plt.plot(t_vals,Tal42_List, "k", label="Structure")
plt.plot(t_vals,Ta42_List, "c-.", label="Foam")
plt.plot(t_vals,Ti42_List, "b", label="Carbon Felt")
plt.plot(t_vals,Tc42_List, "g", label="CubeSat")
plt.title("42 Degrees Temperatures Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()

# Plot temperatures at and anlge 64 degrees from nose as a function of time
plt.figure(64)
plt.plot(t_vals,Ts64_List, "r", label="PICA")
plt.plot(t_vals,Tal64_List, "k", label="Structure")
plt.plot(t_vals,Ta64_List, "c:", label="Foam")
plt.plot(t_vals,Ti64_List, "b", label="Carbon Felt")
plt.plot(t_vals,Tc64_List, "g", label="CubeSat")
plt.title("64 Degrees Temperatures Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()

# Plot temperatures at and anlge 86 degrees from nose as a function of time
plt.figure(86)
plt.plot(t_vals,Ts86_List, "r", label="PICA")
plt.plot(t_vals,Tal86_List, "k", label="Structure")
plt.plot(t_vals,Ta86_List, "c-.", label="Foam")
plt.plot(t_vals,Ti86_List, "b", label="Carbon Felt")
plt.plot(t_vals,Tc86_List, "g", label="CubeSat")
plt.title("86 Degrees Temperatures Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()

# Plot ablation heat flux at nose as a function of time
plt.figure(111)
plt.plot(t_vals, ab_rates )
plt.title("Ablation Heat Flux at Nose Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.show()

# Plot stagnation point convective heat flux as a function of time
plt.figure(1)
plt.plot(t_vals,q_conv_list)
plt.title("Stagnation Point Convective Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.show()

# Plot stagnation point reradiative heat flux as a function of time
plt.figure(2)
plt.plot(t_vals,q_rerad_list)
plt.title("Stagnation Point Reradiative Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.show()



# Plot stagnation point net heat flux as a function of time
plt.figure(6)
plt.plot(t_vals,q_net_list)
plt.title("Stagnation Point Net Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
#plt.legend()
plt.show()


# Plot total heat flux as a function of time
plt.figure(3)
plt.plot(t_vals,Heat_flux_list)
plt.title("Total Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW)')
plt.show()


# Plot ablation rate at nose as a function of convective heat flux
plt.figure(143320094)
plt.plot(q_conv_list, ab_rate1)
plt.xlabel("Convective Heat Flux at Nose (MW/m^2)")
plt.ylabel("Ablation Rate (mm/s)")
plt.title("Nose Ablation Rate Plotted as a Function of Convective Heat Flux")
plt.show()

# Plot stagnation point convective heat flux as a function of altitude
plt.figure(1235234234234)
plt.plot(alt_vals / 1000, q_conv_list)
plt.title("Stagnation Point Convective Heat Flux Plotted as a Function of Altitude")
plt.xlabel('Altitude (km)')
plt.ylabel('Heat Flux (MW/m^2)')
plt.show()

# Print results
print("\n\nMass = {g:0.5f} kg".format(g=mass))
print("PICA thickness = {g:0.3f} mm\nInsulator thickness = {gg:0.3f} mm\nFoam thickness = {ggg:0.3f} mm".format(g = pthc*10**3, gg = sthc*10**3, ggg = fthc*10**3))
print("\nThe maximum temperature of the CubeSat = {g:0.3f} C, {gg:0.3f} K".format(g = Tcube, gg = Tcube+273.15))
print("\nThe maximum temperature of the foam = {g:0.3f} C, {gg:0.3f} K".format(g = altop, gg = altop+273.15))
print("\nThe maximum temperature of the insulation = {g:0.3f} C, {gg:0.3f} K".format(g = Talll, gg = Talll+273.15))
print("\nThe maximum temperature of the structure = {g:0.3f} C, {gg:0.3f} K".format(g = stTalll, gg = stTalll+273.15))
print("\nThe maximum surface temperature of the heat shield = {g:0.3f} C, {gg:0.3f} K".format(g = Hottt, gg = Hottt+273.15))
print("\nThe maximum instantaneous convective heat flux through the shield is {g:0.6f} MW/m^2".format(g = maxcon))
print("\nThe maximum instantaneous net heat flux throught the shield is {g:0.6f} MW/m^2".format(g = QMAX))
print("\nThe total heat load is {g:0.6f} MJ".format(g = HLOAD))
print("\nThe nose ablates {g:0.6f} mm".format(g = ab_len))
print("\nThe nose ablates for {g:0.4f} seconds".format(g = ab_tim))
print("\nTotal mass ablated =", sum(mass_abl), "kg")
print("\nMaximum ablation rate =", max(ab_rate)*10**3, "mm/s")
print("\nTotal heat ablated =", sum(Q_abtot)*10**-6, "MJ")
print("\nMaximum heat ablation rate =", max(Q_abl), "MW/m^2\n\n")


plt.plot(t_vals, alt_vals/1000, 'k', label = "Trajectory")
plt.plot(t_vals[int(round(maxtime, 3)/0.05)], alt_vals[int(round(maxtime, 3)/0.05)]/1000, 'b^', label = "Max Convective Heat Flux")
plt.plot(t_vals[int(round(mgtime, 3)/0.05)], alt_vals[int(round(mgtime, 3)/0.05)]/1000, 'go', label = "Max G Load")
plt.plot(t_vals[int(round(mTtime, 3)/0.05)], alt_vals[int(round(mTtime, 3)/0.05)]/1000, 'rX', label = "Max Temperature")
plt.title("Altitude Plotted as a Function of Time")
plt.ylabel('Altitude (km)')
plt.xlabel('Time (s)')
plt.legend()
plt.show()

plt.plot(t_vals, v_vals/1000, 'k', label = "Trajectory")
plt.plot(t_vals[int(round(maxtime, 3)/0.05)], v_vals[int(round(maxtime, 3)/0.05)]/1000, 'b^', label = "Max Convective Heat Flux")
plt.plot(t_vals[int(round(mgtime, 3)/0.05)], v_vals[int(round(mgtime, 3)/0.05)]/1000, 'go', label = "Max G Load")
plt.plot(t_vals[int(round(mTtime, 3)/0.05)], v_vals[int(round(mTtime, 3)/0.05)]/1000, 'rX', label = "Max Temperature")
plt.title("Velocity Plotted as a Function of Altitude")
plt.ylabel('Velocity (m/s)')
plt.xlabel('Time (s)')
plt.legend()
plt.show()

plt.plot(disp_vals/1000, alt_vals/1000, 'k', label = "Trajectory")
plt.plot(disp_vals[int(round(maxtime, 3)/0.05)]/1000, alt_vals[int(round(maxtime, 3)/0.05)]/1000, 'b^', label = "Max Convective Heat Flux")
plt.plot(disp_vals[int(round(mgtime, 3)/0.05)]/1000, alt_vals[int(round(mgtime, 3)/0.05)]/1000, 'go', label = "Max G Load")
plt.plot(disp_vals[int(round(mTtime, 3)/0.05)]/1000, alt_vals[int(round(mTtime, 3)/0.05)]/1000, 'rX', label = "Max Heat Flux")
plt.title("Altitude Plotted as a Function of Downrange Displacement")
plt.ylabel('Altitude (km)')
plt.xlabel('Downrange Displacement (km)')
plt.legend()
plt.show()


# Plot nose temperatures as a function of time
plt.figure(121213224314324423)
plt.plot(t_vals,np.array(Tal0_List)-273.15, "k", label="Structure")
plt.title("Structure Maximum Temperature Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (C)')
plt.legend()
plt.show()
