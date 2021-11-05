# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 12:00:27 2021
@author: cliff
"""

"""New Thermal Analysis Code - Thermal_Analysis_Final may be bugged"""
import matplotlib.pyplot as plt
import numpy as np
def runbaby(fpa):
    #import scipy.integrate as sp
    import numpy as np
    import matplotlib.pyplot as plt
    
    """Initial Conditions"""
    init_alt = 200E3            #[m]
    init_disp = 0               #[m]
    init_v = 7800               #[m/s]
    
    """Constants"""
    G = 6.67430E-11             #[Nm^2/kg^2]
    Me = 5.97219E24             #[kg]
    Re = 6371E3                 #[m]
    y0 = 6620                   #[m] - Scale Height - AERO4800 Textbook
    beta = 1/y0
    surface_rho = 1.225         #[kg/m^3]
    pica_rho = 280              #[kg/m^3]  (ASSUMED PICA DENSITY - WILL NEED CHANGING)
    char_rho = 280              #[kg/m^3]  (ASSUMED CHAR DENSITY - WILL NEED CHANGING)
    ins_rho = 280               #[kg/m^3] Insualtor density
    al_rho = 2710               #[kg/m^3] Aluminium density
    fo_rho = 96                 #[kg/m^3] Internal foam density
    mass_cube = 1.33            #[kg] CubeSat mass
    mass_struc = 2        #[kg] Aluminium structure mass (both sheets included)
    mass_bolts = 0            #[kg] Assumed bolt mass
    mass_other = 5              #[kg] Assumed mass of parachute and other components 
    
    """Vehicle Properties"""
    
    
    radius = 0.14586
    ang = 87
    #ang = 90 # Angle to in degrees
    scale = 1 # Steps per degree
    step = ((ang * scale) + 1)
    angles = np.linspace(0, ang, step)
    
    
    # Thickness
    pthc = 0.010 # PICA
    sthc = 0.026 # DPCF516 Insulator
    fthc = 0.0135 # Foam
    
    thic_p = []
    thic_s = []
    for i in range(step):
        thicccc = pthc
        thic_p.append(thicccc)
        if i < 30.5:
            thicss = (-6.9921 * (10**-6) * i**2) + (1.9423*i*10**-6) + (0.032418 - sthc)
        else:
            thicss = 0.00
        thic_s.append(thicss + sthc)
        # print(i, pthc, thicss)
        
    #print(thic_p)
    #thic_p = 0.0125
    thic_f = fthc
    #mass_in = 0  # Insulator mass
    #mass = 10 #mass_cube + mass_struc + mass_p + mass_i                  #[kg] Total vehicle mass            (ASSUMED - WILL NEED CHANGING)
    noser = 0.07927                 #[m] Nose radius                    (ASSUMED - WILL NEED CHANGING)
    Cd = 1.4                      #Coefficient of drag                (ASSUMED - WILL NEED CHANGING)
    Cl = 0                      #Coefficient of lift                (ASSUMED - WILL NEED CHANGING)
    S = 0.05106375          #[m^2] - Reference Area
    
    em = 0.9                    #Emissivity                         (ASSUMED - WILL NEED CHANGING)
    sig = 5.67*(10**-8)         #[W / (m^2 k^4)] Stefan-Boltzmann constant
    
            #Cp's:      Carbon = 710 [J/kgK]            Aluminium = 900 [J/kgK]   
             
    cpc = 1925.56  # 900            #[J/kgK] Specific heat of carbon    (ASSUMED - WILL NEED CHANGING)
    cpch = 4605.48                  #[J/kgK] Specific heat of charr 
    cpin = 730                      #[J/kgK] Specific heat of insulator
    cpfo = 1050                     #[J/kgK] Specific heat of foam
    #cpa = 1050            #[m] Heat shield thickness          (ASSUMED - WILL NEED CHANGING) (Potentially do this as a function later on as likely will be thicker at nose)
    #thicc_a = 0.012
    #thic_p = 0.0125
    
    """Loop Properties"""
    dt = 2                  #[s]
    time = 204000                 #[s]
    steps = time/dt + 1
    
    
    ###############################################################################
    ###############################################################################
    
    
    """Area Calculations"""
    
    #Area Slicing Parameters
    
    
    # Area iterative generator - Calculates sliced area based on angle and radius
    def area_slice(r, t0, t1):
        A = 2 * np.pi * r * r * (1 - np.cos(np.radians(t1)) - (1 - np.cos(np.radians(t0))))
        return(A)
    
    def areacone_slice(t1, t2):
        r1 = radius * (1 - np.cos(np.radians(t1)))
        r2 = radius * (1 - np.cos(np.radians(t2)))
        h1 = (r1 / np.tan(np.radians(t1))) - 0.0032
        h2 = (r2 / np.tan(np.radians(t2))) - 0.0032
        A1 = np.pi * r1  * np.sqrt((h1**2) + (r1**2))
        A2 = np.pi * r2  * np.sqrt((h2**2) + (r2**2))
        # Surface area of a cone = π * r * (r + sqrt((h**2) + (r**2)))
        return(A2-A1)
    
    #Area Calculation Sanity Check
    
    ar_it = 0
    vol_it = 0
    vol_s = 0
    
    mass_it = 0
    i = 0
    areas = []
    conestart = 30.502
    
    for i in range(step - 1):
        #print(i)
        if i <= conestart:
            x = area_slice(radius, angles[i], angles[i+1])
            y = x * thic_p[i]
            #print("a", x, y)
            ar_it += x
            vol_it += y
            mass_it += y * pica_rho
            vol_s += x * thic_s[i]
            areas.append(x)
        if i > conestart:
            x = areacone_slice(angles[i], angles[i+1])
            y = x * thic_p[i]
            #print("b", x, y)
            ar_it += x
            vol_it += y
            mass_it += y * pica_rho
            vol_s += x * thic_s[i]
            areas.append(x)
    
    mass_in = vol_s * ins_rho
    ar_ac_cir = 2 * np.pi * noser * noser * (1 - np.cos(np.radians(30.52)))
    
    cone1 = np.pi * radius * np.sqrt((0.0859181**2) + (radius**2))
    cone2 = np.pi * 0.04023 * np.sqrt((0.04023**2) + (0.0208181**2))
    
    ar_ac_con = cone1 - cone2
    ar_ac = ar_ac_cir + ar_ac_con
    
    #print(cone1, cone2, ar_ac_cir, ar_ac_con)
    """
    print("\nSanity Check: Analytical Area = Iterative Area (Exactly)")
    print("Analytical Area =", round(ar_ac, 8), "m^2")
    print("Iterative Area =", round(ar_it, 8), "m^2")
    print("Volume =", round(vol_it , 8), "m^3")
    print("Starting mass of ablative material =", round(mass_it, 8), "kg\n")"""
    ###############################################################################
    ###############################################################################
    
    mass = mass_it + mass_cube + mass_struc + mass_bolts + mass_in + mass_other
    BC = mass/(S*Cd)            #Ballistic coefficient
    
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
    alt_vals,disp_vals,rho_vals,v_vals,a_vals,t_vals,gamma = array_cleaner(fpa)
    # alt_vals1,disp_vals1,rho_vals1,v_vals1,a_vals1,t_vals1,gamma1=array_cleaner(5)
    # alt_vals2,disp_vals2,rho_vals2,v_vals2,a_vals2,t_vals2,gamma2=array_cleaner(10)
    # alt_vals3,disp_vals3,rho_vals3,v_vals3,a_vals3,t_vals3,gamma3=array_cleaner(15)
    
    
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ###############################################################################
    
    ###################################
    
    
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
    #Ts = np.ones((len(v_vals), step))
    #Ts[:,0] = 300
    
    # Pyrolysis starts at 900C
    # Ablation starts at 1100C
    
    
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
    Ti = np.ones((len(v_vals), step))
    Ta = np.ones((len(v_vals), step))
    #Tal = np.ones((len(v_vals), step))
    Ts[:,:] = 300
    Tc[:,:] = 300 # Cubesat
    Ti[:,:] = 300
    Ta[:,:] = 300
    #Tal[:,:] = 300
    #radius_in = noser - thicc
    #pica_mass = pica_rho*((2*np.pi*(noser**3)/3)-(2*np.pi*(radius_in**3)/3))
    #radius_al = noser - thicc - thicc_a
    #pica_mass = pica_rho*((2*np.pi*(radius_in**3)/3)-(2*np.pi*(radius_al**3)/3))
    k_pica = 0.577 #Thermal conductivity for carbon
    # k_glue = 0.55
    k_char = 0.577
    k_ins = 0.13 #0.3385 # Thermal conductivity
    k_foa = 0.05 # Thermal conductiviy
    
    Ti0_List = []
    Tc0_List = []
    Ts0_List = []
    Ta0_List = []
    
    
    Ti225_List = []
    Tc225_List = []
    Ts225_List = []
    Ta225_List = []
    
    Ti45_List = []
    Tc45_List = []
    Ts45_List = []
    Ta45_List = []
    
    Ti675_List = []
    Tc675_List = []
    Ts675_List = []
    Ta675_List = []
    
    Ti90_List = []
    Tc90_List = []
    Ts90_List = []
    Ta90_List = []
    
    #Tal_List = []
    heatrad = []
    heatcond = []
    Q_abl = []
    Q_abtot = []
    mass_abl = []
    line_ab = []
    ab_rate = []
    ab_tot = np.zeros(step)
    ab_mount = np.zeros((len(v_vals), step))
    ab_time = []
    for i in range(len(v_vals)): # i is timestep
        Heat_flux_slices = 0
        Heat_load_slices = 0
        q_conv_stag = 7.455*(10**(-9))*(rho_vals[i]**0.4705)*(v_vals[i]**(3.089))*(noser**(-0.52))*10000
        q_conv_prev = q_conv_stag.copy()
        q_conv_list.append(q_conv_stag/(10**6))
        for j in range(step - 1): # j is angle
            area = areas[j]
            if j < 15.75:
                q_conv_ang = q_conv_stag * (((1.01*(10**-5))*(j**4)) - (0.000293 * (j**3)) + (0.0022571 * (j**2)) - (0.019333 * j) + 1)
                #print("a", q_conv_stag, q_conv_ang)
            if j >= 15.75:
                q_conv_ang = q_conv_stag * ((-0.001289 * j) + 0.75435)
                #print("b", q_conv_stag, q_conv_ang)
            #q_conv_ang = q_conv_stag*(np.cos(np.radians(angles[j+1])))
            q_conv_slice = (q_conv_ang + q_conv_prev)/2
            if i == 0:
                q_rerad = sig*em*((300)**4)
                Ts[i,j] = 300 # MIGHT NOT NEEDS THESE
                Tc[i,j] = 300 
                Ti[i,j] = 300 
                Ta[i,j] = 300
                q_cond = 0
                q_cond_c = 0
                Q_abla = 0
                thic_char = 0
                if j == 0:
                    Q_abl.append(Q_abla)
                    mass_abl.append(0)
                #Tal[i,j] = 300
                #Tal[i+1,j] = 300
            if i != 0:
                Temp = Ts[i-1, j]# + Ts[i-1,j+1])/2
                #"""if Temp > 900:print("t =", i*dt, "angle =", ((j-1)*0.5), "pyrolysis")"""
                #if j == 0: 
                    #print(Ts[i,j])
                if Temp < (1100 + 273.15):
                    Q_abla = 0
                    if j == 0:
                        Q_abl.append(Q_abla)
                    mass_abl.append(0)
                    ab_tot[j] += 0
                if Temp >= (1100 + 273.15):
                    abrate = ( 0.00000001 * Temp - 0.000013732)  # *    0  # m/s
                    ab_rate.append(abrate)
                    abm = abrate * dt
                    ab_mount[i,j] = abm
                    ab_tot[j] += abm
                    Q_abla = pica_rho * abrate * 116299895.55 #232599791.1 #116299895.55
                    
                    if j == 0:
                        ab_time.append(dt)
                        Q_abl.append(Q_abla * 10**-6)
                        line_ab.append(abrate*dt)
                        #print(abrate)
                    Q_abtot.append(area*Q_abla*dt)
                    mass_abl.append(dt * area* Q_abla / 232599791.1)
                    #if j == 0:
                    #    print(Temp, abrate, Q_abla, "\n\n")
                q_rerad = sig*em*((Temp)**4)
                #q_cond = (((Ts[i-1,j]+Ts[i-1,j+1])/2) - ((Ti[i-1,j]+Ti[i-1,j+1])/2))*((k_pica * area)/thic_p[j]) # Conduction from surface to top of insulator
                #q_cond_ins = (((Ti[i-1,j]+Ti[i-1,j+1])/2) - ((Ta[i-1,j]+Ta[i-1,j+1])/2))*((k_ins * area)/thic_s) # Conduction from top of insulator to top of structure
                thic_char = ab_tot[j]
                q_cond = ((Ts[i-1,j]) - (Ti[i-1,j]))*(((k_pica * (1 - thic_char)) + (k_char * thic_char)) * area)/(thic_p[j]) # Conduction from surface to top of insulator
                q_cond_ins = ((Ti[i-1,j]) - (Ta[i-1,j]))*((k_ins * area)/thic_s[j]) # Conduction from top of insulator to top of foam
                q_cond_fo = ((Ta[i-1,j]) - (Tc[i-1,j]))*((k_foa * area)/thic_f) # Conduction from top of foam to top of cubesat
                Ti[i,j] = Ti[i-1,j] + (((q_cond - q_cond_ins) * dt) / (area * pica_rho * thic_p[j] * cpc))
                Ta[i,j] = Ta[i-1,j]  + (((q_cond_ins) * dt) / (area * ins_rho * thic_s[j] * cpin))
                Tc[i,j] = Tc[i-1,j]  + (((q_cond_fo) * dt) / (area * fo_rho * thic_f * cpfo))
                
                if Ti[i,j] < 300:
                    Ti[i,j] = 300
                if Ta[i,j] < 300:
                    Ta[i,j] = 300
                if Tc[i,j] < 300:
                    Tc[i,j] = 300                
                #print(j)
                #if j == 0:
                    #print("Pica Conduction =", round(q_cond, 4))
            """if i > 1:
                q_cond2 = (((Ti[i-1,j]+Ti[i-1,j+1])/2) - ((Tal[i-1,j]+Tal[i-1,j+1])/2)) * ((k_glue * area)/thicc_a)
                Tal[i,j] = Tal[i-1,j] + (((q_cond2) * dt) / (area * al_rho * thicc_a * cpa))
                #Tal[i,j] = Tal[i-1,j] + ((q_cond2*thicc_a*dt)/(0.55 * area))"""
            """if i > 1:
                q_cond2 = (((Ti[i-1,j]+Ti[i-1,j+1])/2) - ((Tal[i-1,j]+Tal[i-1,j+1])/2))*((0.55 * area)/thicc_a)
                Tal[i,j] = Tal[i-1,j] + (((q_cond2) * dt) / (area * al_rho * thicc_a * cpa))
                #Tal[i,j] = Tal[i-1,j] + ((q_cond2*thicc_a*dt)/(0.55 * area))
            """    
                #if j == 0:
                    #print("Aluminium Conduction =", round`(q_cond2, 4))
            Heat_flux_conv = q_conv_slice * area
            Heat_flux_rerad = q_rerad * area
            Heat_flux_cond = q_cond 
            Heat_flux_abla = Q_abla * area
            #if j == 0:
                #print(Ts[i-1,j], q_conv_slice, q_rerad, q_cond, Q_abla, "\n")
            
            if i != 0:
                Ts[i,j] = Ts[i-1,j] + (((Heat_flux_conv - Heat_flux_abla - Heat_flux_rerad - Heat_flux_cond) * dt) / (area * pica_rho * thic_p[j] * cpc))
                if Ts[i,j] < 300:
                    Ts[i,j] = 300
            Heat_Load_in += Heat_flux_conv * dt
            Heat_flux_slices += area * (q_conv_slice - q_rerad - Q_abla)
            Heat_load_slices += area * (q_conv_slice - q_rerad - Q_abla) * dt
            
            q_conv_prev = q_conv_ang.copy()
            if j == 0:    
                q_rerad_list.append(-q_rerad*(10**(-6)))
                q_net_list.append((q_conv_slice-q_rerad-Q_abla)*10**(-6))
                Ti0_List.append(Ti[i,j])
                Tc0_List.append(Tc[i,j])
                Ts0_List.append(Ts[i,j])
                Ta0_List.append(Ta[i,j])
            if j == 44:    
                Ti225_List.append(Ti[i,j])
                Tc225_List.append(Tc[i,j])
                Ts225_List.append(Ts[i,j])
                Ta225_List.append(Ta[i,j])
            if j == 89:    
                Ti45_List.append(Ti[i,j])
                Tc45_List.append(Tc[i,j])
                Ts45_List.append(Ts[i,j])
                Ta45_List.append(Ta[i,j])
            if j == 134:    
                Ti675_List.append(Ti[i,j])
                Tc675_List.append(Tc[i,j])
                Ts675_List.append(Ts[i,j])
                Ta675_List.append(Ta[i,j])
            if j == 179:    
                Ti90_List.append(Ti[i,j])
                Tc90_List.append(Tc[i,j])
                Ts90_List.append(Ts[i,j])
                Ta90_List.append(Ta[i,j])
                
        Heat_flux_list.append(Heat_flux_slices/(10**6))
        Heat_load_list.append(Heat_load_slices)
        heatcond.append(Heat_flux_cond)
    
    surftemp_plot = []
    cubetemp_plot = []
    intemp_plot = []
    altemp_plot = []
    ablation_plot = []
    angle_plot = []
    for n in range(step-1):
        surftemp_plot.append(max(Ts[:,n]))
        cubetemp_plot.append(max(Tc[:,n]))
        intemp_plot.append(max(Ti[:,n]))
        altemp_plot.append(max(Ta[:,n]))
        ablation_plot.append(sum(ab_mount[:,n]) * 10**3)
        angle_plot.append(n)
    
    for n in range(1):
        angle_plot.append(90)
        surftemp_plot.append(300)
        cubetemp_plot.append(300)
        intemp_plot.append(300)
        altemp_plot.append(300)
        ablation_plot.append(0)
        
    
    # plt.figure(453)
    # plt.plot(angle_plot, thic_p)
    # plt.title("PICA Thickness Plotted as a Function of Angle")
    # plt.xlabel('Angle (degree)')
    # plt.ylabel('PICA Thickness (mm)')
    # plt.show()
    
    
    # plt.figure(457)
    # plt.plot(angle_plot, thic_s)
    # plt.title("Mat Thickness Plotted as a Function of Angle")
    # plt.xlabel('Angle (degree)')
    # plt.ylabel('PICA Thickness (mm)')
    # plt.show()
    
    
        
    
    # plt.figure(124242)
    # plt.plot(angle_plot, ablation_plot)
    # plt.title("Linear Ablation Plotted as a Function of Angle")
    # plt.xlabel('Angle (degree)')
    # plt.ylabel('Linear Ablation (mm)')
    # plt.show()
    
    # plt.figure(124243)
    # plt.plot(angle_plot, surftemp_plot)
    # plt.title("Maximum Surface Temperature Plotted as a Function of Angle")
    # plt.xlabel('Angle (degree)')
    # plt.ylabel('Temperature (K)')
    # plt.show()
    
    # plt.figure(124249)
    # plt.plot(angle_plot, intemp_plot)
    # plt.title("Maximum Insulation Temperature Plotted as a Function of Angle")
    # plt.xlabel('Angle (degree)')
    # plt.ylabel('Temperature (K)')
    # plt.show()
    
    # plt.figure(124244)
    # plt.plot(angle_plot, altemp_plot)
    # plt.title("Maximum Foam Temperature Plotted as a Function of Angle")
    # plt.xlabel('Angle (degree)')
    # plt.ylabel('Temperature (K)')
    # plt.show()
    
    # plt.figure(1243)
    # plt.plot(angle_plot, cubetemp_plot)
    # plt.title("Maximum CubeSat Temperature Plotted as a Function of Angle")
    # plt.xlabel('Angle (degree)')
    # plt.ylabel('Temperature (K)')
    # plt.show()
    
    # """
    # ab0 = sum(ab_mount[:,0])
    # ab225 = sum(ab_mount[:,44])
    # ab45 = sum(ab_mount[:,89])
    # ab675 = sum(ab_mount[:,134])
    # ab90 = sum(ab_mount[:,179])
    # ablation_plot = []
    # ablation_plot.append(ab0)
    # ablation_plot.append(ab225)
    # ablation_plot.append(ab45)
    # ablation_plot.append(ab675)
    # ablation_plot.append(ab90)
    # print(ablation_plot)
    # """
    
    # plt.figure(0)
    # plt.plot(t_vals,Ts0_List, "r", label="PICA")
    # #plt.plot(t_vals,Tal_List, "b", label="Tal")
    # plt.plot(t_vals,Tc0_List, "g", label="CubeSat")
    # plt.plot(t_vals,Ti0_List, "b", label="Insulator")
    # plt.plot(t_vals,Ta0_List, "k", label="Aluminium")
    # plt.title("Nose Temperatures Plotted as a Function of Time")
    # plt.xlabel('Time (s)')
    # plt.ylabel('Temperature (K)')
    # plt.legend()
    # plt.show()
    # """
    # plt.figure(22.5)
    # plt.plot(t_vals,Ts225_List, "r", label="PICA")
    # #plt.plot(t_vals,Tal_List, "b", label="Tal")
    # #plt.plot(t_vals,Tc225_List, "g", label="Char")
    # plt.plot(t_vals,Ti225_List, "b", label="Insulator")
    # plt.plot(t_vals,Ta225_List, "k", label="Aluminium")
    # plt.title("22.5 Degrees Temperatures Plotted as a Function of Time")
    # plt.xlabel('Time (s)')
    # plt.ylabel('Temperature (K)')
    # plt.legend()
    # plt.show()
    # plt.figure(45)
    # plt.plot(t_vals,Ts45_List, "r", label="PICA")
    # #plt.plot(t_vals,Tal_List, "b", label="Tal")
    # #plt.plot(t_vals,Tc45_List, "g", label="Char")
    # plt.plot(t_vals,Ti45_List, "b", label="Insulator")
    # plt.plot(t_vals,Ta45_List, "k", label="Aluminium")
    # plt.title("45 Degrees Temperatures Plotted as a Function of Time")
    # plt.xlabel('Time (s)')
    # plt.ylabel('Temperature (K)')
    # plt.legend()
    # plt.show()
    # plt.figure(67.5)
    # plt.plot(t_vals,Ts675_List, "r", label="PICA")
    # #plt.plot(t_vals,Tal_List, "b", label="Tal")
    # #plt.plot(t_vals,Tc675_List, "g", label="Char")
    # plt.plot(t_vals,Ti675_List, "b", label="Insulator")
    # plt.plot(t_vals,Ta675_List, "k", label="Aluminium")
    # plt.title("67.5 Degrees Temperatures Plotted as a Function of Time")
    # plt.xlabel('Time (s)')
    # plt.ylabel('Temperature (K)')
    # plt.legend()
    # plt.show()
    # plt.figure(90)
    # plt.plot(t_vals,Ts90_List, "r", label="PICA")
    # #plt.plot(t_vals,Tal_List, "b", label="Tal")
    # #plt.plot(t_vals,Tc90_List, "g", label="Char")
    # plt.plot(t_vals,Ti90_List, "b", label="Insulator")
    # plt.plot(t_vals,Ta90_List, "k", label="Aluminium")
    # plt.title("90 Degrees Temperatures Plotted as a Function of Time")
    # plt.xlabel('Time (s)')
    # plt.ylabel('Temperature (K)')
    # plt.legend()
    # plt.show()
    # """
    altop = np.max(Ta) - 273.15         
    Talll = np.max(Ti) - 273.15       # Max aluminium temperature (Celsius)
    Hottt = np.max(Ts) - 273.15       # Max surface temperature 
    Hotbox = np.max(Tc) - 273.15      # Max CubeSat temperature 
    # #chartop = max(Tc0_List)
    # maxcon = max(q_conv_list)
    QMAX = max(q_net_list)               # Max instantaneous heat flux (MW/m^2)
    HLOAD = (sum(Heat_flux_list)*dt)     # Total Heat Load (MJ)
    # ab_len = sum(line_ab) * 10**3        # Length of ablation (mm/s)
    # ab_tim = sum(ab_time)                # Time ablating (seconds)
    # ab_rates = Q_abl                     # Ablation heat fluxes at nose (MW/m^2)
    
    # #print("\n\nFor a heat shield thickness of {h:0.5f} m".format(h = thic_p))
    # print("\n\nMass = {g:0.5f} kg".format(g=mass))
    # print("PICA thickness = {g:0.3f} mm\nInsulator thickness = {gg:0.3f} mm\nFoam thickness = {ggg:0.3f} mm".format(g = pthc*10**3, gg = sthc*10**3, ggg = fthc*10**3))
    # print("\nThe maximum temperature of the CubeSat = {g:0.3f} C, {gg:0.3f} K".format(g = Hotbox, gg = Hotbox+273.15))
    # print("\nThe maximum temperature of the foam = {g:0.3f} C, {gg:0.3f} K".format(g = altop, gg = altop+273.15))
    # print("\nThe maximum temperature of the insulation = {g:0.3f} C, {gg:0.3f} K".format(g = Talll, gg = Talll+273.15))
    # #print("\nThe maximum temperature of the char = {c:0.3f} K".format(c = chartop))
    # print("\nThe maximum surface temperature of the heat shield = {g:0.3f} C, {gg:0.3f} K".format(g = Hottt, gg = Hottt+273.15))
    # print("\nThe maximum instantaneous convective heat flux throught the shield is {g:0.6f} MW/m^2".format(g = maxcon))
    # print("\nThe maximum instantaneous net heat flux throught the shield is {g:0.6f} MW/m^2".format(g = QMAX))
    # print("\nThe total heat load is {g:0.6f} MJ".format(g = HLOAD))
    # print("\nThe nose ablates {g:0.6f} mm".format(g = ab_len))
    # print("\nThe nose ablates for {g:0.4f} seconds".format(g = ab_tim))
    # print("\nTotal mass ablated =", sum(mass_abl), "kg")
    # #######print("\nMaximum ablation rate =", max(ab_rate)*10**3, "mm/s")
    # print("\nTotal heat ablated =", sum(Q_abtot)*10**-6, "MJ")
    # print("\nMaximum heat ablation rate =", max(Q_abl), "MW/m^2\n\n")
    
    
    
    # plt.figure(111)
    # plt.plot(t_vals,ab_rates)
    # #plt.plot(t_vals,Tal_List, "b", label="Tal")
    # plt.title("Ablation Heat Flux Plotted as a Function of Time")
    # plt.xlabel('Time (s)')
    # plt.ylabel('Heat Flux (MW/m^2)')
    # #plt.legend()
    # plt.show()
    
    # #Stagnation Point Convective Heat Flux - Time
    # plt.figure(1)
    # plt.plot(t_vals,q_conv_list)
    # plt.title("Stagnation Point Convective Heat Flux Plotted as a Function of Time")
    # plt.xlabel('Time (s)')
    # plt.ylabel('Heat Flux (MW/m^2)')
    # #plt.legend()
    # plt.show()
    
    # #Stagnation Reradiative Heat Flux - Time
    # plt.figure(2)
    # plt.plot(t_vals,q_rerad_list)
    # plt.title("Stagnation Point Reradiative Heat Flux Plotted as a Function of Time")
    # plt.xlabel('Time (s)')
    # plt.ylabel('Heat Flux (MW/m^2)')
    # #plt.legend()
    # plt.show()
    
    
    
    # #Nose Net Heat Flux - Time 
    # plt.figure(6)
    # plt.plot(t_vals,q_net_list)
    # plt.title("Net Heat Flux at Nose Plotted as a Function of Time")
    # plt.xlabel('Time (s)')
    # plt.ylabel('Heat Flux (MW/m^2)')
    # #plt.legend()
    # plt.show()
    
    
    # #Entire Surface Plotting - Time 
    # plt.figure(3)
    # plt.plot(t_vals,Heat_flux_list)
    # plt.title("Total Heat Flux Plotted as a Function of Time")
    # plt.xlabel('Time (s)')
    # plt.ylabel('Heat Flux (MW)')
    # #plt.legend()
    # plt.show()
    
    
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
    return fpa, Ta0_List, Hotbox+273.15, Ts0_List, HLOAD, a_vals, q_net_list,t_vals

fpa1, tempfoam1, tempcube1, tempsurf1,totalheat1,gvals1,heatflux1,t_vals1 = runbaby(8)
fpa2, tempfoam2, tempcube2, tempsurf2,totalheat2,gvals2,heatflux2,t_vals2 = runbaby(10)
fpa3, tempfoam3, tempcube3, tempsurf3,totalheat3,gvals3,heatflux3,t_vals3 = runbaby(20)

plt.figure(figsize=(5.54,5.54 ))
plt.plot(t_vals1,tempfoam1,label = '8 degrees')
plt.plot(t_vals2,tempfoam2,label = '10 degrees')
plt.plot(t_vals3,tempfoam3,label = '20 degrees')
plt.title("Structural Temperature vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.legend()
plt.show()

plt.figure(figsize=(5.54,5.54 ))
plt.plot(t_vals1,gvals1,label = '8 degrees')
plt.plot(t_vals2,gvals2,label = '10 degrees')
plt.plot(t_vals3,gvals3,label = '20 degrees')
plt.title("Deceleration vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Deceleration (g's)")
plt.legend()
plt.show()


plt.figure(figsize=(5.54,5.54 ))
plt.plot(t_vals1,tempsurf1,label = '8 degrees')
plt.plot(t_vals2,tempsurf2,label = '10 degrees')
plt.plot(t_vals3,tempsurf3,label = '20 degrees')
plt.title("Surface Temperature vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.legend()
plt.show()

plt.figure(figsize=(5.54,5.54 ))
plt.plot(t_vals1,heatflux1,label = '8 degrees')
plt.plot(t_vals2,heatflux2,label = '10 degrees')
plt.plot(t_vals3,heatflux3,label = '20 degrees')
plt.title("Heat Flux vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Heat Flux (MW/m^2)")
plt.legend()
plt.show()

#angles = []
#tempfoam = []
#tempcube = []
#tempsurf = []
#totalheat = []
#gvals = []
#heatflux = []
#def checker():
#    for i in np.linspace(0.1,25):
#        a,tf,tc,ts,th,av,hf = runbaby(i) 
#        angles.append(a)
#        tempfoam.append(tf)
#        tempcube.append(tc)
#        tempsurf.append(ts)
#        totalheat.append(th)
#        gvals.append(max(av))
#        heatflux.append(hf)
#checker()
#
#
#
#FoamFail = []
#SurfaceFail = []
#
#def wrecker():
#    for i in enumerate(tempfoam):
#        if tempfoam[i[0]] >= 200 + 273:
#            FoamFail.append(angles[i[0]])
#def mecker():
#    for i in enumerate(tempsurf):
#        if tempsurf[i[0]] >= 3650 + 273:
#            SurfaceFail.append(angles[i[0]])
#wrecker()
#mecker()    
#
#print('Foam Temperature Failure Angles =',FoamFail)
#print('Surface Temperature Failure Angles =',SurfaceFail)  
#    
#plt.plot(angles,tempfoam)
#plt.title("Foam Max Temperature Plotted as a Function of Flight Path Angle")
#plt.xlabel('Angle (degree)')
#plt.ylabel('Temperature (K)')
##plt.legend()
#plt.show()
#
#plt.plot(angles,tempcube)
#plt.title("Cubesat Max Temperature Plotted as a Function of Flight Path Angle")
#plt.xlabel('Angle (degree)')
#plt.ylabel('Temperature (K)')
##plt.legend()
#plt.show()
#
#plt.plot(angles,tempsurf)
#plt.title("Surface Max Temperature Plotted as a Function of Flight Path Angle")
#plt.xlabel('Angle (degree)')
#plt.ylabel('Temperature (K)')
##plt.legend()
#plt.show()
#
#plt.plot(angles,totalheat)
#plt.title("Total Heat Load Plotted as a Function of Flight Path Angle")
#plt.xlabel('Angle (degree)')
#plt.ylabel('Heat (MJ)')
##plt.legend()
#plt.show()
#
#plt.plot(angles,heatflux)
#plt.title("Max Instantaneous Heat Flux Plotted as a Function of\
#          Flight Path Angle")
#plt.xlabel('Angle (degree)')
#plt.ylabel('Heat Flux (MW/m^2)')
##plt.legend()
#plt.show()
#
#plt.plot(angles,gvals)
#plt.title("Max g-load as a Function of Flight Path Angle")
#plt.xlabel('Angle (degree)')
#plt.ylabel('g-load (g)')
##plt.legend()
#plt.show()
