"""
Created on Thu Jun 17 12:00:27 2021
@author: cliff
"""

"""New Thermal Analysis Code - Thermal_Analysis_Final may be bugged"""

#import scipy.integrate as sp
import numpy as np
import matplotlib.pyplot as plt

h_flux = 5.6                  # Test heat flux (MW/m^2)
time = 80                   # Test time (s)
pthic = 18                  # Sample thickness (mm)    
diam = 20.675               # Sample diameter (mm)

pica_rho = 1155.54              #[kg/m^3]  (ASSUMED PICA DENSITY - WILL NEED CHANGING)
char_rho = 280              #[kg/m^3]  (ASSUMED CHAR DENSITY - WILL NEED CHANGING)


"""Vehicle Properties"""


radius = diam/(2 * 1000)
areas = np.pi * radius * radius
# print(angles)

thic_p = pthic / 1000 # PICA
torch_flux = h_flux * 10**6

em = 0.9                    #Emissivity                         (ASSUMED - WILL NEED CHANGING)
sig = 5.67*(10**-8)         #[W / (m^2 k^4)] Stefan-Boltzmann constant

        #Cp's:      Carbon = 710 [J/kgK]            Aluminium = 900 [J/kgK]   
         
cpc = 1925.56  # 900            #[J/kgK] Specific heat of carbon    (ASSUMED - WILL NEED CHANGING)
cpch = 4605.48                  #[J/kgK] Specific heat of charr 


"""Loop Properties"""
dt = 0.005                   #[s]
# time = 80                 #[s]
steps = time/dt + 1
testpoints = int(steps)


# Area iterative generator - Calculates sliced area based on radius
def area_slice(r1, r2):
    A = np.pi * ((r2**2) - (r1**2))
    return(A)


###############################################################################
###############################################################################



Heat_Load_in = []
Heat_flux_list = []
Heat_load_list = []
q_conv_list = []
q_rerad_list = []
q_net_list = []
Ts = np.ones((testpoints + 1))
Ti = np.ones((testpoints + 1))
Ts[:] = 300
Ti[:] = 300
k_pica = 0.577 #Thermal conductivity for carbon
k_char = 0.577

Ti0_List = []
Ts0_List = []

heatrad = []
heatcond = []
Q_abl = []
Q_abtot = []
mass_abl = []
line_ab = []
ab_rate = []
ab_mount = np.zeros((testpoints + 1))
ab_time = []
t_vals = []
ab_tot = 0
Heat_Load_in = 0

for i in range(testpoints + 1): # i is timestep
    t_vals.append(i*dt)
    Heat_flux_slices = 0
    Heat_load_slices = 0
    q_conv_stag = torch_flux
    q_conv_list.append(q_conv_stag/(10**6))
    area = areas
    q_conv_slice = (q_conv_stag)
    if i == 0:
        q_rerad = sig*em*((300)**4)
        Ts[i] = 300 # MIGHT NOT NEEDS THESE
        Ti[i] = 300 
        q_cond = 0
        Q_abla = 0
        thic_char = 0
    
        Q_abl.append(Q_abla)
        mass_abl.append(0)
    if i != 0:        
        # print(i)
        Temp = Ts[i-1]
        if Temp < (1100 + 273.15):
            Q_abla = 0
            Q_abl.append(Q_abla)
            mass_abl.append(0)
            ab_tot += 0
        if Temp >= (1100 + 273.15):
            abrate = ( 0.00000001 * Temp - 0.000013732)  # *    0  # m/s
            ab_rate.append(abrate)
            abm = abrate * dt
            ab_mount[i] = abm
            ab_tot += abm
            Q_abla = pica_rho * abrate * 116299895.55 #232599791.1 #116299895.55
            ab_time.append(dt)
            Q_abl.append(Q_abla * 10**-6)
            line_ab.append(abrate*dt)
            Q_abtot.append(area*Q_abla*dt)
            mass_abl.append(dt * area* Q_abla / 232599791.1)
            #if j == 0:
            #    print(Temp, abrate, Q_abla, "\n\n")
        q_rerad = sig*em*((Temp)**4)
        thic_char = ab_tot
        q_cond = ((Ts[i-1]) - (Ti[i-1]))*(((k_pica * (1 - thic_char)) + (k_char * thic_char)) * area)/(thic_p) # Conduction from surface to top of insulator
        q_cond_ins = sig*em*(((Ti[i-1])-300)**4)
        if q_cond < q_cond_ins:
            print(i*dt, q_cond - q_cond_ins)
        Ti[i] = Ti[i-1] + (((q_cond - q_cond_ins) * dt) / (area * pica_rho * thic_p * cpc))
        
        if Ti[i] < 300:
            Ti[i] = 300  
            
        Heat_flux_conv = q_conv_stag * area
        Heat_flux_rerad = q_rerad * area
        Heat_flux_cond = q_cond 
        Heat_flux_abla = Q_abla * area
        
        if i != 0:
            Ts[i] = Ts[i-1] + (((Heat_flux_conv - Heat_flux_abla - Heat_flux_rerad - Heat_flux_cond) * dt) / (area * pica_rho * thic_p * cpc))
            if Ts[i] < 300:
                Ts[i] = 300
        Heat_Load_in += Heat_flux_conv * dt
        Heat_flux_slices += area * (q_conv_slice - q_rerad - Q_abla)
        Heat_load_slices += area * (q_conv_slice - q_rerad - Q_abla) * dt
       
    q_rerad_list.append(-q_rerad*(10**(-6)))
    q_net_list.append((q_conv_slice-q_rerad-Q_abla)*10**(-6))
    Ti0_List.append(Ti[i])
    Ts0_List.append(Ts[i])
    
    Heat_flux_cond = q_cond        
    Heat_flux_list.append(Heat_flux_slices/(10**6))
    Heat_load_list.append(Heat_load_slices)
    heatcond.append(Heat_flux_cond)


plt.figure(0)
plt.plot(t_vals,Ts0_List, "r", label="PICA")
plt.plot(t_vals,Ti0_List, "b", label="Rear of PICA")
plt.title("Temperatures Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()

plt.figure(0)
plt.plot(t_vals,Ts0_List, "r")
plt.title("PICA Temperatures Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
# plt.legend()
plt.show()


plt.figure(0)
plt.plot(t_vals,Ti0_List, "b")
plt.title("Rear of PICA Temperatures Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
# plt.legend()
plt.show()

Talll = np.max(Ti) - 273.15       # Max aluminium temperature (Celsius)
Hottt = np.max(Ts) - 273.15       # Max surface temperature 
#chartop = max(Tc0_List)
maxcon = max(q_conv_list)
QMAX = max(q_net_list)               # Max instantaneous heat flux (MW/m^2)
HLOAD = (sum(Heat_flux_list)*dt)     # Total Heat Load (MJ)
ab_len = sum(line_ab) * 10**3        # Length of ablation (mm/s)
ab_tim = sum(ab_time)                # Time ablating (seconds)
ab_rates = Q_abl                     # Ablation heat fluxes at nose (MW/m^2)
ab_mmps = ab_len / time
ab_gps = (sum(mass_abl) * 1000/time)


print("\n\nPICA thickness = {g:0.3f} mm".format(g = thic_p*10**3))
print("\nThe maximum temperature of the rear of the PICA = {g:0.3f} C, {gg:0.3f} K".format(g = Talll, gg = Talll+273.15))
print("\nThe maximum surface temperature of the heat shield = {g:0.3f} C, {gg:0.3f} K".format(g = Hottt, gg = Hottt+273.15))
print("\nThe maximum instantaneous convective heat flux throught the shield is {g:0.6f} MW/m^2".format(g = maxcon))
print("\nThe total heat load is {g:0.6f} MJ".format(g = HLOAD))
print("\nThe nose ablates {g:0.6f} mm".format(g = ab_len))
print("\nThe nose ablates for {g:0.4f} seconds".format(g = ab_tim))
print("\nTotal mass ablated =", sum(mass_abl)*1000, "g")
print("\nAblation amount averaged over time =", ab_mmps, "mm/s or", ab_gps, "g/s")
print("\nTotal heat ablated =", sum(Q_abtot)*10**-6, "MJ")
print("\nMaximum heat ablation rate =", max(Q_abl), "MW/m^2\n\n")


plt.figure(111)
plt.plot(t_vals,ab_rates)
#plt.plot(t_vals,Tal_List, "b", label="Tal")
plt.title("Ablation Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
#plt.legend()
plt.show()

#Stagnation Point Convective Heat Flux - Time
plt.figure(1)
plt.plot(t_vals,q_conv_list)
plt.title("Stagnation Point Convective Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
#plt.legend()
plt.show()

#Stagnation Reradiative Heat Flux - Time
plt.figure(2)
plt.plot(t_vals,q_rerad_list)
plt.title("Stagnation Point Reradiative Heat Flux Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
#plt.legend()
plt.show()



#Nose Net Heat Flux - Time 
plt.figure(6)
plt.plot(t_vals,q_net_list)
plt.title("Net Heat Flux at Nose Plotted as a Function of Time")
plt.xlabel('Time (s)')
plt.ylabel('Heat Flux (MW/m^2)')
#plt.legend()
plt.show()


