import numpy as np
import matplotlib.pyplot as plt

"""Initial Conditions"""
init_alt = 200E3 #[m]
init_disp = 0    #[m]
init_v = 7800    #[m/s]
#init_fpa = np.radians(0) #[radians]

"""Constants"""
G = 6.67430E-11  #[Nm^2/kg^2]
Me = 5.97219E24  #[kg]
Re = 6371E3      #[m]
y0 = 6620        #[m] - Scale Height - AERO4800 Textbook
beta = 1/y0
surface_rho = 1.225 #[kg/m^3]

"""Vehicle Properties"""
mass = 10        #[kg]
noser = 0.3      #[m]
Cd = 1.4
Cl = 0
S = 0.05106375 #[m^2] - Reference Area
#BC = mass/(S*Cd)

"""Loop Properties"""
dt = 2 #[s]
time = 204000 #[s]
steps = time/dt + 1


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

def drag(v,rho,Cd):
    
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

def amb_temp(h):
    if h < 11000:
        T = (15) - (6.5/1000)*h
    elif h < 25000:
        T = -56.5
    elif h < 50000:
        T = (-56.5) + (2.16/1000)*(h-25000)
    elif h < 75000:
        T = (-2.5) - (3.6/1000)*(h-50000)
    else:
        T = -92.5
    return T

def Mach(v,h):
    a = np.sqrt(1.4*287*(amb_temp(h)+273.15))
    M = v/a
    return M


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
        
        # if Mach(v,alt) >= 13:
        #     Cd = 1.3
        # elif Mach(v,alt) < 13:
        #     Cd = 1.45
        
        liftval = lift(v,rho)
        dragval = drag(v,rho,Cd)        #Calculates Drag at this step
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
i = np.argmax(a_vals)
print("Maximum Deceleration is",round(max(a_vals),2),"g's, occuring at an altitude of",round(alt_vals[i],2),"m and a velocity of",round(v_vals[i],2),"m/s.")


#alt_vals1,disp_vals1,rho_vals1,v_vals1,a_vals1,t_vals1,gamma1=array_cleaner(5)
#alt_vals2,disp_vals2,rho_vals2,v_vals2,a_vals2,t_vals2,gamma2=array_cleaner(10)
#alt_vals3,disp_vals3,rho_vals3,v_vals3,a_vals3,t_vals3,gamma3=array_cleaner(15)
plotter(t_vals,alt_vals,disp_vals,v_vals,a_vals)
#plot_comparisons(alt_vals,disp_vals,v_vals,a_vals,t_vals,gamma,alt_vals1,\
                      #disp_vals1,v_vals1,a_vals1,t_vals1,gamma1,alt_vals2,\
                       #   disp_vals2,v_vals2,a_vals2,t_vals2,gamma2,alt_vals3,\
                        #      disp_vals3,v_vals3,a_vals3,t_vals3,gamma3)


