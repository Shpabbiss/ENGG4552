import numpy as np
import matplotlib.pyplot as plt

"""Initial Conditions"""
init_alt = 200E3 #[m]
init_disp = 0    #[m]
init_v = 8000    #[m/s]
init_fpa = np.radians(1) #[radians]

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
Cd = 1
Cl = 0
S = np.pi*noser**2 #[m^2] - Reference Area
BC = mass/(S*Cd)

"""Loop Properties"""
dt = 0.05 #[s]
time = 2000 #[s]
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
    tlim [float] Parachute deployment time in s
    a_vals array [floats] acceleration in g's'
    
    output
    Altitude, Velocity, Acceleration vs Time
    """
    
    #Plot Altitude vs Time
    plt.plot(t_vals, alt_vals)
    plt.title("Altitude vs Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Altitude (m)")
    #plt.xlim(0,tlim)
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
        
    return alt_vals,disp_vals,v_vals,a_vals,t_vals

"""Running the Code"""
print(steps)
alt_vals,disp_vals,v_vals,a_vals,t_vals = array_cleaner()
plotter(t_vals,alt_vals,disp_vals,v_vals,a_vals)



