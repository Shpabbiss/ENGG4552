import numpy as np
import matplotlib.pyplot as plt

"""Initial Conditions"""
init_alt = 200E3 #[m]
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
dt = 1 #[s]
time = 1000 #[s]
steps = time/dt + 1


def g_acc(alt): 
    #This function calculates the acceleration due to gravity at the
    #specified altitude
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
    #calculates drag based upon flow conditions
    Fd = 0.5*Cd*(v**2)*rho*S
    return Fd

def wackycalcs():
    t_vals = np.linspace(0,time,int(steps))
    vx = init_v*np.cos(init_fpa)
    vy = init_v*np.sin(init_fpa)
    vx_prev = init_v*np.cos(init_fpa)
    vy_prev = init_v*np.sin(init_fpa)
    fpa = init_fpa
    alt = init_alt
    alt_prev = init_alt
    v_vals = np.zeros(int(steps))
    alt_vals = np.zeros(int(steps))
    i = 0
    for i,t in enumerate(t_vals):

        v = np.sqrt((vx**2)+(vy**2))
        alt_vals[i] = alt
        if t > 485 and t < 486:
            print("t =",t,"alt =",alt,"v = ",v, "fpa = ,",fpa)
        #print("v = ",v)
        #print("t =",t)
        v_vals[i] = v
        rho = density(alt)
        g = g_acc(alt)
        dragval = drag(v,rho)
        Fx = -dragval*np.cos(fpa)
        Fy = mass*g - dragval*np.sin(fpa)
        vx = vx_prev + ((Fx/mass) * dt)
        if vx <= 0:
            vx = 0
        vy = vy_prev + ((Fy/mass) * dt)
        
        if vx == 0:
            fpa = np.radians(90)
        else:
            fpa = np.arctan(vy/vx)
        alt = alt_prev - (vy*dt + 0.5*(Fy/mass)*(dt**2))
        vx_prev = vx
        vy_prev = vy
        alt_prev = alt
        #print(np.degrees(fpa))
        if alt <= 2000:
            print("Parachute Deployed at",alt)
            break 
    return alt_vals,v_vals
    
    
def plotter(alt_vals,v_vals):
    t_vals = np.linspace(0,time,int(steps))
    plt.plot(t_vals, alt_vals)
    #plt.xlim(440,460)
    plt.show()
    plt.plot(t_vals, v_vals)
    #plt.xlim(440,460)
    plt.show()

print(steps)
alt_vals,v_vals = wackycalcs()
plotter(alt_vals,v_vals)
        


