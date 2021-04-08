"""Initial Trajectory Bullshit"""

init_altitude = 400*10**3 #[m]
init_velocity = 8000    #[m/s]
surface_density = 1.225 #[kg/m^3]
gas_constant = 287.05   #[J/kg K]
abs_temp = 0 + 273.15   #[K] - Not 100% sure what this should be??
g = 9.81


def scaleheight():   #Need to check this formula - conflicting information
    b = - (g/gas_constant)/(abs_temp)
    return(b)

    