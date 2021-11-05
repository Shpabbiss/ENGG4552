# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 15:33:49 2021

@author: Tharang Addepalli
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

var = pd.read_excel("book1.xlsx")



x = list(var['x']) #Distance of heatshield 
y = list(var['y']) # Mach 14 Pressure
z = list(var['r']) # Y distance  
q = list(var['p(x) M6']) # Mach 6 Pressure 
w = list(var['Temp M6']) # Mach 6 Temperatures 
u = list(var['Temp M14']) # Mach 14 Temperatures 


# Plotting Mach 14 Pressures 
plt.plot(x,y) # Plot 
plt.xlabel("x (m)") # X axis label 
plt.ylabel("Surface Pressure (Pa)") # Y axis Label 
plt.title("Pressure on Surface of Heatshield") # Plot Title 
plt.savefig("Surface Pressure (Pa)",dpi=300, bbox_inches='tight') # Save Image 
plt.show()

# Comapring Surface Pressures  
plt.plot(x,y) # Plot Mach 14
plt.plot(x,q) # Plot Mach 6
plt.xlabel("x (m)") # X axis label
plt.ylabel("Surface Pressure (Pa)") # Y axis label
plt.title("Comparing Surface Pressure of Various Mach Numbers") # Plot Title
plt.legend([ "Mach 14", "Mach 6"]) # Legend
plt.savefig("Comparing Surface Pressure (Pa)",dpi=300, bbox_inches='tight') # Save Image 
plt.show()


# Comapring Temperature Pressures  
plt.plot(x,u)
plt.plot(x,w)
plt.xlabel("x (m)")
plt.ylabel("Fluid Temperature (K)")
plt.title("Comparing Fluid Temperatures of Various Mach Numbers")
plt.legend([ "Mach 14","Mach 6"])
plt.savefig("Comparing Temperature (Pa)",dpi=300)
plt.show()

plt.plot(x,z)
plt.xlabel("x (m)")
plt.ylabel("r (x)")
plt.title("Shape of Heatshield")
plt.savefig("Shape of Heatshield",dpi=300)
plt.show()




