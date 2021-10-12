# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 22:37:46 2021

@author: natha
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Define function to take data from excel
def DataTaker(file):
    Test_data = pd.read_excel(file, header=None)
    Test_data.columns = ['Time', 'Extension', 'Load']
    Test_time = Test_data['Time'].tolist()
    Test_ext= Test_data['Extension'].tolist()
    Test_load = Test_data['Load'].tolist()
    
    return Test_time, Test_ext, Test_load
    
#Call Desired Files For Data
time1, ext1, load1 = DataTaker('Specimen_RawData_1.xlsx')
time2, ext2, load2 = DataTaker('Specimen_RawData_2.xlsx')

#Plot Desired Graphs
plt.plot(ext1, load1, label = 'Specimen 1')
plt.plot(ext2, load2, label = 'Specimen 2')
plt.title('Load (N) vs Extension')
plt.xlabel('Extension (mm)')
plt.ylabel('Load (N)')
plt.legend()
plt.show()