import os
from os.path import join as pjoin  # tool to make paths
import numpy as np
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import pathlib
from scipy.optimize import curve_fit

datapath = "/Users/hanscastorp/Dropbox/Thing/MIT/nano fluids/Lab Data/Lubrizol Data/MIT fluids August 2021_reformat.xlsx"
fluids = ['Isoparaffin','OS622339H','OS625475','OS625473','OS625474',"Water"]
for fluid in fluids:
    data = pd.read_excel(datapath, comment='#', sheet_name=fluid,engine='openpyxl')

    # print(data)
    def VogelFletcherTammann(x,A,T0,n):

        return np.multiply(n,np.exp(np.divide(A,x-T0)))
    xdata,ydata = data["Temperature"], data["Viscosity"]
    popt, pcov = curve_fit(VogelFletcherTammann, xdata+273.15,ydata,bounds=([0,100,0], [2000, 300, 10**-3]))
    print(fluid,popt)
    plt.scatter(xdata,ydata,s=5,label=fluid)
    plt.plot(xdata, VogelFletcherTammann(xdata+273.15, *popt), label=fluid + '; fit: a=%5.3f, b=%5.6f, c=%5.6f' % tuple(popt))

plt.plot(np.linspace(25,60,100),VogelFletcherTammann(np.linspace(25,60,100)+273.15,657.66,	150.9,	0.0000551),label = "Isoparaffin 2")
plt.legend()
plt.ylabel('Viscosity [Pa•s]')
plt.xlabel('Temperature [°C]')
plt.show()