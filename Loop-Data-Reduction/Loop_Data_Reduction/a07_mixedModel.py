import os
from os.path import join as pjoin  # tool to make paths
import numpy as np
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import pathlib
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from os.path import join as pjoin  # tool to make paths
import itertools


def mixedConvection(c, k, mu, rho, deltaTwall2bulk, heatflux, d_h, CTE):
    h = heatflux / deltaTwall2bulk
    t1 = 0.075 * ((k) / (c * mu)) ** (2 / 3) + 1
    t2 = (d_h ** (3) * h ** (3)) / (k ** (3))
    t3n = 0.18 * (9.81 * CTE * deltaTwall2bulk * rho ** 2 / mu ** 2) ** (3 / 4) * ((c * mu) / (k)) ** (3 / 4)
    t3d = (0.63 * (k / (c * mu)) ** (9 / 16) + 1) ** (4 / 3)
    t4 = rho * d_h * c ** (2 / 3)
    v = (4.7 * mu ** (1 / 3)) * (k * t1 ** (3 / 4) * (t2 - t3n / t3d)) ** (2 / 3) / t4
    w = mu * v ** 2
    FOM = w ** -1
    # print((4.7 * mu ** (1 / 3)) ,(k * t1 ** (3 / 4) ) ,t2 ,t2 - t3n / t3d,(t2 - t3n / t3d) **(2 / 3))
    # print(t1,t2,t3n,t3d,t4,v,w)
    return FOM


datapath = '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/'
data = pd.read_csv(pjoin(datapath, 'measuredFluidPropertiestoFit.csv'))
mixedConvection_vec = np.vectorize(mixedConvection)

if 1==1:
    result = mixedConvection_vec(c=data['cv_nf'],
                                 k=data['k_nf'],
                                 mu=data['Measured μ'],
                                 rho=data['rho_nf'],
                                 deltaTwall2bulk=25,
                                 heatflux=1.175 * 100 ** 2,
                                 d_h=.652 / 100,
                                 CTE=data['CTE'])
    data["FOM"] = result
    data["FOM"] = data["FOM"] /data["FOM"][1]
    print(data)

par_space = [["fluid",[0,1,2]],
    ['deltaTwall2bulk', np.linspace(5,100,10)],
    ['heatflux', np.asarray([1.17]) *100 ** 2],
    ]


parameter = []
possible_vals = []

for par in par_space:
    parameter.append(par[0])
    possible_vals.append(par[1])
par_tests = list(itertools.product(*possible_vals))
par_tests = np.vstack(par_tests)
print(par_tests)

combos_data = pd.DataFrame(data=par_tests, columns=parameter)
fluid_dat_cols = ['cv_nf','k_nf','Measured μ','rho_nf','CTE']
for i in fluid_dat_cols:
    combos_data[i] = [data[i][f] for f in combos_data['fluid'].astype(int)]

combos_data["FOM"] = mixedConvection_vec(c=combos_data['cv_nf'],
                                 k=combos_data['k_nf'],
                                 mu=combos_data['Measured μ'],
                                 rho=combos_data['rho_nf'],
                                 deltaTwall2bulk=combos_data['deltaTwall2bulk'],
                                 heatflux=combos_data['heatflux'],
                                 d_h=.652 / 100,
                                 CTE=combos_data['CTE'])

import seaborn as sns
for i in np.unique(combos_data["deltaTwall2bulk"]):
    set0 = combos_data[combos_data["deltaTwall2bulk"]==i]
    combos_data["FOM"][combos_data["deltaTwall2bulk"]==i] =  set0["FOM"]/np.min(set0["FOM"])


sns.scatterplot('deltaTwall2bulk','FOM', hue='fluid',data=combos_data)
plt.show()
print(combos_data)

