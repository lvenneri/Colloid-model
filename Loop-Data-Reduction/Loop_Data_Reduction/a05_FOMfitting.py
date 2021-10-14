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


def fitFunc(x, a, b, c, d):
    return x[0] ** a * x[1] ** b * x[2] ** c / x[3] ** d


def fitFuncAlphaBeta(x, a, b):
    return x[0] ** (2 * b / a) * x[1] ** 2 * x[2] ** (2 * (1 - b) / a) / x[3] ** ((2 - 2 * b) / a)


def fitFuncCTE(x, a, b, c, d, e):
    return x[0] ** a * x[1] ** b * x[2] ** c * x[4] ** e / x[3] ** d


def mixedConvection(c, k, mu, rho, deltaTwall2bulk, heatflux, d_h,CTE ):
    h = heatflux / deltaTwall2bulk
    t1 = 0.075 * ((k) / (c * mu)) ** (2 / 3) + 1
    t2 = (d_h ** (3) * h ** (3)) / (k ** (3))
    t3n = 0.18 * (9.81*CTE*deltaTwall2bulk*rho**2/mu**2) ** (3 / 4) * ((c * mu) / (k)) ** (3 / 4)
    t3d = (0.63 * (k / (c * mu)) ** (9 / 16) + 1) ** (4 / 3)
    t4 = rho * d_h * c ** (2 / 3)
    v = (4.7 * mu ** (1 / 3)) * (k * t1 ** (3 / 4) * (t2 - t3n / t3d)) **(2 / 3) / t4
    w = mu * v ** 2
    FOM = w ** -1
    # print((4.7 * mu ** (1 / 3)) ,(k * t1 ** (3 / 4) ) ,t2 ,t2 - t3n / t3d,(t2 - t3n / t3d) **(2 / 3))
    # print(t1,t2,t3n,t3d,t4,v,w)
    return FOM


datapath = '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/'
data = pd.read_csv(pjoin(datapath, 'measuredFluidPropertiestoFit.csv'))

if 1 == 0:
    print("Standard")
    # plt.scatter(data["cv_nf"],data["cv_nf"],data["cv_nf"],data["cv_nf"])
    bounds = (.01, 30)

    x_3d = np.array([data["cv_nf"].values,
                     data["rho_nf"].values,
                     data["k_nf"].values,
                     data["Measured μ"].values])

    y = data["Loop FOM"]
    x = [4 / 3, 2, 2, 5 / 3]
    fitParams, fitCovariances = curve_fit(fitFunc, x_3d, y, x, method="trf", bounds=bounds)
    print(fitParams)
    # fit the T vs pump power data
    # popt, pcov = curve_fit(funcFOM, fluid_T_wall_agg, np.log(fluid_pumpPower_agg),
    #                        bounds=bounds)

    data["fitted"] = fitFunc([data["cv_nf"].values,
                              data["rho_nf"].values,
                              data["k_nf"].values,
                              data["Measured μ"].values], *fitParams)

    print(data)

if 1 == 0:
    print('___________________Alpha Beta')
    # plt.scatter(data["cv_nf"],data["cv_nf"],data["cv_nf"],data["cv_nf"])
    bounds = (.1, 2)

    x_3d = np.array([data["cv_nf"].values,
                     data["rho_nf"].values,
                     data["k_nf"].values,
                     data["Measured μ"].values])

    y = data["Loop FOM"]
    x = [.6, .4]
    fitParams, fitCovariances = curve_fit(fitFuncAlphaBeta, x_3d, y, x, method="trf", bounds=bounds)
    print(fitParams)
    # fit the T vs pump power data
    # popt, pcov = curve_fit(funcFOM, fluid_T_wall_agg, np.log(fluid_pumpPower_agg),
    #                        bounds=bounds)

    data["fitted"] = fitFuncAlphaBeta([data["cv_nf"].values,
                                       data["rho_nf"].values,
                                       data["k_nf"].values,
                                       data["Measured μ"].values], *fitParams)

    print(data)

if 1 == 0:
    print('___________________With CTE parameter')
    # plt.scatter(data["cv_nf"],data["cv_nf"],data["cv_nf"],data["cv_nf"])
    bounds = (.1, 30)

    x_3d = np.array([data["cv_nf"].values,
                     data["rho_nf"].values,
                     data["k_nf"].values,
                     data["Measured μ"].values,
                     data["CTE"].values
                     ])

    y = data["Loop FOM"]
    x = [4 / 3, 2, 2, 5 / 3, 2]
    fitParams, fitCovariances = curve_fit(fitFuncCTE, x_3d, y, x, method="trf", bounds=bounds)
    print(fitParams)
    # fit the T vs pump power data
    # popt, pcov = curve_fit(funcFOM, fluid_T_wall_agg, np.log(fluid_pumpPower_agg),
    #                        bounds=bounds)

    data["fitted"] = fitFuncCTE([data["cv_nf"].values,
                                 data["rho_nf"].values,
                                 data["k_nf"].values,
                                 data["Measured μ"].values, data["CTE"].values], *fitParams)

    print(data)

if 1 == 1:
    mixedConvection_vec = np.vectorize(mixedConvection)

    empty = np.ones(data['cv_nf'].size)
    result = mixedConvection(c=data['cv_nf'][0],
                    k=data['k_nf'][0],
                    mu=data['Measured μ'][0],
                    rho=data['rho_nf'][0],
                    deltaTwall2bulk=25,
                    heatflux=1.175*100**2,
                    d_h=.652/100,
                    CTE=data['CTE'][0])
    print(result)
    result = mixedConvection_vec(c=data['cv_nf'],
                    k=data['k_nf'],
                    mu=data['Measured μ'],
                    rho=data['rho_nf'],
                    deltaTwall2bulk=25,
                    heatflux=1.175*100**2,
                    d_h=.652/100,
                    CTE=data['CTE'])
    data["FOM"] = result


    print(data)
