

# quick test to get shear thninning behavior
# have viscoity vs shear rate
# shear rate
#   parallel plates --> v/h
#   tube --> 8v/d
#   set of tubes?
import scipy.io
import os
from os.path import join as pjoin  # tool to make paths
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, fsolve,root

def func_logistic(x, a, b, c,d):
    # use a logistic to model shear thinning as it's well behaved
    # return a /( np.exp(-b * (np.log(x)-c)) + 1)+d
    return (a-d)*(np.tanh(-b*(np.log(x)-c))+1)/2+d

def func_exp(x, a, b, c):
    return a * np.exp(-b * x) + c
dirf = '/Users/hanscastorp/Colloid-model/data/FluidData'
fnames = [f for f in os.listdir(dirf) if f.endswith('.mat')]
print(fnames)
datas = {}


# def func(x,a,b):
#     return [x[0] * np.cos(x[1]) - a,
#             x[1] * x[0] - x[1] - b]
#
# pars = (4,5,1,2,3,1)
# root = fsolve(func, [1, 1],pars)
# print(root)

# np.isclose(func(root), [0.0, 0.0])  # func(root) should be almost 0.0.

def func_velocity(x,pars):
    # x = [velocity, viscosity]
    a, b, c, d, e = pars
    # return [x[1]-a /( np.exp(-b * (np.log(x[0])-c)) + 1)+d,
    #         x[1]-np.divide(e,x[0])]
    return [x[1]-( (a-d)*(np.tanh(-b*(np.log(x[0])-c))+1)/2+d),
            x[1]-np.divide(e,x[0])]


diameter = 1
e_cst = 1
for f in fnames:
    mat = scipy.io.loadmat(pjoin(dirf,f))

    data = np.vstack([[row.flat[0] for row in line] for line in mat[list(mat.keys())[3]]])
    # shear rate 1/s, viscosity in Pa s
    x_dat = data[:,0]#np.multiply(data[:,0],data[:,1])
    x_dat = x_dat*diameter / 8 # velocity for tube
    y_dat = data[:,1]
    plt.scatter(x_dat,y_dat,label=f)
    popt, pcov = curve_fit(func_logistic, x_dat, y_dat,maxfev=4000)
    # get a good fitting function (one that is monotonic
    if 1==1:
        plt.plot(np.linspace(.01,1000,10000), func_logistic(np.linspace(.01,1000,10000), *popt), 'r-',
                 )

#     data['Velocity'] = (Q_out_fixed / (A_s * (T_max - T_out_fixed)) * (
#                 .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(
#             np.multiply((np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])) ** .4,
#                         (np.divide(data['rho_nf'], data['mu_nf']) * D_h) ** .6),
#             data['k_nf']) / D_h) ** -1) ** (5 / 3)
    plt.plot(x_dat,e_cst * x_dat**-1)
    popt = np.append(popt,e_cst)
    print(popt)
    # print(func_velocity([1,1],*popt),'a')
    root0 = fsolve(func_velocity, [50,.02],args = popt)
    root1 = root(func_velocity, [50,.02],args = popt)

    plt.axvline(x=root0[0])


        # func_velocity()
plt.xlabel('Velocity (m/s)')
plt.ylabel('Viscosity (Ps)')
plt.xscale('log')
plt.legend()
plt.ylim([0,.05])






# try out the solver to get conversion
# calc the velocity with initial guess of viscosity, then recalculate viscosity



plt.show()