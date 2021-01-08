

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
from scipy.optimize import curve_fit, fsolve

def func_logistic(x, a, b, c,d):
    # use a logistic to model shear thinning as it's well behaved
    return a /( np.exp(-b * (np.log(x)-c)) + 1)+d

def func_exp(x, a, b, c):
    return a * np.exp(-b * x) + c
dirf = '/Users/hanscastorp/Colloid-model/data/junk'
fnames = [pjoin(dirf, f) for f in os.listdir(dirf) if f.endswith('.mat')]
print(fnames)
datas = {}


def func(x,a,b):
    return [x[0] * np.cos(x[1]) - a,
            x[1] * x[0] - x[1] - b]

pars = (4,5)
root = fsolve(func, [1, 1],pars)
print(root)

np.isclose(func(root), [0.0, 0.0])  # func(root) should be almost 0.0.

def func_velocity(x,a,b,c):

    return [x[0]**-1,x[0]**-1]

for f in fnames:
    mat = scipy.io.loadmat(f)

    data = np.vstack([[row.flat[0] for row in line] for line in mat[list(mat.keys())[3]]])
    # shear rate 1/s, viscosity in Pa s
    x_dat = data[:,0]#np.multiply(data[:,0],data[:,1])
    y_dat = data[:,1]
    plt.scatter(x_dat,y_dat,label=f)
    popt, pcov = curve_fit(func_logistic, x_dat, y_dat,maxfev=1000)
    popt
    # get a good fitting function (one that is monotonic
    if 1==1:
        plt.plot(np.linspace(.01,1000,1000), func_logistic(np.linspace(.01,1000,1000), *popt), 'r-',
                 label='fit'+f)

#     data['Velocity'] = (Q_out_fixed / (A_s * (T_max - T_out_fixed)) * (
#                 .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(
#             np.multiply((np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])) ** .4,
#                         (np.divide(data['rho_nf'], data['mu_nf']) * D_h) ** .6),
#             data['k_nf']) / D_h) ** -1) ** (5 / 3)


    root = fsolve(func_velocity, [1, 1],pars_v)

        # func_velocity()
plt.xlabel('Shear rate')
plt.ylabel('Viscosity')
plt.xscale('log')
plt.legend()






# try out the solver to get conversion
# calc the velocity with initial guess of viscosity, then recalculate viscosity



plt.show()