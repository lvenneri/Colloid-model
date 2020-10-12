import shutil
import struct
import scipy
from utilities.utility_functions import *
from astropy.io import fits
from os.path import join as pjoin
import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter
import xlrd
import pandas as pd
import itertools
from astropy.table import Table, Column, MaskedColumn, vstack, hstack
from astropy.io import ascii
import datetime
from os.path import join as pjoin
from scipy.optimize import fsolve
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import matplotlib.cm as cm
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def compare_Qfixed(T_outs=np.linspace(-20, 30, 100), T_max=70, Q_out_fixed=10, T_air=0, S_t=2.337 / 100):
    # pick best performing nanofluid at phi = .1, and scan the Tout fixed space to get W,Tout,delta coordinates
    # data_loc file produced from /basic_model191203.py

    # data_fluids = ascii.read('/Users/hanscastorp/nanofluids/data/data_1_bestperformers.csv', format='csv')
    data_fluids = ascii.read('/Users/hanscastorp/nanofluids/data/200408_sweep_4.3/simplified_selection.csv', format='csv')

    L_t = 6.510 / 100
    D_t = 1.825 / 100
    D_h = 4*.8666 * S_t ** 2 / np.pi / D_t - D_t
    D_hp = D_h  # (4*(S_t ** 2-(D_t/4)**2*np.pi))/(4*np.pi*D_t+D_h)
    PoD = S_t / D_t
    A_c = S_t ** 2 - (D_t / 2) ** 2 * np.pi
    A_s = np.pi * D_t * L_t
    mu_ratio = 1

    # each fluid requries its own data set
    # so array is fluids x 3 for T_out, W, deltaT

    fluid_names = [data_fluids['np_id'][f] + ':' + data_fluids['bf_id'][f]+ ':' + str(data_fluids['p'][f]) for f in range(data_fluids['np_id'].size)]

    data_WTdT = np.zeros((len(fluid_names), 5, T_outs.size))  # fluids, (W,T,dT), Tout
    for fi, fluid in enumerate(fluid_names):

        data = vstack([data_fluids[fi] for ll in range(T_outs.size)])
        if 1 == 1:
            data['T_out'] = T_outs

            data['psi_c'] = .9090 + .0783 * PoD - .1283 * np.exp(-2.4 * (PoD - 1))

            data['Velocity'] = (Q_out_fixed / (A_s * (T_max - data['T_out'])) * (
                        .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(
                    np.multiply((np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])) ** .4,
                                (np.divide(data['rho_nf'], data['mu_nf']) * D_h) ** .6),
                    data['k_nf']) / D_h) ** -1) ** (5 / 3)

            deltaT = np.divide(Q_out_fixed, (data['c_nf'] * data['rho_nf'] * A_c * data['Velocity']))
            data['deltaT'] = deltaT

            Reynolds = np.divide(data['rho_nf'], data['mu_nf']) * data['Velocity'] * D_h
            Prandtl = np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])

            data['m_flowrate'] = data['rho_nf'] * A_c * data['Velocity']

            data['htc'] = .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(
                np.multiply(Prandtl ** .4, Reynolds ** .6),
                data['k_nf']) / D_h
            data['Q_out'] = data['c_nf'] * data['m_flowrate'] * (deltaT)

            data['Q_out_check'] = A_s * data['htc'] * (T_max - data['T_out'])
            data['W_out'] = 162 * (PoD - 1) ** .435 * (Reynolds) ** -1 * data['rho_nf'] * (L_t / (2 * D_hp)) * data[
                'Velocity'] ** 2 * data[
                                'm_flowrate'] * data['rho_nf'] ** -1
            data['Reynolds'] = Reynolds
            data['Prandtl'] = Prandtl

            data['T_in'] = data['T_out'] - deltaT

            data['FOM Q/W'] = np.divide(data['Q_out'], data['W_out'])
            data['FOM Simplified'] = np.divide(data['c_nf'] ** (4 / 3) * data['rho_nf'] ** (2) * data['k_nf'] ** (2),
                                               data['mu_nf'] ** (5 / 3))

            data['FOM Mouromtseff'] = data['rho_nf'] ** .8 * data['c_nf'] ** .33 * data['k_nf'] ** .67 / data[
                'mu_nf'] ** .47
            data['Cmu : Ck'] = np.divide(np.divide(data['mu_nf'], data['mu_f']), np.divide(data['k_nf'], data['k_f']))
            data['T_outcheck'] = T_max - data['Q_out'] / (A_s * data['htc'])

        # fluids, (W,T,dT), Tout
        data_WTdT[fi, 0, :] = np.array(data['W_out'])
        data_WTdT[fi, 1, :] = np.array(data['T_out'])
        data_WTdT[fi, 2, :] = np.array(data['deltaT'])
        data_WTdT[fi, 3, :] = np.array(data['Reynolds'])
        data_WTdT[fi, 4, :] = np.array(data['m_flowrate'])
    return fluid_names, data_WTdT


Q_out_fixed = 10

import pickle

if 1==1: # data prep
    fluid_names, dataFluidX3 = compare_Qfixed(T_outs=np.power(10,np.linspace(-1,np.log10(38), 100)), T_max=40, Q_out_fixed=Q_out_fixed, T_air=0,
                                              S_t=2.337 / 100)
    with open(str(Q_out_fixed)+'objs.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([fluid_names, dataFluidX3], f)
else:
    with open(str(Q_out_fixed)+'objs.pkl','rb') as f:  # Python 3: open(..., 'rb')
        fluid_names, dataFluidX3 = pickle.load(f)

# plot the data in 3d space to show each fluid is different

if 1==0:
    fig = plt.figure(figsize=(7, 7))
    plt.subplots_adjust(wspace=.5, hspace=.5)
    mpl.rcParams['legend.fontsize'] = 10

    print('plotting...')

    ax = fig.add_subplot(1,1, 1, projection='3d')

    basefluids_unq = np.unique([f.split(':')[1] for f in fluid_names])
    colors = cm.Set1(range(len(basefluids_unq)))

    colordict = {basefluids_unq[f]: colors[f] for f in range(len(basefluids_unq))}
    if 1==1:
        # dataFluidX3[:, 1, :] = np.log10(dataFluidX3[:, 1, :])
        dataFluidX3[:, 0, :] = np.log10(dataFluidX3[:, 0, :]*7000)
        dataFluidX3[:, 2, :] = np.log10(dataFluidX3[:, 2, :])
    else:
        dataFluidX3[:, 0, :] = dataFluidX3[:, 0, :]*7000

    for fi, fluid in enumerate(fluid_names):
        bff = fluid.split(':')[1]
        if fluid.find('-') >= 0:
            ax.plot(dataFluidX3[fi, 0, :] , dataFluidX3[fi, 1, :], dataFluidX3[fi, 2, :], label=fluid.split(':')[1],
                    lw=1, linestyle='--', c=colordict[bff])
            # ax.text(dataFluidX3[fi,0,-1]*7000, dataFluidX3[fi,1,-1], dataFluidX3[fi,2,-1], fluid.split(':')[1], 'y',fontsize=8)
        else:
            if fluid.split(':')[0] == 'ZnO':
                # ax.text(dataFluidX3[fi, 0,-1] , dataFluidX3[fi, 1, -1], dataFluidX3[fi, 2, -1], fluid,'x',fontsize=8)
                if fluid.split(':')[2] == '10':
                    plt.plot(dataFluidX3[fi, 0, :] , dataFluidX3[fi, 1, :], dataFluidX3[fi, 2, :], label=fluid.split(':')[0] +' in '+ bff,
                             lw=2, alpha=1, c=colordict[bff])
                # else:
                #     plt.plot(dataFluidX3[fi, 0, :] , dataFluidX3[fi, 1, :], dataFluidX3[fi, 2, :], label=fluid + ' sphere',
                #              lw=1, linestyle='dashdot', alpha=1, c=colordict[bff])


            elif fluid.split(':')[2] == '10':
                plt.plot(dataFluidX3[fi, 0, :] , dataFluidX3[fi, 1, :], dataFluidX3[fi, 2, :], label='_nolegend_',
                         lw=1, alpha=.2, c=colordict[bff])



    ax.legend(fontsize=6,loc='center right')
    if 1==0:
        ax.set_xlabel('$log(\dot{W}) $ for 7000 cells (W)')
        ax.set_ylabel('$log(T_{out} (C))$')
        ax.set_zlabel('$log(\Delta T_{max})$')
    else:
        ax.set_xlabel('$log(\dot{W} $ for 7000 cells (W))')
        ax.set_ylabel('$T_{out} (C)$')
        ax.set_zlabel('$log(\Delta T_{max})$')
    ax.set_title('$\dot{Q} = ' + str(Q_out_fixed) + '$ and $T_{wall} = 40$°C')

if 1 == 1:
    fig = plt.figure(figsize=(7, 7))
    plt.subplots_adjust(wspace=0, hspace=0)
    mpl.rcParams['legend.fontsize'] = 10

    print('plotting...')

    basefluids_unq = np.unique([f.split(':')[1] for f in fluid_names])
    colors = cm.Set1(range(len(basefluids_unq)))

    colordict = {basefluids_unq[f]: colors[f] for f in range(len(basefluids_unq))}

    dataFluidX3[:, 0, :] = dataFluidX3[:, 0, :] * 7000
    axesnames = ['$\dot{W} $ for 7000 cells (W)','$T_{out} (C)$','$\Delta T_{max}$']
    for fi, fluid in enumerate(fluid_names):
        bff = fluid.split(':')[1]
        against = [[0,2],[2,1],[0,1]]
        for ai, aa in enumerate(against):
            if fluid.find('-') >= 0:
                plt.subplot(2,2,ai+2)
                plt.plot(dataFluidX3[fi, aa[0], :], dataFluidX3[fi, aa[1], :], label=fluid.split(':')[1],
                        lw=1, linestyle='--', c=colordict[bff])
                # ax.text(dataFluidX3[fi,0,-1]*7000, dataFluidX3[fi,1,-1], dataFluidX3[fi,2,-1], fluid.split(':')[1],
                # 'y',fontsize=8)
            else:
                if fluid.split(':')[0] == 'ZnO':
                    # ax.text(dataFluidX3[fi, 0,-1] , dataFluidX3[fi, 1, -1], dataFluidX3[fi, 2, -1], fluid,'x',fontsize=8)
                    plt.subplot(2, 2, ai + 2)
                    if fluid.split(':')[2] == '10':
                        plt.plot(dataFluidX3[fi, aa[0], :], dataFluidX3[fi, aa[1], :], label=fluid + ' rod',
                                 lw=2, alpha=1, c=colordict[bff])
                    else:
                        plt.plot(dataFluidX3[fi, aa[0], :], dataFluidX3[fi, aa[1], :], label=fluid + ' sphere',
                                 lw=1, linestyle='dashdot',alpha=1, c=colordict[bff])


                elif fluid.split(':')[2] == '10':
                    plt.subplot(2, 2, ai + 2)
                    plt.plot(dataFluidX3[fi, aa[0], :], dataFluidX3[fi, aa[1], :], label='_nolegend_',
                            lw=1, alpha=.2, c=colordict[bff])
            if not axesnames[aa[0]] == '$T_{out} (C)$':
                plt.xscale('log')
            if aa[1] == 2:
                plt.yscale('log')
            else:
                plt.yscale('linear')

            plt.tick_params( axis='y', right='True', left='False', labelleft='False', labelright='True')

            if ai ==1:
                plt.tick_params(which="both", axis='y', right='False', left='True', labelleft='True', labelright='False')

            plt.ylabel(axesnames[aa[1]])

            plt.xlabel(axesnames[aa[0]])
            plt.grid(which='major', color='k', linestyle='-', alpha=0.5)
            plt.grid(which='minor', color='k', linestyle='-', alpha=0.2)


    plt.legend(fontsize=6, bbox_to_anchor=(-.3, 1.9))
    # fig.legend( fontsize=6,
    #            loc="center right",  # Position of legend
    #            borderaxespad=0.1,  # Small spacing around legend box
    #            title="Legend Title"  # Title for the legend
    #            )


plt.show()




