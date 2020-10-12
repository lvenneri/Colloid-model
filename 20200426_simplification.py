import shutil
import struct
import scipy
from aLMODEL.utilities.utility_functions import *
from astropy.io import fits
from os.path import join as pjoin
import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter
import xlrd
import pandas as pd
import itertools
from astropy.table import Table, Column, MaskedColumn, vstack
from astropy.io import ascii
import datetime
from os.path import join as pjoin
from scipy.optimize import fsolve
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import matplotlib.cm as cm


def basic2d(data, x, y, save_loc, name='', norm=1, axess=1):
    """
    Function for plotting all three fluids simultaneusly
    :param data:
    :param x:
    :param y:
    :param save_loc:
    :param name:
    :param norm:
    :param axess:
    :return:
    """
    fluids_unq = np.unique(data['bf_id'])
    particles_unq = np.unique(data['np_id'])
    p_unq = np.unique(data['p'])

    colors = cm.tab20(np.linspace(0, 1, len(particles_unq)))

    # plt.figure(figsize=(13.0, 8.0), dpi=400)
    plt.figure(figsize=(13.0, 8.0))
    plt.subplots_adjust(wspace=0.0, hspace=0.0)

    if norm == 1:
        normer = data[data["phi"] == 0][y][0]
    else:
        normer = 1
    fluids_unq = fluids_unq[[0, 1, 3, 2]]
    for bfi, bf in enumerate(fluids_unq):

        plt.subplot(1, 4, bfi + 1)
        bf_data = data[data['bf_id'] == bf]
        for npi, nanop in enumerate(particles_unq):
            np_data = bf_data[bf_data['np_id'] == nanop]

            for pi, pppp in enumerate(p_unq):
                pp_data = np_data[np_data['p'] == pppp]
                if pi == 0:
                    plt.plot(pp_data[x], pp_data[y] / normer, lw=1, c=colors[npi], label=nanop.replace('_', ' '))
                else:
                    plt.plot(pp_data[x], pp_data[y] / normer, lw=1, c=colors[npi], label='_nolegend_')

        #             plt.scatter(np_data[x],np_data[y]/normer,s=2,label=nanop)

        if x == 'phi':
            plt.xlabel('$\phi$')
        else:
            plt.xlabel(x.replace('_', ' '))
        if bfi == 0:
            plt.ylabel(y.replace('_', ' '))

        if bfi == 2 or bfi == 1:
            plt.tick_params(which="both", axis='y', labelleft='False', labelright='False', right='False',
                            left='False', )
            plt.legend()

        if bfi == 3:
            plt.tick_params(which="both", axis='y', right='True', left='False', labelleft='False', labelright='True')
        plt.xlim([0, .35])
        if y == 'deltaT':
            plt.hlines(4,0,.35,colors='r',linestyles='dashed',lw=1)
            plt.text(4, .1, '$\Delta T$', rotation=0)
        if axess == 1:
            plt.ylim([np.nanpercentile(data[y]/normer,0),np.nanpercentile(data[y]/normer,100)])
        #         plt.ylim([.0001,.01])

        plt.yscale('log')
        ax = plt.gca()
        ax.set_yscale('log')
        plt.grid(which='minor', color='k', linestyle='-', alpha=0.2)
        # ax.yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))
        # ax.yaxis.set_major_formatter(
        #     ticker.FuncFormatter(lambda y, pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y), 0)))).format(y)))
        plt.title(bf)

    plt.savefig(pjoin(save_loc, name + '--' + x + '_vs_' + 'FOM' + '_22.png'))
    return

def basic2d_single(data, x, ys, save_loc, name='', norm=1, axess=1,fluid= 0):
    """
    Function for plotting many variables for single fluid
    :param data:
    :param x:
    :param y:
    :param save_loc:
    :param name:
    :param norm:
    :param axess:
    :return:
    """
    fluids_unq = np.unique(data['bf_id'])
    particles_unq = np.unique(data['np_id'])
    p_unq = np.unique(data['p'])

    colors = cm.tab20(np.linspace(0, 1, len(particles_unq)))

    # plt.figure(figsize=(13.0, 8.0), dpi=400)
    plt.figure(figsize=(8.0,13.0))
    plt.subplots_adjust(wspace=.1, hspace=.1)

    if norm == 1:
        normer = data[data["phi"] == 0][y][0]
    else:
        normer = 1
    fluids_unq = fluids_unq[[0]]

    for bfi, bf in enumerate(fluids_unq):

        for yi,y in enumerate(ys):
            plt.subplot(len(ys), 1, bfi + 1)

            bf_data = data[data['bf_id'] == bf]
            for npi, nanop in enumerate(particles_unq):
                np_data = bf_data[bf_data['np_id'] == nanop]

                for pi, pppp in enumerate(p_unq):
                    pp_data = np_data[np_data['p'] == pppp]
                    if pi == 0:
                        plt.plot(pp_data[x], pp_data[y] / normer, lw=1, c=colors[npi], label=nanop.replace('_', ' '))
                    else:
                        plt.plot(pp_data[x], pp_data[y] / normer, lw=1, c=colors[npi], label='_nolegend_')

            #             plt.scatter(np_data[x],np_data[y]/normer,s=2,label=nanop)

            if x == 'phi':
                plt.xlabel('$\phi$')
            else:
                plt.xlabel(x.replace('_', ' '))
            if bfi == 0:
                plt.ylabel(y.replace('_', ' '))

            if bfi == 2 or bfi == 1:
                plt.tick_params(which="both", axis='y', labelleft='False', labelright='False', right='False',
                                left='False', )
                plt.legend()

            if bfi == 3:
                plt.tick_params(which="both", axis='y', right='True', left='False', labelleft='False', labelright='True')
            plt.xlim([0, .35])
            if y == 'deltaT':
                plt.hlines(4,0,.35,colors='r',linestyles='dashed',lw=1)
                plt.text(4, .1, '$\Delta T$', rotation=0)
            if axess == 1:
                plt.ylim([np.nanpercentile(data[y]/normer,0),np.nanpercentile(data[y]/normer,100)])
            #         plt.ylim([.0001,.01])

            plt.yscale('log')
            ax = plt.gca()
            ax.set_yscale('log')
            plt.grid(which='minor', color='k', linestyle='-', alpha=0.2)
            # ax.yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))
            # ax.yaxis.set_major_formatter(
            #     ticker.FuncFormatter(lambda y, pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y), 0)))).format(y)))
            plt.title(bf)

    plt.savefig(pjoin(save_loc, name + '--' + x + '_vs_' + 'FOM' + '_22.png'))
    return



def compare_ToutfixedQfixed(T_out_fixed=20, T_max=70, Q_out_fixed=10, T_air=0, S_t=2.337 / 100, phi_compare=.1):
    """
    This computes FOM by fixing T_out (essentially limited by the air temp and then solving for everything else.
    :param Q_fixed:
    :param T_max:
    :param deltaT:
    :param T_air:
    :param S_t:
    :param phi_compare:
    :return:
    """
    # data_loc file produced from /basic_model191203.py

    data = ascii.read(pjoin(data_loc, 'data_extract.csv'), format='csv')
    phi = data['phi']
    data['mu_nf'] = 1 + 2.5 * phi*9 + 6.2 * (phi*9) ** 2.0
    L_t = 6.510 / 100
    D_t = 1.825 / 100
    D_h = 4*.8666 * S_t ** 2 / np.pi / D_t - D_t
    D_hp = D_h  # (4*(S_t ** 2-(D_t/4)**2*np.pi))/(4*np.pi*D_t+D_h)
    PoD = S_t / D_t
    A_c = S_t ** 2 - (D_t / 2) ** 2 * np.pi
    A_c = np.sqrt(3)/4*(S_t)**2-.5*np.pi*D_t**2/4
    A_s = np.pi * D_t * L_t/2
    data['T_out'] = T_out_fixed
    data['psi_c'] = .9090 + .0783 * PoD - .1283 * np.exp(-2.4 * (PoD - 1))
    mu_ratio = 1
    print(T_out_fixed)

    data['Velocity'] = (Q_out_fixed/(A_s*(T_max-T_out_fixed))*(.128 * mu_ratio ** .14 * data['psi_c'] *np.multiply(np.multiply((np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])) ** .4,
                                                   (np.divide(data['rho_nf'], data['mu_nf']) * D_h) ** .6),
                                                                   data['k_nf']) / D_h)**-1)**(5/3)


    print( data['Velocity'][0:3])

    deltaT = np.divide(Q_out_fixed , (data['c_nf'] * data['rho_nf'] * A_c *data['Velocity']))
    data['deltaT'] = deltaT

    Reynolds = np.divide(data['rho_nf'], data['mu_nf']) * data['Velocity'] * D_h
    Prandtl = np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])

    data['m_flowrate'] = data['rho_nf'] * A_c * data['Velocity']

    data['htc'] = .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(np.multiply(Prandtl ** .4, Reynolds ** .6),
                                                                       data['k_nf']) / D_h
    data['Q_out'] = data['c_nf'] * data['m_flowrate'] * (deltaT)

    data['Q_out_check'] = A_s*data['htc']*(T_max-T_out_fixed)
    data['W_out'] = 162 * (PoD - 1) ** .435 * (Reynolds) ** -1 * data['rho_nf'] * (L_t / (2 * D_hp)) * data['Velocity'] ** 2 * data[
        'm_flowrate'] * data['rho_nf'] ** -1
    data['Reynolds'] = Reynolds
    data['Prandtl'] = Prandtl

    data['T_in'] = data['T_out'] - deltaT
    logmeanT = deltaT * (np.log(np.divide(data['T_out'] - T_air, data['T_in'] - T_air))) ** -1

    data['FOM Q/W'] = np.divide(data['Q_out'], data['W_out'])
    data['FOM Simplified'] = np.divide(data['c_nf']**(4/3)*data['rho_nf']**(2)*data['k_nf']**(2), data['mu_nf']**(5/3))

    data['FOM Mouromtseff'] = data['rho_nf'] ** .8 * data['c_nf'] ** .33 * data['k_nf'] ** .67 / data['mu_nf'] ** .47
    data['Cmu : Ck'] = np.divide(np.divide(data['mu_nf'], data['mu_f']), np.divide(data['k_nf'], data['k_f']))
    data['T_outcheck'] = T_max - data['Q_out'] / (A_s * data['htc'])

    # look at phi = .10
    np.argmin(data['phi'] - phi_compare)

    # for each particle, grab highest performing

    name = 'T_out_fixed=' + str(T_out_fixed) + ', Qoutfixed='+str(Q_out_fixed) + ', T_wall=' + str(T_max) + ', S_t=' + str(
        S_t) + ', T_air=' + str(T_air)

    save_loc = pjoin(data_loc,'Fixed_ToutQout')
    create_dir_safe(save_loc)
    if 1==1:
        ascii.write(data, pjoin(save_loc, 'data_extract'+name+'.csv'), format='csv',
                    fast_writer=True, overwrite=True)
        basic2d(data, 'phi', 'FOM Q/W', save_loc=save_loc, name='FOM -- ' + name, norm=1)
        basic2d(data, 'phi', 'FOM Simplified', save_loc=save_loc, name='FOM simple -- ' + name, norm=1)
        basic2d(data, 'phi', 'c_nf', save_loc=save_loc, name='c_nf -- ' + name, norm=0,axess=1)
        basic2d(data, 'phi', 'mu_nf', save_loc=save_loc, name='mu_nf -- ' + name, norm=0,axess=1)
        basic2d(data, 'phi', 'k_nf', save_loc=save_loc, name='k_nf -- ' + name, norm=0,axess=1)
        basic2d(data, 'phi', 'rho_nf', save_loc=save_loc, name='rho_nf -- ' + name, norm=0,axess=1)
        basic2d(data, 'phi', 'Q_out', save_loc=save_loc, name='Q_out -- ' + name, norm=0,axess=1)
        basic2d(data, 'phi', 'W_out', save_loc=save_loc, name='W_out -- ' + name, norm=0,axess=1)
        basic2d(data, 'phi', 'Velocity', save_loc=save_loc, name='Velocity -- ' + name, norm=0,axess=1)
        basic2d(data, 'phi', 'Reynolds', save_loc=save_loc, name='Reynolds -- ' + name, norm=0,axess=1)
        basic2d(data, 'phi', 'deltaT', save_loc=save_loc, name='deltaT -- ' + name, norm=0,axess=1)
        basic2d(data, 'phi', 'm_flowrate', save_loc=save_loc, name='massflow -- ' + name, norm=0,axess=1)
        basic2d(data, 'phi', 'Cmu : Ck', save_loc=save_loc, name='CmuCk -- ' + name, norm=0,axess=1)

    return data




global data_loc
data_loc = pjoin('/Users/hanscastorp/Dropbox/SandBox2019/aLMODEL/data/200505_sweep_9')


if 1 == 1:  # GOOD -- FOM using T_out fixed, Q' fixed, deltaT loose, solve for W' only
    data = compare_ToutfixedQfixed(T_out_fixed=30, T_max=40, Q_out_fixed=5, T_air=0, S_t=(2.337) / 100)
