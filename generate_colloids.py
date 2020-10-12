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
from astropy.table import Table, Column, MaskedColumn, vstack
from astropy.io import ascii
import datetime
import time
import json
from scipy.optimize import fsolve


def density_nf(phi, rho_p, rho_f):
    """
    Simple volumetric average
    :param phi:
    :param rho_p:
    :param rho_f:
    :return: density of nano fluid
    """
    return np.multiply(phi, (rho_p)) + np.multiply((1 - phi), rho_f)


def sedimentation_rate(rho_p, rho_f, radius, mu):
    """
    Simnple sedimenation rate based on Stoke's law
    :param rho_p: particle density (kg/m^3)
    :param rho_f: fluid density (kg/m^3)
    :param radius: partcle radius (nm)
    :param mu: nanofluid viscosity (PaS)?
    :return: sedimenation rate in (mm/yr)
    """
    return np.divide(2 / 9 * np.multiply(rho_p - rho_f, (radius * 10 ** -9) ** 2) * 9.81, mu) * 1000 * 3600 * 24 * 365


def hc_EMTclassical(k_p, k_f, phi):
    """
    Classical model of heat conduction enhancement in composite

    Effective Medium Theory by Maxwell (1881) for well diseprsed, spherical, non-interacting particles with
    negligible thermal resistance at surface. This works for multiple properties like dielectric constant as well.
    Weird things happen when you have resonances.

    :param k_p: thermal conductivity of particle
    :param k_f: thermal conductivity of base fluid
    :param phi: particle volumetric fraction
    :return: enhancement as fraction of k of nano fluid / k of base fluid

    Example:
        print('Set1,Sample1 EMT Classical',hc_EMTclassical(35 ,  0.61,.01), ' -- From paper: 1.024')

    """
    return np.divide(k_p + 2 * k_f + 2 * np.multiply(phi, (k_p - k_f)),
                     k_p + 2 * k_f - 1 * np.multiply(phi, (k_p - k_f)))


def hc_EMTnan(k_p, k_f, phi, a_33, p, R_bd):
    """
    Nan model of heat conduction enhancement in composite

    Nan et al's enhancement to include geometry and surface resistance.
    https://aip.scitation.org/doi/10.1063/1.3245330
    :param k_p: thermal conductivity of particle
    :param k_f: thermal conductivity of base fluid
    :param phi: particle volumetric fraction
    :param a_33: long axis, ﻿prolate ellipsoids with principal axes a11=a22 < a33
    :param a_11: short axis, ﻿prolate ellipsoids with principal axes a11=a22 < a33
    :param p: a_33/a_11
    :param R_bd: the Kapitza interfacial thermal resistance
    :return: enhancement as fraction of k of nano fluid / k of base fluid

    Example
        print('Set1,Sample1 EMT Nan',hc_EMTnan(k_p=35, k_f=.61, phi=.01, a_33=80, p=8, R_bd=10**-8),
        ' -- From paper: 1.086')

    """
    a_11 = a_33 / p
    if p > 1:
        L_11 = np.divide(p ** 2, 2 * (p ** 2 - 1)) - np.divide(p, 2 * (p ** 2 - 1) ** 1.5) * np.arccosh(p)
        L_33 = 1 - 2 * L_11

        gamma = (2 + 1 / p) * R_bd * k_f / (a_11 / 2)

        k_11 = np.divide(k_p, 1 + gamma * L_11 * k_p / k_f)
        k_33 = np.divide(k_p, 1 + gamma * L_33 * k_p / k_f)
        B_11 = np.divide(k_11 - k_f, k_f + L_11 * (k_11 - k_f))
        B_33 = np.divide(k_33 - k_f, k_f + L_33 * (k_33 - k_f))

        enhancement_fraction = np.divide(3 + phi * (2 * B_11 * (1 - L_11) + B_33 * (1 - L_33)),
                                         3 - phi * (2 * B_11 * L_11 + B_33 * L_33))
    else:
        enhancement_fraction = hc_EMTclassical(k_p, k_f, phi)

    return enhancement_fraction


def shc_nf(phi, cv_p, cv_f, rho_p, rho_f):
    """
    Specific Heat Capacity for nanofluids based on volumetric weighting
    Method II in https: // dx.doi.org / 10.1155 / 2012 / 181079
    test:     print(shc_nf(.1, 1, 4.1814, 1, 1))

    :param phi: particle volumetric fraction
    :param cv_p: nano particle heat capacity
    :param cv_f: fluid heat capacity
    :param rho_p: particle density
    :param rho_f: fluid density
    :return: J/kg/K
    """
    return np.divide(phi * (cv_p * rho_p) + (1 - phi) * (cv_f * rho_f),
                     phi * (rho_p) + (1 - phi) * rho_f)


def v_nf(phi):
    """
    Einstein approx for spherical particles with extra terms from Verberg et al.
    https://journals.aps.org/pre/abstract/10.1103/PhysRevE.55.3143
    lower bound:
    https: // www.ar.ethz.ch / TMPPDF / 27189910425.024 / ApplRheol_20_44582.pdf
    upper bound:
    https://asmedigitalcollection.asme.org/heattransfer/article/128/3/240/477257/Convective-Transport-in-Nanofluids


    :param phi: volumetric fraction
    :return:  viscosity enhancement from Einsten, from Alumina model, from Titania model
    """
    phi = phi * 9
    return [1 + 2.5 * phi + 6.2 * phi ** 2.0, 1 + 39.11 * phi + 533.9 * phi ** 2,
            1 + 5.45 * phi + 108.2 * phi ** 2]





def viscosity_temp_dependence(temp, type='water'):
    """
    Scales mu from 298K to the desired temperature assuming the same behaviour as water or ethanol
    :param temp:
    :param type: assume water or ethanol temp dependence
    :return: scaling factor
    """
    if type == 'water':  # valid 273 - 643
        pars = [1.856 * 10 ** (-11), 4209, 0.04527, 3.376 * 10 ** -5]
    if type == 'ethanol':  # valid 168–516
        pars = [0.00201, 1614, 0.00618, -1.132 * 10 ** -5]
    mu_at298 = pars[0] * np.exp(pars[1] / 298 + pars[2] * 298 - pars[3] * 298 ** 2)
    return pars[0] * np.exp(pars[1] / temp + pars[2] * temp - pars[3] * temp ** 2) / mu_at298


def other_pars():
    # for calculating additional parameters relvant to solution
    return


def htc_basic(k_p, k_f, phi, a_33, p, R_bd, cv_p, cv_f, rho_p, rho_f, mu_f, np_id, bf_id):
    L_t = 6.510 / 100
    S_t = 2.337 / 100
    D_t = 1.825 / 100
    D_h = .8666 * S_t ** 2 / np.pi / D_t - D_t / 4
    D_hp = D_h * 1  # TODO

    PoD = S_t / D_t

    deltaT = 5
    T_max = 40
    T_air = 0
    A_c = .45 * S_t ** 2
    A_s = np.pi * D_t * L_t

    # density
    rho_nf = density_nf(phi, rho_p, rho_f)

    # conductivity - def hc_EMTnan(k_p, k_f, phi, a_33, p, R_bd)
    k_nf = k_f * hc_EMTnan(k_p, k_f, phi, a_33, p, R_bd)

    # heat capacity - def shc_nf(phi, cv_p, cv_f, rho_p, rho_f)
    c_nf = shc_nf(phi, cv_p, cv_f, rho_p, rho_f)

    # viscosity - def v_nf(phi)



    return { "rho_nf": rho_nf,"mu_f": mu_f ,
            "k_nf": k_nf, "c_nf": c_nf}


def par_testing_huge(pars):
    # create parameter testing bank - which is all the combination of parameters, this is very general - can be for
    # any system

    parameter = []
    possible_vals = []
    for par in pars:
        parameter.append(par[0])
        possible_vals.append(par[1])
    par_tests = list(itertools.product(*possible_vals))
    # now run through all the tests
    sim_stack = []
    for ii, test in enumerate(par_tests):
        kwargs = {"k_p": 1, "k_f": 1, "phi": 1, "a_33": 1, "p": 1, "R_bd": 1,
                  "cv_p": 1, "cv_f": 1, "rho_p": 1, "rho_f": 1, "mu_f": 1}

        for p, par in enumerate(parameter):
            kwargs[par] = test[p]

        output = htc_basic(**kwargs)
        # {**dict1, **dict2}

        data_table.add_row()
        this_row = len(data_table) - 1
        data_table['id'][this_row] = ii
        # place the data into the big data_table using a slow lookup
        for kk in kwargs:
            data_table[kk][this_row] = kwargs[kk]
        for kk in output:
            data_table[kk][this_row] = output[kk]

    return


def par_testing_fixed_fluid(pars, fluidData, npData):
    """
    create parameter testing bank - which is all the combination of parameters
    use fluid_data,np_data to define the basic characteristics of particular fluids and np
    the pars are the free parameters, which include fluid and np lookups

    :param pars:
    :param kwargBase: baseline kwargs for things don't change
    :param fluidData:
    :param npData:
    :return:
    """

    parameter = []
    possible_vals = []
    for par in pars:
        parameter.append(par[0])
        possible_vals.append(par[1])
    par_tests = list(itertools.product(*possible_vals))

    print('Tests to run: ', len(par_tests))
    # now run through all the tests
    np_pars = ["k_p", "R_bd", "rho_p", "cv_p", ]
    bf_pars = ["k_f", "cv_f", "rho_f", "mu_f"]
    sim_stack = []
    for ii, test in enumerate(par_tests):
        kwargs = {}

        for p, par in enumerate(parameter):
            kwargs[par] = test[p]
        if kwargs['np_id'] >= 0:
            for p, par in enumerate(np_pars):
                kwargs[par] = npData[par][kwargs['np_id']]
        if kwargs['bf_id'] >= 0:
            for p, par in enumerate(bf_pars):
                kwargs[par] = fluidData[par][kwargs['bf_id']]
        output = htc_basic(**kwargs)
        print(ii, npData['id'][kwargs['np_id']])

        data_table.add_row()
        this_row = len(data_table) - 1
        data_table['id'][this_row] = str(ii)
        # place the data into the big data_table using a slow lookup
        kadjust = [k for k in kwargs if k != 'bf_id' and k != 'np_id']
        for kk in kadjust:
            data_table[kk][this_row] = kwargs[kk]
        for kk in output:
            data_table[kk][this_row] = output[kk]

        data_table['bf_id'][this_row] = fluidData['id'][kwargs['bf_id']]
        data_table['np_id'][this_row] = npData['id'][kwargs['np_id']]

    return len(par_tests)


def grab_table(data, sheet):
    table = Table.from_pandas(data[sheet])
    table.remove_row(0)
    return table


def grab_data(file_path):
    xl = pd.ExcelFile(file_path)
    dfs = {sheet: xl.parse(sheet) for sheet in xl.sheet_names}
    data = pd.read_excel(file_path, index_col=0, comment='#')

    par_space = []
    for i, id in enumerate(dfs["parameter ranges"]['id']):
        if not np.isnan(dfs["parameter ranges"]['samples'][i]):
            par_space.append(
                [id, np.linspace(dfs["parameter ranges"]['min'][i],
                                 dfs["parameter ranges"]['max'][i],
                                 dfs["parameter ranges"]['samples'][i])])
        else:
            par_space.append(
                [id, [dfs["parameter ranges"]['single'][i]]])

    nanoparticle_cand = grab_table(dfs, "nanoparticle")
    basefluid_cand = grab_table(dfs, "basefluid")
    problem_par = grab_table(dfs, "problem")

    return par_space, nanoparticle_cand, basefluid_cand, problem_par


if __name__ == "__main__":
    # save_loc = pjoin('/Users/hanscastorp/nanofluids/data/', datetime.date.today().strftime("%Y%m%d"))
    save_loc = pjoin('/Users/hanscastorp/Dropbox/SandBox2019/aLMODEL/data', '200505_sweep_9')

    global data_table

    # for visualation see basic_model_vis.ipynb
    if 1 == 1:
        # just for given basefluids
        par_space, nanoparticle_cand, basefluid_cand, problem_par = grab_data('data/data_1.xlsx')

        create_dir(save_loc)
        data_table = Table(masked=True,
                           names=["id", "k_p", "k_f", "phi",
                                  "a_33", "p", "R_bd",
                                  "cv_p", "cv_f", "rho_p",
                                  "rho_f", "mu_f",
                                  "rho_nf", "k_nf", "c_nf",
                                   'np_id', 'bf_id',
                                  ],


                           dtype=['S', float, float, float,
                                  float, float, float,
                                  float, float, float,
                                  float, float,
                                  float, float, float,
                                  'S', 'S',
                                  ])

        if 1 == 1:  # testing
            samples = 1
            par_space = [
                # ['k_p', np.linspace(1, 2000, 3)],
                # ['cv_p', np.linspace(.1, 4.183, 3)],
                # ['rho_p', np.linspace(.05, 20, 3)],
                ['R_bd', np.linspace(10 ** -8, 10 ** -8, 1)],
                # ['np_id', [-1]], # -1 for not used
                ['np_id', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]],  # -1 for not used

                ['bf_id', [0, 1, 3, 2, 4]],  # range(len(basefluid_cand['id']))],
                ['p', np.linspace(1, 10, 2)],
                ['phi', np.linspace(0, .25, 20)],

                ['a_33', np.linspace(1, 100, 1)],


                # ['cv_f', np.linspace(.3, .7, samples)],
                # ['k_f', np.linspace(.3, .7, samples)],
                # ['rho_f', np.linspace(.3, .7, samples)],
                # ['mu_f', np.linspace(.3, .7, samples)],


            ]

            par_extent = []
            for i in par_space:
                print([i[0], min(i[1]), max(i[1])])
                par_extent.append([i[0], min(i[1]), max(i[1]), len(i[1])])
            par_extent = np.vstack(par_extent)
            ascii.write(Table(par_extent), pjoin(save_loc, 'sweepspace.csv'), format='csv',
                        fast_writer=True, overwrite=True)

            start = time.time()

            numTests = par_testing_fixed_fluid(pars=par_space, fluidData=basefluid_cand, npData=nanoparticle_cand)
            print((-start + time.time()) / 3600 / numTests * 100000)
            ascii.write(data_table, pjoin(save_loc, 'data_extract.csv'), format='csv',
                        fast_writer=True, overwrite=True)


    plt.show()

