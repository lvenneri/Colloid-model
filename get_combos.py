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
from colloid_combo_fun import *


def grab_data(file_path):
    xl = pd.ExcelFile(file_path)
    dfs = {sheet: xl.parse(sheet) for sheet in xl.sheet_names}

    nanoparticle_cand = grab_table(dfs, "nanoparticle")
    basefluid_cand = grab_table(dfs, "basefluid")
    problem_par = grab_table(dfs, "problem")

    return nanoparticle_cand, basefluid_cand, problem_par


def grab_table(data, sheet):
    table = Table.from_pandas(data[sheet])
    table.remove_row(0)
    return table


if __name__ == "__main__":
    # save_loc = pjoin('/Users/hanscastorp/nanofluids/data/', datetime.date.today().strftime("%Y%m%d"))
    nanoparticle_cand, basefluid_cand, problem_par = grab_data('data/data_1.xlsx')

    samples = 1
    par_space = [

        ['np_id', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]],
        # -1 for not used
        ['bf_id', [0, 1, 3, 2, 4]],  # range(len(basefluid_cand['id']))],
        ['p', np.linspace(1, 10, 2)],
        ['phi', np.linspace(0, .25, 20)],
        ['a_33', np.linspace(1, 100, 1)],
        ['R_bd', np.linspace(10 ** -8, 10 ** -8, 1)],

    ]

    par_extent = []
    for i in par_space:
        print([i[0], min(i[1]), max(i[1])])
        par_extent.append([i[0], min(i[1]), max(i[1]), len(i[1])])
    par_extent = np.vstack(par_extent)
    print(par_extent)

    parameter = []
    possible_vals = []

    # numTests = par_testing_fixed_fluid(pars=par_space, fluidData=basefluid_cand, npData=nanoparticle_cand)

    for par in par_space:
        parameter.append(par[0])
        possible_vals.append(par[1])
    par_tests = list(itertools.product(*possible_vals))
    par_tests = np.vstack(par_tests)
    print(parameter)

    # nanoparticle_cand[par]
    combos_data = pd.DataFrame(data=par_tests, columns=parameter)
    np_pars = ['np_id_s' ,"k_p", "R_bd", "rho_p", "cv_p" ]
    bf_pars = ['bf_id_s', "k_f", "cv_f", "rho_f", "mu_f"]
    for npp in np_pars:
        combos_data[npp] = nanoparticle_cand[npp][combos_data['np_id'].astype(int)]
    for bpp in bf_pars:
        combos_data[bpp] = basefluid_cand[bpp][combos_data['bf_id'].astype(int)]

    # density
    combos_data['rho_nf'] = density_nf(combos_data["phi"], combos_data["rho_p"], combos_data["rho_f"])

    # conductivity - def hc_EMTnan(k_p, k_f, phi, a_33, p, R_bd)
    combos_data['k_nf'] = np.multiply(combos_data["k_f"] , hc_EMTnan(combos_data["k_p"], combos_data["k_f"], combos_data["phi"],
                                          combos_data["a_33"], combos_data["p"], combos_data["R_bd"]))

    # heat capacity - def shc_nf(phi, cv_p, cv_f, rho_p, rho_f)
    combos_data['c_nf'] = shc_nf(combos_data["phi"], combos_data["cv_p"], combos_data["cv_f"],
                                 combos_data["rho_p"], combos_data["rho_f"])

    print(combos_data)

    if 1 == 0:
        print('Tests to run: ', len(par_tests))
        # now run through all the tests

        sim_stack = []
        for ii, test in enumerate(par_tests):
            kwargs = {}

            for p, par in enumerate(parameter):
                kwargs[par] = test[p]
            if kwargs['np_id'] >= 0:
                for p, par in enumerate(np_pars):
                    kwargs[par] = nanoparticle_cand[par][kwargs['np_id']]
            if kwargs['bf_id'] >= 0:
                for p, par in enumerate(bf_pars):
                    kwargs[par] = fluidData[par][kwargs['bf_id']]
            output = htc_basic(**kwargs)
            print(ii, nanoparticle_cand['id'][kwargs['np_id']])

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
            data_table['np_id'][this_row] = nanoparticle_cand['id'][kwargs['np_id']]
