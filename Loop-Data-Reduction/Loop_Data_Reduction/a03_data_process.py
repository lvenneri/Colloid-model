import os
from os.path import join as pjoin  # tool to make paths
import numpy as np
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import pathlib
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d, UnivariateSpline, CubicSpline
from os.path import join as pjoin  # tool to make paths
import streamlit as st
import itertools
import altair as alt
from astropy.stats import sigma_clip
from matplotlib import cm
from datetime import datetime
from viscosity_models import viscosity
from colloid_combo_fun import *
from pandas import ExcelWriter
import inspect


def save_xls(list_dfs, xls_path):
    with ExcelWriter(xls_path) as writer:
        for n, df in enumerate(list_dfs):
            df.to_excel(writer, 'sheet%s' % n)
        writer.save()


@st.cache
def grab_data(file_path):
    nanoparticle_cand = pd.read_excel(file_path, sheet_name="nanoparticle", index_col=0, header=[0, 1], comment='#',
                                      engine='openpyxl')
    basefluid_cand = pd.read_excel(file_path, sheet_name="basefluid", index_col=0, header=[0, 1], comment='#',
                                   engine='openpyxl')
    nanoparticle_cand['id'] = nanoparticle_cand.index
    basefluid_cand['id'] = basefluid_cand.index
    nanoparticle_cand['np_id_s'] = nanoparticle_cand['id'].copy()
    basefluid_cand['bf_id_s'] = basefluid_cand['id'].copy()

    basefluid_cand['mu_f^-1'] = basefluid_cand['mu_f'] ** -1

    return nanoparticle_cand, basefluid_cand


@st.cache
def grab_data2(file_path):
    return pd.read_csv(file_path)


def setup_text_plots(fontsize=12, usetex=True, style='default'):
    """
    e.g. setup_text_plots(fontsize=14, usetex=True,style='default')
    :param fontsize:
    :param usetex:
    :param style:
    :return: setup for nice plotting
    """
    import matplotlib
    matplotlib.style.use(style)
    matplotlib.rcParams['savefig.dpi'] = 200
    matplotlib.rc('legend', fontsize=fontsize, handlelength=3)
    matplotlib.rc('axes', titlesize=fontsize)
    matplotlib.rc('axes', titlesize=fontsize)
    matplotlib.rc('axes', labelsize=fontsize)
    matplotlib.rc('xtick', labelsize=fontsize)
    matplotlib.rc('ytick', labelsize=fontsize)
    matplotlib.rc('text', usetex=usetex)
    matplotlib.rc('font', size=fontsize, family='serif',
                  style='normal', variant='normal',
                  stretch='normal', weight='normal')
    matplotlib.rcParams['grid.color'] = 'k'
    matplotlib.rcParams['axes.grid'] = True
    matplotlib.rcParams['grid.linestyle'] = ':'
    matplotlib.rcParams['grid.linewidth'] = 0.5
    matplotlib.rcParams['figure.figsize'] = [8.0, 8.0]
    matplotlib.rcParams['errorbar.capsize'] = 3
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'
    matplotlib.rcParams['xtick.top'] = True
    matplotlib.rcParams['ytick.right'] = True
    matplotlib.rcParams['image.interpolation'] = 'nearest'
    matplotlib.rcParams['image.resample'] = False
    matplotlib.rcParams['figure.dpi'] = 100
    matplotlib.rcParams['mathtext.fontset'] = 'cm'
    matplotlib.rcParams['mathtext.rm'] = 'serif'
    matplotlib.rc('lines', linewidth=1.0, color='k')
    matplotlib.rc('legend', fancybox=False, shadow=False, framealpha=1, borderpad=0, frameon=True)

setup_text_plots()


def temp_tolerance(T):
    return 0.15 + 0.002 * T


def compare_ToutfixedQfixed(data, T_out_fixed, T_max, Q_out_fixed, S_t,
                            L_t, D_t, D_h, A_c, A_s, option_flow_type):
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

    # data = ascii.read(pjoin(data_loc, 'data_extract.csv'), format='csv')
    D_hp = D_h  # (4*(S_t ** 2-(D_t/4)**2*np.pi))/(4*np.pi*D_t+D_h)
    PoD = (D_t + S_t) / D_t

    data['psi_c'] = .9090 + .0783 * PoD - .1283 * np.exp(-2.4 * (PoD - 1))

    phi = data['phi'].copy()

    a = viscosity()

    data['T_out'] = T_out_fixed
    data['wt fraction'] = np.divide(np.multiply(data['phi'], data['rho_p']), data['rho_nf'])
    Prandtl = np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])

    mu_ratio = 1
    if option_flow_type == "Laminar":
        data['Velocity'] = (Q_out_fixed / (A_s * (T_max - T_out_fixed)) * (
                .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(
            np.multiply((np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])) ** .4,
                        (np.divide(data['rho_nf'], data['mu_nf']) * D_h) ** .6),
            data['k_nf']) / D_h) ** -1) ** (5 / 3)

        Reynolds = np.divide(data['rho_nf'], data['mu_nf']) * data['Velocity'] * D_h

        data['HTC'] = .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(np.multiply(Prandtl ** .4, Reynolds ** .6),
                                                                           data['k_nf']) / D_h

    elif option_flow_type == "Turbulent":
        data['Velocity'] = (Q_out_fixed / (A_s * (T_max - T_out_fixed)) * (
                .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(
            np.multiply((np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])) ** .4,
                        (np.divide(data['rho_nf'], data['mu_nf']) * D_h) ** .8),
            data['k_nf']) / D_h) ** -1) ** (5 / 4)

        Reynolds = np.divide(data['rho_nf'], data['mu_nf']) * data['Velocity'] * D_h

        data['HTC'] = .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(np.multiply(Prandtl ** .4, Reynolds ** .8),
                                                                           data['k_nf']) / D_h
    elif option_flow_type == "Auto":
        v_lam = (Q_out_fixed / (A_s * (T_max - T_out_fixed)) * (
                .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(
            np.multiply((np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])) ** .4,
                        (np.divide(data['rho_nf'], data['mu_nf']) * D_h) ** .6),
            data['k_nf']) / D_h) ** -1) ** (5 / 3)
        Reynolds_lam = np.multiply(np.divide(data['rho_nf'], data['mu_nf']), v_lam) * D_h
        htc_lam = .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(np.multiply(Prandtl ** .4, Reynolds_lam ** .6),
                                                                       data['k_nf']) / D_h
        v_turb = (Q_out_fixed / (A_s * (T_max - T_out_fixed)) * (
                .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(
            np.multiply((np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])) ** .4,
                        (np.divide(data['rho_nf'], data['mu_nf']) * D_h) ** .8),
            data['k_nf']) / D_h) ** -1) ** (5 / 4)

        Reynolds_turb = np.multiply(np.divide(data['rho_nf'], data['mu_nf']), v_turb) * D_h

        htc_turb = .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(np.multiply(Prandtl ** .4, Reynolds_turb ** .8),
                                                                        data['k_nf']) / D_h
        turbulent_idx = Reynolds_lam > 2000
        data['Velocity'] = v_lam.copy()
        data['Velocity'][turbulent_idx] = v_turb[turbulent_idx]
        data['HTC'] = htc_lam.copy()
        data['HTC'][turbulent_idx] = htc_turb[turbulent_idx]

        Reynolds = np.divide(data['rho_nf'], data['mu_nf']) * data['Velocity'] * D_h

    data['m_flowrate'] = data['rho_nf'] * A_c * data['Velocity']

    deltaT = np.divide(Q_out_fixed, (data['c_nf'] * data['rho_nf'] * A_c * data['Velocity']))

    data['deltaT'] = deltaT

    data['Q_out'] = data['c_nf'] * data['m_flowrate'] * (deltaT)

    data['Q_out_check'] = A_s * data['HTC'] * (T_max - T_out_fixed)
    data['W_out'] = 162 * (PoD - 1) ** .435 * (Reynolds) ** -1 * data['rho_nf'] * (L_t / (2 * D_hp)) * data[
        'Velocity'] ** 2 * data['m_flowrate'] * data['rho_nf'] ** -1
    data['Reynolds'] = Reynolds
    data['Prandtl'] = Prandtl

    data['T_in'] = data['T_out'] - deltaT

    data['FOM Standard'] = np.divide(data['Q_out'], data['W_out'])
    data['T_max'] = T_max

    return data

def data_proces(datapath):

    data_wide = pd.read_pickle(pjoin(datapath, 'd_wide.pkl'))
    data_calibration = pd.read_pickle(pjoin(datapath, 'd_cal.pkl'))
    data_calibration_flow = pd.read_pickle(pjoin(datapath, 'd_cal_flow.pkl'))
    data_measurement = pd.read_pickle(pjoin(datapath, 'd_meas.pkl'))
    positions_wall = np.load(pjoin(datapath, 'positions_wall.pkl.npy'))

    data_measurement = data_measurement.sort_values("Mass Flow Rate (kg/hr)")
    data_calibration_flow = data_calibration_flow.sort_values("Mass Flow Rate (kg/hr)")


    # load data from a02 script - csv and pkl files
    data_fluid = pd.read_csv(pjoin(datapath, 'colloid_data.csv'))
    data_fluid_span = data_fluid[0:0]  # empty array to hold the data

    fluid_bulkT = np.mean(data_measurement['T2 TS Out'])

    # data_fluid["mu_nf"][0] = data_fluid["eta"][0]*np.exp(data_fluid["A"][0]/(fluid_bulkT+273.15 - data_fluid["T0"][
    # 0]))
    data_fluid.loc[0,"mu_nf"] = data_fluid["eta"][0]*np.exp(data_fluid["A"][0]/(fluid_bulkT+273.15 - data_fluid["T0"][
    0]))

    if 1 == 1:
        # create prediction using the battery geometry data
        st.sidebar.subheader('Battery Thermal and Geometry')
        Q_out_fixed = st.sidebar.slider('Cell Heat Generation (W)', 1., 100.,
                                        float(np.median(data_measurement['Sapphire Power'])))
        T_max = st.sidebar.slider('Maximum Cell Temperature (C)', 30., 60., 55.)
        T_out_fixed = st.sidebar.slider('Fluid Outlet Temperature (C)', 20., T_max, 40.)
        option_geo = st.sidebar.selectbox('Geometry', ["Cylindrical Cell", "Square Channel"], 1)
        S_t = st.sidebar.slider('Cell Separation (mm)', 0., 6., 2.) / 1000
        L_t = st.sidebar.slider('Cell Length (cm)', 1., 10., 6.510) / 100
        D_t = st.sidebar.slider('Cell Diameter (cm)', 1., 10., 1.825) / 100
        option_flow_type = st.sidebar.selectbox('Laminar or Turbulent Flow', ["Laminar", "Turbulent", "Auto"], 0)
        Wetted_p = np.pi * D_t / 2
        A_c = np.sqrt(3) / 4 * (D_t + S_t) ** 2 - 0.5 * np.pi * (D_t / 2) ** 2
        A_s = np.pi * D_t * L_t / 2
        D_h = 4 * A_c / Wetted_p

        if option_geo == "Square Channel":
            Wetted_p = 4 * D_h
            A_c = D_h ** 2
            A_s = D_h * L_t

        st.subheader('Computed Geometry Parameters')
        "Wetted Perimeter (cm):", '{:.2f}'.format(Wetted_p * 100)
        "Hydraulic Diameter (cm):", '{:.2f}'.format(D_h * 100)
        "Cross Sectional Area (cm^2):", '{:.2f}'.format(A_c * 100 ** 2)
        "Surface Area (cm^2):", '{:.2f}'.format(A_s * 100 ** 2)
        samples = 100
        t_maxes = np.linspace(30, 160, samples)
        for i in range(samples):
            aaaa = compare_ToutfixedQfixed(data_fluid, T_out_fixed, t_maxes[i], Q_out_fixed, S_t,
                                           L_t, D_t, D_h, A_c, A_s, option_flow_type)
            data_fluid_span = data_fluid_span.append(aaaa, ignore_index=True)
    if 1 == 0:
        fig1, (ax12, ax13) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 3]})
        # fig1.subplots_adjust(left=.15, bottom=.15, right=.85, top=.85, wspace=0, hspace=0)
        ax12.scatter(data_fluid_span['m_flowrate'] * 3600, data_fluid_span['T_max'])
        ax13.scatter(data_fluid_span['m_flowrate'] * 3600, data_fluid_span['T_out'] - data_fluid_span['T_in'])
        st.pyplot(fig1)

    fig6, ((ax, ax1), (ax1a, ax1b)) = plt.subplots(2, 2, gridspec_kw={'height_ratios': [3, 3]})
    fig6.set_size_inches(9, 9)

    fig7, (ax2, ax3, ax3t) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [3, 3, 3]})
    fig7.subplots_adjust(left=.15, bottom=.15, right=.85, top=.85, wspace=0, hspace=0)
    # fig7.set_size_inches(5, 10)

    fig7a, (ax4_7a, ax5_7a) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 3]})
    fig7a.subplots_adjust(left=.15, bottom=.15, right=.85, top=.85, wspace=0, hspace=0)
    # fig7.set_size_inches(5, 10)

    deviation = []
    colors = cm.inferno(np.linspace(0, 1, positions_wall.shape[0]))
    warning = "Calibration okay"
    maxTWall = -1000
    minTWall = 1000

    # find calibration polynomial coefficients, and fit the data
    for i in range(positions_wall.shape[0]):
        x, y = data_calibration[positions_wall[i, 0]], data_calibration["T2 TS Out"]
        # by SVD
        N_poly = 6
        A = np.zeros((y.size, N_poly))
        B = y  # B = np.mean(x)
        for m in range(N_poly):
            A[:, m] = x ** m
        NUC_coeff = np.flipud(np.matmul(np.linalg.pinv(A), B))

        # correct all the IR data
        data_wide[positions_wall[i, 0] + '_corr'] = np.polyval(NUC_coeff, data_wide[positions_wall[i, 0]])
        data_calibration[positions_wall[i, 0] + '_corr'] = np.polyval(NUC_coeff, data_calibration[positions_wall[i, 0]])
        data_calibration_flow[positions_wall[i, 0] + '_corr'] = np.polyval(NUC_coeff,
                                                                           data_calibration_flow[positions_wall[i, 0]])
        data_measurement[positions_wall[i, 0] + '_corr'] = np.polyval(NUC_coeff, data_measurement[positions_wall[i, 0]])

        # correct the wall temp measurement for the flow change #TODO
        data_measurement[positions_wall[i, 0] + '_corr'] = data_measurement[
            positions_wall[i, 0] + '_corr']  # + (40-data_calibration_flow[positions_wall[i, 0]+'_corr'])

        # check that all the measurement IR data is within the calibration range
        minRAW = np.min(
            [np.min(data_measurement[positions_wall[i, 0]]), np.min(data_calibration[positions_wall[i, 0]])])
        maxRAW = np.max(
            [np.max(data_measurement[positions_wall[i, 0]]), np.max(data_calibration[positions_wall[i, 0]])])
        if np.min(data_measurement[positions_wall[i, 0]]) < np.min(data_calibration[positions_wall[i, 0]]):
            print('MEASUREMENT LESS THAN CALIBRATION')
            warning = "* Insufficient Calibration Data *"
        if np.max(data_measurement[positions_wall[i, 0]]) > np.max(data_calibration[positions_wall[i, 0]]):
            print('MEASUREMENT GREATER THAN CALIBRATION')
            warning = "* Insufficient Calibration Data *"

        rangeRAW = np.linspace(minRAW, maxRAW, 1000)
        maxTWall = max(maxTWall, np.max(data_calibration[positions_wall[i, 0] + '_corr']))
        minTWall = min(minTWall, np.min(data_calibration[positions_wall[i, 0] + '_corr']))
        ax.plot(y, y - data_calibration[positions_wall[i, 0] + '_corr'], lw=.5, label='Corrected')

        deviation.append(y - np.polyval(NUC_coeff, x))
        ax.set_title('IR LineOut Calibration')
        ax.set_xlabel('TS Outlet Temperature °C')
        ax.set_ylabel('T Residuals')

        ax1a.scatter(data_calibration[positions_wall[i, 0]], data_calibration[positions_wall[i, 0] + '_corr'], c='r',
                     s=5, label='Corrected Data')
        ax1a.plot(rangeRAW, np.polyval(NUC_coeff, rangeRAW), c=colors[i], label='Poly fit over range')
        ax1a.set_title('IR LineOut Calibration')
        ax1a.set_xlabel('Raw IR Values')
        ax1a.set_ylabel('T Wall °C')
        ax1a.scatter(data_measurement[positions_wall[i, 0]], data_measurement[positions_wall[i, 0] + '_corr'], c='k',
                     marker="*", s=5, label='Measurement Data')

        ax1b.scatter(y, data_calibration[positions_wall[i, 0]], c='r', marker="+", s=5, label='Calibration Data')
        ax1b.scatter(data_measurement["T2 TS Out"], data_measurement[positions_wall[i, 0]], c='b', marker="*", s=5,
                     label='Measurement Data')
        ax1b.set_title(warning)
        ax1b.set_xlabel('TS Outlet Temperature °C')
        ax1b.set_ylabel('Counts')

        ax1.scatter(data_calibration_flow["Mass Flow Rate (kg/hr)"],
                    data_calibration_flow[positions_wall[i, 0] + '_corr'], c=colors[i], s=5, label='_nolegend_')
        ax1.set_title('Thin Film Off - Changing Flow at MEASUREMENT SP')
        ax1.set_xlabel('Mass Flow (kg/hr)')
        ax1.set_ylabel('T Wall °C')

        ax2.scatter(data_measurement["Mass Flow Rate (kg/hr)"], data_measurement[positions_wall[i, 0] + '_corr'], s=5,
                    c=colors[i])
        ax2.set_title('Thin Film Measurement %3.3f W' % np.mean(data_measurement['Sapphire Power']))
        ax2.set_ylabel('T Wall °C')

    ax1.scatter(data_calibration_flow["Mass Flow Rate (kg/hr)"], data_calibration_flow["T2 TS Out"], s=5,
                label='T TS Out')
    ax1.scatter(data_calibration_flow["Mass Flow Rate (kg/hr)"], data_calibration_flow["T3 TS In"], s=5,
                label='T TS In')
    ax1.legend()

    ax3.scatter(data_measurement["Mass Flow Rate (kg/hr)"], data_measurement['$\Delta T_corrected$'], s=5, c='k',
                label='Corrected')
    ax3.scatter(data_measurement["Mass Flow Rate (kg/hr)"],
                data_measurement['T2 TS Out'] - data_measurement['T3 TS In'], s=5, c='k', marker="x",
                label='Uncorrected')
    ax3.plot(data_fluid_span['m_flowrate'] * 3600, data_fluid_span['T_out'] - data_fluid_span['T_in'], c='r',
             label='Predicted')
    ax3.set_xlim([0, np.max(data_measurement["Mass Flow Rate (kg/hr)"]) * 1.05])
    # ax3.set_ylim([0,2])
    ax3.set_ylabel('$T_{outlet} - T_{inlet}$')

    ax3t.scatter(data_measurement["Mass Flow Rate (kg/hr)"],
                 np.multiply(data_measurement['Mass Flow Rate (kg/hr)'] / 3600,
                             data_measurement["$\Delta T_corrected$"] * data_fluid["c_nf"][0]), s=5,
                 label='Thin Film Heat Transfer from Energy Conservation')
    ax3t.scatter(data_measurement["Mass Flow Rate (kg/hr)"],
                 data_measurement["Sapphire Power"], s=5, marker="x",
                 label='Thin Film Power Delivered')

    data_measurement["Thin Film Power Supply Power"] = data_measurement["Sapphire Voltage"] * data_measurement[
        "Sapphire Current"]
    ax3t.scatter(data_measurement["Mass Flow Rate (kg/hr)"],
                 data_measurement["Thin Film Power Supply Power"], s=5, marker="+",
                 label='Thin Film Power Delivered + Contact Loss')
    ax3t.set_ylabel('Power')

    ax2.plot(data_fluid_span['m_flowrate'] * 3600, data_fluid_span['T_max'], c='r', label='Predicted')
    ax2.set_xlim([0, np.max(data_measurement["Mass Flow Rate (kg/hr)"]) * 1.05])
    ax2.set_ylim([minTWall * .95, maxTWall * 1.05])

    ax3t.set_xlabel('Mass Flow (kg/hr)')
    ax3.legend()
    ax3t.legend()

    # A_s in meters (A_s / 100 ** 2)
    ax4_7a.scatter(data_measurement["Mass Flow Rate (kg/hr)"], data_measurement["$\dot{Q}_corrected$"], s=5, c='k',
                label='Corrected')
    ax4_7a.scatter(data_calibration_flow["Mass Flow Rate (kg/hr)"], data_calibration_flow["$\dot{Q}$"], s=5, c='r',
                label='Heat Loss')
    ax4_7a.axhline(y=np.mean(data_measurement['Sapphire Power']), linestyle="-.", lw=1,
                label='Heater Power: ' + str(np.mean(data_measurement['Sapphire Power'])) + ' W')
    ax4_7a.set_ylabel('$\dot{Q}$ W')
    ax4_7a.legend()

    poswall = [f + '_corr' for f in positions_wall[:, 0]]

    # Extract represtnative wall temperatures: max, min, median
    data_measurement["T_wall Outlet"] = np.max(data_measurement[poswall], axis=1)
    data_measurement["T_wall Inlet"] = np.min(data_measurement[poswall], axis=1)
    data_measurement["T_wall Midpoint"] = np.median(data_measurement[poswall], axis=1)

    # calculate the HTC and Nu, using the measured Q dump ($\dot{Q}_corrected$) or the power supply measurement (
    # "Sapphire Power")
    A_s_thinfilm = D_h * L_t  # surface area of thin film
    data_measurement['T Bulk'] = (data_measurement['T2 TS Out']+data_measurement['T3 TS In'])/2
    data_measurement["HTC"] = np.divide(data_measurement["Sapphire Power"],
                                        A_s_thinfilm * (data_measurement['T_wall Midpoint'] - data_measurement['T Bulk']))
    data_measurement["HTC_Qcorr"] = np.divide(data_measurement["$\dot{Q}_corrected$"],
                                        A_s_thinfilm * (data_measurement['T_wall Midpoint'] - data_measurement['T Bulk']))

    data_measurement["HTC_err"] = np.power(np.power(.05,2)+
                                  np.power((.5+.1)*(data_measurement['T_wall Midpoint'] - data_measurement['T Bulk'])**-1,2)
                                   ,.5)
    data_measurement["HTC_err"] = np.multiply(data_measurement["HTC_err"],data_measurement["HTC"])

    data_measurement["Nu"] = data_measurement["HTC"] * D_h / data_fluid["k_nf"][0]

    data_measurement["HTC_knownQ"] = np.divide(data_measurement["Sapphire Power"], A_s_thinfilm * (
            data_measurement['T_wall Midpoint'] - data_measurement['T Bulk']))
    data_measurement["Nu_knownQ"] = data_measurement["HTC_knownQ"] * D_h / data_fluid["k_nf"][0]

    data_measurement["velocity_calc"] = data_measurement["Mass Flow Rate (kg/hr)"] / 3600 / data_fluid["rho_nf"][
        0] / A_c

    # Estimate the viscosity and density at the fluid outlet

    # calc Reynolds using viscosity of the bulk at the outlet temeprature
    data_measurement["Re"] = data_measurement["velocity_calc"] * data_fluid["rho_nf"][0] * D_h / data_fluid["mu_nf"][0]

    Prandtl = data_fluid["mu_nf"][0] * data_fluid["c_nf"][0] / data_fluid["k_nf"][0]
    data_measurement["Pr"] = Prandtl
    visc_kinematic = (data_fluid["mu_nf"][0] / data_fluid["rho_nf"][0])

    data_measurement["Gr"] = data_fluid["CTE"][0] * 9.81 * (data_measurement["T_wall Midpoint"] - (
                data_measurement['T2 TS Out'] - data_measurement[
            'T3 TS In']) / 2) * D_h ** 3 / visc_kinematic ** 2  # Grashof
    data_measurement["Ra"] = np.multiply(data_measurement["Gr"], data_measurement["Pr"])  # Rayleigh
    data_measurement["Ri"] = np.divide(data_measurement["Gr"], data_measurement["Re"] ** 2)  # Richardson

    data_measurement["Nu_Mixed_F"] = .464 * np.power(data_measurement["Re"], .5) * Prandtl ** (1 / 3) / (
                1 + (0.0205 / Prandtl) ** (2 / 3)) ** .25
    data_measurement["Nu_Mixed_N"] = .563 * np.power(data_measurement["Ra"], 1 / 4) / (
                1 + (0.437 / Prandtl) ** (9 / 16)) ** (4 / 9)  # Rayleigh

    n_mixed = 3.0
    Prandtl
    n_mixed
    data_measurement["Nu_Mixed"] = np.power(
        data_measurement["Nu_Mixed_F"] ** n_mixed + data_measurement["Nu_Mixed_N"] ** n_mixed, 1 / n_mixed)

    PoD = (D_t + S_t) / D_t

    # calculate the pumping power using the calcualted pressure drop and friction factor
    data_measurement["Pumping Power (corr)[W]"] = 162 * (PoD - 1) ** .435 * (data_measurement["Re"]) ** -1 * data_fluid["rho_nf"][0] * (
            L_t / (2 * D_h)) \
                                  * data_measurement["velocity_calc"] ** 2 * data_measurement["Mass Flow Rate (" \
                                                                                              "kg/s)"] * \
                                  data_fluid["rho_nf"][0] ** -1

    # calucate the pupming power using the measured drop and measured density etc - updated with Rosemount diff
    # pressure sensor output
    data_measurement["Pumping Power (DP)[W]"] = np.multiply(data_measurement["DPTS [inH20]"], data_measurement["Mass Flow Rate (kg/s)"]) / \
                                  data_fluid["rho_nf"][0]*249.09



    # errors for DPTS: .2 inH20, flow .3 kg/hr = 0.000083333333333 kg/s, 0.01 for density
    data_measurement["Pumping Power (DP)[W]_err"] = np.sqrt(
        np.divide(.4, data_measurement["DPTS [inH20]"]) ** 2 + np.divide(0.000083333333333, data_measurement[
            "Mass Flow Rate (kg/s)"]) ** 2 + 0.01 ** 2)  # fractional
    data_measurement["Pumping Power (DP)[W]_err"] = np.multiply(data_measurement["Pumping Power (DP)[W]"], data_measurement[
        "Pumping Power (DP)[W]_err"])  # absolute by multiplying fractional and value
    # which estimate to use for future calcs?
    pumping_power_est_type = "Pumping Power (DP)[W]"

    # now fit a function to the T_wall_max or T_wall_mid vs Pumpin power, and interpolate at a T_wall_max of interest.
    # The function to fit comes from the FOM which relates 1/W' = f(T_max)
    def funcPumpingPower(x, a, b, c):
        """
        :param a:
        :param b:
        :param c: this is the outlet temperatue, can set to the outlet temp or allow fitting
        :param x: temperatures
        :return:
        """
        # return a * np.power(x - b,c)+e*(f*x-d)**-1
        # return np.exp(-a*(x-30) - b)+c
        return (np.power(a, -x + b) + c / 100)

    def funcLine2(x, b, c):
        return b * x + c
    def funcLine(x, a, b, c):
        return a * x ** 2 + b * x + c

    def funcLineExp(x, a, b, c):
        return np.exp(a * x ** 2 + b * x + c)

    def funcLineExpDerivative(x, a, b, c):
        return np.multiply((2*a*x+b),np.exp(a * x ** 2 + b * x + c))

    def funcPupmingPowerEr(x, b, c, d, a):
        return b * x ** .5 + c * x ** 2 + d * x + a

    ydata = data_measurement[pumping_power_est_type]
    use_function_fit = 1
    if 1 == use_function_fit:
        bounds = ([1.1, 0, 0], [4, 100, 100])
        bounds = (-np.inf, np.inf)
        # bounds = ([], [])
        if 1 == 0:
            poptMax, pcovMax = curve_fit(funcPumpingPower, data_measurement["T_wall Outlet"], ydata,
                                         bounds=bounds)
            poptMin, pcovMin = curve_fit(funcPumpingPower, data_measurement["T_wall Inlet"], ydata,
                                         bounds=bounds)
            poptMedian, pcovMedian = curve_fit(funcPumpingPower, data_measurement["T_wall Midpoint"], ydata,
                                             bounds=bounds)
        else:
            poptMax, pcovMax = curve_fit(funcLine, data_measurement["T_wall Outlet"], np.log(ydata),
                                         bounds=bounds)
            poptMin, pcovMin = curve_fit(funcLine, data_measurement["T_wall Inlet"], np.log(ydata),
                                         bounds=bounds)
            poptMedian, pcovMedian = curve_fit(funcLine, data_measurement["T_wall Midpoint"], np.log(ydata),
                                             bounds=bounds)
        popts = [poptMin,poptMedian,poptMax]
        'Fit Values'
        print('Fit Values')
        print(poptMedian)

    else:
        if 1 == 1:
            f_max = interp1d(data_measurement["T_wall Outlet"], ydata, kind='linear', fill_value="extrapolate")
            f_min = interp1d(data_measurement["T_wall Inlet"], ydata, kind='linear', fill_value="extrapolate")
            f_median = interp1d(data_measurement["T_wall Midpoint"], ydata, kind='linear', fill_value="extrapolate")

    ax5_7a.scatter(data_measurement["Mass Flow Rate (kg/hr)"], data_measurement["HTC"], s=5, c='k', label='Corrected')
    ax5_7a.set_ylabel('HTC $W/cm^2K$')
    ax5_7a.set_xlabel('Mass Flow (kg/hr)')
    ax5_7a.legend()

    deviation = np.vstack(deviation)
    st.text('IR Camera: Average Deviation from Fit %3.5f' % np.mean(abs(deviation)))

    fig9, (ax5) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
    fig8, (ax6) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})

    positions_wall_corr = [f + '_corr' for f in positions_wall[:, 0]]
    imaa = ax5.imshow(data_wide[positions_wall_corr].to_numpy().transpose(), aspect='auto', cmap="inferno")
    # cbar = plt.colorbar(imaa, orientation='vertical', aspect=50, fraction=.05, shrink=1.0, pad=0.05,
    #                     label='Counts')
    # imaa = ax6.imshow(data_wide[positions_wall[:,0]].to_numpy().transpose(), aspect='auto', cmap="inferno")

    im2plot = data_wide[positions_wall_corr].to_numpy().transpose()
    for f in range(im2plot.shape[1]):
        im2plot[:, f] -= np.min(im2plot[:, f])

    imaa2 = ax6.imshow(im2plot, vmin=0, vmax=np.max(abs(im2plot)), aspect='auto', cmap="gist_heat")  # RdBu_r") #
    cbar = plt.colorbar(imaa2, orientation='horizontal', aspect=50, fraction=.05, shrink=1.0, pad=0.05,
                        label='Temp Deviation')
    ax5.set_yticks(np.arange(len(positions_wall[:, 1])))
    ax5.set_yticklabels(list(positions_wall[:, 1]))
    # ax5.tick_params(bottom='off')

    ax5.set_xticklabels([])
    ax5.set_xlabel("Time")
    ax5.set_ylabel("Position (cm)")

    ax5.set_xticks([])
    ax5.grid(None)

    ax6.set_yticks(np.arange(len(positions_wall[:, 1][::-1])))
    ax6.set_yticklabels(list(positions_wall[:, 1][::-1]))
    ax6.tick_params(bottom='off')
    ax6.set_xticklabels([])
    ax6.set_xticks([])
    ax6.grid(None)

    fig10, (ax10) = plt.subplots()
    spanRe = np.linspace(10, 1000, 200)
    ax10.scatter(data_measurement["Re"], data_measurement["Nu"], label="From Energy Balance")

    poptReNu, pcov = curve_fit(funcLine2, np.log(data_measurement["Re"]), np.log(data_measurement["Nu"]))
    "Reynolds Fit", poptReNu
    ax10.plot(spanRe, np.exp(funcLine2(np.log(spanRe), *poptReNu)), label="fit")
    ax10.plot(spanRe, np.power(spanRe, .6) * .128 * Prandtl ** .4, label="Laminar")
    ax10.plot(spanRe, np.power(spanRe, .8) * .023 * Prandtl ** .4, label="Turbulent")

    # ax10.scatter(data_measurement["Re"], data_measurement["Nu_Mixed"], label="Mixed Prediction")
    # ax10.scatter(data_me
    # {Ra}^{1/4} + Laminar$")
    ax10.set_ylabel('Nu')
    ax10.set_xlabel('Re')
    ax10.set_yscale('log')
    ax10.set_xscale('log')

    ax10.scatter(data_measurement["Re"], data_measurement["Nu_knownQ"], label="From Power Supply Output")
    poptReNu_knownQ, pcov = curve_fit(funcLine, np.log(data_measurement["Re"]), np.log(data_measurement["Nu_knownQ"]))
    "Reynolds Fit _knownQ", poptReNu_knownQ
    ax10.plot(spanRe, np.exp(funcLine(np.log(spanRe), *poptReNu_knownQ)), label="fit")
    ax10.legend()

    fig10A, (ax10A) = plt.subplots()
    ax10A.scatter(data_measurement["Mass Flow Rate (kg/s)"], data_measurement["Ri"])
    ax10A.set_xlabel('Mass Flow Rate (kg/s)')
    ax10A.set_ylabel('Richardson')
    ax10A.set_yscale('log')

    fig11, (ax12) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
    fig11_13, (ax13) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
    ax12.scatter(data_measurement["Mass Flow Rate (kg/s)"], data_measurement["Pumping Power (corr)[W]"],
                 label="by Calculated Pressure Drop")
    ax12.scatter(data_measurement["Mass Flow Rate (kg/s)"], data_measurement["Pumping Power (DP)[W]"],
                 label="by Measured Pressure Drop")

    ax12.set_ylabel('$Pumping Power \dot{W}$')
    ax12.set_xlabel('Mass Flow Rate (kg/s)')
    ax12.legend()

    # ax12.set_yscale('log')
    # ax12.set_xscale('log')
    ax12.set_title('Pumping Power')
    ax13.errorbar(data_measurement["T_wall Outlet"], data_measurement[pumping_power_est_type], xerr=.5,
                  yerr=data_measurement[pumping_power_est_type + '_err'], ls='none', label="$T_{wall, outlet}$")
    ax13.errorbar(data_measurement["T_wall Inlet"], data_measurement[pumping_power_est_type], xerr=.5,
                  yerr=data_measurement[pumping_power_est_type + '_err'], ls='none', label="$T_{wall, inlet}$")
    ax13.errorbar(data_measurement["T_wall Midpoint"], data_measurement[pumping_power_est_type], xerr=.5,
                  yerr=data_measurement[pumping_power_est_type + '_err'], ls='none', label="$T_{wall, midpoint}$")
    xdata_span = np.linspace(np.min(data_measurement["T_wall Inlet"]), np.max(data_measurement["T_wall Outlet"]), 300)

    if 1 == use_function_fit:
        ax13.plot(xdata_span, funcLineExp(xdata_span, *poptMax), label="Fit Outlet")
        ax13.plot(xdata_span, funcLineExp(xdata_span, *poptMin), label="Fit Inlet")
        ax13.plot(xdata_span, funcLineExp(xdata_span, *poptMedian), label="Fit Midpoint")
    else:
        ax13.plot(xdata_span, f_max(xdata_span), label="Fit Outlet")
        ax13.plot(xdata_span, f_min(xdata_span), label="Fit Inlet")
        ax13.plot(xdata_span, f_median(xdata_span), label="Fit Midpoint")
    ax13.set_title('Pumping Power vs Wall Temp')
    ax13.set_xlabel('$T_{wall}$ [°C]')

    # interpolated pupming power at target temp
    targetTemp = 55
    if 1 == use_function_fit:
        interpolatePumpingPowerMax = funcLineExp(targetTemp, *poptMax)
        interpolatePumpingPowerMin = funcLineExp(targetTemp, *poptMin)
        interpolatePumpingPowerMediean = funcLineExp(targetTemp, *poptMedian)
    else:
        interpolatePumpingPowerMax = f_max(targetTemp)
        interpolatePumpingPowerMin = f_min(targetTemp)
        interpolatePumpingPowerMediean = f_median(targetTemp)
    print("interpolatePumpingPowerMax", interpolatePumpingPowerMax)
    print("interpolatePumpingPowerMin", interpolatePumpingPowerMin)
    print("interpolatePumpingPowerMediean", interpolatePumpingPowerMediean)

    pumpingPower3Points = np.asarray([interpolatePumpingPowerMin, interpolatePumpingPowerMediean, interpolatePumpingPowerMax])
    IRCAMERA_ERR = .5  # deg C

    # pumping power error fit function
    popt_err, pcov_err = curve_fit(funcPupmingPowerEr, data_measurement["Pumping Power (DP)[W]"], data_measurement["Pumping Power (DP)[W]_err"],
                                   bounds=bounds)
    pumpingPower3Points_err = []
    for ff in range(pumpingPower3Points.size):

        pumpingPower3Points_err.append(np.power((funcPupmingPowerEr(funcLineExp(targetTemp, *popts[ff]), *popt_err) ** 2
                                     + IRCAMERA_ERR * funcLineExpDerivative(targetTemp,*popts[ff]) ** 2),
                                    .5))


    np.save(pjoin(datapath, 'pumpingPower3Points.npy'), pumpingPower3Points)
    np.save(pjoin(datapath, 'pumpingPower3Points_err.npy'), np.hstack(pumpingPower3Points_err) )

    ax13.axvline(x=targetTemp, color='r', linestyle='--')
    ax13.axhline(y=interpolatePumpingPowerMax, color='r', linestyle='--')
    ax13.axhline(y=interpolatePumpingPowerMin, color='r', linestyle='--')
    ax13.axhline(y=interpolatePumpingPowerMediean, color='r', linestyle='--')
    ax13.legend()
    # ax13.set_yscale('log')
    # ax13.set_xscale('log')

    for i in range(positions_wall.shape[0]):
        data_measurement[positions_wall[i, 0] + "HTC"] = np.divide(data_measurement["Sapphire Power"],
                                                                   A_s_thinfilm * (
                                                                           data_measurement[
                                                                               positions_wall[i, 0] + '_corr'] -
                                                                           data_measurement['T Bulk']))
        data_measurement[positions_wall[i, 0] + "Nu"] = data_measurement[positions_wall[i, 0] + "HTC"] * D_h / \
                                                        data_fluid["k_nf"][0]

    fig12, (ax12A) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
    colors_mflow = cm.inferno(np.linspace(0, 1, data_measurement["Mass Flow Rate (kg/s)"].shape[0]))

    for i, mass_flow_rate in enumerate(data_measurement["Mass Flow Rate (kg/s)"]):
        d_local = data_measurement[data_measurement["Mass Flow Rate (kg/s)"] == mass_flow_rate]
        colnames = [f + "Nu" for f in positions_wall[:, 0]]
        # print(len(list(positions_wall[:, 1])))
        # print(data_measurement[colnames].shape)
        ax12A.plot(np.asarray(list(positions_wall[:, 1][::-1])).astype(float), d_local[colnames].values.ravel(),
                   c=colors_mflow[i],
                   label=f'{mass_flow_rate * 3600:.2f} kg/hr')
    data_measurement["Nu Midpoint"] = np.median(data_measurement[colnames], axis=1)
    data_measurement["Nu Outlet"] = np.min(data_measurement[colnames], axis=1)
    data_measurement["Nu Inlet"] = np.max(data_measurement[colnames], axis=1)
    ax12A.set_xlabel('Wall Position (cm)')
    ax12A.set_title("Developing Nusselt")
    ax12A.set_ylabel('Nusselt')
    ax12A.legend()
    st.pyplot(fig12)
    fig12.savefig(pjoin(datapath, '1_Nusselt_wall.png'))


    fig12HTC, (ax12AHTC) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
    colors_mflow = cm.inferno(np.linspace(0, 1, data_measurement["Mass Flow Rate (kg/s)"].shape[0]))

    for i, mass_flow_rate in enumerate(data_measurement["Mass Flow Rate (kg/s)"]):
        d_local = data_measurement[data_measurement["Mass Flow Rate (kg/s)"] == mass_flow_rate]
        colnames = [f + "HTC" for f in positions_wall[:, 0]]
        # print(len(list(positions_wall[:, 1])))
        # print(data_measurement[colnames].shape)
        ax12AHTC.plot(np.asarray(list(positions_wall[:, 1][::-1])).astype(float), d_local[colnames].values.ravel(),
                   c=colors_mflow[i],
                   label=f'{mass_flow_rate * 3600:.1f} kg/hr, Re:{np.mean(d_local["Re"]):.1f}')

    data_measurement["HTC Midpoint"] = np.median(data_measurement[colnames], axis=1)
    data_measurement["HTC Outlet"] = np.min(data_measurement[colnames], axis=1)
    data_measurement["HTC Inlet"] = np.max(data_measurement[colnames], axis=1)


    ax12AHTC.set_xlabel('Wall Position (cm)')
    ax12AHTC.set_title("Developing HTC")
    ax12AHTC.set_ylabel('HTC $[W/m^2K]$')
    ax12AHTC.legend()
    st.pyplot(fig12)
    fig12HTC.savefig(pjoin(datapath, '1_HTC_wall.png'))



    st.pyplot(fig6)
    st.pyplot(fig7)
    st.pyplot(fig7a)
    st.pyplot(fig8)
    st.pyplot(fig10)
    st.pyplot(fig10A)
    st.pyplot(fig11)
    st.pyplot(fig11_13)
    fig6.savefig(pjoin(datapath, '1_calibration.png'))
    fig7.savefig(pjoin(datapath, '1_measurement.png'))
    fig7a.savefig(pjoin(datapath, '1_HTC_v_pumping.png'))
    fig8.savefig(pjoin(datapath, '1_lineout_time.png'))
    fig10.savefig(pjoin(datapath, '1_Re_Nu.png'))
    ax10.set_yscale('linear')
    ax10.set_xscale('linear')
    fig10.savefig(pjoin(datapath, '1_Re_Nu_linear.png'))
    fig10A.savefig(pjoin(datapath, '1_massFlow_Ri.png'))
    fig11.savefig(pjoin(datapath, '1_Mflow_PumpingPower_knownQ.png'))
    fig11_13.savefig(pjoin(datapath, '1_T_max_PumpingPower_knownQ_linear.png'))
    ax13.set_yscale('log')
    # ax13.set_xscale('log')
    fig11_13.savefig(pjoin(datapath, '1_T_max_PumpingPower_knownQ.png'))

    save_xls([data_measurement, data_calibration, data_calibration_flow], pjoin(datapath, 'data_excl_a03.xlsx'))
    data_wide.to_pickle(pjoin(datapath, 'd_wide_a03.pkl'))
    data_calibration.to_pickle(pjoin(datapath, 'd_cal_a03.pkl'))
    data_calibration_flow.to_pickle(pjoin(datapath, 'd_cal_flow_a03.pkl'))
    data_measurement.to_pickle(pjoin(datapath, 'd_meas.pkl_a03'))


if __name__ == "__main__":
    datapath = "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210325_water"
    datapath = "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210801_isoparaffin_5W"
    datapath ="/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210726_isoparaffin_5W"
    # datapath ="/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210729_isoparaffin_5W"

    data_proces(datapath)
    # plt.show()
