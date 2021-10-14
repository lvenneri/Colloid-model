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
import streamlit as st
import itertools
import altair as alt
from astropy.stats import sigma_clip
from matplotlib import cm
from datetime import datetime
from pandas import ExcelWriter
import matplotlib.colors

"""

Collects data from all the specified datapaths and does analysis. Uses the folder name as YYYYMMDD_FLUID to parse. 
The plotting and analysis is interspersed and disorganized. See a06_multifluidB for easier to understand analysis.


"""


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

datapaths1 = ["/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210726_isoparaffin_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210727_isoparaffin_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210728_isoparaffin_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210729_isoparaffin_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210730_isoparaffin_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210801_isoparaffin_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210802_isoparaffin_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210804_isoparaffin_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210817_isoparaffin_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210824_isoparaffin_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210826_isoparaffin_5W",

              ]

datapaths2 = ["/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210805_OS622339H_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210806_OS622339H_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210808_OS622339H_5W", ]

datapaths3 = ["/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210809_OS625475_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210810_OS625475_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210811_OS625475_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210812_OS625475_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210813_OS625475_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210814_OS625475_5W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210815_OS625475_5W", ]
datapaths4 = ["/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210924_Water_10W",
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210925_Water_10W",
              ]
datapaths5 = [ "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210930_Water_5W", ]
fig4, (ax4) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]}, figsize=(6, 6))


def forward(x):
    return 5 / (x - 30) / (4.25 / 10 ** 4)


def inverse(x):
    return 5 / x / (4.25 / 10 ** 4) + 30


def analysis(name="median", idx_ref=1):
    sets = [datapaths1, datapaths2, datapaths3,datapaths4,datapaths5]
    compiled_pumpingPower3Points = []
    fluids = []

    pumping_power_est_type = "Pumping Power (DP)[W]"
    if 1 == 1:  # plot pumping power versus wall temperature
        colors_a = plt.cm.tab10(range(len(sets)))
        print(colors_a)
        # colors_a = ['r', 'b', 'g','k','orange']
        linestyle = {"Midpoint": '-', "Outlet": "--", "Inlet": "-.-"}

        fig1, (ax1) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})

        fig8, (ax8) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
        fig6, (ax6, ax7) = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw={'width_ratios': [3, 3]})
        fig6.subplots_adjust(left=.08, bottom=.2, right=.92, top=.88, wspace=0, hspace=0)
        fig4.subplots_adjust(bottom=.2)
        fig1.subplots_adjust(bottom=.2)

        fig5, ax5 = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]}, figsize=(5, 5))
        fig5.subplots_adjust(left=.08, bottom=.2, right=.92, top=.88, wspace=0, hspace=0)

        fig2, (ax2_1, ax2_2, ax2_3) = plt.subplots(1, 3, figsize=(12, 4), gridspec_kw={'width_ratios': [3, 3, 3]})
        fig2.subplots_adjust(left=.08, bottom=.2, right=.92, top=.88, wspace=0, hspace=0)

        sets_len = [len(f) for f in sets]
        print(sets_len)
        fig3, (ax3_0) = plt.subplots(1, len(sets), figsize=(12, 4), gridspec_kw={'width_ratios': sets_len})
        fig3.subplots_adjust(left=.08, bottom=.2, right=.92, top=.88)
        fig3B, (ax3B) = plt.subplots(1, 1, figsize=(12, 4), gridspec_kw={'width_ratios': [1]})
        fig3B.subplots_adjust(left=.08, bottom=.2, right=.92, top=.88)
        c = []
        interpolatePumpingPower_agg = []
        interpolatePumpingPowerErr_agg = []
        fluids_agg = []

        popt_fluids = []
        for iset, datapaths in enumerate(sets):
            dataSamples = len(sets)
            fluid_T_wall_agg = []
            fluid_HTC_wall_agg = []
            fluid_HTC_wallErr_agg = []
            fluid_Re_agg = []
            fluid_volFlow_agg = []
            fluid_pumpPower_agg = []
            fluid_pumpPowerErr_agg = []
            fluid_heatingAggNoSaph = []
            fluid_Richard = []
            fluid_heatingAggSaph = []
            dates = []
            pump_powerMidpoint = []
            pump_powerMidpointErr = []
            fluid_HTC_wall = []
            for i, datapath in enumerate(datapaths):
                print(datapath)
                test_name = datapath.split('/Lubrizol Loop/')[1]
                date = datapath.split('/Lubrizol Loop/')[1].split('_')[0]
                dates.append(date)
                fluid_name = datapath.split('_')[1]
                fluid_power = datapath.split('_')[2]
                fluids.append(fluid_name + ' ' + str(i))
                print(fluid_name)
                pumpingPower3Points = np.load(pjoin(datapath, "pumpingPower3Points.npy"))
                pumpingPower3Points_err = np.load(pjoin(datapath, "pumpingPower3Points_err.npy"))
                data_measurement = pd.read_pickle(pjoin(datapath, 'd_meas.pkl_a03'))
                data_calibration_flow = pd.read_pickle(pjoin(datapath, 'd_cal_flow_a03.pkl'))
                positions_wall = np.load(pjoin(datapath, 'positions_wall.pkl.npy'))

                colnames = [f + "HTC" for f in positions_wall[:, 0]]

                # wall position, Re, HTC or Nu
                for i, mass_flow_rate in enumerate(data_measurement["Re"]):
                    d_local = data_measurement[data_measurement["Mass Flow Rate (kg/s)"] == mass_flow_rate]
                    np.asarray(list(positions_wall[:, 1][::-1])).astype(float), d_local[colnames].values.ravel()

                fluid_HTC_wall = []
                fluid_T_wall_agg.append(data_measurement["T_wall " + name].values)
                fluid_HTC_wall_agg.append(data_measurement["HTC " + name].values)
                fluid_HTC_wallErr_agg.append(data_measurement["HTC_err"].values)
                fluid_Re_agg.append(data_measurement["Re"].values)
                fluid_heatingAggNoSaph.append(data_measurement["PreHeat Power"].values)
                fluid_Richard.append(data_measurement["Ri"].values)
                fluid_heatingAggSaph.append(data_calibration_flow["PreHeat Power"].values)
                fluid_pumpPower_agg.append(data_measurement[pumping_power_est_type].values)
                fluid_volFlow_agg.append(data_measurement["Volume Flow Rate (L/min)"].values)
                fluid_pumpPowerErr_agg.append(data_measurement['Pumping Power (DP)[W]' + '_err'].values)
                compiled_pumpingPower3Points.append(
                    [date, fluid_name, float(pumpingPower3Points[0]), float(pumpingPower3Points[1]),
                     float(pumpingPower3Points[2]), float(pumpingPower3Points_err[1])])
                pump_powerMidpoint.append(pumpingPower3Points[1])
                pump_powerMidpointErr.append(pumpingPower3Points_err[1])
            dates = np.hstack(dates)
            pump_powerMidpoint = np.hstack(pump_powerMidpoint)
            fluid_heatingAggNoSaph = np.hstack(fluid_heatingAggNoSaph)
            fluid_Richard = np.hstack(fluid_Richard)
            fluid_heatingAggSaph = np.hstack(fluid_heatingAggSaph)
            pump_powerMidpointErr = np.hstack(pump_powerMidpointErr)
            fluid_T_wall_agg = np.hstack(fluid_T_wall_agg)
            fluid_HTC_wall_agg = np.hstack(fluid_HTC_wall_agg)
            fluid_HTC_wallErr_agg = np.hstack(fluid_HTC_wallErr_agg)
            fluid_Re_agg = np.hstack(fluid_Re_agg)
            fluid_volFlow_agg = np.hstack(fluid_pumpPower_agg)
            fluid_pumpPower_agg = np.hstack(fluid_pumpPower_agg)
            fluid_pumpPowerErr_agg = np.hstack(fluid_pumpPowerErr_agg)

            # plot for looking at the midpoint FOM for each data capture.
            if 1 == 1 and name == "Midpoint":
                x_pos = np.arange(len(dates))

                FOMs = 5 / pump_powerMidpoint  # 5W
                ax3 = ax3_0[iset]
                ax3.bar(x_pos, FOMs, yerr=np.divide(pump_powerMidpointErr, pump_powerMidpoint) * FOMs, align='center',
                        ecolor='blue',
                        capsize=10, color="grey")

                ax3_0[0].set_ylabel('FOM at ' + name)
                ax3.set_xticks(x_pos)
                ax3.set_xticklabels(dates, rotation=45)
                ax3.set_title(fluid_name)
                ax3.yaxis.grid(True)
                ax3.xaxis.grid(False)

            # interpolate to all the data collected

            def funcLine(x, a, b, c):

                return a * x ** 2 + b * x + c

            def funcLinear(x, a, b, c):

                return b * np.log(x) + c

            def funcLineExp(x, a, b, c):
                return np.exp(a * x ** 2 + b * x + c)

            def funcLineExpDerivative(x, a, b, c):
                return np.multiply((2 * a * x + b), np.exp(a * x ** 2 + b * x + c))

            def funcPupmingPowerEr(x, b, c, d, a):
                return b * x ** .5 + c * x ** 2 + d * x + a

            def funcInterpError(xdata_span, popt, popt_err, IRCAMERA_ERR, samples):
                return np.power((funcPupmingPowerEr(funcLineExp(xdata_span, *popt),
                                                    *popt_err) ** 2 + IRCAMERA_ERR * funcLineExpDerivative(xdata_span,
                                                                                                           *popt) ** 2),
                                .5) / np.sqrt(samples)

            bounds = (-np.inf, np.inf)

            # fit the T vs pump power data
            popt, pcov = curve_fit(funcLine, fluid_T_wall_agg, np.log(fluid_pumpPower_agg),
                                   bounds=bounds)
            # fit the pump power vs pump power error data
            popt_err, pcov_err = curve_fit(funcPupmingPowerEr, fluid_pumpPower_agg, fluid_pumpPowerErr_agg,
                                           bounds=bounds)
            popt_fluids.append(popt)
            xdata_span = np.linspace(np.min(fluid_T_wall_agg), np.max(fluid_T_wall_agg), 300)
            y_data = funcLineExp(xdata_span, *popt)
            ax1.plot(xdata_span, y_data, c= colors_a[iset],label=fluid_name + ' ' +fluid_power+ ' Fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
            ax1.errorbar(fluid_T_wall_agg, fluid_pumpPower_agg, xerr=.5,
                         yerr=fluid_pumpPowerErr_agg, ls='none', alpha=1, lw=1, label='_nolegend_',c= colors_a[iset])
            secax1 = ax1.secondary_xaxis('bottom', functions=(forward, inverse))
            secax1.set_frame_on(True)
            secax1.patch.set_visible(False)
            secax1.xaxis.set_label_position('bottom')
            secax1.spines['bottom'].set_position(('outward', 40))
            secax1.set_xlabel('HTC $[W/m^2K]$')

            ydata_heat = fluid_Richard
            # popt_preheat, pcov_preheat = curve_fit(funcLinear, fluid_pumpPower_agg,ydata_heat)
            xdata_pspan = np.linspace(np.min(fluid_pumpPower_agg), np.max(fluid_pumpPower_agg), 300)
            # ax6.plot(xdata_pspan, funcLinear(xdata_pspan,*popt_preheat), c=colors_a[iset],label=fluid_name)
            ax6.scatter(fluid_T_wall_agg, ydata_heat, s=15, c=colors_a[iset],
                        label=fluid_name)
            # label=fluid_name + " Mean: %5.2f , std: %5.2f" % tuple([np.mean(ydata_heat),np.std(ydata_heat)]))

            ax7.scatter(fluid_pumpPower_agg, ydata_heat, s=15, c=colors_a[iset],
                        label=fluid_name)
            ax7.set_xlabel('Pumping Power $\dot{W}$ [W]')
            # ax7.set_ylabel('Richardson Number')
            ax7.set_xscale('log')
            # ax7.set_title('Richardson Number')
            ax7.legend(fontsize=10)

            IRCAMERA_ERR = .5  # deg C
            interpErrorTotal = funcInterpError(xdata_span, popt, popt_err, IRCAMERA_ERR, dataSamples)
            ax1.fill_between(xdata_span, y_data - interpErrorTotal, y_data + interpErrorTotal, facecolor='lightgrey',
                             alpha=.3)

            if name == "Midpoint":
                dataframe_fluid = pd.DataFrame(
                    {"Wall Temp": xdata_span, "Pumping Power": y_data, "Pumping Power Error": interpErrorTotal})
                dataframe_fluid.to_csv('/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/04_Compiled_data_temp_pumpingpower' + fluid_name + '.csv')

            ax2_1.errorbar(fluid_pumpPower_agg, fluid_HTC_wall_agg, xerr=fluid_pumpPowerErr_agg,
                           yerr=fluid_HTC_wallErr_agg, ls='none', alpha=1, lw=.5, label=fluid_name)
            ax2_2.errorbar(fluid_volFlow_agg, fluid_HTC_wall_agg, xerr=0,
                           yerr=fluid_HTC_wallErr_agg, ls='none', alpha=1, lw=.5, label=fluid_name)
            ax2_3.errorbar(fluid_Re_agg, fluid_HTC_wall_agg, xerr=0,
                           yerr=fluid_HTC_wallErr_agg, ls='none', alpha=1, lw=.5, label=fluid_name)

            # interpolate HTC vs Reynolds
            def funcHTC(x, a, b):
                return np.exp(a * np.log(x) + b)

            popt_ReHTC, pcov = curve_fit(funcHTC, fluid_Re_agg, fluid_HTC_wall_agg,
                                         # bounds=([0,-np.inf], [np.inf,np.inf])
                                         )
            ax8.errorbar(fluid_Re_agg, fluid_HTC_wall_agg, xerr=0,
                         yerr=fluid_HTC_wallErr_agg, ls='none', c=colors_a[iset], alpha=1, lw=.5, label=fluid_name)
            xrange_re = np.linspace(1, 1000, 1000)
            print('_____________________ FLUID HTC vs Re fit ____ ', fluid_name)
            print(popt_ReHTC)

            ax8.plot(fluid_Re_agg, funcHTC(fluid_Re_agg, *popt_ReHTC), c=colors_a[iset],
                     label=fluid_name + ' Fit: a=%5.3f, b=%5.3f' % tuple(popt_ReHTC))
            targetReynolds = 200
            interpHTC = funcHTC(targetReynolds, *popt_ReHTC)
            print(f"HTC at ", targetReynolds, ": ", interpHTC)

            targetTemp = 55

            interpolatePumpingPower = funcLineExp(targetTemp, *popt)
            interpolatePumpingPowerErr = funcInterpError(targetTemp, popt, popt_err, IRCAMERA_ERR, dataSamples)
            interpolatePumpingPower_agg.append(interpolatePumpingPower)
            interpolatePumpingPowerErr_agg.append(interpolatePumpingPowerErr)
            fluids_agg.append(fluid_name)
            ax1.axhline(y=interpolatePumpingPower, color='r', linestyle='--')

        uniqueFluids = np.unique(fluids_agg)
        colors = {}
        for iic, ic in enumerate(uniqueFluids):
            colors[ic] = plt.cm.tab10(iic)

        labels = list(colors.keys())

        ax1.axvline(x=targetTemp, color='r', linestyle='--')

        interpolatePumpingPower_agg = np.hstack(interpolatePumpingPower_agg)
        interpolatePumpingPowerErr_agg = np.hstack(interpolatePumpingPowerErr_agg)
        interpolatePumpingPowerErr_agg = np.divide(interpolatePumpingPowerErr_agg, interpolatePumpingPower_agg)
        FOM = interpolatePumpingPower_agg ** -1
        FOM_normed = FOM / FOM[idx_ref]
        FOM_normed_err = np.multiply(FOM_normed, interpolatePumpingPowerErr_agg)
        i = 0
        left, right = ax1.get_xlim()
        # xdata_span = np.linspace(45, 65 , 300)

        for datapaths in sets:
            # ax1.text(x=60,y=interpolatePumpingPower_agg[i],s=fluids_agg[i]+ " FOM: %10.1E" %FOM_normed[i])

            ax1.text(x=60, y=interpolatePumpingPower_agg[i],
                     s=fluids_agg[i] + " FOM: %10.1f $\pm$%10.1f" % (FOM_normed[i], FOM_normed_err[i]))
            ax1.text(x=right + (right - left) * .01, y=interpolatePumpingPower_agg[i], c="r",
                     s="%10.2E W" % interpolatePumpingPower_agg[i], fontsize=10)

            i += 1

        if 1 == 1:
            x_pos = np.arange(len(fluids_agg))
            p1 = ax5.bar(x_pos, np.round(FOM_normed, 1), yerr=FOM_normed_err, align='center',
                         ecolor='blue',
                         capsize=10, facecolor="white", edgecolor="black", fill=True)

            ax5.set_ylabel('FOM')
            ax5.set_xlabel('Fluid')
            ax5.set_xticks(x_pos)
            ax5.set_xticklabels(fluids_agg, rotation=45)
            ax5.set_title("Midpoint FOM Norm to Basefluid")
            ax5.yaxis.grid(True)
            ax5.xaxis.grid(False)
            ax5.bar_label(p1, label_type='center', padding=-5)
            ax5.bar_label(p1, labels=['±%.1f' % e for e in FOM_normed_err],
                          padding=2, color='b')

        i = 0
        if name != "Inlet":
            for datapaths in sets:
                scale_factor = np.divide(funcLineExp(xdata_span, *popt_fluids[i]) ** -1,
                                         funcLineExp(xdata_span, *popt_fluids[idx_ref]) ** -1)
                ax4.plot(xdata_span, scale_factor, c=colors_a[i], linestyle=linestyle[name],
                         label=fluids_agg[i] + ' at Wall ' + name)

                i += 1
            ax4.set_ylabel('Normalized FOM ')
            ax4.set_xlabel('$T_{wall}$ [°C]')
            # ax4.set_yscale('log')
            ax4.set_title('FOM Sensitivity Normaalzing to ' + fluids_agg[idx_ref])

            secax = ax4.secondary_xaxis('bottom', functions=(forward, inverse))
            secax.set_frame_on(True)
            secax.patch.set_visible(False)
            secax.xaxis.set_label_position('bottom')
            secax.spines['bottom'].set_position(('outward', 40))
            secax.set_xlabel('HTC $[W/m^2K]$')
            ax4.legend(fontsize=10)

        ax1.text(x=48, y=.002, s="$y = e^{a x^2 + b x +c}$", fontsize=18,
                 bbox=dict(boxstyle="square",
                           ec=(0, 0, 0),
                           fc=(1., 1, 1),
                           ))

        ax1.set_ylabel('Pumping Power $\dot{W}$ [W]')
        ax1.set_xlabel('$T_{wall}$ [°C]')
        ax1.set_yscale('log')
        ax1.set_ylim([10**-5,1.1])
        ax1.set_title('Pumping Power vs Wall Temp at: ' + name)

        ax8.set_ylabel('HTC $[W/m^2K]$')
        ax8.set_xlabel('Re')
        ax8.yaxis.grid(True)
        ax8.xaxis.grid(True)

        ax8.set_yscale('log')
        ax8.set_xscale('log')
        ax8.set_title('Reynolds vs HTC at: ' + name)
        ax8.legend(fontsize=10)

        ax6.set_xlabel('$T_{wall}$ [°C]')
        ax6.set_ylabel('Richardson Number')
        # ax6.set_yscale('log')
        # ax6.set_title('Richardson Number')
        ax6.legend(fontsize=10)
        ax7.tick_params(labelleft=False, labelright=True)

        ax2_1.set_xlabel('Test Section Pumping Power $\dot{W}$ [W]')
        ax2_1.set_ylabel('HTC $[W/m^2K]$')
        ax2_1.set_xscale('log')
        # ax2_2.set_xscale('log')
        # ax2_3.set_xscale('log')
        ax2_1.set_title('Pumping Power')
        ax2_2.set_title('$Volume Flow Rate $')
        ax2_3.set_title('Re')
        ax2_2.tick_params(labelleft=False)
        ax2_3.tick_params(labelleft=False)
        ax2_3.tick_params(labelright=True)
        ax2_1.set_xlabel('Pumping Power $\dot{W}$ [W]')
        ax2_2.set_xlabel('Volume Flow Rate [L/min]')
        # ax2_2.set_xlabel('$T_{wall}$ [C]')
        ax2_3.set_xlabel('Reynolds')
        ax2_1.legend(fontsize=10)

        fig2.suptitle('Wall Location: ' + name, fontsize=16)
        fig6.suptitle('Richardson at Wall Location: ' + name, fontsize=16)
        compiled_pumpingPower3Points = np.vstack(compiled_pumpingPower3Points)
        dataframe = pd.DataFrame(compiled_pumpingPower3Points, index=fluids,
                                 columns=["Date", "Fluid", "Pumping Power Inlet", "Pumping Power Midpoint",
                                          "Pumping Power Outlet", "Pumping Power Midpoint Error"])
        for ffff in ["Pumping Power Inlet", "Pumping Power Midpoint",
                     "Pumping Power Outlet", "Pumping Power Midpoint Error"]:
            dataframe[ffff] = pd.to_numeric(dataframe[ffff], downcast="float")

        dataframe = dataframe.sort_values("Date")
        dataframe["Fluid"] = [f.split(' ')[0] for f in dataframe["Fluid"]]

        if 1 == 1:

            x_pos = np.arange(len(dataframe["Date"]))
            # with pd.option_context('display.max_rows', None, 'display.max_columns',
            #                        None):  # more options can be specified also
            #     print(dataframe)
            FOMs = 5 / dataframe["Pumping Power Midpoint"].values  # 5W
            bars = ax3B.bar(x_pos, FOMs, yerr=np.divide(dataframe["Pumping Power Midpoint Error"],
                                                        dataframe["Pumping Power Midpoint"]) * FOMs, align='center',
                            ecolor='blue',
                            capsize=10, color="grey")

            ax3B.set_ylabel('FOM at ' + name)
            ax3B.set_xticks(x_pos)
            ax3B.set_xticklabels(dataframe["Date"], rotation=45)
            ax3B.set_title("FOM Measurements")
            ax3B.yaxis.grid(True)
            ax3B.xaxis.grid(False)
            for it, test in enumerate(dataframe["Fluid"]):
                bars[it].set_color(colors[test])

            handles = [plt.Rectangle((0, 0), 1, 1, color=colors[label]) for label in labels]
            ax3B.legend(handles, labels)

        dataframe.to_csv(
            '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210827_Compiled_data.csv')
        fig1.savefig(
            '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/4_pumppower_vs_Twall_' + name +
            '.png')
        fig8.savefig(
            '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/4_Rey_vs_HTC_' + name +
            '.png')
        fig6.savefig(
            '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/4_Ri_vs_WallTemp' + '_' + name +
            '.png')
        fig5.savefig(
            '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/4_fluids_comp_' + name + '.png')
        fig2.savefig(
            '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/4_pumppower_vs_HTC_' + name +
            '.png')
        fig3.savefig(
            '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/4_fluid_Time' + '_' + fluid_name
            + '_' + name + '.png')
        fig3B.savefig(
            '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/4_fluid_TimeB' + '_' +
            fluid_name + '_' + name + '.png')


if __name__ == "__main__":
    analysis(name="Midpoint")
    analysis(name="Outlet")
    analysis(name="Inlet")
    fig4.savefig(
        '/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/4_pumppower_vs_Twall_scale.png')
