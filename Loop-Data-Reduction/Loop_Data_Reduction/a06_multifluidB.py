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


def mixedConvection(c, k, mu, rho, deltaTwall2bulk, heatflux, d_h, ):
    h = heatflux / deltaTwall2bulk
    t1 = 0.075 * ((k) / (c * mu)) ** (2 / 3) + 1
    t2 = (d_h ** (3) * h ** (3)) / (k ** (3))
    t3n = 0.18 * 9.81 ** (3 / 4) * ((c * mu) / (k)) ** (3 / 4)
    t3d = (0.63 * (k / (c * mu)) ** (9 / 16) + 1) ** (4 / 3)
    t4 = rho * d_h * c ** (2 / 3)
    v = (4.7 * mu ** (1 / 3)) * (k * t1 ** (3 / 4) * (t2 - t3n / t3d)) ** (2 / 3) / (t4)
    w = mu * v ** 2
    FOM = w ** -1
    return v, w, FOM


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

"""
Collects data from all the specified datapaths and does analysis. Uses the folder name as YYYYMMDD_FLUID to parse. 
Similar to a04_multifluid except the data is easier to work with 

"""
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
              "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210930_Water_5W", ]


def forward(x):
    return 5 / (x - 30) / (4.25 / 10 ** 4)


def inverse(x):
    return 5 / x / (4.25 / 10 ** 4) + 30


import seaborn as sns
import sys


# np.set_printoptions(threshold=sys.maxsize)
#
# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)

def analysis(dataFrame_agg_t, xdata_new):
    dataFrame_agg = dataFrame_agg_t[dataFrame_agg_t["Set Type"] == "Measurement"]

    newWalls = ["Wall Position " + str(ii) for ii in np.round(xdata_new, 2)]

    print(dataFrame_agg.columns.tolist())
    dataFrame_agg["HTC Mixed Predicted"] = np.divide(np.multiply(dataFrame_agg["Nu Mixed"], dataFrame_agg["k nf"]),
                                                     dataFrame_agg["D h"])
    dataFrame_agg["Nu laminar"] = np.multiply(np.power(dataFrame_agg["Re"], .6) * .128, dataFrame_agg["Pr"] ** .4)
    dataFrame_agg["HTC Laminar Predicted"] = np.divide(np.multiply(dataFrame_agg["Nu laminar"], dataFrame_agg["k nf"]),
                                                       dataFrame_agg["D h"])
    uniqueFluids = pd.unique(dataFrame_agg["Fluid"])
    print(uniqueFluids)

    # sns.scatterplot('Re', 'HTC', data=dataFrame_agg, hue='Fluid')
    if 1 == 1:
        plt.figure()

    if 1 == 1:
        fig, (ax1) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
        fig11b, (ax1b) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
        fig11, (ax11) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
        import matplotlib.colors
        colors = {}
        for iic, ic in enumerate(uniqueFluids):
            colors[ic] = plt.cm.tab10(iic)
        popt_ReHTC_s = []
        laminaryvals = []
        for i in range(len(uniqueFluids)):
            dataFrame_agg_f = dataFrame_agg[dataFrame_agg["Fluid"] == uniqueFluids[i]]

            # interpolate HTC vs Reynolds
            def funcHTC(x, a, b):
                return np.exp(a * np.log(x) + b)

            popt_ReHTC, pcov = curve_fit(funcHTC, dataFrame_agg_f['Re'], dataFrame_agg_f['HTC'],
                                         # bounds=([0,-np.inf], [np.inf,np.inf])
                                         )
            popt_ReHTC_s.append(popt_ReHTC)
            ax1.errorbar(dataFrame_agg_f['Re'], dataFrame_agg_f['HTC'], xerr=0,
                         yerr=dataFrame_agg_f['HTC err'], ls='none', c=colors[uniqueFluids[i]], alpha=1, lw=.5,
                         label=uniqueFluids[i])
            ax1b.scatter(dataFrame_agg_f['Re'], dataFrame_agg_f['Nu'],  c=colors[uniqueFluids[i]], alpha=1,
                         label=uniqueFluids[i])

            xrange_re = np.power(10,np.linspace(0, 5, 100))
            # ax1.plot(dataFrame_agg_f['Re'].sort_values(), (
            #             dataFrame_agg_f['Re'] ** .6 * dataFrame_agg_f['Pr'] ** .4 * .128 * dataFrame_agg_f['k nf'] *
            #             dataFrame_agg_f['D h'] ** -1).sort_values()
            #          , c=colors[uniqueFluids[i]], linestyle="--", label=uniqueFluids[i] + " Laminar")

            yvals = xrange_re ** .6 * np.mean(dataFrame_agg_f['Pr']) ** .4 * .128 * np.mean(
                dataFrame_agg_f['k nf']) * np.mean(dataFrame_agg_f['D h']) ** -1

            laminaryvals.append(yvals)
            ax1.plot(xrange_re, yvals, c=colors[uniqueFluids[i]], linestyle="--", label=uniqueFluids[i] + " Laminar Hewitt")
            ax1b.plot(xrange_re, xrange_re ** .6 * np.mean(dataFrame_agg_f['Pr']) ** .4 * .128, c=colors[uniqueFluids[i]], linestyle="--", label=uniqueFluids[i] + " Laminar Hewitt")

            # Perkins et al from Shah text for 4 sided heated square channel in thermally developing hydrodynamically
            # developed flow
            Prandtl = np.mean(dataFrame_agg_f['Pr'])
            print('------------ ', np.mean(dataFrame_agg_f['Pr']))
            xrange_xstar = ((4 / .651) * xrange_re ** -1 * Prandtl ** -1)
            yvals_Perkins = np.divide(1, .277 - .152 *
                                      np.power(2.71828, -38.6 * xrange_xstar))
            ax1.plot(xrange_re, yvals_Perkins* np.mean(
                dataFrame_agg_f['k nf']) * np.mean(dataFrame_agg_f['D h']) ** -1, c=colors[uniqueFluids[i]], linestyle="-.",
                     label=uniqueFluids[i] + " Laminar Perkins")


            ax1b.plot(xrange_re,yvals_Perkins, c=colors[uniqueFluids[i]], linestyle="-.",
                     label=uniqueFluids[i] + " Laminar Perkins")

            yvals_Churchill = np.power(1+np.power(220/np.pi*xrange_xstar,-10/9),3/10)*5.365-1
            ax1.plot(xrange_re, yvals_Churchill * np.mean(
                dataFrame_agg_f['k nf']) * np.mean(dataFrame_agg_f['D h']) ** -1, c=colors[uniqueFluids[i]],
                     linestyle="-.",
                     label=uniqueFluids[i] + " Laminar Churchill")
            ax1b.plot(xrange_re,yvals_Churchill, c=colors[uniqueFluids[i]], linestyle=":",
                     label=uniqueFluids[i] + " Laminar Churchill") # churchill oand Ozoe from Shah p128 for thermally developing hydrodynamically developed specified heat flux in cirucular tube
            # ax1.plot(dataFrame_agg_f['Re'].sort_values(), dataFrame_agg_f['HTC Mixed Predicted'].sort_values(),
            # c=colors[uniqueFluids[i]],
            #          linestyle=":", label=uniqueFluids[i] + " Mixed")

            print('_____________________ FLUID HTC vs Re fit ____ ', uniqueFluids[i])
            print(popt_ReHTC)
            ax1.plot(xrange_re, funcHTC(xrange_re, *popt_ReHTC), c=colors[uniqueFluids[i]],
                     label=uniqueFluids[i] + ' Fit: a=%5.3f, b=%5.3f' % tuple(popt_ReHTC))

            ax1.set_ylim([100, 10000])

            targetReynolds = 200
            print()
            interpHTC = funcHTC(targetReynolds, *popt_ReHTC)
            print(f"HTC at ", targetReynolds, ": ", interpHTC)

        ax1.set_xlabel('Re')
        ax1.set_ylabel('HTC $[W/m^2K]$')
        ax1.set_title("HTC Estimates")
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1b.set_xlabel('Re')
        ax1b.set_ylabel('Nu')
        ax1b.set_title("Nu")
        ax1b.set_xscale('log')
        ax1b.set_yscale('log')
        fig11b.subplots_adjust(left=.15, bottom=.1, right=.7, top=.95, wspace=0, hspace=0)
        fig11b.legend(loc='center right', borderaxespad=0.2, fontsize=8, frameon=False)
        fig.subplots_adjust(left=.15, bottom=.1, right=.7, top=.95, wspace=0, hspace=0)
        fig.legend(loc='center right', borderaxespad=0.2, fontsize=8, frameon=False)

        idx_compare = 1
        for i in range(len(uniqueFluids)):
            yvals = np.divide(funcHTC(xrange_re, *popt_ReHTC_s[i]), funcHTC(xrange_re, *popt_ReHTC_s[idx_compare]))
            ax11.plot(xrange_re, yvals, c=colors[uniqueFluids[i]],
                      label=uniqueFluids[i] + " Meas. Data Fit")
            ax11.plot(xrange_re, np.divide(laminaryvals[i], laminaryvals[idx_compare]), c=colors[uniqueFluids[i]],
                      label=uniqueFluids[i] + " Laminar Meas. Prop.", linestyle='--')

        ax11.set_xlabel('Re')
        ax11.set_ylabel('Normalized HTC')
        ax11.set_title("HTC Estimates Normalized to Basefluid")
        plt.xscale('log')
        # plt.yscale('log')
        ax11.legend()
        fig11.savefig(pjoin(dir_path,"06_HTC_Re_fits.png"))
        fig.savefig(pjoin(dir_path,"06_HTC_Re.png"))
        fig11b.savefig(pjoin(dir_path,"06_Nu_Re.png"))

    if 1 == 1:
        fig4, (ax4) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
        dtemp_t = dataFrame_agg_t[dataFrame_agg_t["Wall Position"] == newWalls[-1]]
        temp_names = ['T1 PreHeat', 'T3 TS In', 'T2 TS Out', 'T4 HX In', 'T5 HX Out']
        idvars = ["Fluid", "Mass Flow Rate (kg/s)", "Set Type", 'T wall Midpoint', "HTC", "HTC knownQ", "Re",
                  "$\dot{Q} corrected$",
                  "Pumping Power (DP)[W]", "Pumping Power (corr)[W]"]

        dtemp_t = pd.melt(dtemp_t, id_vars=idvars, var_name="Temp Type", value_vars=temp_names,
                          value_name="Temp [°C]")
        # 'T1 PreHeat', 'T2 TS Out',
        #        'T3 TS In', 'T2B TS Out', 'T4 HX In', 'T5 HX Out'
        sns.violinplot(x="Temp [°C]", y="Temp Type", hue="Set Type", data=dtemp_t, inner=None, color=".8", height=6,
                       split=True)

        sns.stripplot(x="Temp [°C]", y="Temp Type", hue="Fluid", data=dtemp_t, dodge=True)
        plt.title('Loop Temperatures')
        fig4.savefig(pjoin(dir_path,"06_loopTemps.png"))

    if 1 == 1:
        fig2, (ax2) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
        dtemp = dataFrame_agg_t[dataFrame_agg_t["Wall Position"] == newWalls[-1]]
        # dataFrame_i["mu_nf"]
        dtemp["Pressure Drop Correrlation"] = dtemp["Pumping Power (corr)[W]"] / np.max(
            dataFrame_agg["Pumping Power (corr)[W]"])
        dtemp["Measured Pressure Drop"] = dtemp["Pumping Power (DP)[W]"] / np.max(dtemp["Pumping Power (DP)[W]"])
        dtemp = pd.melt(dtemp, id_vars=["Fluid", "Mass Flow Rate (kg/s)", "Set Type", 'T wall Midpoint'],
                        var_name="Pumping Power Estimate",
                        value_vars=["Pressure Drop Correrlation", "Measured Pressure Drop"],
                        value_name="Pumping Power [W]")
        sns.scatterplot('T wall Midpoint', "Pumping Power [W]",
                        data=dtemp, hue='Fluid', style="Pumping Power Estimate")

        # plt.xscale('log')
        plt.yscale('log')
        plt.title("Pumping Power Estimate Comparison")
        fig2.savefig(pjoin(dir_path,"06_pumpingPowers2.png"))

    if 1 == 1:
        height = 1  # m

        fig5, (ax5) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
        dtemp = dataFrame_agg[dataFrame_agg["Wall Position"] == newWalls[-1]]
        dtemp["$\Delta P_{fric}$"] = dtemp["Pumping Power (corr)[W]"] * 10 / dtemp["Mass Flow Rate (kg/s)"] * dtemp[
            "rho nf"]
        # 'T1 PreHeat', 'T2 TS ut', 'T3 TS In',
        dtemp["$\Delta P_{grav}$"] = -height * 9.81 * np.multiply(dtemp["T1 PreHeat"] - dtemp["T2 TS Out"],
                                                                  np.multiply(dtemp["CTE"], dtemp["rho nf"]))
        dtemp["$\Delta P_{measure}$"] = dtemp["DPTS [inH20]"] * 249.09
        press_names = ['$\Delta P_{fric}$', '$\Delta P_{grav}$', '$\Delta P_{measure}$']
        dtemp1 = pd.melt(dtemp, id_vars=idvars, var_name="Pressure Type",
                         value_vars=press_names,
                         value_name="Pressure [Pa]")
        sns.violinplot(x="Pressure [Pa]", y="Pressure Type", data=dtemp1, inner=None, color=".8", height=6)

        sns.stripplot(x="Pressure [Pa]", y="Pressure Type", hue="Fluid", data=dtemp1, dodge=True)

        plt.title('Pressure Drop Estimates')
        fig5.savefig(pjoin(dir_path,"06_delPs.png"))

        fig6, (ax6) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
        fig6.subplots_adjust(bottom=.2)
        dtemp["$\Delta P_{grav} / \Delta P_{measure}$"] = np.divide(dtemp["$\Delta P_{grav}$"],
                                                                    dtemp['$\Delta P_{measure}$'])
        sns.scatterplot('T wall Midpoint', "$\Delta P_{grav} / \Delta P_{measure}$", hue="Fluid", data=dtemp, ax=ax6)
        plt.yscale('log')
        plt.title('Calculated Pressure Drop Components')
        secax1 = ax6.secondary_xaxis('bottom', functions=(forward, inverse))
        secax1.set_frame_on(True)
        secax1.patch.set_visible(False)
        secax1.xaxis.set_label_position('bottom')
        secax1.spines['bottom'].set_position(('outward', 40))
        secax1.set_xlabel('HTC $[W/m^2K]$')
        fig6.savefig(pjoin(dir_path,"06_delPs_ratioTwall.png"))

        fig7, (ax7) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
        fig7.subplots_adjust(bottom=.2)
        dtemp["$\Delta T_{outlet-inlet}$"] = dtemp['T2 TS Out'] - dtemp['T3 TS In']
        sns.scatterplot('T wall Midpoint', "$\Delta T_{outlet-inlet}$", hue="Fluid", data=dtemp, ax=ax7)
        plt.yscale('log')
        plt.title('$\Delta T_{outlet-inlet}$')
        secax1 = ax7.secondary_xaxis('bottom', functions=(forward, inverse))
        secax1.set_frame_on(True)
        secax1.patch.set_visible(False)
        secax1.xaxis.set_label_position('bottom')
        secax1.spines['bottom'].set_position(('outward', 40))
        secax1.set_xlabel('HTC $[W/m^2K]$')
        fig7.savefig(pjoin(dir_path,"06_delTs_Twall.png"))

        plt.figure()
        sns.scatterplot('Re', "$\Delta P_{grav} / \Delta P_{measure}$", hue="Fluid", data=dtemp)
        plt.yscale('log')
        plt.xscale('log')
        plt.title('Calculated Pressure Drop Components')
        plt.savefig(pjoin(dir_path,"06_delPs_ratioRe.png"))

        plt.figure()
        sns.scatterplot('HTC knownQ', "$\Delta P_{grav} / \Delta P_{measure}$", hue="Fluid", data=dtemp)
        plt.yscale('log')
        # plt.xscale('log')
        plt.xlabel("HTC $[W/m^2K]$")
        plt.title('Calculated Pressure Drop Components')
        plt.savefig(pjoin(dir_path,"06_delPs_ratioHTC.png"))

        plt.figure()
        sns.scatterplot('Re', "$\dot{Q} corrected$", hue="Fluid", data=dtemp)
        # plt.yscale('log')
        plt.xscale('log')
        plt.axhline(y=5, c='r', linestyle='--')
        plt.title('Bulk $\dot{Q}$ [W]')
        plt.ylabel('$\dot{Q}$ [W]')
        plt.savefig(pjoin(dir_path,"06_QCorrectedRe.png"))

        plt.figure()
        sns.scatterplot('Pumping Power (DP)[W]', "Pumping Power (corr)[W]", hue="Fluid", data=dtemp)
        # plt.yscale('log')
        # plt.xscale('log')
        plt.title('')

        plt.savefig(pjoin(dir_path,"06_pumpingpowerEstimates.png"))

    if 1 == 1:
        fig3 = plt.figure()
        import matplotlib.colors
        colors = {}
        for iic, ic in enumerate(uniqueFluids):
            colors[ic] = plt.cm.tab10(iic)

        labels = list(colors.keys())
        ax = fig3.gca(projection='3d')
        for i in range(len(uniqueFluids)):
            dataFrame_agg_f = dataFrame_agg[dataFrame_agg["Fluid"] == uniqueFluids[i]]
            idvars = dataFrame_agg_f.columns
            idvars = [f for f in idvars if not f in newWalls]

            ax.scatter3D(dataFrame_agg_f['Wall Position [cm]'], dataFrame_agg_f['Pumping Power (DP)[W]'],
                         dataFrame_agg_f['HTC Local'],
                         color=colors[uniqueFluids[i]], linewidth=0.01, label=uniqueFluids[i])

            ax.legend()

        ax.set_xlabel('Wall Position [cm]')
        ax.set_ylabel('Pumping Power [W]')
        ax.set_zlabel('HTC Local $[W/m^2K]$')

    return


if __name__ == "__main__":
    n_points = 20
    xdata_new = np.linspace(0, 6.5, n_points)
    newWalls = ["Wall Position " + str(ii) for ii in np.round(xdata_new, 2)]

    dir_path = "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop"
    sets = [datapaths1, datapaths2, datapaths3,datapaths4]
    # parse the data and create melted dataframe with row for each wall position which is easier to work
    # with. only run this if statement if you have changed something in the original data.
    if 1 == 0:
        dataFrame_agg = pd.DataFrame()
        for iset, datapaths in enumerate(sets):
            for i, datapath in enumerate(datapaths):
                test_name = datapath.split('/Lubrizol Loop/')[1]
                date = datapath.split('/Lubrizol Loop/')[1].split('_')[0]
                fluid_name = datapath.split('_')[1]
                print('---',date,fluid_name)
                data_measurement = pd.read_pickle(pjoin(datapath, 'd_meas.pkl_a03'))
                data_calibration_flow = pd.read_pickle(pjoin(datapath, 'd_cal_flow_a03.pkl'))
                positions_wall = np.load(pjoin(datapath, 'positions_wall.pkl.npy'))
                data_fluid = pd.read_csv(pjoin(datapath, 'colloid_data.csv'))

                data_measurement["Set_Type"] = "Measurement"
                data_calibration_flow["Set_Type"] = "Calibration Flow"
                fluid_bulkT = np.mean(data_measurement['T2 TS Out'])

                dataFrame_i = pd.concat([data_measurement, data_calibration_flow], axis=0, ignore_index=False)
                # dataFrame_i = data_measurement
                dataFrame_i["Fluid"] = fluid_name
                dataFrame_i["Date"] = date
                dataFrame_i["k_nf"] = data_fluid["k_nf"].values[0]
                dataFrame_i["c_nf"] = data_fluid["c_nf"].values[0]
                dataFrame_i["CTE"] = data_fluid["CTE"].values[0]

                dataFrame_i["mu_nf"] = data_fluid["eta"][0] * np.exp(
                    data_fluid["A"][0] / (fluid_bulkT + 273.15 - data_fluid["T0"][
                        0]))
                dataFrame_i["rho_nf"] = data_fluid["rho_nf"].values[0]
                D_h = .00652581426
                L_t = 6.510 / 100
                dataFrame_i["D_h"] = D_h
                dataFrame_i["A_s_thinfilm"] = D_h * L_t

                poswall = [f + '_corr' for f in positions_wall[:, 0]]


                def func_entrance(x, a, b, c, d):
                    return a - b * (x / 100 / D_h + c) ** (-2 / 3) + d * x ** 2
                    return a * x ** 3 + a * x ** 3 + a * x ** 3 + a * x ** 3 + a * x ** 3


                # resample wall temperatures to same grid defined above as xdata_new
                xdata = positions_wall[:, 1].astype(float)[::-1]
                for iii, newWall in enumerate(newWalls):
                    dataFrame_i[newWall] = 0

                for i, mass_flow_rate in enumerate(dataFrame_i["Mass Flow Rate (kg/s)"]):
                    d_local = dataFrame_i[dataFrame_i["Mass Flow Rate (kg/s)"] == mass_flow_rate]
                    ydata = d_local[poswall].values.ravel()

                    print(xdata.shape)
                    print(ydata.shape)
                    if 1 == 0:
                        popt, pcov = curve_fit(func_entrance, xdata, ydata)
                        plt.plot(xdata, func_entrance(xdata, *popt))

                    else:
                        z = np.polyfit(xdata, ydata, 6)
                        p = np.poly1d(z)
                        # plt.plot(xdata, p(xdata))

                    # plt.scatter(xdata,ydata)
                    for iii, newWall in enumerate(newWalls):
                        dataFrame_i[newWall][dataFrame_i["Mass Flow Rate (kg/s)"] == mass_flow_rate] = p(xdata_new[iii])

                dataFrame_agg = pd.concat([dataFrame_agg, dataFrame_i], axis=0, ignore_index=False)
                print(dataFrame_agg.columns)

        print(dataFrame_agg)
        idvars = dataFrame_agg.columns
        idvars = [f for f in idvars if not f in newWalls]
        dataFrame_agg = pd.melt(dataFrame_agg, id_vars=idvars, var_name="Wall Position", value_vars=newWalls,
                                value_name="Wall Temp [C] local")
        dataFrame_agg["Wall Position [cm]"] = [np.float(f.split("Wall Position")[1]) for f in
                                               dataFrame_agg["Wall Position"]]

        dataFrame_agg["HTC Local"] = np.divide(dataFrame_agg["$\dot{Q}_corrected$"],
                                               dataFrame_agg["A_s_thinfilm"] * (
                                                       dataFrame_agg['Wall Temp [C] local'] - dataFrame_agg[
                                                   'T Bulk']))
        idvars = dataFrame_agg.columns
        idvars = {f: f.replace('_', ' ') for f in idvars}
        dataFrame_agg = dataFrame_agg.rename(columns=idvars)

        dataFrame_agg.to_pickle(pjoin(dir_path, '06_agg.pkl'))

    # this will use the pickled and melted dataframe to produce some plots and fits.
    if 1 == 1:
        dataFrame_agg = pd.read_pickle(pjoin(dir_path, '06_agg.pkl'))

        analysis(dataFrame_agg, xdata_new)
        plt.show()
