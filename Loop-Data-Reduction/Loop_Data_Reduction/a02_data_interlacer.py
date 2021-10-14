import os
from os.path import join as pjoin  # tool to make paths
import numpy as np
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import pathlib
from scipy import interpolate
from os.path import join as pjoin  # tool to make paths
import streamlit as st
import itertools
import altair as alt
from astropy.stats import sigma_clip
from matplotlib import cm
from datetime import datetime


def setup_text_plots(fontsize=8, usetex=True, style='default'):
    """
    e.g. setup_text_plots(fontsize=14, usetex=True,style='default')
    :param fontsize:
    :param usetex:
    :param style:
    :return: setup for nice plotting
    """
    import matplotlib
    matplotlib.style.use(style)
    matplotlib.rcParams['savefig.dpi'] = 150
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
    matplotlib.rcParams['figure.figsize'] = [6.0, 6.0]
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
    matplotlib.rcParams['axes.xmargin'] = 0.0
    matplotlib.rcParams['axes.ymargin'] = 0.0


setup_text_plots()


def fast_nearest_interp(xi, x, y):
    """Assumes that x is monotonically increasing!!."""
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]


@st.cache
def importdata_interlace(datapath, sampling_IR=[100, 5], change_time=[60], ):
    # sampling_IR = [ # frames, time between in minutes]

    csv_loc = \
    [f for f in os.listdir(datapath) if f.endswith('.csv') and not f.startswith('Rec-') and not f.startswith('coll')][0]
    IR_csv_loc = [f for f in os.listdir(datapath) if f.endswith('.csv') and f.startswith('Rec-')][0]
    data_fluid = pd.read_csv(pjoin(datapath, 'colloid_data.csv'))

    # import csv data file
    data = pd.read_csv(pjoin(datapath, csv_loc))[::1]  # resample

    data.columns = [f.strip() for f in data.columns]  # removing leading space
    print(data.columns)
    # "%d:%H:%M:%S%3u"

    try:
        data["Time"] = pd.to_datetime(data['Time'], format="%H:%M:%S.%f")  # DD:HH:MM:SS.fff
    except:
        # data["Time"] = [f[0:2]+':' + f[2:] for f in data["Time"]] # temp fix
        print(data)
        data["Time"] = pd.to_datetime(data['Time'], format="%d:%H:%M:%S.%f")  # DD:HH:MM:SS.fff

        print('exception')
        print(data)
    data["Time Delta"] = data['Time'].diff().astype('timedelta64[ms]')
    # data["Time Delta"][0] = 0
    data.loc[0,"Time Delta"] = 0

    # data["Time Delta"][data["Time Delta"] < 0] = data["Time Delta"][data["Time Delta"] < 0] + 3600 * 24 * 1000  #
    data.loc[data["Time Delta"] < 0,"Time Delta"] = data["Time Delta"][data["Time Delta"] < 0] + 3600 * 24 * 1000  #

    # correction for data straddling two days
    data["Time (ms) old"] = data["Time (ms)"]
    data["Time (ms)"] = data["Time Delta"].cumsum()

    startTime = pd.to_datetime(pd.DataFrame({'year': [csv_loc.split('_')[0]],
                                             'month': [csv_loc.split('_')[1]],
                                             'day': [csv_loc.split('_')[2]],
                                             'hour': [csv_loc.split('_')[3]],
                                             'minute': [csv_loc.split('_')[4]],
                                             'second': [csv_loc.split('_')[5].split('.')[0]]
                                             }))
    print("START TIME", startTime[0])
    data["Time Stamp"] = pd.to_datetime(data['Time (ms)'], unit='ms', origin=startTime[0])  # DD:HH:MM:SS.fff
    data["Time (hr)"] = data["Time (ms)"] / 3600 / 1000
    data["$\Delta T$"] = data["T2 TS Out"] - data["T3 TS In"]
    data["$\Delta T_{setpoint}$"] = data["T2 TS Out"] - data["T Preheat Setpoint"]
    data["Mass Flow Rate (kg/hr)"] = data["Mass Flow Rate (kg/s)"] * 3600
    data["Volume Flow Rate (L/min)"] = data["Mass Flow Rate (kg/s)"] * 60 / data_fluid["rho_nf"][0] * 1000

    data["$\dot{Q}$"] = np.multiply(data["Mass Flow Rate (kg/hr)"] / 3600, data["T2 TS Out"] - data["T3 TS In"]) * \
                        data_fluid['c_nf'][0]
    data = data.rename(columns={"DPTS": 'DPTS [inH20]'})

    # import corresponding fits image file
    ## This
    positions = []

    if 1 == 1:
        data_IR = pd.read_csv(pjoin(datapath, IR_csv_loc), index_col=0)
        IR_z_pos = data_IR.columns
        data_IR.index = pd.to_datetime(data_IR.index)

        # for each slice in the main data file, we need to assign IR data, which we then use to calibrate
        for i, pos_i in enumerate(IR_z_pos):
            # getting the
            data["T_{wall} Raw P" + str(i)] = fast_nearest_interp(data["Time Stamp"].values.astype(int) / 10 ** 9,
                                                                  data_IR.index.values.astype(int) / 10 ** 9,
                                                                  data_IR[pos_i].values)
            positions.append(["T_{wall} Raw P" + str(i), pos_i])

    return data, np.vstack(positions)


def plot_set(data_wide, idx_calibration, idx_change_1, idx_change_0, datapath, y_var="$\Delta T$", title=""):
    fig, (ax) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})
    plt.subplots_adjust(left=.15, bottom=.333, right=.85, top=.95, wspace=0, hspace=0)

    capture_start = .75
    capture_range = .2
    idx_look = np.argwhere(idx_calibration == 1).ravel()
    colors = cm.inferno(np.linspace(0, 1, idx_look.size + 2))
    data_a = data_wide[0:0]  # empty array to hold the data
    for ii, i in enumerate(idx_look):
        data_s = data_wide.iloc[idx_change_0[i]:idx_change_1[i]]

        idx_range = idx_change_1[i] - idx_change_0[i]
        idx_a, idx_b = int(capture_start * idx_range), int(capture_start * idx_range + capture_range * idx_range)
        # local_dataframe= data_s.iloc[idx_a:idx_b].agg(['mean'])
        local_dataframe = data_s.iloc[idx_a:idx_b].mean(numeric_only=None)
        data_a = data_a.append(local_dataframe, ignore_index=True)
        label = "%1.0i:  Temp: %3.1f, Mass Flow: %3.2f" % (i,
                                                           local_dataframe["T Preheat Setpoint"],
                                                           local_dataframe["Mass Flow Rate (kg/hr)"])
        ax.plot(data_s["Time (hr)"] - data_s["Time (hr)"].iloc[0], data_s[y_var], label=label, c=colors[ii])

        ax.axhline(y=local_dataframe[y_var], linestyle="-.", lw=1, c=colors[ii], label='_nolegend_')
    ax.set_ylabel(y_var)  # we already handled the x-label with ax1
    ax.set_xlabel('Time (hr)')  # we already handled the x-label with ax1
    ax.set_title(title)  # we already handled the x-label with ax1
    # ax.xaxis.tick_top()
    # ax.xaxis.set_label_position("top")
    # ax.set_xlim([0,np.max(data_s["Time (hr)"])])

    fig.legend(loc='lower center', borderaxespad=0.2, fontsize=8, frameon=False)
    fig.savefig(pjoin(datapath, '0_' + title.replace(' ', '') + '.png'))

    return fig, data_a


def plot_setim(data_wide, positions):
    fig, (ax1) = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3]})

    if 1 == 1:
        imaa = ax1.imshow(data_wide[positions[:, 0]].to_numpy().transpose(), aspect='auto', cmap="inferno")
        cbar = plt.colorbar(imaa, orientation='horizontal', aspect=50, fraction=.05, shrink=1.0, pad=0.05,
                            label='Counts')
        ax1.set_yticks(np.arange(len(positions[:, 1])))
        ax1.set_yticklabels(list(positions[:, 1]))
        ax1.tick_params(bottom='off')
        ax1.set_xticklabels([])
        ax1.set_xticks([])
        ax1.grid(None)
        ax1.margins(x=0)

    return fig


def interlace_main(datapath):
    data_fluid = pd.read_csv(pjoin(datapath, 'colloid_data.csv'))

    data_wide, positions_wall = importdata_interlace(
        datapath=datapath)

    # mark the time zones of interest, anytime the setpoint temperature or flow rate is changed
    change_saph = abs(data_wide["Sapphire Status"].diff()) > 0
    change = abs(data_wide["T Preheat Setpoint"].diff() + data_wide["Sapphire Status"].diff()) > 0

    idx_change_1 = np.asarray([i for i in range(change.size) if change.iloc[i] == 1])
    idx_change_0 = np.concatenate((np.asarray([0]), idx_change_1[0:-1]))

    # initialize defaulting to 0
    idx_calibration = np.zeros(idx_change_1.size)
    idx_measurement = np.zeros(idx_change_1.size)
    idx_calibration_flow = np.zeros(idx_change_1.size)
    for i, idx_i in enumerate(idx_change_1):
        # print(idx_i,data_wide["Time (hr)"][idx_change_1])

        # set temp changes are cal. We already ignore the final temp change because there is no index at the end
        if i == 0:
            idx_calibration[i] = 1
        elif abs(data_wide["T Preheat Setpoint"].diff()[idx_change_1[i - 1]]) > 0:
            idx_calibration[i] = 1
        if i < idx_change_1.size - 1 and data_wide["Sapphire Status"].diff()[idx_change_1[i + 1]] < 0:
            idx_calibration_flow[i] = 1
        # find the sapphire on points
        if data_wide["Sapphire Status"].diff()[idx_i] < 0:  # we're looking at the end, so when it turns off
            idx_measurement[i] = 1

    data_wide["Change"] = abs(data_wide["T Preheat Setpoint"].diff() + data_wide["Sapphire Status"].diff())

    # big picture
    fig, (ax, ax3, ax3a, ax4) = plt.subplots(4, 1, gridspec_kw={'height_ratios': [3, 2, 2, 1]})
    plt.subplots_adjust(left=.15, bottom=.15, right=.85, top=.85, wspace=0, hspace=0)

    ax.plot(data_wide["Time (hr)"], data_wide["T2 TS Out"], c='b', label="TS Out")
    ax.plot(data_wide["Time (hr)"], data_wide["T Preheat Setpoint"], '--', c='b', label="TS Out Setpoint ")
    ax.set_ylabel('Temp TS Out (°C)')  # we already handled the x-label with ax1
    ax.xaxis.set_label_position("top")
    ax.set_xlabel('Time (hr)')  # we already handled the x-label with ax1
    ax.xaxis.tick_top()

    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('Mass Flow (kg/hr)', color=color)
    ax2.plot(data_wide["Time (hr)"], data_wide["Mass Flow Rate (kg/hr)"], color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    ax3.plot(data_wide["Time (hr)"], data_wide["Pump Pressure Out"])
    ax3.plot(data_wide["Time (hr)"], data_wide["Pump Pressure In"])

    ax3.set_ylabel('Pressure')
    data_wide['Pumping Power Over TS'] = np.multiply(data_wide["DPTS [inH20]"], data_wide["Mass Flow Rate (kg/s)"]) / \
                                         data_fluid["rho_nf"][0]
    ax3a.plot(data_wide["Time (hr)"], data_wide["Pumping Power Over TS"])
    # ax3a.plot(data_wide["Time (hr)"], data_wide["PreHeat Power"])
    ax3a.set_ylabel('Pumping Power (W)')
    ax3a.yaxis.set_label_position("right")
    ax3a.yaxis.tick_right()

    indicators = ['Sapphire Status', 'Preheat Status', 'Temp Trip', 'Pressure Trip', 'Flow Trip']
    swithc_im = data_wide[indicators].transpose()
    ax4.imshow(swithc_im, aspect='auto', cmap="binary")
    ax4.set_yticks(np.arange(len(indicators)))
    ax4.set_yticklabels(list(indicators))
    ax4.tick_params(bottom='off')
    ax4.set_xticklabels([])
    ax4.set_xticks([])
    ax4.grid(None)
    ax.margins(x=0)
    ax2.margins(x=0)
    ax3a.margins(x=0)
    ax3.margins(x=0)

    st.pyplot(fig)
    fig.savefig(pjoin(datapath, '0_summary.png'))

    if 1 == 1:
        fig0, (ax02) = plt.subplots(1, gridspec_kw={'height_ratios': [3]})
        # fig1.subplots_adjust(left=.15, bottom=.15, right=.85, top=.85, wspace=0, hspace=0)
        x, y = data_wide['Mass Flow Rate (kg/hr)'], data_wide['DPTS [inH20]']
        ax02.scatter(x, y, c=data_wide['T2 TS Out'], s=5)
        x0 = np.linspace(0, np.max(x) * 1.2, 100)

        ax02.set_ylabel('Pressure Drop TS (inH20)')
        ax02.set_xlabel('Mass Flow Rate (kg/s)')
        ax02.set_title('Pressure Drop over Test Section')

        st.pyplot(fig0)
        fig0.savefig(pjoin(datapath, '0_flow_delP_check.png'))

    fig2, data_a = plot_set(data_wide, idx_calibration, idx_change_1, idx_change_0, datapath,
                            y_var="$\Delta T_{setpoint}$",
                            title="Pre Heater Performance")
    st.pyplot(fig2)
    fig3, data_calibration = plot_set(data_wide, idx_calibration, idx_change_1, idx_change_0, datapath,
                                      y_var="$\Delta T$",
                                      title="Calibration Data")
    st.pyplot(fig3)

    fig4, data_measurement = plot_set(data_wide, idx_measurement, idx_change_1, idx_change_0, datapath,
                                      y_var="$\Delta T$",
                                      title="Flow Change With Heat")
    st.pyplot(fig4)
    fig5, data_calibration_flow = plot_set(data_wide, idx_calibration_flow, idx_change_1, idx_change_0, datapath,
                                           y_var="$\Delta T$", title="Flow Change Without Heat")

    fig7, throwaway = plot_set(data_wide, idx_calibration, idx_change_1, idx_change_0, datapath,
                               y_var="PreHeat Power",
                               title="Calibration Preheat Power [W]")
    st.pyplot(fig7)
    fig7a, throwaway = plot_set(data_wide, idx_measurement, idx_change_1, idx_change_0, datapath,
                                y_var="PreHeat Power",
                                title="Measurement Preheat Power [W]")
    st.pyplot(fig7a)
    fig7b, throwaway = plot_set(data_wide, idx_calibration_flow, idx_change_1, idx_change_0, datapath,
                                y_var="PreHeat Power",
                                title="Calibration No Flow Preheat Power [W]")
    st.pyplot(fig7b)

    fig8, throwaway = plot_set(data_wide, idx_measurement, idx_change_1, idx_change_0, datapath,
                               y_var="Sapphire Power",
                               title="Sapphire Power [W]")
    st.pyplot(fig8)

    st.pyplot(fig5)
    fig6 = plot_setim(data_wide, positions_wall)
    st.pyplot(fig6)

    calibration_setpoint_delT_correction = \
    data_calibration[data_calibration['T Preheat Setpoint'] == data_measurement['T Preheat Setpoint'][0]][
        "$\Delta T$"].values
    calibration_setpoint_delT_correction

    # corrction from the flow rate change, and correction for the particular temperature
    data_measurement["$\Delta T_corrected$"] = data_measurement["$\Delta T$"] - data_calibration_flow["$\Delta T$"]
    data_measurement["$\dot{Q}_corrected$"] = data_measurement["$\dot{Q}$"] - data_calibration_flow["$\dot{Q}$"]

    data_calibration
    data_calibration_flow
    data_measurement

    data_wide.to_pickle(pjoin(datapath, 'd_wide.pkl'))
    data_calibration.to_pickle(pjoin(datapath, 'd_cal.pkl'))
    data_calibration_flow.to_pickle(pjoin(datapath, 'd_cal_flow.pkl'))
    data_measurement.to_pickle(pjoin(datapath, 'd_meas.pkl'))
    np.save(pjoin(datapath, 'positions_wall.pkl'), positions_wall)


if __name__ == "__main__":
    datapath = "/Users/hanscastorp/Dropbox/MIT-LUBRIZOL/LabData/Lubrizol Loop/20210726_isoparaffin_5W"
    interlace_main(datapath)
    plt.show()
