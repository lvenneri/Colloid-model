import streamlit as st
import pandas as pd
import itertools
import altair as alt
import matplotlib.pyplot as plt
from altair import datum
from os.path import join as pjoin
import matplotlib.cm as cm
import inspect
from viscosity_models import viscosity
from colloid_combo_fun import *
import base64
import os
"""
# Colloid Model for MIT-Lubrizol Electric Vehicle Battery Cooling Application"""

"""
This script allows the user to look at the colloid design space and compare colloids in different application 
environments.
"""


def basic2d(data, x, y, norm=1, axess=1):
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
    fluids_unq = np.unique(data['bf_id_s'])
    particles_unq = np.unique(data['np_id_s'])
    p_unq = np.unique(data['p'])

    colors = cm.tab10(np.linspace(0, 1, len(particles_unq)))

    # plt.figure(figsize=(13.0, 8.0), dpi=400)
    fig = plt.figure()
    plt.subplots_adjust(wspace=0.0, hspace=0.0)

    if norm == 1:
        normer = data[data["phi"] == 0][y][0]
    else:
        normer = 1
    # fluids_unq = fluids_unq[[0, 1, 3, 2]]

    numfluids = len(fluids_unq)
    for bfi, bf in enumerate(fluids_unq):

        plt.subplot(1, numfluids, bfi + 1)
        bf_data = data[data['bf_id_s'] == bf]
        for npi, nanop in enumerate(particles_unq):
            np_data = bf_data[bf_data['np_id_s'] == nanop]

            for pi, pppp in enumerate(p_unq):
                pp_data = np_data[np_data['p'] == pppp]
                if pi == 0:
                    plt.plot(pp_data[x], pp_data[y] / normer, lw=2, c=colors[npi], label=nanop.replace('_', ' '))
                else:
                    plt.plot(pp_data[x], pp_data[y] / normer, lw=2, c=colors[npi], label='_nolegend_')

        #             plt.scatter(np_data[x],np_data[y]/normer,s=2,label=nanop)

        if x == 'phi':
            plt.xlabel('$\phi$')
        else:
            plt.xlabel(x.replace('_', ' '))
        if bfi == 0:
            plt.ylabel(y.replace('_', ' '))

        if bfi > 0 and bfi < numfluids - 1:
            plt.tick_params(which="both", axis='y', labelleft='False', labelright='False', right='False',
                            left='False', )
            plt.legend()

        if bfi == numfluids - 1:
            plt.tick_params(which="both", axis='y', right='True', left='False', labelleft='False', labelright='True')
        plt.xlim([np.min(pp_data[x]), np.max(pp_data[x])])
        if y == 'deltaT':
            plt.hlines(4, 0, .35, colors='r', linestyles='dashed', lw=1)
            plt.text(4, .1, '$\Delta T$', rotation=0)
        if axess == 1:
            plt.ylim([np.nanpercentile(data[y] / normer, 0), np.nanpercentile(data[y] / normer, 100)])
        #         plt.ylim([.0001,.01])

        plt.yscale('log')
        ax = plt.gca()
        ax.set_yscale('log')
        plt.grid(which='minor', color='k', linestyle='-', alpha=0.2)
        # ax.yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))
        # ax.yaxis.set_major_formatter(
        #     ticker.FuncFormatter(lambda y, pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y), 0)))).format(y)))
        plt.title(bf, size=8)

    return fig


def basic2d_old(data, x, y, save_loc, name='', norm=1, axess=1):
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
    fluids_unq = np.unique(data['bf_id_s'])
    particles_unq = np.unique(data['np_id_s'])
    p_unq = np.unique(data['p'])

    colors = cm.tab10(np.linspace(0, 1, len(particles_unq)))

    # plt.figure(figsize=(13.0, 8.0), dpi=400)
    fig = plt.figure(figsize=(13.0, 8.0))
    plt.subplots_adjust(wspace=0.0, hspace=0.0)

    if norm == 1:
        normer = data[data["phi"] == 0][y][0]
    else:
        normer = 1
    # fluids_unq = fluids_unq[[0, 1, 3, 2]]

    numfluids = len(fluids_unq)
    for bfi, bf in enumerate(fluids_unq):

        plt.subplot(1, numfluids, bfi + 1)
        bf_data = data[data['bf_id_s'] == bf]
        for npi, nanop in enumerate(particles_unq):
            np_data = bf_data[bf_data['np_id_s'] == nanop]

            for pi, pppp in enumerate(p_unq):
                pp_data = np_data[np_data['p'] == pppp]
                if pi == 0:
                    plt.plot(pp_data[x], pp_data[y] / normer, lw=2, c=colors[npi], label=nanop.replace('_', ' '))
                else:
                    plt.plot(pp_data[x], pp_data[y] / normer, lw=2, c=colors[npi], label='_nolegend_')

        #             plt.scatter(np_data[x],np_data[y]/normer,s=2,label=nanop)

        if x == 'phi':
            plt.xlabel('$\phi$')
        else:
            plt.xlabel(x.replace('_', ' '))
        if bfi == 0:
            plt.ylabel(y.replace('_', ' '))

        if bfi > 0 and bfi < numfluids - 1:
            plt.tick_params(which="both", axis='y', labelleft='False', labelright='False', right='False',
                            left='False', )
            plt.legend()

        if bfi == numfluids - 1:
            plt.tick_params(which="both", axis='y', right='True', left='False', labelleft='False', labelright='True')
        plt.xlim([np.min(pp_data[x]), np.max(pp_data[x])])
        if y == 'deltaT':
            plt.hlines(4, 0, .35, colors='r', linestyles='dashed', lw=1)
            plt.text(4, .1, '$\Delta T$', rotation=0)
        if axess == 1:
            plt.ylim([np.nanpercentile(data[y] / normer, 0), np.nanpercentile(data[y] / normer, 100)])
        #         plt.ylim([.0001,.01])

        plt.yscale('log')
        ax = plt.gca()
        ax.set_yscale('log')
        plt.grid(which='minor', color='k', linestyle='-', alpha=0.2)
        # ax.yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))
        # ax.yaxis.set_major_formatter(
        #     ticker.FuncFormatter(lambda y, pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y), 0)))).format(y)))
        plt.title(bf)

    # plt.savefig(pjoin(save_loc, name + '--' + x + '_vs_' + 'FOM' + '_22.png'))
    return fig


def compare_ToutfixedQfixed(data, T_out_fixed, T_max, Q_out_fixed, S_t,
                            L_t, D_t, D_h, A_c, A_s, visctype, visctype_var, save_loc, option_flow_type):
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

    visc_vars['phi'] = phi.copy()
    visc_vars['p'] = data['p'].copy()
    visc_vars['d_p'] = data['a_33'].copy()
    data['mu_nf_enhancement'] = getattr(a, option_visc)(**visc_vars)
    data['mu_nf'] = np.multiply(data['mu_f'], data['mu_nf_enhancement'])

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
        turbulent_idx = Reynolds_lam>2000
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
        'Velocity'] ** 2 * data[
                        'm_flowrate'] * data['rho_nf'] ** -1
    data['Reynolds'] = Reynolds
    data['Prandtl'] = Prandtl

    data['T_in'] = data['T_out'] - deltaT

    data['FOM Standard'] = np.divide(data['Q_out'], data['W_out'])
    data['FOM User'] = np.divide(
        np.multiply(data['c_nf'] ** (user_c_coeff), np.multiply(data['rho_nf'] ** (user_rho_coeff), data['k_nf'] ** (user_k_coeff))),
        data['mu_nf'] ** (user_mu_coeff))


    data['FOM Mouromtseff'] = data['rho_nf'] ** .8 * data['c_nf'] ** .33 * data['k_nf'] ** .67 / data['mu_nf'] ** .47
    data['FOM k/mu enhancement'] = np.divide( np.divide(data['k_nf'], data['k_f']),np.divide(data['mu_nf'], data['mu_f']))
    data['T_outcheck'] = T_max - data['Q_out'] / (A_s * data['HTC'])
    data['Shear Rate'] = 8 * np.divide(data['Velocity'], D_h)

    # look at phi = .10
    FOMS = ['FOM Standard','FOM Mouromtseff','FOM k/mu enhancement','FOM User']
    if option_FOM_norm:
        for iFOM in FOMS:
            data[iFOM]= data[iFOM] / data[data["phi"] == 0][iFOM][0]

    data['FOM Standard normed'] = data['FOM Standard'] / data[data["phi"] == 0]['FOM Standard'][0]
    # for each particle, grab highest performing

    name = 'T_out_fixed=' + str(T_out_fixed) + ', Qoutfixed=' + str(Q_out_fixed) + ', T_wall=' + str(
        T_max) + ', S_t=' + str(
        S_t)

    return data


@st.cache
def grab_data(file_path):
    nanoparticle_cand = pd.read_excel(file_path, sheet_name="nanoparticle", index_col=0, header=[0, 1], comment='#')
    basefluid_cand = pd.read_excel(file_path, sheet_name="basefluid", index_col=0, header=[0, 1], comment='#')
    nanoparticle_cand['id'] = nanoparticle_cand.index
    basefluid_cand['id'] = basefluid_cand.index
    nanoparticle_cand['np_id_s'] = nanoparticle_cand['id'].copy()
    basefluid_cand['bf_id_s'] = basefluid_cand['id'].copy()

    basefluid_cand['mu_f^-1'] = basefluid_cand['mu_f'] ** -1

    return nanoparticle_cand, basefluid_cand


@st.cache
def grab_data2(file_path):
    return pd.read_csv(file_path)


def get_table_download_link(df, filename):
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}.csv">Download csv file</a>'
    return href


def generate_colloids(nanoparticle_cand, basefluid_cand, par_space):
    parameter = []
    possible_vals = []

    for par in par_space:
        parameter.append(par[0])
        possible_vals.append(par[1])
    par_tests = list(itertools.product(*possible_vals))
    par_tests = np.vstack(par_tests)

    # nanoparticle_cand[par]
    combos_data = pd.DataFrame(data=par_tests, columns=parameter)
    np_pars = ['np_id_s', "k_p", "R_bd", "rho_p", "cv_p"]
    bf_pars = ['bf_id_s', "k_f", "cv_f", "rho_f", "mu_f"]
    for npp in np_pars:
        combos_data[npp] = nanoparticle_cand[npp].values[combos_data['np_id'].astype(int)]
    for bpp in bf_pars:
        combos_data[bpp] = basefluid_cand[bpp].values[combos_data['bf_id'].astype(int)]

    # density
    combos_data['rho_nf'] = density_nf(combos_data["phi"], combos_data["rho_p"], combos_data["rho_f"])

    # conductivity - def hc_EMTnan(k_p, k_f, phi, a_33, p, R_bd)
    combos_data['k_nf'] = np.multiply(combos_data["k_f"], hc_EMTnan(combos_data["k_p"], combos_data["k_f"],
                                                                    combos_data["phi"], combos_data["a_33"],
                                                                    combos_data["p"], combos_data["R_bd"]))

    # heat capacity - def shc_nf(phi, cv_p, cv_f, rho_p, rho_f)
    combos_data['c_nf'] = shc_nf(combos_data["phi"], combos_data["cv_p"], combos_data["cv_f"],
                                 combos_data["rho_p"], combos_data["rho_f"])
    return combos_data

st.sidebar.title('Colloid Space Exploration')
password = st.sidebar.text_input("Password:", value="")

if password=='notarealpwd':
    # Create a text element and let the reader know the data is loading.
    data_load_state = st.text('Loading data...')
    # Load 10,000 rows of data into the dataframe.
    nanoparticle_cand, basefluid_cand = grab_data('data/data_1.xlsx')
    # Notify the reader that the data was successfully loaded.
    data_load_state.text("Done! (using st.cache)")

    st.header('Available Nanoparticles and Basefluids')


# df = pd.DataFrame(np.random.randn(200, 3),
# columns=['a', 'b', 'c'])
# cc = alt.Chart(df).mark_circle().encode(x='a', y='b', size='c', color='c', tooltip=['a', 'b', 'c'])
# st.altair_chart(cc)
def review():
    st.subheader('Nanoparticle Candidates')
    "Table shows property values. Plot shows values normalized to the maximum."
    nanoparticle_cand

    parallelchart_np = alt.Chart(nanoparticle_cand.reset_index()).transform_window(
        index='count()'
    ).transform_fold(
        ['k_p', 'cv_p', 'rho_p']
    ).transform_joinaggregate(
        min='min(value)',
        max='max(value)',
        groupby=['key']
    ).transform_calculate(
        minmax_value=(datum.value - datum.min) / (datum.max - datum.min),
        mid=(datum.min + datum.max) / 2
    ).mark_line(size=3).encode(
        x='key:N',
        y='minmax_value:Q',
        tooltip=['id', 'k_p', 'cv_p', 'rho_p'],
        color='id:N',
        detail='index:N',
        opacity=alt.value(.5)
    ).properties(height=400,width=500).interactive()

    st.write(parallelchart_np)

    st.subheader('Basefluid Candidates')
    "Table shows property values. Plot shows values normalized to the maximum."
    basefluid_cand

    parallelchart_bf = alt.Chart(basefluid_cand.reset_index()).transform_window(
        index='count()'
    ).transform_fold(
        ['k_f', 'cv_f', 'rho_f', 'mu_f^-1']
    ).transform_joinaggregate(
        min='min(value)',
        max='max(value)',
        groupby=['key']
    ).transform_calculate(
        minmax_value=(datum.value - datum.min) / (datum.max - datum.min),
        mid=(datum.min + datum.max) / 2
    ).mark_line(size=3).encode(
        x='key:N',
        y='minmax_value:Q',
        tooltip=['id', 'k_f', 'cv_f', 'rho_f', 'mu_f^-1'],
        color='id:N',
        detail='index:N',
        opacity=alt.value(1)
    ).properties(height=200,width=500).interactive()
    st.write(parallelchart_bf)
    return


# aLMODEL/data/data_1.xlsx

if password=='notarealpwd':
    review()

# ____________________________________________________
# import pairings or generate them






if 1 == 0:
    st.sidebar.subheader('Configuration')
    config_options = ['None - define yourself', 'Preset 1', 'Preset 2']
    config_picked = st.sidebar.selectbox('Config will pick viscosity and FOM models', config_options)
if password=='notarealpwd':

    st.sidebar.subheader('Define Colloid Space')
    st.sidebar.text("Indicate the range for each parameter. A colloid space will be generated that covers that range.")

    option_phi = st.sidebar.slider('Volume fraction (%)', 0., 30., (0., 15.))

    option_aspectratio = st.sidebar.slider('Particle aspect ratio', 1., 50., (1., 10.))
    option_size = st.sidebar.slider('Particle size (nm)', 1., 100., (10., 100.))

    samp_res_phi = st.sidebar.slider('Volume fraction samples', 1, 50, 20)
    samp_res_ar = st.sidebar.slider('Particle aspect ratio samples', 1, 10, 3)
    samp_res_size = st.sidebar.slider('Particle size samples', 1, 10, 1)

    st.sidebar.text("Select basefluids to include")
    bfs = []
    default_bf = [True, True, True, True, False, False]
    for i, ii in enumerate(basefluid_cand['id'].unique()):
        bfs.append(st.sidebar.checkbox(ii, default_bf[i], key='bf' + ii))

    st.sidebar.text("Select nanoparticles to include")
    nps = []
    default_np = np.repeat(False, len(nanoparticle_cand['id'].unique()))
    default_np[3] = True
    default_np[7] = True
    for i, ii in enumerate(nanoparticle_cand['id'].unique()):
        nps.append(st.sidebar.checkbox(ii, default_np[i], key='np' + ii))

    bfs = [i for i, x in enumerate(bfs) if x]
    nps = [i for i, x in enumerate(nps) if x]

    par_space = [
        ['np_id', nps],
        ['bf_id', bfs],
        ['p', np.linspace(float(option_aspectratio[0]), float(option_aspectratio[1]), samp_res_ar)],
        ['phi', np.linspace(float(option_phi[0]) / 100., float(option_phi[1]) / 100., samp_res_phi)],
        ['a_33', np.linspace(float(option_size[0]), float(option_size[1]), samp_res_size)],
        ['R_bd', np.linspace(10 ** -8, 10 ** -8, 1)]]

    data_loc = pjoin('/Users/hanscastorp/Dropbox/SandBox2019/aLMODEL/data/200505_sweep_9')

    st.header('Generated Colloid Candidates')
    """
    Use the sidebar options to set the colloid parameter space to analyze. Ranges and sampling should be set 
    appropriately. The output colloid space is shown below.
    """
    if 1 == 1:
        colloid_combos = generate_colloids(nanoparticle_cand=nanoparticle_cand,
                                           basefluid_cand=basefluid_cand,
                                           par_space=par_space)
    else:
        colloid_combos = grab_data2(pjoin(data_loc, 'data_extract.csv'))
        colloid_combos['bf_id_s'] = colloid_combos['bf_id']
        colloid_combos['np_id_s'] = colloid_combos['np_id']
    colloid_combos

    colloid_combos['combo_name'] = colloid_combos['bf_id_s'] + '-' + colloid_combos['np_id_s']

    st.sidebar.subheader('Thermophysical Parameters')

    object_methods_visc = [method_name for method_name in dir(viscosity)
                           if callable(getattr(viscosity, method_name)) and not method_name.startswith('__')]

    option_visc = st.sidebar.selectbox('Viscosity Model', object_methods_visc, 1)

    var_assign = [f for f in inspect.getfullargspec(getattr(viscosity, option_visc))[0] if
                  not f == 'self' and not f == 'phi' and not f == 'p' and not f=='d_p']
    visc_vars = {}

    st.header('Thermophysical Property Models')

    """
    A variety of models are available to predict the thermophysical properties of nanofluids. The key properties to 
    consider are viscosity, thermal conducvity, specific heat capacity, and density. The latter two are straightforward. 
    Viscosity and thermal conductivity have a variety of model expressions reflecting various physical mechanisms and 
    empirical measurements. You can pick different models and see how various input parameters affect the properties as 
    well as the final Figure of Merit. For a reasonable range of input values, see (1) and (2). The FOM can be calcualted 
    in 
    various ways depending on the users particular optimization goals. We have collected a few FOMs including a standard 
    one for minimizing pumping power in laminar or turbulent flows. 
    
    
    """
    st.subheader('Intrinsic viscosities and maximum volume fraction from (1)')
    sample_data = grab_data2('data/kandlikar.csv')
    sample_data

    dlt_sldr_vals = {'phi_e': ['Effective Volume Factor', 'Factor', 0., 10., 4.],
                     'nu': ['Intrinsic Viscosity', 'Factor', 0., 50., .65],
                     'k_h': ['Higgins Constant', 'Factor', 0., 50., 1.],
                     'D': ['Fractal Index', 'D', 0., 3., 2.5],
                     'd_a': ['Agglomerate Dimension', 'nm', 1., 1000., 100.],
                     'phi_max': ['Max Volume Fraction', 'Fraction', 0., 1., .35],
                     }

    for d in var_assign:
        visc_vars[d] = float(st.sidebar.slider(dlt_sldr_vals[d][0] + ' (' + dlt_sldr_vals[d][1] + ')',
                                               dlt_sldr_vals[d][2],
                                               dlt_sldr_vals[d][3],
                                               dlt_sldr_vals[d][4]))
    a = viscosity()

    st.subheader('Figures of Merit')
    "We include a figure of merit based on the pumping. The unconstrained parameters in this FOM are the fluid bulk " \
    "temperature difference, and the pumping power. The comparison is useful because the heat load is fixed at the " \
    "relevant value, and the outlet temperature is fixed, as might be expected in a typical system design where the " \
    "ambient or chiller temps are dictated by other system considerations. Temperature is unconstrained because it " \
    "matters less than we think and it is relatively easy to achieve acceptable values in an axial flow geometry. " \
    "Acceptable solutions are those that achieve temperature change less than 4 degrees, and better solutions simply have " \
    "" \
    "" \
    "a lower pumping power. The analytical form for the FOM simplifies to a simple proportionality when dividing out " \
    "constant " \
    "quantities for comparison."
    st.latex(r''' FOM  = \frac{\dot{Q}}{\dot{W}} 
      = \frac{\dot{Q} 2 D_{h} D_{hp}}{162 \cdot  (\frac{P}{D_t} -1)^{.435} L_{T} A_{c}} \left( \frac{A_s (T_{max} - T_{
      out}) 0.128 \phi^.14 \psi }{\dot{Q} D_{hp}}\right)^\frac{10}{3} \left[ \frac{c^{\frac{4}{3}} \rho^2 k^2}{\mu^{
      \frac{5}{3}}} \right] ''')

    st.latex(r' \propto  \left[ \frac{c^{\frac{4}{3}} \rho^2 k^2}{\mu^{\frac{5}{3}}} \right]')


    "The literature includes several figures of merit like Mouromtseff's number or the ratio of conductivity to viscosity " \
    "" \
    "enhancement."
    st.latex(r'''M_{0} = \frac{\rho^{.8} k^{.67} c^{.33}}{\mu^{.47}}''')
    st.latex(r'''\frac{C_{k} }{C_{\mu}}>4''')

    "You can create a user specified figure of merit by changing the exponents on the thermophysical parameters in the sidebar."



    st.sidebar.text(inspect.getsource(getattr(a, option_visc)))

    option_cond = st.sidebar.selectbox('Conductivity Model', ["EMT-Nan"])

    st.sidebar.subheader('Figure of Merit')
    object_methods_FOM = ['FOM Standard', 'FOM User', 'FOM Mouromtseff', 'FOM k/mu enhancement']
    option_FOM = st.sidebar.selectbox('Pick a FOM to compare colloids by.', object_methods_FOM, 0)
    option_FOM_norm = st.sidebar.checkbox('Normalize to first basefluid', True)


    if option_FOM == 'FOM User':
        st.sidebar.text('Select FOM exponents for user specified FOM.')
        user_c_coeff = st.sidebar.number_input('Specific Heat Power', value=4 / 3)
        user_rho_coeff = st.sidebar.number_input('Density Power', value=2.)
        user_k_coeff = st.sidebar.number_input('Conductivity Power', value=2.)
        user_mu_coeff = st.sidebar.number_input('Viscosity Power', value=5 / 3)


        string_toshowFOMuser = r"FOM_{user} \propto \left[ \frac{c^{" + str(round(user_c_coeff, 2)) + r"} \rho^{" + str(
            round(user_rho_coeff, 2)) + r"} k^{" + str(round(
            user_k_coeff, 2)) + r"}}{\mu^{" + str(round(user_mu_coeff, 2)) + r"}} \right]"
    else:
        string_toshowFOMuser = 'NA'
        user_c_coeff = 4./3.
        user_rho_coeff = 2.
        user_k_coeff = 2.
        user_mu_coeff = 5./3
    FOM_strings = {'FOM Standard':r'FOM \propto  \left[ \frac{c^{\frac{4}{3}} \rho^2 k^2}{\mu^{\frac{5}{3}}} \right]',
                   'FOM User':string_toshowFOMuser,
                   'FOM Mouromtseff':r'''FOM = M_{0} = \frac{\rho^{.8} k^{.67} c^{.33}}{\mu^{.47}}''',
                   'FOM k/mu enhancement':r'''FOM = \frac{C_{k} }{C_{\mu}}>4'''}

    st.sidebar.latex(r'' + FOM_strings[option_FOM])

    st.sidebar.subheader('Battery Thermal and Geometry')
    Q_out_fixed = st.sidebar.slider('Cell Heat Generation (W)', 1., 10., 5.)
    T_max = st.sidebar.slider('Maximum Cell Temperature (C)', 30., 50., 40.)
    T_out_fixed = st.sidebar.slider('Fluid Outlet Temperature (C)', 20., T_max, 30.)
    option_geo = st.sidebar.selectbox('Geometry', ["Cylindrical Cell", "Square Channel"])
    S_t = st.sidebar.slider('Cell Separation (mm)', 0., 6., 2.) / 1000
    L_t = st.sidebar.slider('Cell Length (cm)', 1., 10., 6.510) / 100
    D_t = st.sidebar.slider('Cell Diameter (cm)', 1., 10., 1.825) / 100

    option_flow_type = st.sidebar.selectbox('Laminar or Turbulent Flow', ["Laminar", "Turbulent", "Auto"], 0)
    # TODO ADD THIS!!!

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

    st.header("Model Outputs")
    autocompute_toggle = st.checkbox('Auto compute?', False)
    plotToggle = st.checkbox('Plot Results?', True)
    newplot_toggle = st.checkbox('Additional Plotting?', False)

    if st.button('Press to Compute') or autocompute_toggle:

        data = compare_ToutfixedQfixed(data=colloid_combos,
                                       T_out_fixed=T_out_fixed, T_max=T_max,
                                       Q_out_fixed=Q_out_fixed,
                                       S_t=S_t, L_t=L_t, D_t=D_t, D_h=D_h, A_c=A_c, A_s=A_s, save_loc=data_loc,
                                       visctype=option_visc, visctype_var=visc_vars, option_flow_type=option_flow_type)
        st.subheader('Computed Performance')
        data = data.round({'phi': 3, 'a_33': 3, 'p': 3})
        data
        st.markdown(get_table_download_link(data, 'colloid_data'), unsafe_allow_html=True)

        unq_bf = np.unique(data['bf_id_s'])
        unq_np = np.unique(data['np_id_s'])
        unq_p = np.unique(data['p'])
        unq_a33 = np.unique(data['a_33'])

        if plotToggle:
            # fig1 = basic2d(data, 'phi', option_FOM, norm=1)
            # st.write(fig1)



            if 1 == 1:
                # peak perfromance
                data_peak = data[data['p'] == np.max(unq_p)]
                data_peak = data_peak[data_peak['a_33'] == np.max(unq_a33)]

                base_peak = alt.Chart(data_peak).mark_line(size=4).encode(
                    x=alt.X('phi', type='ordinal'),
                    y=alt.Y(option_FOM, type='quantitative',
                            scale=alt.Scale(type='log')),
                    facet=alt.Facet('bf_id_s:N', columns=len(unq_bf), spacing=-.5),
                    color='np_id_s:N',
                    # tooltip=['colloid name','p']
                ).properties(title="Max Performance (max aspect ratio and max size)", width=180).interactive()

                st.write(base_peak)

            if 1==1:
                base = alt.Chart(data).mark_circle(size=60).encode(
                    x=alt.X('phi', type='ordinal'),
                    y=alt.Y(option_FOM, type='quantitative',
                            scale=alt.Scale(type='log')),
                    facet=alt.Facet('bf_id_s:N', columns=len(unq_bf), spacing=-.5),
                    color='np_id_s:N',
                    size='a_33:O',

                ).properties(width=180).interactive()

                p_radio = alt.binding_radio(options=unq_p.tolist())
                p_select = alt.selection_single(fields=['p'], bind=p_radio, name="Aspect ratio, p")
                p_color_condition = alt.condition(p_select,
                                                  alt.Color('p:N', legend=None),
                                                  alt.value('lightgray'))

                np_dropdown = alt.binding_select(options=unq_np.tolist())
                np_select = alt.selection_single(fields=['np_id_s'], bind=np_dropdown, name="Nanoparticle")

                radio_p = base.add_selection(p_select
                                             ).encode(color=p_color_condition,
                                                      ).add_selection(np_select
                                                                      ).transform_filter(np_select
                                                                                         ).properties(
                    title="Select Aspect Ratio (p) and Nanoparticle")
                st.write(radio_p)


        if newplot_toggle:
            st.subheader('Extra Plots')
            "Pick variables to plot"
            plot_variable_y = st.selectbox('Y Variable', data.columns, 35)
            if 1==1:
                base = alt.Chart(data).mark_circle(size=60).encode(
                    x=alt.X('phi', type='ordinal'),
                    y=alt.Y(plot_variable_y, type='quantitative',
                            scale=alt.Scale(type='log')),
                    facet=alt.Facet('bf_id_s:N', columns=len(unq_bf), spacing=-.5),
                    color='np_id_s:N',
                    size='a_33:O',

                ).properties(width=180).interactive()

                p_radio = alt.binding_radio(options=unq_p.tolist())
                p_select = alt.selection_single(fields=['p'], bind=p_radio, name="Aspect ratio, p")
                p_color_condition = alt.condition(p_select,
                                                  alt.Color('p:N', legend=None),
                                                  alt.value('lightgray'))

                np_dropdown = alt.binding_select(options=unq_np.tolist())
                np_select = alt.selection_single(fields=['np_id_s'], bind=np_dropdown, name="Nanoparticle")

                radio_p = base.add_selection(p_select
                                             ).encode(color=p_color_condition,
                                                      ).add_selection(np_select
                                                                      ).transform_filter(np_select
                                                                                         ).properties(
                    title="Select Aspect Ratio (p) and Nanoparticle")
                st.write(radio_p)
            else:
                fig5 = basic2d(data, "phi", plot_variable_y, norm=st.checkbox('Normalize', False))
                st.write(fig5)


        # st.markdown(
        #     "Laminar Conditions (Re<2000) for all colloids? **" + str((data['Reynolds'].values < 2000).all()) + "**")
        data.to_csv('data/data_out.csv')
        # unique combos
        colloid_combos_uq = np.unique(data['combo_name'])

        uq_data = []
        for cc in colloid_combos_uq:
            loc_data = data[data['combo_name'] == cc].reset_index()
            idxmax = np.argmax(loc_data[option_FOM])
            dict1 = {}
            # get input row in dictionary format
            # key = col_name
            dict1.update({'Colloid':cc, option_FOM:np.around(loc_data[option_FOM][idxmax], 2),
                            'Relative FOM normed':np.around(loc_data['FOM Standard normed'][idxmax], 2), 'phi':np.around(loc_data['phi'][idxmax], 3),
                            'p':loc_data['p'][idxmax], 'a_33':loc_data['a_33'][idxmax]})

            uq_data.append(dict1)

        summary = pd.DataFrame(uq_data)

        summary = summary.sort_values(option_FOM,ascending=False)
        st.subheader('Colloid Max Performance Summary')
        summary
        st.markdown(get_table_download_link(summary, 'summary'), unsafe_allow_html=True)

        # TODO add calculated LOOP crap
        # st.subheader('Lubrizol Flow Loop Parameters')



    st.header("References")
    """
    1) Intrinsic viscosity and max volume fraction --> Ye, X.; Kandlikar, S. G.; Li, C. Viscosity of Nanofluids 
    Containing Anisotropic Particles: A Critical Review and a 
    Comprehensive Model. Eur. Phys. J. E 2019, 42 (12), 60–65. https://doi.org/10.1140/epje/i2019-11923-7. \n
    2) Higgins constant --> Pamies, R.; Hernández Cifre, J. G.; del Carmen López Martínez, M.; García de la Torre, 
    J. Determination of 
    Intrinsic Viscosities of Macromolecules and Nanoparticles. Comparison of Single-Point and Dilution Procedures. 
    Colloid Polym. Sci. 2008, 286 (11), 1223–1231. https://doi.org/10.1007/s00396-008-1902-2.
    """



    '-- Created at the Massachusetts Institute of Technology --'
    ' '
    ' '
    ' '
    'http://accessibility.mit.edu'

