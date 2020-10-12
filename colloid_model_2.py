import streamlit as st
import numpy as np
import pandas as pd
import time
import itertools
import altair as alt
import matplotlib.pyplot as plt
from altair import datum
from os.path import join as pjoin
import matplotlib.cm as cm
import inspect
from viscosity_models import viscosity
from colloid_combo_fun import *

"""
# Colloid Model for Lubrizol Electric Vehicle Battery Cooling Application
This script allows the user to look at the colloid design space and compare colloids in different application 
environments.
"""


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
    fluids_unq = fluids_unq[[0, 1, 3, 2]]

    numfluids = len(fluids_unq)
    for bfi, bf in enumerate(fluids_unq):

        plt.subplot(1, numfluids, bfi + 1)
        bf_data = data[data['bf_id_s'] == bf]
        for npi, nanop in enumerate(particles_unq):
            np_data = bf_data[bf_data['np_id_s'] == nanop]

            for pi, pppp in enumerate(p_unq):
                pp_data = np_data[np_data['p'] == pppp]
                if pi == 0:
                    plt.plot(pp_data[x], pp_data[y] / normer, lw=2, c=colors[npi], label=nanop.replace('_',' '))
                else:
                    plt.plot(pp_data[x], pp_data[y] / normer, lw=2, c=colors[npi], label='_nolegend_')

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


def compare_ToutfixedQfixed(data, T_out_fixed, T_max, Q_out_fixed, S_t, phi_compare,
                            L_t, D_t, D_h, A_c, A_s, visctype, visctype_var, save_loc):
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

    phi = data['phi']

    a = viscosity()

    visc_vars['phi'] = phi.copy()
    visc_vars['p'] = data['p'].copy()
    visc_vars['d_p'] = data['a_33'].copy()
    data['mu_nf'] = getattr(a, option_visc)(**visc_vars)

    # data['mu_nf'] = 1 + 2.5 * phi * 9 + 6.2 * (phi * 9) ** 2.0

    data['T_out'] = T_out_fixed
    mu_ratio = 1

    data['Velocity'] = (Q_out_fixed / (A_s * (T_max - T_out_fixed)) * (
            .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(
        np.multiply((np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])) ** .4,
                    (np.divide(data['rho_nf'], data['mu_nf']) * D_h) ** .6),
        data['k_nf']) / D_h) ** -1) ** (5 / 3)

    deltaT = np.divide(Q_out_fixed, (data['c_nf'] * data['rho_nf'] * A_c * data['Velocity']))
    data['deltaT'] = deltaT

    Reynolds = np.divide(data['rho_nf'], data['mu_nf']) * data['Velocity'] * D_h
    Prandtl = np.divide(np.multiply(data['c_nf'], data['mu_nf']), data['k_nf'])

    data['m_flowrate'] = data['rho_nf'] * A_c * data['Velocity']

    data['htc'] = .128 * mu_ratio ** .14 * data['psi_c'] * np.multiply(np.multiply(Prandtl ** .4, Reynolds ** .6),
                                                                       data['k_nf']) / D_h
    data['Q_out'] = data['c_nf'] * data['m_flowrate'] * (deltaT)

    data['Q_out_check'] = A_s * data['htc'] * (T_max - T_out_fixed)
    data['W_out'] = 162 * (PoD - 1) ** .435 * (Reynolds) ** -1 * data['rho_nf'] * (L_t / (2 * D_hp)) * data[
        'Velocity'] ** 2 * data[
                        'm_flowrate'] * data['rho_nf'] ** -1
    data['Reynolds'] = Reynolds
    data['Prandtl'] = Prandtl

    data['T_in'] = data['T_out'] - deltaT

    data['FOM Q/W'] = np.divide(data['Q_out'], data['W_out'])
    data['FOM Simplified'] = np.divide(data['c_nf'] ** (4 / 3) * data['rho_nf'] ** (2) * data['k_nf'] ** (2),
                                       data['mu_nf'] ** (5 / 3))

    data['FOM Mouromtseff'] = data['rho_nf'] ** .8 * data['c_nf'] ** .33 * data['k_nf'] ** .67 / data['mu_nf'] ** .47
    data['Cmu : Ck'] = np.divide(np.divide(data['mu_nf'], data['mu_f']), np.divide(data['k_nf'], data['k_f']))
    data['T_outcheck'] = T_max - data['Q_out'] / (A_s * data['htc'])
    data['Shear Rate'] = 8 * np.divide(data['Velocity'], D_h)

    # look at phi = .10
    np.argmin(data['phi'] - phi_compare)

    # for each particle, grab highest performing

    name = 'T_out_fixed=' + str(T_out_fixed) + ', Qoutfixed=' + str(Q_out_fixed) + ', T_wall=' + str(
        T_max) + ', S_t=' + str(
        S_t)

    fig = basic2d(data, 'phi', 'FOM Q/W', save_loc=save_loc, name='FOM -- ' + name, norm=1)

    return data, fig


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
    combos_data['k_nf'] = np.multiply(combos_data["k_f"],hc_EMTnan(combos_data["k_p"], combos_data["k_f"],
                                                                   combos_data["phi"], combos_data["a_33"],
                                                                   combos_data["p"], combos_data["R_bd"]))

    # heat capacity - def shc_nf(phi, cv_p, cv_f, rho_p, rho_f)
    combos_data['c_nf'] = shc_nf(combos_data["phi"], combos_data["cv_p"], combos_data["cv_f"],
                                 combos_data["rho_p"], combos_data["rho_f"])
    combos_data['combo_name'] = combos_data['bf_id_s'] +'-' +combos_data['np_id_s']
    return combos_data


# Create a text element and let the reader know the data is loading.
data_load_state = st.text('Loading data...')
# Load 10,000 rows of data into the dataframe.
nanoparticle_cand, basefluid_cand = grab_data('data/data_1.xlsx')
# Notify the reader that the data was successfully loaded.
data_load_state.text("Done! (using st.cache)")

# df = pd.DataFrame(np.random.randn(200, 3),
# columns=['a', 'b', 'c'])
# cc = alt.Chart(df).mark_circle().encode(x='a', y='b', size='c', color='c', tooltip=['a', 'b', 'c'])
# st.altair_chart(cc)
def review():
    st.subheader('Nanoparticle Candidates')
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
    ).mark_line().encode(
        x='key:N',
        y='minmax_value:Q',
        tooltip=['id', 'k_p', 'cv_p', 'rho_p'],
        color='id:N',
        detail='index:N',
        opacity=alt.value(.5)
    ).properties(height=600).interactive()

    st.altair_chart(parallelchart_np)

    st.subheader('Basefluid Candidates')
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
    ).mark_line().encode(
        x='key:N',
        y='minmax_value:Q',
        tooltip=['id', 'k_f', 'cv_f', 'rho_f', 'mu_f^-1'],
        color='id:N',
        detail='index:N',
        opacity=alt.value(1)
    ).properties().interactive()

    st.altair_chart(parallelchart_bf)
    return


review()

# ____________________________________________________
# import pairings or generate them

# aLMODEL/data/data_1.xlsx
st.sidebar.title('Colloid Space Exploration')
# fixed space to look at
# if 1==0:
#     colloid_cand = pd.read_csv("/Users/hanscastorp/Dropbox/SandBox2019/aLMODEL/data/200505_sweep_9/data_extract.csv")
#     colloid_cand
#
#     if modality == 'Single Colloid':
#         st.sidebar.subheader('Create Colloid')
#         option_bf = st.sidebar.selectbox(
#             'Basefluid',
#             colloid_cand['bf_id'].unique())
#
#         option_np = st.sidebar.selectbox(
#             'Nanoparticle',
#             colloid_cand['np_id'].unique())
#
#         option_aspectratio = st.sidebar.slider('Particle aspect ratio', 1., 50., 10.)
#         option_size = st.sidebar.slider('Particle size (nm)', 1., 100., 50.)
#         option_phi = st.sidebar.slider('Volume fraction (%)', 0., 30., 10.) / 100

st.subheader('Generated Colloid Candidates')

st.sidebar.subheader('Define Colloid Space')
st.sidebar.text("Indicate the range for each parameter")

samp_res = st.sidebar.slider('Sampling Resolution', 1, 5, 3)
option_phi = st.sidebar.slider('Volume fraction (%)', 0., 30., (0., 15.))
option_aspectratio = st.sidebar.slider('Particle aspect ratio', 1., 50., (1., 10.))
option_size = st.sidebar.slider('Particle size (nm)', 1., 100., (50., 100.))

st.sidebar.text("Select basefluids to include")
bfs = []
default_bf = [True,True,True,True,False,False]
for i,ii in enumerate(basefluid_cand['id'].unique()):
    bfs.append(st.sidebar.checkbox(ii,default_bf[i]))

st.sidebar.text("Select nanofluids to include")
nps = []
default_np = np.repeat(False,len(nanoparticle_cand['id'].unique()))
default_np[3] = True
default_np[7] = True
for i,ii in enumerate(nanoparticle_cand['id'].unique()):
    nps.append(st.sidebar.checkbox(ii,default_np[i]))

bfs = [i for i, x in enumerate(bfs) if x]
nps = [i for i, x in enumerate(nps) if x]

par_space = [
    ['np_id', nps],
    ['bf_id', bfs],
    ['p', np.linspace(float(option_aspectratio[0]), float(option_aspectratio[1]), samp_res)],
    ['phi', np.linspace(float(option_phi[0])/100., float(option_phi[1])/100., 20)],
    ['a_33', np.linspace(float(option_size[0]), float(option_size[1]), 1)],
    ['R_bd', np.linspace(10 ** -8, 10 ** -8, 1)]]

colloid_combos = generate_colloids(nanoparticle_cand=nanoparticle_cand,
                                   basefluid_cand=basefluid_cand,
                                   par_space=par_space)

data_loc = pjoin('/Users/hanscastorp/Dropbox/SandBox2019/aLMODEL/data/200505_sweep_9')

@st.cache
def grab_data2(file_path):
    return pd.read_csv(file_path)


# colloid_combos = grab_data2(pjoin(data_loc, 'data_extract.csv'))

st.sidebar.header('Model Parameters')
st.sidebar.subheader('Fluid Models')

object_methods = [method_name for method_name in dir(viscosity)
                  if callable(getattr(viscosity, method_name)) and not method_name.startswith('__')]

option_visc = st.sidebar.selectbox('Visocsity Model', object_methods)

var_assign = [f for f in inspect.getfullargspec(getattr(viscosity, option_visc))[0] if
              not f == 'self' and not f == 'phi' and not f == 'p']
visc_vars = {}

"""
For reasonable range of values and decent starting values, see \n
1) Intrinsic viscosity and max volume fraction --> Ye, X.; Kandlikar, S. G.; Li, C. Viscosity of Nanofluids 
Containing Anisotropic Particles: A Critical Review and a 
Comprehensive Model. Eur. Phys. J. E 2019, 42 (12), 60–65. https://doi.org/10.1140/epje/i2019-11923-7. \n
2) Higgins --> Pamies, R.; Hernández Cifre, J. G.; del Carmen López Martínez, M.; García de la Torre, 
J. Determination of 
Intrinsic Viscosities of Macromolecules and Nanoparticles. Comparison of Single-Point and Dilution Procedures. 
Colloid Polym. Sci. 2008, 286 (11), 1223–1231. https://doi.org/10.1007/s00396-008-1902-2.
"""
dlt_sldr_vals = {'phi_e': ['Effective Volume Factor', 'Factor', 0., 10., 1.],
                 'nu': ['Intrinsic Viscosity', 'Factor', 0., 50., .65],
                 'k_h': ['Higgins Constant', 'Factor', 0., 50., 1.],
                 'D': ['Fractal Index', 'D', 0., 3., 2.5],
                 'd_a': ['Agglomerate Dimension', 'nm', 1., 1000., 100.],
                 'd_p': ['Particle Dimension', 'nm', 1., 100., 50.],
                 'phi_max': ['Max Volume Fraction', 'Fraction', 0., 1., .35],
                 }

for d in var_assign:
    visc_vars[d] = float(st.sidebar.slider(dlt_sldr_vals[d][0] + ' (' + dlt_sldr_vals[d][1] + ')',
                                           dlt_sldr_vals[d][2],
                                           dlt_sldr_vals[d][3],
                                           dlt_sldr_vals[d][4]))
a = viscosity()
st.sidebar.text(inspect.getsource(getattr(a, option_visc)))

option_cond = st.sidebar.selectbox('Conductivity Model', ["EMT-Nan", "EMT-Nan Cluster", "EMT"])

st.sidebar.subheader('Battery Thermal and Geometry')
Q_out_fixed = st.sidebar.slider('Cell Heat Generation (W)', 1., 10., 5.)
T_max = st.sidebar.slider('Maximum Cell Temperature (C)', 30., 50., 40.)
T_out_fixed = st.sidebar.slider('Fluid Outlet Temperature (C)', 20., T_max, 30.)
option_geo = st.sidebar.selectbox('Geometry', ["Cylindrical Cell", "Square Channel"])
S_t = st.sidebar.slider('Cell Separation (mm)', 0., 6., 2.) / 1000
L_t = st.sidebar.slider('Cell Length (cm)', 1., 10., 6.510) / 100
D_t = st.sidebar.slider('Cell Diameter (cm)', 1., 10., 1.825) / 100

Wetted_p = np.pi * D_t / 2
A_c = np.sqrt(3) / 4 * (D_t + S_t) ** 2 - 0.5 * np.pi * (D_t / 2) ** 2
A_s = np.pi * D_t * L_t / 2
D_h = 4 * A_c / Wetted_p

if option_geo == "Square Channel":
    Wetted_p = 4 * D_h
    A_c = D_h ** 2
    A_s = D_h * L_t



st.subheader('Computer Geometry Parameters')
"Wetted Perimeter (cm):", '{:.2f}'.format(Wetted_p * 100)
"Hydraulic Diameter (cm):", '{:.2f}'.format(D_h * 100)
"Cross Sectional Area (cm^2):", '{:.2f}'.format(A_c * 100 ** 2)
"Surface Area (cm^2):", '{:.2f}'.format(A_s * 100 ** 2)

if 1 == 1:
    data, fig = compare_ToutfixedQfixed(data=colloid_combos,
                                        T_out_fixed=T_out_fixed, T_max=T_max,
                                        Q_out_fixed=Q_out_fixed,
                                        phi_compare=.1,
                                        S_t=S_t, L_t=L_t, D_t=D_t, D_h=D_h, A_c=A_c, A_s=A_s, save_loc=data_loc,
                                        visctype=option_visc, visctype_var=visc_vars)
    st.subheader('Computed Performance')
    data
    st.write(fig)

# unique combos
colloid_combos_uq =np.unique(data['combo_name'])
uq_data = []
for cc in colloid_combos_uq:
    loc_data = colloid_combos[data['combo_name']==cc].reset_index()

    idxmax = np.argmax(loc_data['FOM Simplified'])


    uq_data.append([cc,np.around(loc_data['FOM Q/W'][idxmax],2),np.around(loc_data['phi'][idxmax],3),loc_data['p'][idxmax],loc_data['a_33'][idxmax]])


summary = pd.DataFrame(data=np.vstack(uq_data), columns=['Colloid','FOM','phi','p','a_33'])

st.subheader('Colloid Summary Max Performance')
summary



# TODO add a laminar vs Turbulent check

# TODO add calculated LOOP crap
# st.subheader('Lubrizol Flow Loop Parameters')
