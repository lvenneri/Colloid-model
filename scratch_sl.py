import altair as alt
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
from viscosity_models import viscosity, figure_of_merit
from colloid_combo_fun import *
import base64


@st.cache
def grab_data2(file_path):
    colloid_combos = pd.read_csv(file_path)
    # colloid_combos['bf_id_s'] = colloid_combos['bf_id']
    # colloid_combos['np_id_s'] = colloid_combos['np_id']
    colloid_combos['phi']= np.around( colloid_combos['phi'],3)
    colloid_combos['a_33']= np.around( colloid_combos['a_33'],3)
    colloid_combos['p']= np.around( colloid_combos['p'],3)

    return colloid_combos


data_loc = pjoin('/Users/hanscastorp/Downloads/colloid_data (1).csv')

colloid_combos = grab_data2(pjoin(data_loc))
unq_bf = np.unique(colloid_combos['bf_id_s'])
unq_np = np.unique(colloid_combos['np_id_s'])
unq_p = np.unique(colloid_combos['p'])
unq_a33 = np.unique(colloid_combos['a_33'])



if 1==1:
    # peak perfromance
    colloid_combos_peak = colloid_combos[colloid_combos['p']==np.max(unq_p)]
    colloid_combos_peak = colloid_combos_peak[colloid_combos_peak['a_33']==np.max(unq_a33)]

    base_peak = alt.Chart(colloid_combos_peak).mark_line(size=4).encode(
        x=alt.X('phi', type='ordinal'),
        y=alt.Y('FOM Standard', type='quantitative',
                scale=alt.Scale(type='log')),
        facet=alt.Facet('bf_id_s:N', columns=len(unq_bf),spacing=-.5),
        color='np_id_s:N',
        # tooltip=['colloid name','p']
    ).properties(title="Max Performance (max aspect ratio and max size)",width=180).interactive()

    st.write(base_peak)

if 1==1:

    base = alt.Chart(colloid_combos).mark_circle(size=60).encode(
        x=alt.Y('phi',type='ordinal'),
        y=alt.Y('FOM Standard',type='quantitative',
          scale=alt.Scale(type='log')),
        facet=alt.Facet('bf_id_s:N', columns=len(unq_bf),spacing=-.5),
        color='np_id_s:N',
        size = 'a_33:O',

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
    ).add_selection( np_select
    ).transform_filter(np_select
    ).properties(title="Select Aspect Ratio (p) and Nanoparticle")
    st.write(radio_p)


if 1==0:

    source = colloid_combos.copy()

    # Create a selection that chooses the nearest point & selects based on x-value
    nearest = alt.selection(type='single', nearest=True, on='mouseover',
                            fields=['phi'], empty='none')

    # The basic line
    line = alt.Chart(source).mark_line(interpolate='basis').encode(
        x='phi:Q',
        y='FOM Standard:Q',
        color='category:N'
    )

    # Transparent selectors across the chart. This is what tells us
    # the x-value of the cursor
    selectors = alt.Chart(source).mark_point().encode(
        x='x:Q',
        opacity=alt.value(0),
    ).add_selection(
        nearest
    )

    # Draw points on the line, and highlight based on selection
    points = line.mark_point().encode(
        opacity=alt.condition(nearest, alt.value(1), alt.value(0))
    )

    # Draw text labels near the points, and highlight based on selection
    text = line.mark_text(align='left', dx=5, dy=-5).encode(
        text=alt.condition(nearest, 'y:Q', alt.value(' '))
    )

    # Draw a rule at the location of the selection
    rules = alt.Chart(source).mark_rule(color='gray').encode(
        x='x:Q',
    ).transform_filter(
        nearest
    )

    # Put the five layers into a chart and bind the data
    aa = alt.layer(
        line, selectors, points, rules, text
    ).properties(
        width=600, height=300
    )
    st.write(aa)

if 1==0:
    # # A slider filter
    p_slider = alt.binding_range(min=0, max=10, step=1)
    slider_selection = alt.selection_single(bind=p_slider, fields=['p'], name="p_")

    filter_p = base.add_selection(
        slider_selection
    ).transform_filter(
        slider_selection
    ).properties(title="Slider Filtering")
    st.write(filter_p)
