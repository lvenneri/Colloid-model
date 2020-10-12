
import streamlit as st
# To make things easier later, we're also importing numpy and pandas for
# working with sample data.
import numpy as np
import pandas as pd
import time
import itertools
from astropy.table import Table, Column, MaskedColumn, vstack
from astropy.io import ascii
import altair as alt
import matplotlib.pyplot as plt


def v_nf(phi, type):
    """
    Several viscosity models are available. See - Ye, X.; Kandlikar, S. G.; Li, C. Viscosity of Nanofluids Containing
    Anisotropic Particles: A Critical Review and a Comprehensive Model. Eur. Phys. J. E 2019, 42 (12),
    60–65. https://doi.org/10.1140/epje/i2019-11923-7.


    :param phi: volumetric fraction
    :return:  viscosity enhancement
    """

    if type == "batchelor":
        visc = 1 + 2.5 * phi + 6.2 * phi ** 2.0

    return visc



print(v_nf)
