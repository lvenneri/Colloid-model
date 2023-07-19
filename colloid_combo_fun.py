import numpy as np

"""
functions for use in colloid_model.py

"""

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


    L_11 = np.divide(p ** 2, 2 * (p ** 2 - 1)) - np.divide(p, 2 * (p ** 2 - 1) ** 1.5) * np.arccosh(p)
    L_33 = 1 - 2 * L_11

    gamma = (2 + 1 / p) * R_bd * k_f / (a_11 / 2)

    k_11 = np.divide(k_p, 1 + gamma * L_11 * k_p / k_f)
    k_33 = np.divide(k_p, 1 + gamma * L_33 * k_p / k_f)
    B_11 = np.divide(k_11 - k_f, k_f + L_11 * (k_11 - k_f))
    B_33 = np.divide(k_33 - k_f, k_f + L_33 * (k_33 - k_f))

    enhancement_fraction = np.divide(3 + phi * (2 * B_11 * (1 - L_11) + B_33 * (1 - L_33)),
                                     3 - phi * (2 * B_11 * L_11 + B_33 * L_33))

    p_lt_1 = np.argwhere(p.values <= 1).flatten()
    enhancement_fraction[p_lt_1] = hc_EMTclassical(k_p[p_lt_1], k_f[p_lt_1], phi[p_lt_1])

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


