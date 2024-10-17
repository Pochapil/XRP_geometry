import numpy as np
import scipy.integrate
from scipy import interpolate

import config


def create_ds_for_integral(surface):
    # это dl phi * dl theta но без d theta d phi !!! S с волной ~S в моем дипломе
    # для интеграла по simpson
    # dS единичная площадка при интегрировании
    # ds = tilda_s dtheta dphi
    dl = surface.surf_R_e * ((3 * np.cos(surface.theta_range) ** 2 + 1) ** (1 / 2)) * np.sin(surface.theta_range)
    dphi = surface.surf_R_e * (np.sin(surface.theta_range) ** 3)
    tilda_s = dl * dphi
    return tilda_s



def calc_L(surface, T_eff, cos_tensor):
    tilda_s = create_ds_for_integral(surface)
    L = np.abs(4 * config.sigm_Stf_Bolc * scipy.integrate.simps(
        scipy.integrate.simps(T_eff ** 4 * cos_tensor * tilda_s, surface.theta_range), surface.phi_range))
    return L


def calc_L_nu():
    ...


# save


# shadows
# calc shadowed_matrix (ns + columns)
# tau
# calc tau_matrix

# УЧЕСТЬ tau в отражаемой точке!!!!!!!!!

'''
np.array([np.sin(theta_obs) * np.cos(phi),
          np.sin(theta_obs) * np.sin(phi),
          np.cos(theta_obs)]) - типо куда луч стреляем. нужно для нахождения tau1
'''


def calc_scatter_L():
    ...


# save
# calc nu L nu

def calc_scatter_L_nu():
    ...
# save
# calc nu L nu

# tij -> ti -> t
# scipy.integrate(axis=-1) --- можно интегрировать по тензору.
