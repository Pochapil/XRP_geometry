import numpy as np
import scipy.integrate
from scipy import interpolate

import config
import newService


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


def calc_scatter_L(surface, L_x, cos_tensor, tau_scatter_matrix=None):
    tilda_s = create_ds_for_integral(surface)
    d = surface.surf_R_e * np.sin(surface.theta_range) ** 2  # distance to area
    coeff = 1 / (4 * np.pi * d ** 2)
    if tau_scatter_matrix is not None:
        coeff = coeff * (1 - np.exp(-tau_scatter_matrix))  # УЧЕСТЬ tau в отражаемой точке!!!!!!!!! np.abs ???
        # z = (1 - np.exp(-tau_scatter_matrix)) tau = 45-50 ????
    L = L_x * np.abs(scipy.integrate.simps(scipy.integrate.simps(cos_tensor * tilda_s * coeff, surface.theta_range),
                                           surface.phi_range))
    return L


def calculate_total_luminosity(surface, T_eff):
    # emission power
    tilda_s = create_ds_for_integral(surface)
    tensor = np.ones((config.N_phi_accretion, config.N_theta_accretion))

    emission_power = (4 * config.sigm_Stf_Bolc * scipy.integrate.simps(
        scipy.integrate.simps(T_eff ** 4 * tensor * tilda_s, surface.theta_range), surface.phi_range))
    return emission_power


def calc_L_nu(surface, T_eff, cos_tensor):
    ''' распределение L_nu от фазы на какой-то энергии излучения '''

    L_nu = np.empty(config.N_energy, dtype=object)

    tilda_s = create_ds_for_integral(surface)
    for i, energy in enumerate(config.energy_arr):
        freq = newService.get_frequency_from_energy(energy)
        plank_func = newService.plank_energy_on_frequency(freq, T_eff)

        L_nu[i] = 4 * np.pi * np.abs(scipy.integrate.simps(
            scipy.integrate.simps(plank_func * cos_tensor * tilda_s, surface.theta_range), surface.phi_range))

    return L_nu


# for energy_index in range(config.N_energy - 1):
#     current_energy_min = config.energy_arr[energy_index]
#     current_energy_max = config.energy_arr[energy_index + 1]


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


# def calc_scatter_L():
#     ...


# save
# calc nu L nu

def calc_scatter_L_nu():
    ...
# save
# calc nu L nu

# tij -> ti -> t
# scipy.integrate(axis=-1) --- можно интегрировать по тензору.
