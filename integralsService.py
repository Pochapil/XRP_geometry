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


def calculate_total_luminosity(surface, T_eff):
    # emission power = integral without cos
    tensor = np.ones((config.N_phi_accretion, config.N_theta_accretion))
    return calc_L(surface, T_eff, tensor)

    # emission_intensity_inner = calc_L(emission_column.inner_surface, T_eff, tensor)
    # emission_intensity_outer = calc_L(emission_column.outer_surface, T_eff, tensor)
    # emission_intensity = 1 / 2 * (emission_intensity_outer + emission_intensity_inner)

    # tilda_s = create_ds_for_integral(surface)
    # emission_power = (4 * config.sigm_Stf_Bolc * scipy.integrate.simps(
    #     scipy.integrate.simps(T_eff ** 4 * tensor * tilda_s, surface.theta_range), surface.phi_range))
    #
    # return emission_power


def calc_L(surface, T_eff, cos_tensor):
    tilda_s = create_ds_for_integral(surface)
    # (пи входит в сигма) sigm_Stf_Bolc = integral pi B_nu
    L = np.abs(4 * config.sigm_Stf_Bolc * scipy.integrate.simps(
        scipy.integrate.simps(T_eff ** 4 * cos_tensor * tilda_s, surface.theta_range), surface.phi_range))
    return L


def calc_L_nu(surface, T_eff, cos_tensor):
    ''' распределение L_nu от фазы на какой-то энергии излучения '''

    L_nu = np.empty((config.N_energy, config.N_phase))

    tilda_s = create_ds_for_integral(surface)
    for i, energy in enumerate(config.energy_arr):
        freq = newService.get_frequency_from_energy(energy)
        plank_func = newService.plank_energy_on_frequency(freq, T_eff)

        L_nu[i] = 4 * np.pi * np.abs(scipy.integrate.simps(
            scipy.integrate.simps(plank_func * cos_tensor * tilda_s, surface.theta_range), surface.phi_range))

    return L_nu


def calc_scatter_L(surface, L_x, cos_tensor, tau_scatter_matrix=None):
    # расчет отраженной светимости
    tilda_s = create_ds_for_integral(surface)
    d = surface.surf_R_e * np.sin(surface.theta_range) ** 2  # distance to area
    coeff = 1 / (4 * np.pi * d ** 2)
    if tau_scatter_matrix is not None:
        coeff = coeff * (1 - np.exp(-np.abs(tau_scatter_matrix)))  # УЧЕСТЬ tau в отражаемой точке!!!!!!!!!
        # np.abs ??? хотя вроде бы косинус между нормалью и на НЗ не может быть отрицательным
        # z = (1 - np.exp(-tau_scatter_matrix)) tau = 45-50 ????
    L = L_x * np.abs(scipy.integrate.simps(scipy.integrate.simps(cos_tensor * tilda_s * coeff, surface.theta_range),
                                           surface.phi_range))
    return L


def calc_scatter_L_nu(magnet_surface, emission_surface, T_eff, cos_tensor, tau_scatter_matrix=None):
    ''' распределение L_nu от фазы на какой-то энергии излучения
    1/2 - так как при расчете L_nu ы брали 4 pi чтобы из потока F_nu получить светимость L_nu (предполагая
    изотропность источника). Теперь же надо посчитать поток F_nu для рассеяния и отражения - поэтому так как
    площадка излучает в полуплоскость надо взять 1/2 от интеграла

    поток, который должен быть отражен на самом деле в 2 раза меньше, так как излучает только в полуплоскость,
    а мы считаем L_nu как изотропное и поэтому умножали на 4 pi. надо на 2 pi'''

    tensor = np.ones((config.N_phi_accretion, config.N_theta_accretion))
    emission_intensity = calc_L_nu(emission_surface, T_eff, tensor)

    # emission_intensity_inner = calc_L_nu(emission_column.inner_surface, T_eff, tensor)
    # emission_intensity_outer = calc_L_nu(emission_column.outer_surface, T_eff, tensor)
    # emission_intensity = 1 / 2 * (emission_intensity_outer + emission_intensity_inner)

    tilda_s = create_ds_for_integral(magnet_surface)
    d = magnet_surface.surf_R_e * np.sin(magnet_surface.theta_range) ** 2  # distance to area
    coeff = 1 / (4 * np.pi * d ** 2)
    if tau_scatter_matrix is not None:
        coeff = coeff * (1 - np.exp(-np.abs(tau_scatter_matrix)))  # УЧЕСТЬ tau в отражаемой точке!!!!!!!!!
        # np.abs ??? хотя вроде бы косинус между нормалью и на НЗ не может быть отрицательным
        # z = (1 - np.exp(-tau_scatter_matrix)) tau = 45-50 ????
    L_nu = emission_intensity * np.abs(
        scipy.integrate.simps(scipy.integrate.simps(cos_tensor * tilda_s * coeff, magnet_surface.theta_range),
                              magnet_surface.phi_range))
    return L_nu
