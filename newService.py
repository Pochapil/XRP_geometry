import numpy as np

import config


def get_delta_distance_at_surface_NS(R_e):
    # R=R_e * sin_theta ** 2
    theta = get_theta_accretion_begin(R_e)
    # print(theta)
    return get_delta_distance(theta, R_e)


def get_A_normal_at_surface_NS(R_e, a_portion):
    theta = get_theta_accretion_begin(R_e)
    # print(theta)
    # a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ
    return get_A_normal(theta, R_e, a_portion)


# формула 2 в статье
# the distance between the field lines δ
# ширина аккреции на поверхности м
def get_delta_distance(theta, R_e):
    # R=R_e * sin_theta ** 2
    return R_e * np.sin(theta) ** 3 / (1 + 3 * np.cos(theta) ** 2) ** (1 / 2) * config.dRe_div_Re


# формула 3 в статье
# площадь поперечного сечения аккреционного потока
# a cross-section
def get_A_normal(theta, R_e, a_portion):
    # a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ
    return 2 * get_delta_distance(theta, R_e) * 2 * np.pi * a_portion * R_e * np.sin(theta) ** 3


# из усл силовой линии МП : r = R_e sin**2
def get_theta_accretion_begin(R_e):
    return np.arcsin((config.R_ns / R_e) ** (1 / 2))


def get_free_fall_velocity(theta, R_e):
    r = R_e * np.sin(theta) ** 2
    return (2 * config.G * config.M_ns / r) ** (1 / 2)


def get_e_obs(i_angle, phi_angle):
    e_obs = np.array([np.sin(i_angle) * np.cos(phi_angle),
                      np.sin(i_angle) * np.sin(phi_angle),
                      np.cos(i_angle)])
    return e_obs


def get_pulsed_fraction(arr):
    min_value = min(arr)
    max_value = max(arr)
    if max_value == 0:
        return 0
    return (max_value - min_value) / (max_value + min_value)


def get_PF(L):
    # берем по последней размерности - это фаза (ожидаем что поступают данные в таком виде)
    min_val = np.min(L, axis=-1)
    max_val = np.max(L, axis=-1)
    return (max_val - min_val) / (max_val + min_val)


# нужно ли на пи?
def plank_energy_on_wavelength(wavelength, T):
    return 2 * config.h_plank_ergs * config.c ** 2 / wavelength ** 5 \
           * 1 / (np.e ** (config.h_plank_ergs * config.c / (wavelength * config.k_bolc * T)) - 1)


def plank_energy_on_frequency(frequency, T):
    # erg / s / sr / cm**2 / hz
    return 2 * config.h_plank_ergs * frequency ** 3 / config.c ** 2 \
           * 1 / (np.e ** (config.h_plank_ergs * frequency / (config.k_bolc * T)) - 1)


def get_frequency_from_energy(energy):
    coefficient = 1000  # КэВ а не эВ
    frequency = coefficient * energy / config.h_plank_evs  # E = h f
    return frequency  # гц!


def extend_arr_for_plot(arr):
    # нужно расширить массивы, чтобы покрыть фазу [0,2]
    append_index = config.N_phase_for_plot - config.N_phase
    array_to_plot = np.append(arr, arr[0:append_index])
    return array_to_plot
