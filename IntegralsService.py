import numpy as np

import config


def create_ds_for_integral(surf):
    # это dl phi * dl theta но без d theta d phi !!! S с волной ~S в моем дипломе
    # для интеграла по simpson
    # dS единичная площадка при интегрировании
    # ds = tilda_s dtheta dphi
    dl = surf.surf_R_e * ((3 * np.cos(surf.theta_range) ** 2 + 1) ** (1 / 2)) * np.sin(surf.theta_range)
    dphi = surf.surf_R_e * (np.sin(surf.theta_range) ** 3)
    tilda_s = dl * dphi
    return tilda_s


# shadows
# calc shadowed_matrix (ns + columns)
# tau
# calc tau_matrix

# для разных матриц можем посчитать L и посмотреть какой вклад будет.

def calc_L():
    ...


# save
# calc PF
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
