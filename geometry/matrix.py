import numpy as np
from numba import njit

import config


# матрицы поворота
def newRx(f):
    return np.array([[1, 0, 0], [0, np.cos(f), np.sin(f)], [0, -np.sin(f), np.cos(f)]])  # поворот против часовой
    # return np.matrix([[1, 0, 0], [0, np.cos(f), -np.sin(f)], [0, np.sin(f), np.cos(f)]]) # поворот по часовой


def newRy(f):
    return np.array([[np.cos(f), 0, -np.sin(f)], [0, 1, 0], [np.sin(f), 0, np.cos(f)]])
    # return np.matrix([[np.cos(f), 0, np.sin(f)], [0, 1, 0], [-np.sin(f), 0, np.cos(f)]])


def newRz(f):
    return np.array([[np.cos(f), np.sin(f), 0], [-np.sin(f), np.cos(f), 0], [0, 0, 1]])
    # return np.matrix([[np.cos(f), -np.sin(f), 0], [np.sin(f), np.cos(f), 0], [0, 0, 1]])


# аналитическая матрица поворота из двойной системы в СК магнитную
# R_y(betta_mu) @ R_z(phi_mu) @ R_y(betta_rotate) @ A_z(phi_rotate)
# A_matrix_calc = matrix.newRy(betta_mu) @ matrix.newRz(fi_mu) @ matrix.newRy(betta_rotate) \
#                 @ matrix.newRz(fi_rotate)

def A_matrix(fi_rotate, betta_rotate, fi_mu, betta_mu):
    return newRy(betta_mu) @ newRz(fi_mu) @ newRy(betta_rotate) @ newRz(fi_rotate)


def newMatrixAnalytic(fi_rotate, betta_rotate, fi_mu, betta_mu):
    a_11 = np.cos(fi_rotate) * (
            np.cos(betta_rotate) * np.cos(betta_mu) * np.cos(fi_mu) - np.sin(betta_rotate) * np.sin(betta_mu)) \
           - np.cos(betta_mu) * np.sin(fi_rotate) * np.sin(fi_mu)

    a_12 = np.sin(fi_rotate) * (
            np.cos(betta_rotate) * np.cos(betta_mu) * np.cos(fi_mu) - np.sin(betta_rotate) * np.sin(betta_mu)) \
           + np.cos(fi_rotate) * np.cos(betta_mu) * np.sin(fi_mu)

    a_13 = - np.cos(betta_mu) * np.cos(fi_mu) * np.sin(betta_rotate) - np.cos(betta_rotate) * np.sin(betta_mu)

    a_21 = - np.cos(fi_mu) * np.sin(fi_rotate) - np.cos(betta_rotate) * np.cos(fi_rotate) * np.sin(fi_mu)

    a_22 = np.cos(fi_rotate) * np.cos(fi_mu) - np.cos(betta_rotate) * np.sin(fi_rotate) * np.sin(fi_mu)

    a_23 = np.sin(betta_rotate) * np.sin(fi_mu)

    a_31 = np.cos(fi_rotate) * (
            np.cos(betta_mu) * np.sin(betta_rotate) + np.cos(betta_rotate) * np.cos(fi_mu) * np.sin(betta_mu)) \
           - np.sin(fi_rotate) * np.sin(betta_mu) * np.sin(fi_mu)

    a_32 = np.sin(fi_rotate) * (
            np.cos(betta_mu) * np.sin(betta_rotate) + np.cos(betta_rotate) * np.cos(fi_mu) * np.sin(betta_mu)) \
           + np.cos(fi_rotate) * np.sin(betta_mu) * np.sin(fi_mu)

    a_33 = np.cos(betta_rotate) * np.cos(betta_mu) - np.cos(fi_mu) * np.sin(betta_rotate) * np.sin(betta_mu)

    # return np.array([[a_11, a_12, a_13], [a_21, a_22, a_23], [a_31, a_32, a_33]])
    return np.array([[a_11, a_12, a_13], [a_21, a_22, a_23], [a_31, a_32, a_33]])


# базисы в сферической СК
def newE_r(fi_sphere, theta_sphere):
    return np.array(
        [np.sin(theta_sphere) * np.cos(fi_sphere), np.sin(theta_sphere) * np.sin(fi_sphere), np.cos(theta_sphere)])


def newE_theta(fi_sphere, theta_sphere):
    return np.array(
        [np.cos(theta_sphere) * np.cos(fi_sphere), np.cos(theta_sphere) * np.sin(fi_sphere), -np.sin(theta_sphere)])


def newE_fi(fi_sphere):
    return np.array([-np.sin(fi_sphere), np.cos(fi_sphere), 0])


# единичный вектор вдоль силовых линий
def newE_l(fi_sphere, theta_sphere):
    return (2 * np.cos(theta_sphere) * newE_r(fi_sphere, theta_sphere) + np.sin(theta_sphere)
            * newE_theta(fi_sphere, theta_sphere)) / ((3 * np.cos(theta_sphere) ** 2 + 1) ** (1 / 2))


# нормаль к силовым линиям
def newE_n(fi_sphere, theta_sphere):
    return np.cross(newE_l(fi_sphere, theta_sphere), newE_fi(fi_sphere))


# ============================================
# numpy way:
#@njit
def newE_r_n(phi_sphere, theta_sphere):
    # res = тензор размером phi x theta x 3 (x,y,z)
    # x,y,z = тензор размером  phi x theta
    x = np.sin(theta_sphere)[np.newaxis, :] * np.cos(phi_sphere)[:, np.newaxis]
    y = np.sin(theta_sphere)[np.newaxis, :] * np.sin(phi_sphere)[:, np.newaxis]
    z = np.cos(theta_sphere)[np.newaxis, :] * np.ones_like(phi_sphere)[:, np.newaxis]
    res = np.hstack((x, y, z)).reshape((phi_sphere.shape[0], theta_sphere.shape[0], 3), order='F')  # phi x theta x 3
    return res


#@njit
def newE_theta_n(phi_sphere, theta_sphere):
    # res = тензор размером phi x theta x 3 (x,y,z)
    # x,y,z = тензор размером  phi x theta
    x = np.cos(theta_sphere)[np.newaxis, :] * np.cos(phi_sphere)[:, np.newaxis]
    y = np.cos(theta_sphere)[np.newaxis, :] * np.sin(phi_sphere)[:, np.newaxis]
    z = -np.sin(theta_sphere) * np.ones_like(phi_sphere)[:, np.newaxis]
    res = np.hstack((x, y, z)).reshape((phi_sphere.shape[0], theta_sphere.shape[0], 3), order='F')
    return res


#@njit
def newE_fi_n(fi_sphere):
    # res = тензор размером phi x 3 (x,y,z)
    # x,y,z = тензор размером  phi
    x = np.cos(fi_sphere)
    y = -np.sin(fi_sphere)
    z = np.zeros(fi_sphere.shape[0])
    res = np.hstack((x, y, z)).reshape((-1, 3), order='F')
    return res


#@njit
def newE_l_n(phi_sphere, theta_sphere):
    # тензор размером phi x theta x 3 (x,y,z)
    # 4.4 формула
    res = 2 * np.cos(theta_sphere)[np.newaxis, :, np.newaxis] * newE_r_n(phi_sphere, theta_sphere)
    res += np.sin(theta_sphere)[np.newaxis, :, np.newaxis] * newE_theta_n(phi_sphere, theta_sphere)
    res /= ((3 * np.cos(theta_sphere)[np.newaxis, :, np.newaxis] ** 2 + 1) ** (1 / 2))
    return res


#@njit
def newE_n_n(phi_sphere, theta_sphere):
    # тензор размером phi x theta x 3 (x,y,z)
    # вместо np.cross(newE_l(fi_sphere, theta_sphere), newE_fi(fi_sphere)) беру формулу для n
    # был вывод в дипломе 4.12 формула
    res = -2 * np.cos(theta_sphere)[np.newaxis, :, np.newaxis] * newE_theta_n(phi_sphere, theta_sphere)
    res += np.sin(theta_sphere)[np.newaxis, :, np.newaxis] * newE_r_n(phi_sphere, theta_sphere)
    res /= ((3 * np.cos(theta_sphere)[np.newaxis, :, np.newaxis] ** 2 + 1) ** (1 / 2))
    return res


# scatter_point_coord_matrix = np.dstack(np.meshgrid(x, y, z))
#@njit
def get_xyz_coord(surface, normalize=False):
    r = surface.surf_R_e / config.R_ns * np.sin(surface.theta_range) ** 2
    if normalize:
        r = 1
    x = r * np.sin(surface.theta_range)[np.newaxis, :] * np.cos(surface.phi_range)[:, np.newaxis]
    y = r * np.sin(surface.theta_range)[np.newaxis, :] * np.sin(surface.phi_range)[:, np.newaxis]
    z = r * np.cos(surface.theta_range)[np.newaxis, :] * np.ones_like(surface.phi_range)[:, np.newaxis]
    # order='F' = чтобы было N_phi_accretion x N_theta_accretion x (xyz) иначе он бы взял x потом y потом z и решейп
    res = np.hstack((x, y, z)).reshape((config.N_phi_accretion, config.N_theta_accretion, 3), order='F')
    return res


@njit
def vec_to_angles(vector):
    x = vector[0]
    y = vector[1]
    z = vector[2]
    r = np.linalg.norm(vector)
    # r = (x ** 2 + y ** 2 + z ** 2) ** (1 / 2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return phi, theta

#@njit
def vec_to_coord(vector):
    x = vector[0]
    y = vector[1]
    z = vector[2]
    return x, y, z

#@njit
def get_cartesian_from_spherical(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

#@njit
def get_xyz_coord_angles(r, phi_arr, theta_arr):
    x = r * np.sin(theta_arr)[np.newaxis, :] * np.cos(phi_arr)[:, np.newaxis]
    y = r * np.sin(theta_arr)[np.newaxis, :] * np.sin(phi_arr)[:, np.newaxis]
    z = r * np.cos(theta_arr) * np.ones_like(phi_arr)[:, np.newaxis]
    res = np.hstack((x, y, z)).reshape((phi_arr.shape[0], theta_arr.shape[0], 3), order='F')
    return res
