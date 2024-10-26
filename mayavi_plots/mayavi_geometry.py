import numpy as np

import geometry.matrix


def get_projection_on_vector(input_vector, project_vector):
    # cos * project_vector
    return np.dot(input_vector, project_vector) * project_vector / (np.linalg.norm(project_vector)) ** 2


def get_projection_on_surface(input_vector, surface_norm_vector):
    # на нормаль к плоскости и вычитаем из вектора
    projection_on_norm = get_projection_on_vector(input_vector, surface_norm_vector)
    return input_vector - projection_on_norm


def get_roll_angle(first_vector, second_vector, phase, theta_obs, beta_mu):
    if phase % 180 == 0:
        roll_angle = 0
        if theta_obs < beta_mu and phase % 360 == 0:
            roll_angle = 180
    else:
        r_first_vector = np.linalg.norm(first_vector)
        r_second_vector = np.linalg.norm(second_vector)
        if r_first_vector == 0 or r_second_vector == 0:
            roll_angle = 0
        else:
            roll_angle = -np.arccos(np.dot(first_vector, second_vector) / (r_first_vector * r_second_vector))
            roll_angle = np.rad2deg(roll_angle)

    if beta_mu == 0:
        roll_angle = 0

    if 180 < phase < 360 or phase > 540:
        roll_angle = - roll_angle

    return roll_angle


def calculate_roll_angle(theta_obs, beta_mu, e_obs_mu, phase):
    # phase in deg!
    '''Метод view init задает координаты картинной плоскости. Поэтому для того, чтобы перейти из магнитной СК обратно
    в СК, связанную с omega, надо достать угол между проекциями mu и omega на картинную плоскость и повернуть картинную плоскость
    на этот угол. Для этого нахожу проекцию на плоскость векторов mu и omega  и нахожу угол между ними через косинус'''
    view_plane_normal = list(geometry.matrix.vec_to_coord(e_obs_mu))  # x,y,z
    omega_vector = [np.sin(np.deg2rad(-beta_mu)) * np.cos(0),
                    np.sin(np.deg2rad(-beta_mu)) * np.sin(0),
                    np.cos(np.deg2rad(-beta_mu))]

    view_plane_normal = np.array(view_plane_normal)
    omega_vector = np.array(omega_vector)

    # беру проекцию на картинную плоскость, оттуда достаю угол.
    omega_projection_on_view_plane = get_projection_on_surface(omega_vector, view_plane_normal)
    mu_projection_on_view_plane = get_projection_on_surface(np.array([0, 0, 1]), view_plane_normal)

    roll_angle = get_roll_angle(omega_projection_on_view_plane, mu_projection_on_view_plane, phase, theta_obs, beta_mu)
    # print(roll_angle)
    return roll_angle
