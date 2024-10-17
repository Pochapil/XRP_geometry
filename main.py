import numpy as np
import time

import accretingNS
import config
from geometricTask import matrix
import pathService
import shadows
import integralsService
import plot_package.plot_scripts

mu = 0.1e31
beta_mu = 40
mc2 = 100
a_portion = 0.44
phi_0 = 0

curr_configuration = accretingNS.AccretingPulsarConfiguration(mu, beta_mu, mc2, a_portion, phi_0)

theta_obs = 60

cur_path = pathService.PathSaver(mu, theta_obs, beta_mu, mc2, a_portion, phi_0).get_path()
print(cur_path)

theta_obs_rad = np.deg2rad(theta_obs)
beta_mu_rad = np.deg2rad(beta_mu)

x = np.sin(theta_obs_rad) * np.cos(0)
y = np.sin(theta_obs_rad) * np.sin(0)
z = np.cos(theta_obs_rad)
e_obs = np.array([x, y, z])

# obs_matrix = np.zeros((config.N_theta_accretion, config.N_phi_accretion, 3))
obs_matrix = np.zeros((config.N_phase, 3))

# rotate и сохранить матрицу векторов наблюдателя
for phase_index in range(config.N_phase):
    phi_mu = config.phi_mu_0 + config.omega_ns_rad * phase_index

    A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu, beta_mu_rad)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)

    obs_matrix[phase_index] = e_obs_mu

# чтобы соответствовать порядку в старом
# surfaces = {0: top_column.outer_surface, 1: top_column.inner_surface,
#             2: bot_column.outer_surface, 3: bot_column.inner_surface}
accr_col_surfs = [curr_configuration.top_column.outer_surface, curr_configuration.top_column.inner_surface,
                  curr_configuration.bot_column.outer_surface, curr_configuration.bot_column.inner_surface]

T_eff = curr_configuration.top_column.T_eff
print(integralsService.calculate_total_luminosity(curr_configuration.top_column.inner_surface, T_eff))
print(curr_configuration.top_column.L_x)

L = np.empty(4, dtype=object)
L_nu = np.empty(4, dtype=object)

for i, surface in enumerate(accr_col_surfs):
    # тензор косинусов между нормалью и направлением на наблюдателя размером phase x phi x theta
    # умножаем скалярно phi x theta x 3 на phase x 3 (по последнему индексу) и делаем reshape.
    cos_psi_rotation_matrix = np.einsum('ijl,tl->tij', surface.array_normal, obs_matrix)
    # возможно хранить будем все матрицы.
    # для разных матриц можем посчитать L и посмотреть какой вклад будет.
    tensor_shadows_NS = np.ones_like(cos_psi_rotation_matrix)
    tensor_shadows_columns = np.ones_like(cos_psi_rotation_matrix)
    tensor_tau = np.ones_like(cos_psi_rotation_matrix)

    # old_cos_psi_range = cos_psi_rotation_matrix.copy()
    new_cos_psi_range = cos_psi_rotation_matrix.copy()

    for phase_index in range(config.N_phase):
        for phi_index in range(config.N_phi_accretion):
            for theta_index in range(config.N_theta_accretion):
                if new_cos_psi_range[phase_index, phi_index, theta_index] > 0:
                    origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]
                    if shadows.intersection_with_sphere(surface, origin_phi, origin_theta, obs_matrix[phase_index]):
                        tensor_shadows_NS[phase_index, phi_index, theta_index] = 0
                    else:
                        solutions = shadows.get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta,
                                                                                  obs_matrix[phase_index])
                        tensor_shadows_columns[phase_index, phi_index, theta_index] = \
                            shadows.check_shadow_with_dipole(surface, phi_index, theta_index,
                                                             obs_matrix[phase_index], solutions,
                                                             curr_configuration.top_column.inner_surface,
                                                             curr_configuration.bot_column.inner_surface)

                        tensor_tau[phase_index, phi_index, theta_index] = \
                            shadows.get_tau_with_dipole(surface, phi_index, theta_index, obs_matrix[phase_index],
                                                        solutions, surface.surf_R_e, beta_mu,
                                                        curr_configuration.M_accretion_rate,
                                                        curr_configuration.top_column.inner_surface,
                                                        curr_configuration.bot_column.inner_surface,
                                                        a_portion)
                else:
                    new_cos_psi_range[phase_index, phi_index, theta_index] = 0
    new_cos_psi_range = new_cos_psi_range * tensor_shadows_NS * tensor_shadows_columns * tensor_tau

    L[i] = integralsService.calc_L(surface, T_eff, new_cos_psi_range)
    # save
    L_nu[i] = integralsService.calc_L_nu(surface, T_eff, new_cos_psi_range)
    # save
    # calc PF

print(L)
print(L_nu)
plot_package.plot_scripts.plot_L(L)

magnet_line_surfs = []
for surf in magnet_line_surfs:
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

    # save
    # calc nu L nu

    # calc_scatter_L_nu()

    # save
    # calc nu L nu

    # tij -> ti -> t
    # scipy.integrate(axis=-1) --- можно интегрировать по тензору.

# False - т.к. мне нужны внутренние нормали
#         normal_array = self.outer_surface.create_array_normal(self.magnet_lines_phi_range,
#                                                               self.magnet_lines_theta_range, False)
#         # sum_intense изотропная светимость ( * 4 pi еще надо)
#         # для интеграла по simpson
#         cos_psi_range_final = []
#         for t in range(config.t_max):
#             cos_psi_range = np.empty([config.N_phi_accretion, config.N_theta_accretion])
#             # поворот
#             phi_mu = config.phi_mu_0 + config.omega_ns * config.grad_to_rad * t
#             # расчет матрицы поворота в магнитную СК и вектора на наблюдателя
#             A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu,
#                                                          updated_betta_mu)
#             e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК
#             for i in range(config.N_phi_accretion):
#                 for j in range(config.N_theta_accretion):
#                     # умножать на N_theta
#                     # cos_psi_range[i, j] = np.dot(e_obs_mu, self.array_normal[i * config.N_theta_accretion + j])
#                     cos_psi_range[i, j] = np.dot(e_obs_mu, normal_array[i][j])
#                     if self.magnet_lines_mask_array[i, j]:
#                         cos_psi_range[i, j] = 0
#                     elif cos_psi_range[i, j] > 0:
#                         # проверка на пересечения
#                         r = self.outer_surface.R_e / config.R_ns * np.sin(self.magnet_lines_theta_range[j]) ** 2
#
#                         origin_x = np.sin(self.magnet_lines_theta_range[j]) * np.cos(self.magnet_lines_phi_range[i]) * r
#                         origin_y = np.sin(self.magnet_lines_theta_range[j]) * np.sin(self.magnet_lines_phi_range[i]) * r
#                         origin_z = np.cos(self.magnet_lines_theta_range[j]) * r
#
#                         direction_x = e_obs_mu[0, 0]
#                         direction_y = e_obs_mu[0, 1]
#                         direction_z = e_obs_mu[0, 2]
#
#                         if accretionColumnService.intersection_with_sphere(origin_x, origin_y, origin_z, direction_x,
#                                                                            direction_y, direction_z):
#                             cos_psi_range[i, j] = 0
#                         else:
#                             cos_psi_range[i, j] *= accretionColumnService.get_vals(self.magnet_lines_phi_range[i],
#                                                                                    self.magnet_lines_theta_range[j],
#                                                                                    e_obs_mu, origin_x, origin_y,
#                                                                                    origin_z, direction_x, direction_y,
#                                                                                    direction_z,
#                                                                                    self.outer_surface.ksi_shock,
#                                                                                    top_column_phi_range,
#                                                                                    bot_column_phi_range,
#                                                                                    top_column_theta_range,
#                                                                                    bot_column_theta_range,
#                                                                                    self.outer_surface.R_e,
#                                                                                    r, updated_betta_mu)
#                     else:
#                         cos_psi_range[i, j] = 0
#             cos_psi_range_final.append(cos_psi_range)
#         self.magnet_lines_cos_psi_range = cos_psi_range_final
