import numpy as np
import time

import accretingNS
import config
from geometry import matrix
import pathService
import shadows
import integralsService
import plot_package.plot_scripts


def calc_shadows_and_tau(curr_configuration, surface, obs_matrix, mask_flag=False):
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

    if mask_flag:
        mask = np.zeros_like(new_cos_psi_range).astype(bool)
        mask += surface.mask_array
        new_cos_psi_range[mask] = 0

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
                                                        solutions, surface.surf_R_e, curr_configuration.beta_mu,
                                                        curr_configuration.M_accretion_rate,
                                                        curr_configuration.top_column.inner_surface,
                                                        curr_configuration.bot_column.inner_surface,
                                                        curr_configuration.a_portion)
                else:
                    new_cos_psi_range[phase_index, phi_index, theta_index] = 0
    new_cos_psi_range = new_cos_psi_range * tensor_shadows_NS * tensor_shadows_columns * tensor_tau
    if tensor_tau[tensor_tau != 1].shape[0] > 0:
        print(np.max(tensor_tau[tensor_tau != 1]))
    return new_cos_psi_range


if __name__ == '__main__':

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
    # ------------------------------------------------ L_calc ----------------------------------------------------
    for i, surface in enumerate(accr_col_surfs):
        new_cos_psi_range = calc_shadows_and_tau(curr_configuration, surface, obs_matrix)

        L[i] = integralsService.calc_L(surface, T_eff, new_cos_psi_range)
        # save
        L_nu[i] = integralsService.calc_L_nu(surface, T_eff, new_cos_psi_range)
        # save
        # calc PF

    print(L)
    print(L_nu)
    plot_package.plot_scripts.plot_L(L)

    magnet_line_surfs = [curr_configuration.top_magnet_lines, curr_configuration.bot_magnet_lines]

    L = np.empty(2, dtype=object)
    L_nu = np.empty(2, dtype=object)

    # ------------------------------------------------ L_scatter ----------------------------------------------------
    for i, surface in enumerate(magnet_line_surfs):
        new_cos_psi_range = calc_shadows_and_tau(curr_configuration, surface, obs_matrix, True)
        print(new_cos_psi_range)

        # УЧЕСТЬ угол падения от центра в отражаемой точке!!!!!!!!!
        xyz_magnet_line = matrix.get_xyz_coord(surface, normalize=True)  # вектор на отражающую площадку
        # -xyz потому что берем угол от нормали к НЗ
        # cos_alpha_matrix - это также матрица необходимая для расчета tau - как накрест лежащие
        cos_alpha_matrix = np.einsum('ijl,ijl->ij', surface.array_normal, -xyz_magnet_line)
        # cos_alpha_1_matrix = np.einsum('ijl,ijl->ij', -surface.array_normal, xyz_magnet_line)
        new_cos_psi_range *= cos_alpha_matrix

        # УЧЕСТЬ tau в отражаемой точке!!!!!!!!!
        # tau_scatter_matrix = np.ones_like(cos_alpha_matrix)
        tau_scatter_matrix = shadows.get_tau_for_scatter_with_cos(surface.theta_range, surface.surf_R_e,
                                                                  curr_configuration.M_accretion_rate, a_portion,
                                                                  cos_alpha_matrix)

        L[i] = integralsService.calc_scatter_L(surface, curr_configuration.top_column.L_x, new_cos_psi_range,
                                               tau_scatter_matrix)

        # save

        # calc_scatter_L_nu()

        # save

    plot_package.plot_scripts.plot_L(L)
