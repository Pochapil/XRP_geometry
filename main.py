import numpy as np
import time
from itertools import repeat
import multiprocessing as mp

import accretingNS
import config
import newService
from geometry import matrix
import pathService
import shadows
import integralsService
import plot_package.plot_scripts
import save


def calc_shadows_and_tau(curr_configuration, surface, obs_matrix, mask_flag=False):
    # тензор косинусов между нормалью и направлением на наблюдателя размером phase x phi x theta
    # умножаем скалярно phi x theta x 3 на phase x 3 (по последнему индексу) и делаем reshape.
    cos_psi_rotation_matrix = np.einsum('ijl,tl->tij', surface.array_normal, obs_matrix)

    # print(f'obs_matrix = {np.max(np.linalg.norm(obs_matrix, axis=1))}')
    # print(f'array_normal_max = {np.max(np.linalg.norm(surface.array_normal, axis=2))}')

    # возможно хранить будем все матрицы.
    # для разных матриц можем посчитать L и посмотреть какой вклад будет.
    tensor_shadows_NS = np.ones_like(cos_psi_rotation_matrix)
    tensor_shadows_columns = np.ones_like(cos_psi_rotation_matrix)
    tensor_tau = np.ones_like(cos_psi_rotation_matrix)
    new_cos_psi_range = cos_psi_rotation_matrix.copy()

    if mask_flag:
        mask = np.zeros_like(new_cos_psi_range).astype(bool)
        mask += surface.mask_array
        new_cos_psi_range[mask] = 0
    # print(surface)
    for phase_index in range(config.N_phase):
        for phi_index in range(config.N_phi_accretion):
            for theta_index in range(config.N_theta_accretion):
                if new_cos_psi_range[phase_index, phi_index, theta_index] > 0:
                    origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]
                    # сначала проверяем на затмение НЗ
                    if shadows.intersection_with_sphere(surface, origin_phi, origin_theta, obs_matrix[phase_index]):
                        tensor_shadows_NS[phase_index, phi_index, theta_index] = 0
                    # иначе тяжелые вычисления
                    else:
                        # расчет для полинома - находим корни для пересечения с внутр поверхностью на магн линии!
                        solutions = shadows.get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta,
                                                                                  obs_matrix[phase_index])
                        # расчитываем затмение колонкой
                        tensor_shadows_columns[phase_index, phi_index, theta_index] = \
                            shadows.check_shadow_with_dipole(surface, phi_index, theta_index,
                                                             obs_matrix[phase_index], solutions,
                                                             curr_configuration.top_column.inner_surface,
                                                             curr_configuration.bot_column.inner_surface)
                        if tensor_shadows_columns[phase_index, phi_index, theta_index] > 0:
                            # если затмения нет то считаем ослабление тау
                            tensor_tau[phase_index, phi_index, theta_index] = \
                                shadows.get_tau_with_dipole(surface, phi_index, theta_index, obs_matrix[phase_index],
                                                            solutions, curr_configuration.top_column.R_e,
                                                            curr_configuration.beta_mu,
                                                            curr_configuration.M_accretion_rate,
                                                            curr_configuration.top_column.inner_surface,
                                                            curr_configuration.bot_column.inner_surface,
                                                            curr_configuration.a_portion)
                else:
                    # если косинус < 0 -> поверхность излучает от наблюдателя и мы не получим вклад в интеграл
                    new_cos_psi_range[phase_index, phi_index, theta_index] = 0
    new_cos_psi_range = new_cos_psi_range * tensor_shadows_NS * tensor_shadows_columns * tensor_tau
    # if tensor_tau[tensor_tau != 1].shape[0] > 0:
    #     print(np.max(tensor_tau[tensor_tau != 1]))
    return new_cos_psi_range


def calc_number_pow(num):
    pow = 0
    while num > 10:
        num = num / 10
        pow += 1
    return num, pow


if __name__ == '__main__':
    mu = 0.1e31
    beta_mu = 40
    mc2 = 100
    a_portion = 0.44
    phi_0 = 0

    curr_configuration = accretingNS.AccretingPulsarConfiguration(mu, beta_mu, mc2, a_portion, phi_0)

    theta_obs = 60

    cur_dir_saved = pathService.PathSaver(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    cur_path = cur_dir_saved.get_path()
    print(cur_path)
    save.create_file_path(cur_path)

    file_name = 'surfaces_T_eff.txt'
    save.save_arr_as_txt(curr_configuration.top_column.T_eff, cur_path, file_name)

    file_name = "save_phi_range.txt"
    save.save_arr_as_txt(curr_configuration.top_column.outer_surface.phi_range, cur_path, file_name)

    file_name = "save_theta_range.txt"
    save.save_arr_as_txt(curr_configuration.top_column.outer_surface.theta_range, cur_path, file_name)

    file_name = 'save_values.txt'
    with open(cur_path / file_name, 'w') as f:
        f.write(f'R_e = {curr_configuration.top_column.R_e / config.R_ns:.6f}\n')
        f.write(f'ksi_shock = {curr_configuration.top_column.ksi_shock:.6f}\n')
        f.write(f'beta = {curr_configuration.top_column.beta:.6f}\n')

        number, power = calc_number_pow(curr_configuration.top_column.L_x)
        f.write(f'total L_x = {number:.6f} * 10**{power}\n')

        L_calc = integralsService.calculate_total_luminosity(curr_configuration.top_column.inner_surface,
                                                             curr_configuration.top_column.T_eff)

        f.write(f'difference L_x / L_calc - 1 : {(curr_configuration.top_column.L_x / L_calc - 1) * 100:.6f} %\n')

        number, power = calc_number_pow(L_calc / 4)
        f.write(f'calculated total L_x of single surface = {number:.6f} * 10**{power}\n')

    theta_obs_rad = np.deg2rad(theta_obs)
    beta_mu_rad = np.deg2rad(beta_mu)

    e_obs = matrix.get_cartesian_from_spherical(1, theta_obs_rad, 0)

    # obs_matrix = np.zeros((config.N_theta_accretion, config.N_phi_accretion, 3))
    obs_matrix = np.zeros((config.N_phase, 3))

    # rotate и сохранить матрицу векторов наблюдателя
    for phase_index in range(config.N_phase):
        phi_mu = config.phi_mu_0 + config.omega_ns_rad * phase_index

        A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu, beta_mu_rad)
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)

        obs_matrix[phase_index] = e_obs_mu

    T_eff = curr_configuration.top_column.T_eff
    # print(integralsService.calculate_total_luminosity(curr_configuration.top_column.inner_surface, T_eff))
    # print(curr_configuration.top_column.L_x)

    # чтобы соответствовать порядку в старом
    # surfaces = {0: top_column.outer_surface, 1: top_column.inner_surface,
    #             2: bot_column.outer_surface, 3: bot_column.inner_surface}
    accr_col_surfs = [curr_configuration.top_column.outer_surface, curr_configuration.top_column.inner_surface,
                      curr_configuration.bot_column.outer_surface, curr_configuration.bot_column.inner_surface]

    L_surfs = np.empty((4, config.N_phase))
    L_nu_surfs = np.empty((4, config.N_energy, config.N_phase))

    # PF_L = np.empty(4)
    # PF_L_nu = np.empty((4, config.N_energy))
    # print(mp.cpu_count())
    t1 = time.perf_counter()
    with mp.Pool(processes=4) as pool:
        new_cos_psi_range_async = pool.starmap(calc_shadows_and_tau,
                                               zip(repeat(curr_configuration), accr_col_surfs, repeat(obs_matrix)))
    t2 = time.perf_counter()
    print(f'{t2 - t1} seconds')
    # ------------------------------------------------ L_calc ----------------------------------------------------
    for i, surface in enumerate(accr_col_surfs):
        # new_cos_psi_range = calc_shadows_and_tau(curr_configuration, surface, obs_matrix)
        # print(new_cos_psi_range_async[i][np.abs(new_cos_psi_range_async[i] - new_cos_psi_range) > 1e-3].shape)
        new_cos_psi_range = new_cos_psi_range_async[i]
        L_surfs[i] = integralsService.calc_L(surface, T_eff, new_cos_psi_range)
        L_nu_surfs[i] = integralsService.calc_L_nu(surface, T_eff, new_cos_psi_range)

    PF_L_surf = integralsService.get_PF(np.sum(L_surfs, axis=0))
    PF_L_nu_surf = integralsService.get_PF(np.sum(L_nu_surfs, axis=0))

    plot_package.plot_scripts.plot_L(L_surfs)

    file_name = 'save_values.txt'
    with open(cur_path / file_name, 'a') as f:
        avg_L_on_phase = np.mean(np.sum(L_surfs, axis=0))
        number, power = calc_number_pow(avg_L_on_phase)
        f.write(f'avg_L_on_phase = {number:.6f} * 10**{power}\n')
        f.write(f'avg_L_on_phase / L_x = {avg_L_on_phase / curr_configuration.top_column.L_x:.3}\n')

        d0 = newService.get_delta_distance(curr_configuration.top_column.inner_surface.theta_range[0],
                                           curr_configuration.top_column.inner_surface.surf_R_e)
        l0 = newService.get_A_normal(curr_configuration.top_column.inner_surface.theta_range[0],
                                     curr_configuration.top_column.inner_surface.surf_R_e, a_portion) / d0

        f.write(f'width / length = {(d0 / l0):.6}\n')
        f.write(f'A_perp / 2 d0**2 = {(l0 / (2 * d0)):.6}\n')
        f.write('L_data = dummy\n')
        # f.write(f'L_data = {number:.6f} * 10**{power}\n')
    # save

    magnet_line_surfs = [curr_configuration.top_magnet_lines, curr_configuration.bot_magnet_lines]
    L_scatter = np.empty((2, config.N_phase))
    L_nu_scatter = np.empty((2, config.N_energy, config.N_phase))

    t1 = time.perf_counter()
    with mp.Pool(processes=2) as pool:
        new_cos_psi_range_async = pool.starmap(calc_shadows_and_tau,
                                               zip(repeat(curr_configuration), magnet_line_surfs, repeat(obs_matrix),
                                                   repeat(True)))
    t2 = time.perf_counter()
    print(f'{t2 - t1} seconds')
    # ------------------------------------------------ L_scatter ----------------------------------------------------
    for i, magnet_surface in enumerate(magnet_line_surfs):
        # new_cos_psi_range = calc_shadows_and_tau(curr_configuration, magnet_surface, obs_matrix, True)
        new_cos_psi_range = new_cos_psi_range_async[i]

        # УЧЕСТЬ угол падения от центра в отражаемой точке!!!!!!!!!
        xyz_magnet_line = matrix.get_xyz_coord(magnet_surface, normalize=True)  # вектор на отражающую площадку
        # -xyz потому что берем угол от нормали к НЗ
        # cos_alpha_matrix - это также матрица необходимая для расчета tau - как накрест лежащие
        cos_alpha_matrix = np.einsum('ijl,ijl->ij', magnet_surface.array_normal, -xyz_magnet_line)
        # cos_alpha_1_matrix = np.einsum('ijl,ijl->ij', -surface.array_normal, xyz_magnet_line)
        new_cos_psi_range *= cos_alpha_matrix

        # УЧЕСТЬ tau в отражаемой точке!!!!!!!!!
        # tau_scatter_matrix = np.ones_like(cos_alpha_matrix)
        tau_scatter_matrix = shadows.get_tau_for_scatter_with_cos(magnet_surface.theta_range, magnet_surface.surf_R_e,
                                                                  curr_configuration.M_accretion_rate, a_portion,
                                                                  cos_alpha_matrix)

        L_x = integralsService.calculate_total_luminosity(curr_configuration.top_column.inner_surface, T_eff)
        # надо брать не curr_configuration.top_column.L_x а посчитанный через интеграл! хоть они и отличаются на сильно
        L_scatter[i] = integralsService.calc_scatter_L(magnet_surface, L_x, new_cos_psi_range, tau_scatter_matrix)

        L_nu_scatter[i] = integralsService.calc_scatter_L_nu(magnet_surface,
                                                             curr_configuration.top_column.inner_surface, T_eff,
                                                             new_cos_psi_range, tau_scatter_matrix)

    plot_package.plot_scripts.plot_L(L_scatter)

    PF_L_surf = integralsService.get_PF(np.sum(L_scatter, axis=0))
    PF_L_nu_surf = integralsService.get_PF(np.sum(L_nu_scatter, axis=0))

    file_name = 'save_values.txt'
    with open(cur_path / file_name, 'a') as f:
        number, power = calc_number_pow(avg_L_on_phase + np.mean(np.sum(L_scatter, axis=0)))
        f.write(f'avg L with scatter = {number:.6f} * 10**{power}\n')
        number = (avg_L_on_phase + np.mean(np.sum(L_scatter, axis=0))) / curr_configuration.top_column.L_x
        f.write(f'avg L with scatter / L_x = {number:.6f}\n')
