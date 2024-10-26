from itertools import repeat, cycle
import multiprocessing as mp
import pathlib
import time
import numpy as np

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
    '''
    функция расчитывает матрицу косинусов на каждой фазе для направлений наблюдателя

    mask_flag --- для расчета по магнитным линиям - их надо обрезать
    '''
    # тензор косинусов между нормалью на площадке и направлением на наблюдателя. размер phase x phi x theta
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
        mask += surface.mask_array  # происходит broadcast в матрицу размера phase x phi x theta
        new_cos_psi_range[mask] = 0
    # print(surface)
    # for phase_index in range(config.N_phase):
    # для возможности распараллелить еще больше чем на 4 поверхности (разрезая наблюдателя на чанки)
    for phase_index in range(obs_matrix.shape[0]):
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


def make_save_values_file(curr_configuration, L_surfs, L_scatter, cur_path_data, cur_dir_saved):
    file_name = 'save_values.txt'

    if config.old_path_flag:
        old_path = cur_dir_saved.get_old_path()
        cur_path = old_path
    else:
        cur_path = cur_path_data

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

        avg_L_on_phase = np.mean(np.sum(L_surfs, axis=0))
        number, power = calc_number_pow(avg_L_on_phase)
        f.write(f'avg_L_on_phase = {number:.6f} * 10**{power}\n')
        f.write(f'avg_L_on_phase / L_x = {avg_L_on_phase / curr_configuration.top_column.L_x:.3}\n')

        d0 = newService.get_delta_distance(curr_configuration.top_column.inner_surface.theta_range[0],
                                           curr_configuration.top_column.inner_surface.surf_R_e)
        l0 = newService.get_A_normal(curr_configuration.top_column.inner_surface.theta_range[0],
                                     curr_configuration.top_column.inner_surface.surf_R_e,
                                     curr_configuration.a_portion) / d0

        f.write(f'width / length = {(d0 / l0):.6}\n')
        f.write(f'A_perp / 2 d0**2 = {(l0 / (2 * d0)):.6}\n')
        f.write('L_data = dummy\n')
        # f.write(f'L_data = {number:.6f} * 10**{power}\n')

        number, power = calc_number_pow(avg_L_on_phase + np.mean(np.sum(L_scatter, axis=0)))
        f.write(f'avg L with scatter = {number:.6f} * 10**{power}\n')
        number = (avg_L_on_phase + np.mean(np.sum(L_scatter, axis=0))) / curr_configuration.top_column.L_x
        f.write(f'avg L with scatter / L_x = {number:.6f}\n')


def save_some_files(curr_configuration, obs_matrix, L_surfs, L_scatter, L_nu_surfs, L_nu_scatter, PF_L_nu_surfs,
                    cur_path_data, cur_dir_saved):
    if config.old_path_flag:
        cur_path = cur_dir_saved.get_old_path()
    else:
        cur_path = cur_path_data

    file_name = 'surfaces_T_eff.txt'
    save.save_arr_as_txt(curr_configuration.top_column.T_eff, cur_path, file_name)

    file_name = "save_phi_range.txt"
    save.save_arr_as_txt(curr_configuration.top_column.outer_surface.phi_range, cur_path, file_name)

    file_name = "save_theta_range.txt"
    save.save_arr_as_txt(curr_configuration.top_column.outer_surface.theta_range, cur_path, file_name)

    tau_array = shadows.get_tau_for_opacity_old(curr_configuration.top_magnet_lines.theta_range,
                                                curr_configuration.R_e, curr_configuration.M_accretion_rate,
                                                curr_configuration.a_portion)
    file_name = "save_tau_range.txt"
    save.save_arr_as_txt(tau_array, cur_path, file_name)

    ans = np.apply_along_axis(matrix.vec_to_angles, axis=1, arr=obs_matrix)
    observer_phi = ans[:, 0]
    observer_theta = ans[:, 1]
    file_name = 'observer_angles.txt'
    np.savetxt(cur_path / file_name, np.vstack((observer_phi, observer_theta)), delimiter=',')
    # save.save_arr_as_txt(np.vstack((observer_phi, observer_theta)), cur_path, file_name)

    file_name = 'total_luminosity_of_surfaces.txt'
    save.save_arr_as_txt(np.append(L_surfs, np.sum(L_surfs, axis=0)[np.newaxis, :], axis=0), cur_path, file_name)

    file_name = "energy.txt"
    save.save_arr_as_txt(config.energy_arr, cur_path, file_name)

    # --------------------------------------------folders------------------------------------------------

    if config.old_path_flag:

        old_path = cur_dir_saved.get_old_path()

        for i, energy in enumerate(config.energy_arr):
            file_name = f"L_nu_of_energy_{energy:.2f}_KeV_of_surfaces.txt"
            save.save_arr_as_txt(np.sum(L_nu_surfs, axis=0)[i], old_path / ('L_nu/' + 'txt/'), file_name)

            file_name = f"nu_L_nu_of_energy_{energy:.2f}_KeV_of_surfaces.txt"
            freq = newService.get_frequency_from_energy(energy)
            save.save_arr_as_txt(np.sum(L_nu_surfs, axis=0)[i] * freq, old_path / ('nu_L_nu/' + 'txt/'), file_name)

        file_name = "PF.txt"
        save.save_arr_as_txt(PF_L_nu_surfs, old_path / 'L_nu/', file_name)
        save.save_arr_as_txt(PF_L_nu_surfs, old_path / 'nu_L_nu/', file_name)

        file_name = "scattered_energy_top.txt"
        save.save_arr_as_txt(L_scatter[0], old_path / 'scattered_on_magnet_lines/', file_name)

        file_name = "scattered_energy_bot.txt"
        save.save_arr_as_txt(L_scatter[1], old_path / 'scattered_on_magnet_lines/', file_name)

        file_name = "top_column_scatter_L_nu.txt"
        save.save_arr_as_txt(L_nu_scatter[0], old_path / ('scattered_on_magnet_lines/' + 'L_nu/'), file_name)

        file_name = "bot_column_scatter_L_nu.txt"
        save.save_arr_as_txt(L_nu_scatter[1], old_path / ('scattered_on_magnet_lines/' + 'L_nu/'), file_name)


def save_new_way(L_surfs, L_scatter, L_nu_surfs, L_nu_scatter, cur_path_data):
    cur_path = cur_path_data
    file_name = "L_surfs.txt"
    save.save_arr_as_txt(L_surfs, cur_path / 'surfs/', file_name)

    file_name = "L_scatter.txt"
    save.save_arr_as_txt(L_scatter, cur_path / 'scatter/', file_name)

    for i, energy in enumerate(config.energy_arr):
        file_name = f"L_nu_surfs_of_energy_{energy:.2f}_KeV.txt"
        save.save_arr_as_txt(L_nu_surfs[:, i], cur_path / 'surfs/', file_name)

        file_name = f"L_nu_scatter_of_energy_{energy:.2f}_KeV.txt"
        save.save_arr_as_txt(L_nu_scatter[:, i], cur_path / 'scatter/', file_name)


def calc_async_with_split(curr_configuration, obs_matrix, surfs_arr, mask_flag):
    '''для возможности распаралелить еще больше чем на 4 поверхности
    тут мы разбиваем массив наблюдателя на чанки и передаем на расчет, потом сливаем их вместе
    '''

    obs_matrix_new = np.split(obs_matrix, [obs_matrix.shape[0] // 2])
    obs_matrix_new_to_func = []
    for i in range(len(surfs_arr)):
        obs_matrix_new_to_func.append(obs_matrix_new[0])
    for i in range(len(surfs_arr)):
        obs_matrix_new_to_func.append(obs_matrix_new[1])

    t1 = time.perf_counter()
    with mp.Pool(processes=2 * len(surfs_arr)) as pool:
        new_cos_psi_range_async_8 = pool.starmap(calc_shadows_and_tau,
                                                 zip(repeat(curr_configuration), cycle(surfs_arr),
                                                     obs_matrix_new_to_func, repeat(mask_flag)))
    t2 = time.perf_counter()
    print(f'{t2 - t1} seconds')

    new_cos_psi_range_async = np.empty(
        (len(surfs_arr), config.N_phase, config.N_phi_accretion, config.N_theta_accretion))

    for i in range(len(surfs_arr)):
        new_cos_psi_range_async[i] = np.vstack(
            (new_cos_psi_range_async_8[i], new_cos_psi_range_async_8[i + len(surfs_arr)]))
        # print(new_cos_psi_range_async_4[i][np.abs(new_cos_psi_range_async_4[i] - new_cos_psi_range) > 1e-3].shape)
    return new_cos_psi_range_async


def calc_cos_psi(curr_configuration, obs_matrix, surfs_arr, mask_flag):
    '''объединил все методы расчета в 1 функцию'''
    # 4 или 8 процессов будет запущено
    if config.ASYNC_FLAG:
        if config.N_cpus == 8:
            new_cos_psi_range = calc_async_with_split(curr_configuration, obs_matrix, surfs_arr, mask_flag)
        else:
            t1 = time.perf_counter()
            with mp.Pool(processes=len(surfs_arr)) as pool:
                new_cos_psi_range = pool.starmap(calc_shadows_and_tau,
                                                 zip(repeat(curr_configuration), surfs_arr,
                                                     repeat(obs_matrix), repeat(mask_flag)))
            t2 = time.perf_counter()
            print(f'{t2 - t1} seconds')
    else:
        new_cos_psi_range = np.empty(4, dtype=object)
        for i, surface in enumerate(surfs_arr):
            new_cos_psi_range[i] = calc_shadows_and_tau(curr_configuration, surface, obs_matrix, mask_flag)

    return new_cos_psi_range


def calc_and_save_for_configuration(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, figs_flag=False):
    curr_configuration = accretingNS.AccretingPulsarConfiguration(mu, beta_mu, mc2, a_portion, phi_0)

    cur_dir_saved = pathService.PathSaver(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    folder = 'txt/'
    # cur_path = path to save txt !!!
    cur_path = cur_dir_saved.get_path()
    cur_path_data = cur_path / folder
    print(cur_path_data)
    save.create_file_path(cur_path_data)
    print(cur_dir_saved.save_dir)
    theta_obs_rad = np.deg2rad(theta_obs)
    beta_mu_rad = np.deg2rad(beta_mu)

    e_obs = matrix.get_cartesian_from_spherical(1, theta_obs_rad, 0)

    obs_matrix = np.empty((config.N_phase, 3))
    # rotate и сохранить матрицу векторов наблюдателя
    # удобнее всего работать в системе координат связанной с магнитной осью - там конфигурация колонок постоянна
    # поэтому переходим от вращения звезды к вращению наблюдателя
    for phase_index in range(config.N_phase):
        # поворот
        phi_mu = config.phi_mu_0 + config.omega_ns_rad * phase_index
        # расчет матрицы поворота в магнитную СК для вектора на наблюдателя
        # (изначально задавали наблюдателя в системе координат оси вращения)
        A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu, beta_mu_rad)
        # переход в магнитную СК
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)
        obs_matrix[phase_index] = e_obs_mu

    # ------------------------------------------------ L_calc -------------------------------------------------------
    accr_col_surfs = [curr_configuration.top_column.outer_surface, curr_configuration.top_column.inner_surface,
                      curr_configuration.bot_column.outer_surface, curr_configuration.bot_column.inner_surface]
    # чтобы соответствовать порядку в старом
    # surfaces = {0: top_column.outer_surface, 1: top_column.inner_surface,
    #             2: bot_column.outer_surface, 3: bot_column.inner_surface}

    print('start calc surfs')
    new_cos_psi_range_surfs = calc_cos_psi(curr_configuration, obs_matrix, accr_col_surfs, False)

    L_surfs = np.empty((4, config.N_phase))
    L_nu_surfs = np.empty((4, config.N_energy, config.N_phase))
    for i, surface in enumerate(accr_col_surfs):
        # new_cos_psi_range = calc_shadows_and_tau(curr_configuration, surface, obs_matrix)
        # print(new_cos_psi_range_async[i][np.abs(new_cos_psi_range_async[i] - new_cos_psi_range) > 1e-3].shape)
        new_cos_psi_range = new_cos_psi_range_surfs[i]
        L_surfs[i] = integralsService.calc_L(surface, curr_configuration.top_column.T_eff, new_cos_psi_range)
        L_nu_surfs[i] = integralsService.calc_L_nu(surface, curr_configuration.top_column.T_eff, new_cos_psi_range)

    PF_L_surfs = newService.get_PF(np.sum(L_surfs, axis=0))
    PF_L_nu_surfs = newService.get_PF(np.sum(L_nu_surfs, axis=0))
    print('finish calc surfs')
    # plot_package.plot_scripts.plot_L(L_surfs)

    # ----------------------------------------------- Scatter -----------------------------------------------------
    print('start calc scatter')
    magnet_line_surfs = [curr_configuration.top_magnet_lines, curr_configuration.bot_magnet_lines]

    new_cos_psi_range_surfs = calc_cos_psi(curr_configuration, obs_matrix, magnet_line_surfs, True)

    L_scatter = np.empty((2, config.N_phase))
    L_nu_scatter = np.empty((2, config.N_energy, config.N_phase))
    # ------------------------------------------------ L_scatter ----------------------------------------------------
    for i, magnet_surface in enumerate(magnet_line_surfs):
        # new_cos_psi_range = calc_shadows_and_tau(curr_configuration, magnet_surface, obs_matrix, True)
        new_cos_psi_range = new_cos_psi_range_surfs[i]

        # УЧЕСТЬ угол падения от центра в отражаемой точке!!!!!!!!!
        xyz_magnet_line = matrix.get_xyz_coord(magnet_surface, normalize=True)  # вектор на отражающую площадку
        # -xyz потому что берем угол от нормали к НЗ
        # cos_alpha_matrix - это также матрица необходимая для расчета tau - как накрест лежащие
        cos_alpha_matrix = np.einsum('ijl,ijl->ij', magnet_surface.array_normal, -xyz_magnet_line)
        # cos_alpha_1_matrix = np.einsum('ijl,ijl->ij', -surface.array_normal, xyz_magnet_line)
        new_cos_psi_range *= cos_alpha_matrix  # учитываю угол падения от центра в отражаемой точке!!!!!!!!!

        # УЧЕСТЬ tau в отражаемой точке!!!!!!!!!
        # tau_scatter_matrix = np.ones_like(cos_alpha_matrix)
        tau_scatter_matrix = shadows.get_tau_for_scatter_with_cos(magnet_surface.theta_range, magnet_surface.surf_R_e,
                                                                  curr_configuration.M_accretion_rate, a_portion,
                                                                  cos_alpha_matrix)

        L_x = integralsService.calculate_total_luminosity(curr_configuration.top_column.inner_surface,
                                                          curr_configuration.top_column.T_eff)
        # надо брать не curr_configuration.top_column.L_x а посчитанный через интеграл! хоть они и отличаются на сильно
        L_scatter[i] = integralsService.calc_scatter_L(magnet_surface, L_x, new_cos_psi_range, tau_scatter_matrix)

        L_nu_scatter[i] = integralsService.calc_scatter_L_nu(magnet_surface,
                                                             curr_configuration.top_column.inner_surface,
                                                             curr_configuration.top_column.T_eff,
                                                             new_cos_psi_range, tau_scatter_matrix)

    # plot_package.plot_scripts.plot_L(L_scatter)

    PF_L_scatter = newService.get_PF(np.sum(L_scatter, axis=0))
    PF_L_nu_scatter = newService.get_PF(np.sum(L_nu_scatter, axis=0))
    print('finish calc scatter')
    # ------------------------------------------------- save txt -----------------------------------------------------
    make_save_values_file(curr_configuration, L_surfs, L_scatter, cur_path_data, cur_dir_saved)
    print('save_values')
    save_some_files(curr_configuration, obs_matrix, L_surfs, L_scatter, L_nu_surfs, L_nu_scatter, PF_L_nu_surfs,
                    cur_path_data, cur_dir_saved)
    print('some_files')
    save_new_way(L_surfs, L_scatter, L_nu_surfs, L_nu_scatter, cur_path_data)
    print('save_new_way')
    if figs_flag:
        print('start plot')
        t1 = time.perf_counter()
        # ------------------------------------------------- save figs -------------------------------------------------
        cur_path_fig = cur_path / 'fig'
        save.create_file_path(cur_path_fig)
        plot_package.plot_scripts.plot_total_luminosity_of_surfaces(L_surfs, cur_path_fig)

        ans = np.apply_along_axis(matrix.vec_to_angles, axis=1, arr=obs_matrix)
        observer_phi = ans[:, 0]
        observer_theta = ans[:, 1]
        plot_package.plot_scripts.plot_observer_angles(observer_phi, observer_theta, cur_path_fig)
        plot_package.plot_scripts.plot_Teff_to_ksi(curr_configuration.R_e, curr_configuration.top_column.T_eff,
                                                   curr_configuration.top_column.inner_surface.theta_range,
                                                   cur_path_fig)

        # -------------------------------------------------- L_nu --------------------------------------------------
        plot_package.plot_scripts.plot_PF_to_energy(L_nu_surfs, cur_path_fig)
        plot_package.plot_scripts.plot_L_nu(L_nu_surfs, cur_path_fig)
        plot_package.plot_scripts.plot_L_nu_all_in_one(L_nu_surfs, cur_path_fig)
        plot_package.plot_scripts.plot_L_nu_on_phase(L_nu_surfs, cur_path_fig)
        plot_package.plot_scripts.plot_L_nu_avg(L_nu_surfs, cur_path_fig)
        plot_package.plot_scripts.plot_L_nu_with_bb(L_nu_surfs, curr_configuration.top_column.T_eff, cur_path_fig)
        # ------------------------------------------------- L_scatter -------------------------------------------------
        plot_package.plot_scripts.plot_scatter_L(L_surfs, L_scatter, cur_path_fig)
        plot_package.plot_scripts.plot_PF_to_energy_with_scatter(L_nu_surfs, L_nu_scatter, cur_path_fig)
        plot_package.plot_scripts.plot_scatter_L_nu(L_nu_surfs, L_nu_scatter, cur_path_fig)
        # -------------------------------------------------------------------------------------------------------------
        t2 = time.perf_counter()
        print(f'{t2 - t1} seconds')
        print('finish plot')


if __name__ == '__main__':
    # ------------------------------------------------- start -------------------------------------------------------
    mu = 0.1e31
    beta_mu = 40
    mc2 = 100
    a_portion = 0.44
    phi_0 = 0

    theta_obs = 60

    calc_and_save_for_configuration(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, True)
