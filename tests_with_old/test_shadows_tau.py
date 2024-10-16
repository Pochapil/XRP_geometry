import numpy as np

import accretingNS
import config
from geometricTask import matrix
import pathService
import shadows

mu = 0.1e31
beta_mu = 40
mc2 = 100
a_portion = 0.44
phi_0 = 0

'''переделать phi_0 !!!!'''

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

# возможно хранить будем все матрицы.
accr_col_surfs = [curr_configuration.top_column.inner_surface, curr_configuration.top_column.outer_surface,
                  curr_configuration.bot_column.inner_surface, curr_configuration.bot_column.outer_surface]

for surface in accr_col_surfs:
    # тензор косинусов между нормалью и направлением на наблюдателя размером phase x phi x theta
    # умножаем скалярно phi x theta x 3 на phase x 3 (по последнему индексу) и делаем reshape.
    cos_psi_rotation_matrix = np.einsum('ijl,tl->tij', surface.array_normal, obs_matrix)

    tensor_shadows = np.zeros_like(cos_psi_rotation_matrix)
    tensor_tau = np.zeros_like(cos_psi_rotation_matrix)

    old_cos_psi_range = cos_psi_rotation_matrix.copy()
    new_cos_psi_range = cos_psi_rotation_matrix.copy()

    for phase_index in range(config.N_phase):
        for phi_index in range(config.N_phi_accretion):
            for theta_index in range(config.N_theta_accretion):
                if old_cos_psi_range[phase_index, phi_index, theta_index] > 0:

                    origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]
                    solutions = shadows.get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta,
                                                                              obs_matrix[phase_index])

                    tensor_shadows[phase_index, phi_index, theta_index] = shadows.check_shadow_with_dipole(surface,
                                                                                                           phi_index,
                                                                                                           theta_index,
                                                                                                           obs_matrix[
                                                                                                               phase_index],
                                                                                                           solutions,
                                                                                                           curr_configuration.top_column.inner_surface,
                                                                                                           curr_configuration.bot_column.inner_surface)
                    tensor_tau[phase_index, phi_index, theta_index] = shadows.get_tau_with_dipole(surface, phi_index,
                                                                                                  theta_index,
                                                                                                  obs_matrix[
                                                                                                      phase_index],
                                                                                                  solutions,
                                                                                                  surface.surf_R_e,
                                                                                                  beta_mu,
                                                                                                  curr_configuration.M_accretion_rate,
                                                                                                  curr_configuration.top_column.inner_surface,
                                                                                                  curr_configuration.bot_column.inner_surface,
                                                                                                  a_portion)

                    new_cos_psi_range[phase_index, phi_index, theta_index] = new_cos_psi_range[
                                                                                 phase_index, phi_index, theta_index] * \
                                                                             tensor_shadows[
                                                                                 phase_index, phi_index, theta_index] * \
                                                                             tensor_tau[
                                                                                 phase_index, phi_index, theta_index]

                    top_column_phi_range = curr_configuration.top_column.inner_surface.phi_range
                    top_column_theta_range = curr_configuration.top_column.inner_surface.theta_range
                    bot_column_phi_range = curr_configuration.bot_column.inner_surface.phi_range
                    bot_column_theta_range = curr_configuration.bot_column.inner_surface.theta_range

                    old_cos_psi_range[phase_index, phi_index, theta_index] *= shadows.get_vals(surface, phi_index,
                                                                                               theta_index,
                                                                                               obs_matrix[phase_index],
                                                                                               top_column_phi_range,
                                                                                               bot_column_phi_range,
                                                                                               top_column_theta_range,
                                                                                               bot_column_theta_range,
                                                                                               beta_mu,
                                                                                               curr_configuration.M_accretion_rate,
                                                                                               a_portion)
                else:
                    old_cos_psi_range[phase_index, phi_index, theta_index] = 0
                    new_cos_psi_range[phase_index, phi_index, theta_index] = 0

    result = old_cos_psi_range - new_cos_psi_range
    print(result[result > 0])
    print(result[result > 0].shape)
    # print(np.max(result[result > 0]))
print('done')

folder = 'D:\MyProgs\Python\Diplom\\new_magnet_lines\\new_condition\\new_data\mu=0.1e31\\tau\data\cos\i=60 betta_mu=40\mc2=100\\a=0.44 fi_0=280\\top_inner\\'
print(folder)

file_name = 'save_cos_' + 'top_inner' + f'_{phase_index}_phase' + '.txt'
print(file_name)

z1 = np.loadtxt(folder + file_name)
z1_old = z1 - old_cos_psi_range
z1_new = z1 - new_cos_psi_range
print(z1_old[z1_old > 0])
print(z1_old[z1_old > 0].shape)

print(z1_new[z1_new > 0])
print(z1_new[z1_new > 0].shape)

#    for phase_index in range(config.N_phase):
#         for phi_index in range(config.N_phi_accretion):
#             for theta_index in range(config.N_theta_accretion):
#                 if old_cos_psi_range[phase_index, phi_index, theta_index] > 0:
#
#                     origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]
#                     solutions = shadows.get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta,
#                                                                               obs_matrix[phase_index])
#
#                     if shadows.intersection_with_sphere(surface, origin_phi, origin_theta, obs_matrix[phase_index]):
#                         tensor_shadows[phase_index, phi_index, theta_index] = 0
#                     else:
#                         tensor_shadows[phase_index, phi_index, theta_index] = shadows.check_shadow_with_dipole(surface,
#                                                                                                                phi_index,
#                                                                                                                theta_index,
#                                                                                                                obs_matrix[
#                                                                                                                    phase_index],
#                                                                                                                solutions,
#                                                                                                                curr_configuration.top_column.inner_surface,
#                                                                                                                curr_configuration.bot_column.inner_surface)
#
#                     tensor_tau[phase_index, phi_index, theta_index] = shadows.get_tau_with_dipole(surface, phi_index,
#                                                                                                   theta_index,
#                                                                                                   obs_matrix[
#                                                                                                       phase_index],
#                                                                                                   solutions,
#                                                                                                   surface.surf_R_e,
#                                                                                                   beta_mu,
#                                                                                                   curr_configuration.M_accretion_rate,
#                                                                                                   curr_configuration.top_column.inner_surface,
#                                                                                                   curr_configuration.bot_column.inner_surface,
#                                                                                                   a_portion)
#
#                     new_cos_psi_range[phase_index, phi_index, theta_index] = new_cos_psi_range[
#                                                                                  phase_index, phi_index, theta_index] * \
#                                                                              tensor_shadows[
#                                                                                  phase_index, phi_index, theta_index] * \
#                                                                              tensor_tau[
#                                                                                  phase_index, phi_index, theta_index]
#
#                     top_column_phi_range = curr_configuration.top_column.inner_surface.phi_range
#                     top_column_theta_range = curr_configuration.top_column.inner_surface.theta_range
#                     bot_column_phi_range = curr_configuration.bot_column.inner_surface.phi_range
#                     bot_column_theta_range = curr_configuration.bot_column.inner_surface.theta_range
#
#                     old_cos_psi_range[phase_index, phi_index, theta_index] *= shadows.get_vals(surface, phi_index,
#                                                                                                theta_index,
#                                                                                                obs_matrix[phase_index],
#                                                                                                top_column_phi_range,
#                                                                                                bot_column_phi_range,
#                                                                                                top_column_theta_range,
#                                                                                                bot_column_theta_range,
#                                                                                                beta_mu,
#                                                                                                curr_configuration.M_accretion_rate,
#                                                                                                a_portion)
#                 else:
#                     old_cos_psi_range[phase_index, phi_index, theta_index] = 0
#                     new_cos_psi_range[phase_index, phi_index, theta_index] = 0
#
#         folder = 'D:\MyProgs\Python\Diplom\\new_magnet_lines\\new_condition\\new_data\mu=0.1e31\\tau\data\cos\i=60 betta_mu=40\mc2=100\\a=0.44 fi_0=280\\top_inner\\'
#         print(folder)
#
#         file_name = 'save_cos_' + 'top_inner' + f'_{phase_index}_phase' + '.txt'
#         print(file_name)
#
#         z1 = np.loadtxt(folder + file_name)
#         z1_old = z1[phase_index] - old_cos_psi_range[phase_index]
#         z1_new = z1[phase_index] - new_cos_psi_range[phase_index]
#         print(z1_old[z1_old > 0])
#         print(z1_old[z1_old > 0].shape)
#
#         print(z1_new[z1_new > 0])
#         print(z1_new[z1_new > 0].shape)
#         if z1_new[z1_new > 0].shape[0] > 0:
#             print(np.max(z1_new[z1_new > 0]))
#
#     result = old_cos_psi_range - new_cos_psi_range
#     print(result[result > 0])
#     print(result[result > 0].shape)
#
#     print(new_cos_psi_range)


''' 
                                                        функция 
def get_vals(surface, phi_index, theta_index, direction_vector,
             top_column_phi_range, bot_column_phi_range, top_column_theta_range, bot_column_theta_range,
             betta_mu, M_accretion_rate, a_portion):
    origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]
    r = surface.surf_R_e / config.R_ns * np.sin(origin_theta) ** 2

    direction_x, direction_y, direction_z = vec_to_coord(direction_vector)
    origin_x, origin_y, origin_z = get_cartesian_from_spherical(r, origin_theta, origin_phi)

    solutions = get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta, direction_vector)

    # мб добавить условие на минимум t.real -- and solution.real > 1e-2
    for solution in solutions:
        if solution.real > 0 and solution.imag == 0:
            direction_t = solution.real * r
            intersect_point = np.array([origin_x, origin_y, origin_z]) + direction_t * np.array(
                [direction_x, direction_y, direction_z])

            intersect_phi, intersect_theta = vec_to_angles(intersect_point)
            if intersect_phi < 0:
                intersect_phi += 2 * np.pi

            intersect_r = (intersect_point[0] ** 2 + intersect_point[1] ** 2 + intersect_point[2] ** 2) ** (1 / 2)

            # if not (abs(R_e / config.R_ns * np.sin(intersect_theta) ** 2 - intersect_r) < 0.001):
            #     print(abs(R_e / config.R_ns * np.sin(intersect_theta) ** 2 - intersect_r))

            # curr = abs(R_e / config.R_ns * np.sin(intersect_theta) ** 2 - intersect_r)
            # новые линии магнитосферы
            # theta_end = np.pi / 2 - betta_mu * np.cos(intersect_phi) - ОШИБКА !!
            theta_end = np.pi / 2 - np.arctan(np.tan(np.deg2rad(betta_mu)) * np.cos(intersect_phi))

            # для верхней колонки:
            top_column_intersect_phi_correct = (top_column_phi_range[0] <= intersect_phi <= top_column_phi_range[
                -1]) or (0 <= intersect_phi <= top_column_phi_range[-1] - 2 * np.pi)
            top_column_intersect_theta_correct = intersect_theta < theta_end

            # для нижней колонки:
            bot_column_intersect_phi_correct = (bot_column_phi_range[0] <= intersect_phi <= bot_column_phi_range[
                -1]) or (0 <= intersect_phi <= bot_column_phi_range[-1] - 2 * np.pi)
            bot_column_intersect_theta_correct = intersect_theta > theta_end

            # проверяем на пересечение с колонками
            # intersect_r_correct = intersect_r > ksi_shock
            # if not intersect_r_correct and (top_column_intersect_phi_correct or bot_column_intersect_phi_correct):
            #     return 0

            if (intersect_theta < top_column_theta_range[-1] and top_column_intersect_phi_correct) or (
                    intersect_theta > bot_column_theta_range[-1] and bot_column_intersect_phi_correct):
                return 0

            intersection_condition = (top_column_intersect_phi_correct and top_column_intersect_theta_correct) or (
                    bot_column_intersect_phi_correct and bot_column_intersect_theta_correct)

            # if not (top_column_intersect_phi_correct and top_column_intersect_theta_correct) and (
            #         bot_column_intersect_phi_correct and bot_column_intersect_theta_correct):
            #     print('пересекли bot в точке')
            #     print(intersect_phi / config.grad_to_rad, intersect_theta / config.grad_to_rad)
            #     print(f'theta_end={theta_end / config.grad_to_rad}')

            if config.tau_flag:
                if intersection_condition:
                    tau = get_tau_for_opacity(intersect_theta, surface.surf_R_e, M_accretion_rate, a_portion)
                    if tau > config.tau_cutoff:
                        return np.exp(-1 * tau)
                    else:
                        return 1

            elif config.opacity_above_shock > 0:
                if intersection_condition:
                    return 1 - config.opacity_above_shock

            else:
                return 1
    return 1
    '''

'''
эта версия дает много 0!!
   cos_psi_rotation_matrix = np.einsum('ijl,tl->tij', surface.array_normal, obs_matrix)

    tensor_shadows = np.ones_like(cos_psi_rotation_matrix)
    tensor_tau = np.ones_like(cos_psi_rotation_matrix)

    old_cos_psi_range = cos_psi_rotation_matrix.copy()
    new_cos_psi_range = cos_psi_rotation_matrix.copy()

    for phase_index in range(config.N_phase):
        for phi_index in range(config.N_phi_accretion):
            for theta_index in range(config.N_theta_accretion):
                if old_cos_psi_range[phase_index, phi_index, theta_index] > 0:

                    origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]

                    if shadows.intersection_with_sphere(surface, origin_phi, origin_theta, obs_matrix[phase_index]):
                        tensor_shadows[phase_index, phi_index, theta_index] = 0
                    else:
                        solutions = shadows.get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta,
                                                                                  obs_matrix[phase_index])
                        tensor_shadows[phase_index, phi_index, theta_index] = \
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

                    new_cos_psi_range[phase_index, phi_index, theta_index] = \
                        new_cos_psi_range[phase_index, phi_index, theta_index] * \
                        tensor_shadows[phase_index, phi_index, theta_index] * \
                        tensor_tau[phase_index, phi_index, theta_index]

                    top_column_phi_range = curr_configuration.top_column.inner_surface.phi_range
                    top_column_theta_range = curr_configuration.top_column.inner_surface.theta_range
                    bot_column_phi_range = curr_configuration.bot_column.inner_surface.phi_range
                    bot_column_theta_range = curr_configuration.bot_column.inner_surface.theta_range

                    old_cos_psi_range[phase_index, phi_index, theta_index] *= \
                        shadows.get_vals(surface, phi_index, theta_index, obs_matrix[phase_index],
                                         top_column_phi_range, bot_column_phi_range,
                                         top_column_theta_range, bot_column_theta_range,
                                         beta_mu, curr_configuration.M_accretion_rate, a_portion)
                else:
                    new_cos_psi_range[phase_index, phi_index, theta_index] = 0
        folder = 'D:\MyProgs\Python\Diplom\\new_magnet_lines\\new_condition\\new_data\mu=0.1e31\\tau\data\cos\i=60 betta_mu=40\mc2=100\\a=0.44 fi_0=280\\top_inner\\'
        print(folder)

        file_name = 'save_cos_' + 'top_inner' + f'_{phase_index}_phase' + '.txt'
        print(file_name)

        z1 = np.loadtxt(folder + file_name)
        z1_old = z1[phase_index] - old_cos_psi_range[phase_index]
        z1_new = z1[phase_index] - new_cos_psi_range[phase_index]
        # print(z1_old[z1_old > 0].shape)

        print(z1_new[z1_new > 0].shape)
        if z1_new[z1_new > 0].shape[0] > 0:
            print(np.max(z1_new[z1_new > 0]))
    print('done')

'''
