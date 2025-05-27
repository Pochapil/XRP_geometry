import numpy as np

import config


def check_adjacent(arr, i, j, di_arr, dj_arr, val=0):
    for di, dj in zip(di_arr, dj_arr):
        try:
            if arr[i + di, j + dj] > val:
                return False
        except:
            pass
    return True


# cos_psi_rotation_matrix = np.ones((config.N_phase, config.N_phi_accretion, config.N_theta_accretion))

N_phi_accretion = 10
N_theta_accretion = 10

cos_psi_rotation_matrix = np.ones((2, N_phi_accretion, N_theta_accretion))
tensor_shadows_NS = np.ones_like(cos_psi_rotation_matrix)
tensor_shadows_columns = np.ones_like(cos_psi_rotation_matrix)
tensor_tau = np.ones_like(cos_psi_rotation_matrix)

phase_index = 1
trial = cos_psi_rotation_matrix[phase_index]

# fill first
for phi_index in range(0, N_phi_accretion - 1, 2):
    for theta_index in range(0, N_theta_accretion - 1, 2):
        trial[phi_index, theta_index] = round(np.random.rand(), 2)

di_arr = [-1, -1, 1, 1]
dj_arr = [-1, 1, -1, 1]

for phi_index in range(1, N_phi_accretion - 1, 2):
    for theta_index in range(1, N_theta_accretion - 1, 2):

        if check_adjacent(trial, phi_index, theta_index, di_arr, dj_arr, 0.7):
            trial[phi_index, theta_index] = 3
        else:
            trial[phi_index, theta_index] = -1

di_arr = [-1, 1, 0, 0]
dj_arr = [0, 0, -1, 1]

for phi_index in range(0, N_phi_accretion - 1, 2):
    for theta_index in range(1, N_theta_accretion - 1, 2):
        if check_adjacent(trial, phi_index, theta_index, di_arr, dj_arr, 0.8):
            trial[phi_index, theta_index] = 4
        else:
            trial[phi_index, theta_index] = -2

for phi_index in range(1, N_phi_accretion - 1, 2):
    for theta_index in range(0, N_theta_accretion - 1, 2):
        # print(phi_index, theta_index)
        if check_adjacent(trial, phi_index, theta_index, di_arr, dj_arr, 0.8):
            trial[phi_index, theta_index] = 4
        else:
            trial[phi_index, theta_index] = -2

print(trial)


a = np.linspace(0,1,10)

def f(a):
    a[0] = -100

f(a)
print(a)

# TODO last rows columns

#
# def calc_shadows_and_tau_new(curr_configuration: accretingNS.AccretingPulsarConfiguration, surface, obs_matrix,
#                              mask_flag=False):
#     '''
#     использую в расчетах!!!! по сути самая главная функция
#
#     функция расчитывает матрицу косинусов на каждой фазе для направлений наблюдателя
#
#     mask_flag --- для расчета по магнитным линиям - их надо обрезать
#     '''
#
#     def check_adjacent(arr, i, j, di_arr, dj_arr, val=0):
#         for di, dj in zip(di_arr, dj_arr):
#             try:
#                 if arr[i + di, j + dj] > val:
#                     return False
#             except:
#                 pass
#         return True
#
#     # тензор косинусов между нормалью на площадке и направлением на наблюдателя. размер phase x phi x theta
#     # умножаем скалярно phi x theta x 3 на phase x 3 (по последнему индексу) и делаем reshape.
#     # phase x phi x theta
#     # i,j = phi, theta; l=x,y,z
#     cos_psi_rotation_matrix = np.einsum('ijl,tl->tij', surface.array_normal, obs_matrix)
#     # print(f'obs_matrix = {np.max(np.linalg.norm(obs_matrix, axis=1))}')
#     # print(f'array_normal_max = {np.max(np.linalg.norm(surface.array_normal, axis=2))}')
#
#     # возможно хранить будем все матрицы. чтобы расчитывать вклады
#     # для разных матриц можем посчитать L и посмотреть какой вклад будет.
#     tensor_shadows_NS = np.ones_like(cos_psi_rotation_matrix)
#     tensor_shadows_columns = np.ones_like(cos_psi_rotation_matrix)
#     tensor_shadows_optimized = np.ones_like(cos_psi_rotation_matrix)
#
#     tensor_tau = np.ones_like(cos_psi_rotation_matrix)
#     new_cos_psi_range = cos_psi_rotation_matrix.copy()
#
#     tensor_tau_save = np.empty_like(cos_psi_rotation_matrix)
#     tensor_alpha_save = np.empty_like(cos_psi_rotation_matrix)
#
#     if mask_flag:
#         mask = np.zeros_like(new_cos_psi_range).astype(bool)
#         mask += surface.mask_array  # происходит broadcast в матрицу размера phase x phi x theta
#         new_cos_psi_range[mask] = 0
#     # print(surface)
#     # for phase_index in range(config.N_phase):
#     # для возможности распараллелить еще больше чем на 4 поверхности (разрезая наблюдателя на чанки)
#     for phase_index in range(obs_matrix.shape[0]):
#
#         N_phi_accretion = config.N_phi_accretion
#         N_theta_accretion = config.N_theta_accretion
#
#         # trial = cos_psi_rotation_matrix[phase_index]
#
#         def calc(new_cos_psi_range, tensor_shadows_NS, tensor_shadows_columns,
#                  phase_index, phi_index, theta_index,
#                  surface, obs_matrix, curr_configuration,
#                  tensor_tau, tensor_tau_save, tensor_alpha_save):
#
#             '''тензоры передаются по ссылке. так что изменяемые объекты'''
#
#             if new_cos_psi_range[phase_index, phi_index, theta_index] > 0:
#                 origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]
#                 # сначала проверяем на затмение НЗ
#                 if config.NS_shadow_flag and shadows.intersection_with_sphere(surface, origin_phi, origin_theta,
#                                                                               obs_matrix[phase_index]):
#                     tensor_shadows_NS[phase_index, phi_index, theta_index] = 0
#                 # иначе тяжелые вычисления
#                 else:
#                     # расчет для полинома - находим корни для пересечения с внутр поверхностью на магн линии!
#                     solutions = shadows.get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta,
#                                                                               obs_matrix[phase_index])
#                     # сортируем по действительной части, чтобы
#                     solutions = sorted(list(solutions), key=lambda x: x.real)  # [::-1]
#                     # расчитываем затмение колонкой
#                     tensor_shadows_columns[phase_index, phi_index, theta_index] = \
#                         shadows.check_shadow_with_dipole(surface, phi_index, theta_index,
#                                                          obs_matrix[phase_index], solutions,
#                                                          curr_configuration.top_column_for_calc.inner_surface,
#                                                          curr_configuration.bot_column_for_calc.inner_surface)
#                     if tensor_shadows_columns[phase_index, phi_index, theta_index] > 0:
#                         if config.flag_attenuation_above_shock:
#                             # если затмения нет то считаем ослабление тау с внутренними магнитными линиями!!
#                             # здесь нужно запомнить тау и альфа!
#
#                             buf = shadows.get_tau_with_dipole(surface, phi_index, theta_index,
#                                                               obs_matrix[phase_index],
#                                                               solutions, curr_configuration)
#
#                             tensor_tau[phase_index, phi_index, theta_index] = buf[0]
#                             tensor_tau_save[phase_index, phi_index, theta_index] = buf[1]
#                             tensor_alpha_save[phase_index, phi_index, theta_index] = buf[2]
#             else:
#                 # если косинус < 0 -> поверхность излучает от наблюдателя и мы не получим вклад в интеграл
#                 new_cos_psi_range[phase_index, phi_index, theta_index] = 0
#
#         # fill first
#         for phi_index in range(0, N_phi_accretion - 1, 2):
#             for theta_index in range(0, N_theta_accretion - 1, 2):
#                 calc(new_cos_psi_range, tensor_shadows_NS, tensor_shadows_columns,
#                      phase_index, phi_index, theta_index,
#                      surface, obs_matrix, curr_configuration,
#                      tensor_tau, tensor_tau_save, tensor_alpha_save)
#
#         trial = new_cos_psi_range[phase_index] * tensor_shadows_NS[phase_index] \
#                 * tensor_shadows_columns[phase_index] * tensor_tau[phase_index]
#
#         di_arr = [-1, -1, 1, 1]
#         dj_arr = [-1, 1, -1, 1]
#         # fill x
#         for phi_index in range(1, N_phi_accretion - 1, 2):
#             for theta_index in range(1, N_theta_accretion - 1, 2):
#                 if check_adjacent(trial, phi_index, theta_index, di_arr, dj_arr, 0):
#                     trial[phi_index, theta_index] = 0
#                     tensor_shadows_optimized[phase_index, phi_index, theta_index] = 0
#                 else:
#                     trial[phi_index, theta_index] = 10
#
#         di_arr = [-1, 1, 0, 0]
#         dj_arr = [0, 0, -1, 1]
#         # fill o
#         for phi_index in range(0, N_phi_accretion - 1, 2):
#             for theta_index in range(1, N_theta_accretion - 1, 2):
#                 if check_adjacent(trial, phi_index, theta_index, di_arr, dj_arr, 0):
#                     trial[phi_index, theta_index] = 0
#                     tensor_shadows_optimized[phase_index, phi_index, theta_index] = 0
#                 else:
#                     trial[phi_index, theta_index] = 10
#         # fill o
#         for phi_index in range(1, N_phi_accretion - 1, 2):
#             for theta_index in range(0, N_theta_accretion - 1, 2):
#                 # print(phi_index, theta_index)
#                 if check_adjacent(trial, phi_index, theta_index, di_arr, dj_arr, 0):
#                     trial[phi_index, theta_index] = 0
#                     tensor_shadows_optimized[phase_index, phi_index, theta_index] = 0
#                 else:
#                     trial[phi_index, theta_index] = 10
#
#         # fill not
#         idx = np.argwhere(trial > 0)
#         for phi_index, theta_index in idx:
#             calc(new_cos_psi_range, tensor_shadows_NS, tensor_shadows_columns,
#                  phase_index, phi_index, theta_index,
#                  surface, obs_matrix, curr_configuration,
#                  tensor_tau, tensor_tau_save, tensor_alpha_save)
#         # last rows columns
#         for phi_index in range(N_phi_accretion):
#             theta_index = -1
#             calc(new_cos_psi_range, tensor_shadows_NS, tensor_shadows_columns,
#                  phase_index, phi_index, theta_index,
#                  surface, obs_matrix, curr_configuration,
#                  tensor_tau, tensor_tau_save, tensor_alpha_save)
#
#         for theta_index in range(N_theta_accretion):
#             phi_index = -1
#             calc(new_cos_psi_range, tensor_shadows_NS, tensor_shadows_columns,
#                  phase_index, phi_index, theta_index,
#                  surface, obs_matrix, curr_configuration,
#                  tensor_tau, tensor_tau_save, tensor_alpha_save)
#
#     new_cos_psi_range = new_cos_psi_range * tensor_shadows_NS * tensor_shadows_columns * tensor_tau * tensor_shadows_optimized
#     # if tensor_tau[tensor_tau != 1].shape[0] > 0:
#     #     print(np.max(tensor_tau[tensor_tau != 1]))
#     return new_cos_psi_range, tensor_tau_save, tensor_alpha_save