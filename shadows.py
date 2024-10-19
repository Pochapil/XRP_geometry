import numpy as np
import time

import config
import newService
from geometry import matrix


def intersection_with_sphere(surface, origin_phi, origin_theta, direction_vector):
    '''
        проверка на пересечение (затмение) НЗ
    '''
    # sphere x**2 + y**2 + z**2 == 1
    r = surface.surf_R_e / config.R_ns * np.sin(origin_theta) ** 2
    direction_x, direction_y, direction_z = matrix.vec_to_coord(direction_vector)
    origin_x, origin_y, origin_z = matrix.get_cartesian_from_spherical(r, origin_theta, origin_phi)

    def find_intersect_solution(a, b, c):
        # решение кв. уравнения; мы не падаем в комплексные числа! вернем -1
        if b ** 2 - 4 * a * c >= 0:
            t_1 = (-b + (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
            t_2 = (-b - (b ** 2 - 4 * a * c) ** (1 / 2)) / (2 * a)
            return t_1, t_2
        else:
            return -1, -1

    # коэф-ты получаются из уравнений рейтрейсинга
    a_sphere = direction_x ** 2 + direction_y ** 2 + direction_z ** 2
    b_sphere = 2 * (origin_x * direction_x + origin_y * direction_y + origin_z * direction_z)
    c_sphere = origin_x ** 2 + origin_y ** 2 + origin_z ** 2 - 1

    # solution for sphere
    t_sphere = find_intersect_solution(a_sphere, b_sphere, c_sphere)

    # print("t_sphere1 = %f,t_sphere2 = %f" % (t_sphere[0], t_sphere[1]))
    # если корни положительны то было пересечение. если отрицательны то пересечения нет в сторону наблюдателя
    if t_sphere[0] > 0 or t_sphere[1] > 0:
        return True

    return False


def get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta, direction_vector):
    '''
        очень затратная операция - необходимо параллелить

        есть аналитическое уравнение для пересечения с дипольной линией - полином 5 степени
        находим уравнение в сферических координатах.

        достаем корни, ищем положительные
    '''
    # вывод формулы был для 0 угла наблюдателя по фи (его смещали в выводе). поэтому находим phi_delta
    direction_phi, direction_theta = matrix.vec_to_angles(direction_vector)
    # вспомогательные переменные, были введены для упрощения аналитического вывода
    phi_delta = origin_phi - direction_phi
    eta = np.sin(direction_theta) / np.sin(origin_theta)
    cos_alpha = np.sin(origin_theta) * np.cos(phi_delta) * np.sin(direction_theta) + np.cos(origin_theta) * np.cos(
        direction_theta)

    c_x_5 = 1
    c_x_4 = 6 * cos_alpha
    c_x_3 = 3 + 12 * cos_alpha ** 2 - eta ** 4
    c_x_2 = 12 * cos_alpha + 8 * cos_alpha ** 3 - 4 * np.cos(phi_delta) * eta ** 3
    c_x_1 = 3 + 12 * cos_alpha ** 2 - 2 * eta ** 2 - 4 * np.cos(phi_delta) ** 2 * eta ** 2
    c_x_0 = 6 * cos_alpha - 4 * np.cos(phi_delta) * eta

    coefficients = [c_x_0, c_x_1, c_x_2, c_x_3, c_x_4, c_x_5]
    solutions = np.polynomial.polynomial.polyroots(coefficients)
    # coefficients = [c_x_5, c_x_4, c_x_3, c_x_2, c_x_1, c_x_0]
    # np.roots
    return solutions


def get_intersection_from_solution(r, origin_phi, origin_theta, obs_vector, solution):
    if solution.real > 0 and solution.imag == 0:
        direction_x, direction_y, direction_z = matrix.vec_to_coord(obs_vector)
        origin_x, origin_y, origin_z = matrix.get_cartesian_from_spherical(r, origin_theta, origin_phi)

        direction_t = solution.real * r  # ???? wtf = это был вывод (4.2 в дипломе) t = R_0 * x
        intersect_point = np.array([origin_x, origin_y, origin_z]) + direction_t * np.array(
            [direction_x, direction_y, direction_z])

        intersect_phi, intersect_theta = matrix.vec_to_angles(intersect_point)
    else:
        intersect_phi, intersect_theta = None, None
    return intersect_phi, intersect_theta


def get_intersection_phi_with_column(surface, intersect_phi):
    # top_column_intersect_phi_correct = (top_column_phi_range[0] <= intersect_phi <= top_column_phi_range[
    #     -1]) or (top_column_phi_range[0] <= intersect_phi + 2 * np.pi <= top_column_phi_range[-1])
    column_phi_range_begin, column_phi_range_end = surface.phi_range[0], surface.phi_range[-1]
    return (column_phi_range_begin <= intersect_phi <= column_phi_range_end) \
           or (column_phi_range_begin <= intersect_phi + 2 * np.pi <= column_phi_range_end) \
           or (column_phi_range_begin <= intersect_phi - 2 * np.pi <= column_phi_range_end)


def check_shadow_with_dipole(surface, phi_index, theta_index, obs_vector, solutions, top_column, bot_column):
    # checks for 1 element area!
    # obs_vector = direction_vector
    '''
       есть аналитическое уравнение для полинома дипольной линии 5 степени
       находим уравнение в сферических координатах.

       достаем корни
       ищем положительные
       находим пересечение с колонками
           если пересекается то истина
           если нет то ложь

       корни в магнитной СК - betta_mu
       '''
    origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]
    r = surface.surf_R_e / config.R_ns * np.sin(origin_theta) ** 2

    for solution in solutions:
        intersect_phi, intersect_theta = get_intersection_from_solution(r, origin_phi, origin_theta, obs_vector,
                                                                        solution)
        if intersect_phi is not None:
            # для верхней колонки:
            top_column_intersect_phi = get_intersection_phi_with_column(top_column, intersect_phi)

            # для нижней колонки:
            bot_column_intersect_phi = get_intersection_phi_with_column(bot_column, intersect_phi)

            if (intersect_theta < top_column.theta_range[-1] and top_column_intersect_phi) \
                    or (intersect_theta > bot_column.theta_range[-1] and bot_column_intersect_phi):
                return 0
    return 1


def get_tau_for_opacity_old(theta, R_e, M_accretion_rate, a_portion):
    tau = config.k * M_accretion_rate * newService.get_delta_distance(theta, R_e) / (
            newService.get_A_normal(theta, R_e, a_portion) * newService.get_free_fall_velocity(theta, R_e))
    return tau


def get_tau_for_opacity(phi, theta, R_e_of_atenuation_surf, M_accretion_rate, a_portion, obs_vector):
    '''очень маленькие cos = 0.005??'''
    # phi, theta - в точке пересечения; R_e =
    normal = matrix.newE_n(phi, theta)
    cos_alpha = np.dot(normal, obs_vector) / (np.linalg.norm(normal) * np.linalg.norm(obs_vector))
    tau = config.k * M_accretion_rate * newService.get_delta_distance(theta, R_e_of_atenuation_surf)
    tau /= (newService.get_A_normal(theta, R_e_of_atenuation_surf, a_portion)
            * newService.get_free_fall_velocity(theta, R_e_of_atenuation_surf))
    # print(f'{tau=}')
    tau /= cos_alpha
    # print(f'{cos_alpha=}')
    return tau


def get_tau_for_scatter_with_cos(theta_range, R_e_emission_surf, M_accretion_rate, a_portion, cos_alpha):
    tau = config.k * M_accretion_rate * newService.get_delta_distance(theta_range, R_e_emission_surf)
    tau /= (newService.get_A_normal(theta_range, R_e_emission_surf, a_portion)
            * newService.get_free_fall_velocity(theta_range, R_e_emission_surf))
    tau = tau[np.newaxis, :] / cos_alpha
    # print(f'cos_alpha_min = {np.min(cos_alpha)}, cos_alpha_max = {np.max(cos_alpha)}')
    # print(f'tau_min = {np.min(tau)}, max = {np.max(tau)}')
    return tau


def get_tau_with_dipole(surface, phi_index, theta_index, obs_vector, solutions, R_e_of_atenuation_surf, beta_mu,
                        M_accretion_rate, top_column, bot_column, a_portion):
    '''R_e_of_atenuation_surf = inner surf!!!!! я беру пересечения только с внутренней (пока что)'''
    origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]
    # здесь радиус основания, радиус излучающей так как r = для положения излучающей
    r = surface.surf_R_e / config.R_ns * np.sin(origin_theta) ** 2

    for solution in solutions:
        intersect_phi, intersect_theta = get_intersection_from_solution(r, origin_phi, origin_theta, obs_vector,
                                                                        solution)
        if intersect_phi is not None:
            # из геометрических выводов = конец магнитосферной линии
            theta_end = np.pi / 2 - np.arctan(np.tan(np.deg2rad(beta_mu)) * np.cos(intersect_phi))

            top_column_intersect_theta_correct = intersect_theta < theta_end
            bot_column_intersect_theta_correct = intersect_theta > theta_end

            top_column_intersect_phi = get_intersection_phi_with_column(top_column, intersect_phi)
            bot_column_intersect_phi = get_intersection_phi_with_column(bot_column, intersect_phi)

            intersection_condition = (top_column_intersect_phi and top_column_intersect_theta_correct) or (
                    bot_column_intersect_phi and bot_column_intersect_theta_correct)

            if intersection_condition:
                tau = get_tau_for_opacity(intersect_phi, intersect_theta, R_e_of_atenuation_surf, M_accretion_rate,
                                          a_portion, obs_vector)
                if tau > config.tau_cutoff:
                    return np.exp(-1 * tau)
                else:
                    return 1
    return 1


if __name__ == '__main__':
    a_portion = 0.44
    phi_range = np.linspace(-np.pi * a_portion, np.pi * a_portion, config.N_phi_accretion)
    theta_range = np.linspace(0.4, np.pi / 2, config.N_theta_accretion)

    theta_obs = 60
    theta_obs_rad = np.deg2rad(theta_obs)

    x = np.sin(theta_obs_rad) * np.cos(0)
    y = np.sin(theta_obs_rad) * np.sin(0)
    z = np.cos(theta_obs_rad)
    direction_vector = np.array([x, y, z])

    t1 = time.perf_counter()
    for i in range(45):
        for origin_phi in phi_range:
            for origin_theta in theta_range:
                get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta, direction_vector)
    t2 = time.perf_counter()
    print(t2 - t1)
