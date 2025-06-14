import numpy as np
import time
from numba import njit
from scipy.optimize import newton

import accretingNS
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


@njit(cache=True, fastmath=True)
def get_solutions_for_dipole_magnet_lines(origin_phi, origin_theta, direction_vector):
    '''
        очень затратная операция - необходимо параллелить

        есть аналитическое уравнение для пересечения с дипольной линией - полином 5 степени
        находим уравнение в сферических координатах.
        достаем корни, ищем положительные

        вывод - диплом 4.3.2 Затмения, вызванные аккреционной колонкой

        по идее для пересечений верхней колонки с верхними линиями тоже пойдет.
        так как сокращается в расчете Re или Re+delta

        НО если искать пересечения между внешней и внутр или внутр и внеш то нужно менять !!!!
    '''
    direction_phi, direction_theta = matrix.vec_to_angles(direction_vector)
    # вспомогательные переменные, были введены для упрощения аналитического вывода
    # вывод формулы был для 0 угла наблюдателя по фи (его смещали в выводе). поэтому находим phi_delta
    phi_delta = origin_phi - direction_phi
    eta = np.sin(direction_theta) / np.sin(origin_theta)

    cos_phi = np.cos(phi_delta)
    cos_alpha = np.sin(origin_theta) * cos_phi * np.sin(direction_theta) + np.cos(origin_theta) * np.cos(
        direction_theta)

    # коэффициенты в полиноме
    c_x_5 = 1
    c_x_4 = 6 * cos_alpha
    c_x_3 = 3 + 12 * cos_alpha ** 2 - eta ** 4
    c_x_2 = 12 * cos_alpha + 8 * cos_alpha ** 3 - 4 * cos_phi * eta ** 3
    c_x_1 = 3 + 12 * cos_alpha ** 2 - 2 * eta ** 2 - 4 * cos_phi ** 2 * eta ** 2
    c_x_0 = 6 * cos_alpha - 4 * cos_phi * eta

    # теперь коэффициенты надо укзывать в порядке с 0 до 5 степени
    # coefficients = np.array([c_x_0, c_x_1, c_x_2, c_x_3, c_x_4, c_x_5])
    # solutions = np.polynomial.polynomial.polyroots(coefficients)

    # для numba
    # coefficients = np.array([c_x_5, c_x_4, c_x_3, c_x_2, c_x_1, c_x_0])
    coefficients = np.array([c_x_5, c_x_4, c_x_3, c_x_2, c_x_1, c_x_0], dtype=np.complex64)
    solutions = np.roots(coefficients)
    return solutions


# def get_solutions_for_dipole_magnet_lines_prev_sol(origin_phi, origin_theta, direction_vector, prev_solution=None):
#     '''
#         попытка в оптимизацию, но стало хуже
#     '''
#     direction_phi, direction_theta = matrix.vec_to_angles(direction_vector)
#     # вспомогательные переменные, были введены для упрощения аналитического вывода
#     # вывод формулы был для 0 угла наблюдателя по фи (его смещали в выводе). поэтому находим phi_delta
#     phi_delta = origin_phi - direction_phi
#     eta = np.sin(direction_theta) / np.sin(origin_theta)
#
#     cos_phi = np.cos(phi_delta)
#     cos_alpha = np.sin(origin_theta) * cos_phi * np.sin(direction_theta) + np.cos(origin_theta) * np.cos(
#         direction_theta)
#
#     # коэффициенты в полиноме
#     c_x_5 = 1
#     c_x_4 = 6 * cos_alpha
#     c_x_3 = 3 + 12 * cos_alpha ** 2 - eta ** 4
#     c_x_2 = 12 * cos_alpha + 8 * cos_alpha ** 3 - 4 * cos_phi * eta ** 3
#     c_x_1 = 3 + 12 * cos_alpha ** 2 - 2 * eta ** 2 - 4 * cos_phi ** 2 * eta ** 2
#     c_x_0 = 6 * cos_alpha - 4 * cos_phi * eta
#
#     # теперь коэффициенты надо укзывать в порядке с 0 до 5 степени
#     # coefficients = np.array([c_x_0, c_x_1, c_x_2, c_x_3, c_x_4, c_x_5])
#     # solutions = np.polynomial.polynomial.polyroots(coefficients)
#
#     # для numba
#     # coefficients = np.array([c_x_5, c_x_4, c_x_3, c_x_2, c_x_1, c_x_0])
#     coefficients = np.array([c_x_5, c_x_4, c_x_3, c_x_2, c_x_1, c_x_0], dtype=np.complex64)
#     if prev_solution is None:
#         return np.roots(coefficients)
#     else:
#         # Уточнение предыдущих корней
#         roots_refined = [newton(lambda x: np.polyval(coefficients[::-1], x), x0=r) for r in prev_solution]
#         return np.array(roots_refined)


# @njit(cache=True, fastmath=True)
# @njit
def real_roots_numba(origin_phi, origin_theta, direction_vector):
    # попытка deepseek
    '''работает с такой же скоростью'''
    direction_phi, direction_theta = matrix.vec_to_angles(direction_vector)
    # вспомогательные переменные, были введены для упрощения аналитического вывода
    # вывод формулы был для 0 угла наблюдателя по фи (его смещали в выводе). поэтому находим phi_delta
    phi_delta = origin_phi - direction_phi
    eta = np.sin(direction_theta) / np.sin(origin_theta)

    cos_phi = np.cos(phi_delta)
    cos_alpha = np.sin(origin_theta) * cos_phi * np.sin(direction_theta) + np.cos(origin_theta) * np.cos(
        direction_theta)

    # коэффициенты в полиноме
    c_x_5 = 1
    c_x_4 = 6 * cos_alpha
    c_x_3 = 3 + 12 * cos_alpha ** 2 - eta ** 4
    c_x_2 = 12 * cos_alpha + 8 * cos_alpha ** 3 - 4 * cos_phi * eta ** 3
    c_x_1 = 3 + 12 * cos_alpha ** 2 - 2 * eta ** 2 - 4 * cos_phi ** 2 * eta ** 2
    c_x_0 = 6 * cos_alpha - 4 * cos_phi * eta

    coeffs = [c_x_5, c_x_4, c_x_3, c_x_2, c_x_1, c_x_0]  # [::-1]

    # Приводим коэффициенты к вещественному типу
    coeffs = np.array(coeffs, dtype=np.float64)

    # # Создаем companion matrix
    # n = len(coeffs) - 1
    # mat = np.zeros((n, n), dtype=np.float64)
    # mat[0, :] = -coeffs[1:] / coeffs[0]
    # mat[1:, :-1] = np.eye(n - 1)
    #
    # # Находим собственные значения (корни)
    # roots = np.linalg.eigvals(mat)

    n = len(coeffs) - 1

    # Создаем companion matrix (правильная форма)
    mat = np.zeros((n, n), dtype=np.float64)
    mat[0, :] = -coeffs[1:] / coeffs[0]  # Первая строка

    # https://ru.stackoverflow.com/questions/419456/%D0%A1-%D0%BF%D0%BE%D0%BC%D0%BE%D1%89%D1%8C%D1%8E-numpy-%D0%B7%D0%B0%D0%BF%D0%BE%D0%BB%D0%BD%D0%B8%D1%82%D1%8C-%D1%81%D0%BE%D1%81%D0%B5%D0%B4%D0%BD%D0%B8%D0%B5-%D1%81-%D0%B3%D0%BB%D0%B0%D0%B2%D0%BD%D0%BE%D0%B9-%D0%B4%D0%B8%D0%B0%D0%B3%D0%BE%D0%BD%D0%B0%D0%BB%D0%B8
    k = -1
    rows, cols = np.indices(mat.shape)
    row_values = np.diag(rows, k=k)
    col_values = np.diag(cols, k=k)
    mat[row_values, col_values] = 1

    # for i in range(1, n):
    #     mat[i, i - 1] = 1.0  # Единицы на поддиагонали

    # Вычисляем собственные значения (корни)
    roots = np.linalg.eigvals(mat)

    # Фильтруем вещественные корни
    return roots[np.abs(roots.imag) < 1e-7].real


def get_intersection_from_solution(r, origin_phi, origin_theta, obs_vector, solution):
    # ищем положительные корни и при этом вещественные
    # if solution.real > 0 and solution.imag == 0:
    # return -pi/2, pi/2 angles
    if solution.real > 0 and np.abs(solution.imag) < 1e-6:
        # direction - направление на наблюдателя
        direction_x, direction_y, direction_z = matrix.vec_to_coord(obs_vector)
        # origin - откуда бьет луч. площадка на поверхности
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

    # смотрим чтобы intersect_phi лежал внутри диапазона по фи. + двигаем на +- 2pi
    column_phi_range_begin, column_phi_range_end = surface.phi_range[0], surface.phi_range[-1]
    return (column_phi_range_begin <= intersect_phi <= column_phi_range_end) \
        or (column_phi_range_begin <= intersect_phi + 2 * np.pi <= column_phi_range_end) \
        or (column_phi_range_begin <= intersect_phi - 2 * np.pi <= column_phi_range_end) \
        or (column_phi_range_begin <= intersect_phi + 4 * np.pi <= column_phi_range_end) \
        or (column_phi_range_begin <= intersect_phi - 4 * np.pi <= column_phi_range_end)


def check_shadow_with_dipole(surface, phi_index, theta_index, obs_vector, solutions, top_column, bot_column):
    # checks for 1 element area!
    # obs_vector = direction_vector
    '''
        принимает на вход массив решений для затмений с колонкой
        находим пересечение с колонками
           если пересекается то истина
           если нет то ложь

        корни в магнитной СК - betta_mu
    '''

    # origin - откуда бьет луч. площадка на поверхности
    origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]
    # работаем в нормировке на радиус звезды
    r = surface.surf_R_e / config.R_ns * np.sin(origin_theta) ** 2

    # ищем положительные корни и при этом вещественные
    for solution in solutions:
        # ищем корни которые нам подходят, значит произошло затмение
        intersect_phi, intersect_theta = get_intersection_from_solution(r, origin_phi, origin_theta, obs_vector,
                                                                        solution)
        if intersect_phi is not None:
            # условие по phi для верхней колонки:
            top_column_intersect_phi = get_intersection_phi_with_column(top_column, intersect_phi)

            # условие по phi для нижней колонки:
            bot_column_intersect_phi = get_intersection_phi_with_column(bot_column, intersect_phi)

            # если есть пересечение то вернем 0 так как затмение
            if (intersect_theta < top_column.theta_range[-1] and top_column_intersect_phi) \
                    or (intersect_theta > bot_column.theta_range[-1] and bot_column_intersect_phi):
                return 0
    # не нашли затмений - вернем 1
    return 1


def get_tau_for_opacity_old(theta, R_e, R_e_for_delta, M_accretion_rate, a_portion, dRe_div_Re):
    # R_e_delta
    tau = config.k * M_accretion_rate * newService.get_delta_distance(theta, R_e_for_delta, dRe_div_Re) / (
            newService.get_A_normal(theta, R_e_for_delta, a_portion, dRe_div_Re)
            * newService.get_free_fall_velocity(theta, R_e))
    return tau


def get_tau_for_opacity(phi, theta, R_e_of_atenuation_surf, R_e_for_delta, M_accretion_rate, a_portion, obs_vector,
                        dRe_div_Re):
    '''
    считаем коэффициент ослабления

    + теперь учитываем угол alpha - угол между нормалью в точке пересечения и направлением на наблюдателя.
    так как нормаль != пути вдоль вещества. примерно это прямоугольный треугольник. надо делить на косинус.

    мы уже поняли что будет пересечение. отрицательный угол может получиться поэтому вернем модуль

    сейчас почти всегда >>1 и ослабление сразу в 0. почти как затмение

    очень маленькие cos = 0.005?? - хз вроде уже норм

    + save tau and alpha here in tensors

    '''
    # phi, theta - в точке пересечения; R_e =
    normal = matrix.newE_n(phi, theta)
    cos_alpha = np.dot(normal, obs_vector) / (np.linalg.norm(normal) * np.linalg.norm(obs_vector))
    # формула в статье есть
    tau = config.k * M_accretion_rate * newService.get_delta_distance(theta, R_e_for_delta, dRe_div_Re)
    tau /= (newService.get_A_normal(theta, R_e_for_delta, a_portion, dRe_div_Re)
            * newService.get_free_fall_velocity(theta, R_e_of_atenuation_surf))
    # print(f'{tau=}')
    # tau /= cos_alpha
    # print(f'{cos_alpha=}')
    # мы уже поняли что будет пересечение. отрицательный угол может получиться поэтому берем модуль
    return np.abs(tau / cos_alpha), tau, cos_alpha


def get_tau_for_scatter_with_cos(theta_range, R_e_emission_surf, M_accretion_rate, a_portion, cos_alpha, dRe_div_Re):
    tau = config.k * M_accretion_rate * newService.get_delta_distance(theta_range, R_e_emission_surf, dRe_div_Re)
    tau /= (newService.get_A_normal(theta_range, R_e_emission_surf, a_portion, dRe_div_Re)
            * newService.get_free_fall_velocity(theta_range, R_e_emission_surf))
    tau = tau[np.newaxis, :] / cos_alpha
    # print(f'cos_alpha_min = {np.min(cos_alpha)}, cos_alpha_max = {np.max(cos_alpha)}')
    # print(f'tau_min = {np.min(tau)}, max = {np.max(tau)}')
    return tau  # np.abs?


def get_tau_with_dipole(surface, phi_index, theta_index, obs_vector, solutions,
                        curr_configuration: accretingNS.AccretingPulsarConfiguration):
    # срочно ввести учет маски магнитных линий !!!!
    '''
        считаем коэффициент ослабления tau
        R_e_of_atenuation_surf = inner surf!!!!!
        я беру пересечения только с внутренней (пока что)
    '''
    R_e_for_delta = curr_configuration.top_column_for_calc.R_e_for_delta
    R_e_of_atenuation_surf = curr_configuration.top_column_for_calc.R_e
    beta_mu = curr_configuration.beta_mu
    M_accretion_rate = curr_configuration.M_accretion_rate
    a_portion = curr_configuration.a_portion
    dRe_div_Re = curr_configuration.dRe_div_Re
    top_column = curr_configuration.top_column_for_calc.inner_surface
    bot_column = curr_configuration.bot_column_for_calc.inner_surface

    origin_phi, origin_theta = surface.phi_range[phi_index], surface.theta_range[theta_index]
    # здесь радиус основания, радиус излучающей так как r = для положения излучающей
    r = surface.surf_R_e / config.R_ns * np.sin(origin_theta) ** 2

    for solution in solutions:
        # смотрим по тем же решениям что и для колонки. так как они для дипольной линии
        intersect_phi, intersect_theta = get_intersection_from_solution(r, origin_phi, origin_theta, obs_vector,
                                                                        solution)
        if intersect_phi is not None:
            # из геометрических выводов = конец магнитосферной линии. был вывод в тетрадке этой формулы
            theta_end = np.pi / 2 - np.arctan(np.tan(np.deg2rad(beta_mu)) * np.cos(intersect_phi))
            # учет маски!!
            spherical_R_condition = curr_configuration.top_magnet_lines_for_calc.surf_R_e * np.sin(
                intersect_theta) ** 2 < (1 + dRe_div_Re) * curr_configuration.R_disk
            spherical_R_condition = True  # мы же убрали сферическое обрезание!

            # условия пересечений тета
            top_column_intersect_theta_correct = intersect_theta < theta_end
            bot_column_intersect_theta_correct = intersect_theta > theta_end  # wtf ????????????

            # условия пересечений фи - такие же как и у колонки. можно их вызвать
            top_column_intersect_phi = get_intersection_phi_with_column(top_column, intersect_phi)
            bot_column_intersect_phi = get_intersection_phi_with_column(bot_column, intersect_phi)

            # условия пересечений с верхней или нижней линией магнитной
            intersection_condition = ((top_column_intersect_phi and top_column_intersect_theta_correct) or (
                    bot_column_intersect_phi and bot_column_intersect_theta_correct)) and spherical_R_condition

            # if intersection_condition:
            #     # учет маски!!
            #
            #     if top_column_intersect_phi:
            #         phi_range = curr_configuration.top_magnet_lines.phi_range
            #     else:  # bot_column_intersect_phi
            #         phi_range = curr_configuration.bot_magnet_lines.phi_range
            #
            #     phi_range_plus = phi_range + 2 * np.pi
            #     phi_range_minus = phi_range - 2 * np.pi
            #     phi_range = np.hstack([phi_range, phi_range_plus, phi_range_minus])
            #
            #     phi_range = np.abs(phi_range - intersect_phi)
            #
            #     closest_id = np.argmin(phi_range) % curr_configuration.top_magnet_lines.phi_range.shape[0]
            #     # mask_array - общий ???? да
            #     # ~ потому что тут mask = True -- значит нет линий
            #     # берем top_mag_lines так как у bot np.pi - theta; и все ок
            #
            #     if top_column_intersect_phi:
            #         theta_arr = curr_configuration.top_magnet_lines.theta_range[
            #             ~ curr_configuration.top_magnet_lines.mask_array[closest_id]]
            #     else:
            #         theta_arr = curr_configuration.bot_magnet_lines.theta_range[
            #             ~ curr_configuration.top_magnet_lines.mask_array[closest_id]]
            #
            #     if theta_arr.shape[0] > 0:
            #         max_theta = theta_arr[-1]
            #         # if theta_arr.shape[0] != config.N_theta_accretion:
            #         #     next_theta = curr_configuration.top_magnet_lines.theta_range[theta_arr.shape[0]]
            #         #     if next_theta > theta_end:
            #         #         max_theta = theta_end
            #
            #         top_column_intersect_theta_correct = intersect_theta < max_theta
            #         bot_column_intersect_theta_correct = intersect_theta > max_theta  # np.pi - max_theta ????
            #
            #         intersection_condition = (top_column_intersect_phi and top_column_intersect_theta_correct) or (
            #                 bot_column_intersect_phi and bot_column_intersect_theta_correct)
            #     else:
            #         intersection_condition = False

            # если есть пересечение считаем коэффициент ослабления
            if intersection_condition:
                tau_cos, tau, cos_alpha = get_tau_for_opacity(intersect_phi, intersect_theta, R_e_of_atenuation_surf,
                                                              R_e_for_delta, M_accretion_rate,
                                                              a_portion, obs_vector, dRe_div_Re)
                # tau = get_tau_for_opacity(intersect_phi, intersect_theta, R_e_of_atenuation_surf, M_accretion_rate,
                #                           a_portion, obs_vector, dRe_div_Re)
                if tau_cos > config.tau_cutoff:
                    return np.exp(-1 * tau_cos), tau, np.abs(cos_alpha)
                else:
                    # pass
                    return 1, 0, -1
    # если не нашли пересечений, то ослабления нет
    return 1, 0, -1


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
