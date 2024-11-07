import numpy as np
import scipy.integrate
from scipy import interpolate

import newService
import config
import BS_approach as get_T_eff
import geometry.matrix as matrix

# surfs - polar = outer; equatorial = inner
column_surf_types = {'bot': 'bot', 'top': 'top'}
surface_surf_types = {'outer': 'outer', 'inner': 'inner'}


class AccretingPulsarConfiguration:
    '''решил сделать класс который хранит текущую конфигурацию.
    некоторые поля также являются классами

    top_column - верхняя колонка
    bot_column - нижняя колонка

    top_magnet_lines - верхние магнитные линии
    bot_magnet_lines - нижние
    '''

    def __init__(self, mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
        self.mu = mu
        self.theta_obs = theta_obs
        self.beta_mu = beta_mu
        self.mc2 = mc2
        self.a_portion = a_portion
        self.phi_0 = phi_0

        self.M_accretion_rate = mc2 * config.L_edd / config.c ** 2
        self.R_alfven = (mu ** 2 / (2 * self.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
        self.R_e = config.ksi_param * self.R_alfven

        self.top_column = AccretionColumn(self.R_e, mu, self.beta_mu, mc2, self.a_portion, self.phi_0,
                                          column_type=column_surf_types['top'])
        # Accret_col : equatorial, polar surfaces
        self.bot_column = AccretionColumn(self.R_e, mu, self.beta_mu, mc2, self.a_portion, self.phi_0,
                                          column_type=column_surf_types['bot'])

        self.top_magnet_lines = MagnetLine(self.top_column.R_e_inner_surface, self.beta_mu,
                                           self.top_column.inner_surface.phi_range,
                                           self.top_column.inner_surface.theta_range[-1], self.top_column.column_type)

        self.bot_magnet_lines = MagnetLine(self.top_column.R_e_inner_surface, self.beta_mu,
                                           self.top_column.inner_surface.phi_range,
                                           self.top_column.inner_surface.theta_range[-1], self.bot_column.column_type)


class AccretionColumn:
    '''класс колонки. хранит  себе значния конфигурации + знает свои поверхности

    outer_surface - внешняя поверхность/полярная ?
    inner_surface - внутрення поверхность/экваториальная
    '''

    def __init__(self, R_e, mu, beta_mu, mc2, a_portion, phi_0, column_type):
        self.R_e = R_e
        self.column_type = column_type
        self.a_portion = a_portion
        self.phi_0 = phi_0

        # попробую что
        self.R_e_outer_surface, self.R_e_inner_surface = (1 + config.dRe_div_Re) * R_e, R_e
        if config.FLAG_R_E_OLD:
            # если допущение что толщина = 0
            self.R_e_outer_surface, self.R_e_inner_surface = R_e, R_e

        M_accretion_rate = mc2 * config.L_edd / config.c ** 2
        # вызов БС чтобы расчитать все распределения. МБ стоит поднять выше так как расчет только 1 раз нужен
        self.T_eff, self.ksi_shock, self.L_x, self.beta = get_T_eff.get_Teff_distribution(
            newService.get_delta_distance_at_surface_NS(self.R_e),
            newService.get_A_normal_at_surface_NS(self.R_e, self.a_portion),
            mu, M_accretion_rate
        )

        self.outer_surface = Surface(self.R_e_outer_surface, self.R_e, beta_mu, self.a_portion, self.phi_0,
                                     self.ksi_shock, self.column_type, surface_type=surface_surf_types['outer'])
        self.inner_surface = Surface(self.R_e_inner_surface, self.R_e, beta_mu, self.a_portion, self.phi_0,
                                     self.ksi_shock, self.column_type, surface_type=surface_surf_types['inner'])
        # вызов подправки распределения Т
        self.correct_T_eff()

    def correct_T_eff(self):
        # так как нахожу распределение Teff по отрезку кси, то нужно привести к виду по сетке theta !!
        # нужно проинтерполировать внутри отрезка

        # создание сетки по которой было найдено распределение T_eff
        ksi_range = np.linspace(1, self.ksi_shock, config.N_theta_accretion)
        if self.ksi_shock > self.R_e / config.R_ns:
            ksi_range = np.linspace(1, self.R_e / config.R_ns, config.N_theta_accretion)

        # интерполяция
        x = ksi_range
        y = self.T_eff
        f = interpolate.interp1d(x, y, kind='cubic')

        # новый диапазон по х - начинаем с 1 индекса так как 0 элемент совпадает - x[0] = 1
        x_new = self.R_e / config.R_ns * np.sin(self.inner_surface.theta_range[1:-1]) ** 2
        y_new = f(x_new)

        # выравнивание
        for i in range(0, (len(y_new))):
            self.T_eff[i + 1] = y_new[i]


class Surface:
    '''класс поверхности. можно сказать основной - будет участвовать во всех расчетах.
    знает распределение сетки по theta phi

    наверно стоит добавить промежуточный слой между поверхностью и магнитными линиями чтобы учесть обрезание колонки

    тоже теперь через маску - для обрезания колонок. будет верно так как мы не будем учитывать магнитные линии - они
    будут больше 90 градусов и обрезаны. а если нет обрезания колонки тогда есть шанс для магнитных линий

    phi_0 - центр колонки!!!!!!!!!!!!
    '''

    def __init__(self, surf_R_e, col_R_e, beta_mu, a_portion, phi_0, ksi_shock, column_type, surface_type):
        self.surf_R_e = surf_R_e  # R_e на котором сидит - R_e либо R_e + \delta R_e
        self.surface_type = surface_type

        self.theta_accretion_begin = newService.get_theta_accretion_begin(self.surf_R_e)

        if config.outer_R_e_ksi_flag:
            # обрезать большую колонку по кси. (старый способ)
            if (config.R_ns * ksi_shock / self.surf_R_e) >= 1:
                self.theta_accretion_end = np.pi / 2
            else:
                self.theta_accretion_end = np.arcsin((config.R_ns * ksi_shock / self.surf_R_e) ** (1 / 2))
        else:
            # новый способ
            # беру col_R_e - чтобы обрезать большую колонку по тета а не кси.
            if (config.R_ns * ksi_shock / col_R_e) >= 1:
                # есть набор параметров при которых модель не работает и ударная волна дальше магнитосферы, берем 90
                '''вопрос - мб нужна формула с arctan... скорее всего нужно сначала phi потом расчитывать тета
                с np.pi/2 есть проблема --- 
                '''
                # когда наклонен от pi/2 - beta до pi/2 + beta в зависимости от phi. хотим регулярную а отрежем маской
                # поэтому берем +
                self.theta_accretion_end = np.pi / 2 + np.deg2rad(beta_mu)
            else:
                # из усл силовой линии МП : r = R_e sin**2; end: ksi_shock = R_e sin**2
                self.theta_accretion_end = np.arcsin((config.R_ns * ksi_shock / col_R_e) ** (1 / 2))

        self.theta_range = np.linspace(self.theta_accretion_begin, self.theta_accretion_end, config.N_theta_accretion)

        phi_0_rad = np.deg2rad(phi_0)
        phi_delta = 0
        # phi_0 - центр колонки!!!!!!!!!!!!
        self.phi_range = np.linspace(-np.pi * a_portion, np.pi * a_portion, config.N_phi_accretion) + phi_0_rad

        self.mask_array = np.zeros((config.N_phi_accretion, config.N_theta_accretion)).astype(bool)
        # mask = np.zeros_like(x).astype(bool)
        for i, phi in enumerate(self.phi_range):
            for j, theta in enumerate(self.theta_range):
                theta_end = np.pi / 2 - np.arctan(np.tan(np.deg2rad(beta_mu)) * np.cos(phi))
                # ограничиваю диском True - значит не надо использовать эту площадку
                if theta > theta_end:
                    self.mask_array[i][j] = True

        # для нижней сместить углы. маска будет той же
        if column_type == column_surf_types['bot']:
            phi_delta = np.pi
            self.phi_range = self.phi_range + phi_delta
            self.theta_range = (np.pi - self.theta_range)

        if config.FLAG_PHI_0_OLD:
            # - в терминах старого phi
            fi_0_old = np.deg2rad(phi_0 + config.fi_0_dict[a_portion])
            self.phi_range = np.linspace(0, 2 * np.pi * a_portion, config.N_phi_accretion) + fi_0_old + phi_delta

        self.array_normal = self.create_array_normal(self.phi_range, self.theta_range, self.surface_type)

    def create_array_normal(self, phi_range, theta_range, surface_type):
        '''
        array_normal - матрица нормалей ==== тензор размером phi x theta x 3 (x,y,z)

        работало раньше, когда R_e были одинаковые!!!!!!!!!
        мол outer surf (pol) - стреляет наружу; iner (eq) - внутрь
        по факту нужно переделать:
        для всех углов считать array_normal = matrix.newE_n_n(phi_range, theta_range)
        а потом разбираться с пересечениями и ослаблениями

        хотя вроде ок. так как внешняя наружу, а внутрення внутрь
        (погрешности если и есть то лишь на небольшом наборе углов и причем небольшие)
        '''
        # coefficient = -1 значит вектор смотрит вовнутрь (посмотри на векторное умножение)! - для внутренней поверхности
        coefficient = -1
        if surface_type == surface_surf_types['outer']:  # True - внешняя поверхность, False - внутренняя
            coefficient = 1
        array_normal = coefficient * matrix.newE_n_n(phi_range, theta_range)
        return array_normal


class MagnetLine:
    '''класс магнитных поверхностей (вещество над ударной волной) для рассеяния.
        знает распределение сетки по theta phi. начинается над колонкой.'''

    def __init__(self, R_e, beta_mu, top_column_phi_range, top_column_theta_end, column_type):

        '''
        top_column_theta_end = для ВНУТРЕННЕЙ поверхности пока что!!!!

        смог реализовать только с помощью маски

        маска - чтобы не портить сетку по teta и phi и оставить ее регулярной
        мы просто смотрим что ячейки сетки лежат над диском (для верхних магнитных линий)
        для нижней колонки mask_array будет такой же

        в маске значение Ture значит не надо использовать эту площадку в расчетах - она ниже диска (для верхней колонки)
        '''
        self.column_type = column_type
        self.surf_R_e = R_e
        # так как диск связан с осью вращения максимальный угол = np.pi / 2 + beta_mu - угол на котором диск от оси омега
        theta_range_end = np.pi / 2 + np.deg2rad(beta_mu)
        # ограничиваю колонкой
        # theta_range_end = min((np.pi - top_column_theta_end), theta_range_end)
        theta_range_begin = top_column_theta_end

        self.theta_range = np.linspace(theta_range_begin, theta_range_end, config.N_theta_accretion)
        self.phi_range = top_column_phi_range

        self.mask_array = np.zeros((config.N_phi_accretion, config.N_theta_accretion)).astype(bool)
        # mask = np.zeros_like(x).astype(bool)
        for i, phi in enumerate(self.phi_range):
            for j, theta in enumerate(self.theta_range):
                theta_end = np.pi / 2 - np.arctan(np.tan(np.deg2rad(beta_mu)) * np.cos(phi))
                # ограничиваю диском True - значит не надо использовать эту площадку
                if theta > theta_end:
                    self.mask_array[i][j] = True

        if column_type == column_surf_types['bot']:
            # для нижней колонки меняем phi, theta. mask_array будет такой же
            self.theta_range = (np.pi - self.theta_range)  # [::-1]?
            self.phi_range = self.phi_range + np.pi

        self.array_normal = self.create_array_normal(self.phi_range, self.theta_range)

    def create_array_normal(self, phi_range, theta_range):
        # беру только внутренние нормали, так как для расчета отраженного необходимо как раз нормали вовнутрь
        array_normal = -1 * matrix.newE_n_n(phi_range, theta_range)
        return array_normal
