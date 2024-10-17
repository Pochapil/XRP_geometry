import numpy as np
import scipy.integrate
from scipy import interpolate

import newService
import config
import BS_approach as get_T_eff
import geometricTask.matrix as matrix

# surfs - polar = outer; equatorial = inner
column_surf_types = {'bot': 'bot', 'top': 'top'}
surface_surf_types = {'outer': 'outer', 'inner': 'inner'}


class AccretingPulsarConfiguration:

    def __init__(self, mu, beta_mu, mc2, a_portion, phi_0):
        self.mu = mu
        self.beta_mu = beta_mu
        self.mc2 = mc2
        self.a_portion = a_portion
        self.phi_0 = phi_0

        self.M_accretion_rate = mc2 * config.L_edd / config.c ** 2
        self.R_alfven = (mu ** 2 / (2 * self.M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
        self.R_e = config.ksi_param * self.R_alfven

        self.top_column = AccretionColumn(self.R_e, mu, mc2, self.a_portion, self.phi_0,
                                          column_type=column_surf_types['top'])
        # Accret_col : equatorial, polar surfaces
        self.bot_column = AccretionColumn(self.R_e, mu, mc2, self.a_portion, self.phi_0,
                                          column_type=column_surf_types['bot'])

        self.top_magnet_lines = MagnetLine(self.beta_mu, self.top_column.inner_surface.phi_range,
                                           self.top_column.inner_surface.theta_range[-1], self.top_column.column_type)
        self.bot_magnet_lines = MagnetLine(self.beta_mu, self.top_column.inner_surface.phi_range,
                                           self.top_column.inner_surface.theta_range[-1], self.bot_column.column_type)


class AccretionColumn:

    def __init__(self, R_e, mu, mc2, a_portion, phi_0, column_type):
        self.R_e = R_e
        self.column_type = column_type
        self.a_portion = a_portion
        self.phi_0 = phi_0

        # попробую что
        # R_e_outer_surface, R_e_inner_surface = self.R_e, self.R_e
        self.R_e_outer_surface, self.R_e_inner_surface = (1 + config.dRe_div_Re) * R_e, R_e
        if config.FLAG_R_E_OLD:
            self.R_e_outer_surface, self.R_e_inner_surface = R_e, R_e  # допущение что толщина = 0

        M_accretion_rate = mc2 * config.L_edd / config.c ** 2
        self.T_eff, self.ksi_shock, self.L_x, self.beta = get_T_eff.get_Teff_distribution(
            newService.get_delta_distance_at_surface_NS(self.R_e),
            newService.get_A_normal_at_surface_NS(self.R_e, self.a_portion),
            mu, M_accretion_rate
        )

        self.outer_surface = Surface(self.R_e_outer_surface, self.a_portion, self.phi_0, self.ksi_shock,
                                     self.column_type, surface_type=surface_surf_types['outer'])
        self.inner_surface = Surface(self.R_e_inner_surface, self.a_portion, self.phi_0, self.ksi_shock,
                                     self.column_type, surface_type=surface_surf_types['inner'])

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
    def __init__(self, surf_R_e, a_portion, phi_0, ksi_shock, column_type, surface_type):
        self.surf_R_e = surf_R_e  # R_e на котором сидит - R_e либо R_e + \delta R_e
        self.surface_type = surface_type

        self.theta_accretion_begin = newService.get_theta_accretion_begin(self.surf_R_e)
        if (config.R_ns * ksi_shock / self.surf_R_e) >= 1:
            # есть набор параметров при которых модель не работает и ударная волна дальше магнитосферы, берем 90
            '''вопрос - мб нужна формула с arctan...'''
            self.theta_accretion_end = np.pi / 2
        else:
            # из усл силовой линии МП : r = R_e sin**2; end: ksi_shock = R_e sin**2
            self.theta_accretion_end = np.arcsin((config.R_ns * ksi_shock / self.surf_R_e) ** (1 / 2))

        # для нижней сместить углы
        if column_type == column_surf_types['bot']:
            self.theta_accretion_begin = np.pi - self.theta_accretion_begin
            self.theta_accretion_end = np.pi - self.theta_accretion_end
        self.theta_range = np.linspace(self.theta_accretion_begin, self.theta_accretion_end, config.N_theta_accretion)

        phi_0_rad = np.deg2rad(phi_0)
        phi_delta = 0
        if column_type == column_surf_types['bot']:
            # для нижней сместить углы
            phi_delta = np.pi
        # phi_0 - центр колонки!!!!!!!!!!!!
        self.phi_range = np.linspace(-np.pi * a_portion, np.pi * a_portion, config.N_phi_accretion)
        self.phi_range = self.phi_range + phi_0_rad + phi_delta
        if config.FLAG_PHI_0_OLD:
            # - в терминах старого phi
            fi_0_old = np.deg2rad(phi_0 + config.fi_0_dict[a_portion])
            self.phi_range = np.linspace(0, 2 * np.pi * a_portion, config.N_phi_accretion) + fi_0_old + phi_delta

        self.array_normal = self.create_array_normal(self.phi_range, self.theta_range, self.surface_type)

    def create_array_normal(self, phi_range, theta_range, surface_type):
        '''работало раньше, когда R_e были одинаковые!!!!!!!!!
        мол outer surf (pol) - стреляет наружу
            iner (eq) - внутрь
        по факту нужно переделать:
        для всех углов считать
        array_normal = matrix.newE_n_n(phi_range, theta_range)
        а потом разбираться с пересечениями и ослаблениями
        '''
        # array_normal - матрица нормалей
        # тензор размером phi x theta x 3 (x,y,z)
        coefficient = -1
        if surface_type == surface_surf_types['outer']:  # True - внешняя поверхность, False - внутренняя
            coefficient = 1
        array_normal = coefficient * matrix.newE_n_n(phi_range, theta_range)
        return array_normal


class MagnetLine:
    def __init__(self, beta_mu, top_column_phi_range, top_column_theta_end, column_type):
        self.column_type = column_type

        theta_range_end = np.pi / 2 + beta_mu
        # ограничиваю колонкой
        theta_range_end = min((np.pi - top_column_theta_end), theta_range_end)
        theta_range_begin = top_column_theta_end

        self.theta_range = np.linspace(theta_range_begin, theta_range_end, config.N_theta_accretion)
        self.phi_range = top_column_phi_range

        self.mask_array = np.zeros((config.N_phi_accretion, config.N_theta_accretion)).astype(bool)
        # mask = np.zeros_like(x).astype(bool)

        for i, phi in enumerate(self.phi_range):
            for j, theta in enumerate(self.theta_range):
                theta_end = np.pi / 2 - np.arctan(np.tan(np.deg2rad(beta_mu)) * np.cos(phi))
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
