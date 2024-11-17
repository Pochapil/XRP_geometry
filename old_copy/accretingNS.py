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
        self.R_e_outer_surface, self.R_e_inner_surface = (1 + config.dRe_div_Re) * self.R_e, self.R_e
        if config.FLAG_R_E_OLD:
            # если допущение что толщина = 0
            self.R_e_outer_surface, self.R_e_inner_surface = self.R_e, self.R_e

        M_accretion_rate = mc2 * config.L_edd / config.c ** 2
        # вызов БС чтобы расчитать все распределения. МБ стоит поднять выше так как расчет только 1 раз нужен
        self.T_eff, self.ksi_shock, self.L_x, self.beta = get_T_eff.get_Teff_distribution(
            newService.get_delta_distance_at_surface_NS(self.R_e),
            newService.get_A_normal_at_surface_NS(self.R_e, self.a_portion),
            mu, M_accretion_rate
        )

        phi_range = np.linspace(-np.pi * a_portion, np.pi * a_portion, config.N_phi_accretion) + np.deg2rad(phi_0)
        min_angle = np.min(np.pi / 2 - np.arctan(np.tan(np.deg2rad(beta_mu)) * np.cos(phi_range)))
        # скорее всего здесь необходимо пересчитать чтобы было ок ?
        if self.ksi_shock > self.R_e / config.R_ns * np.sin(min_angle) ** 2:
            self.R_e = self.ksi_shock / (np.sin(min_angle) ** 2) * config.R_ns
            self.R_e_outer_surface, self.R_e_inner_surface = (1 + config.dRe_div_Re) * self.R_e, self.R_e

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

        # print(y_new.shape)
        self.T_eff[1:-1] = y_new
        # for i in range(0, (len(y_new))):
        #     self.T_eff[i + 1] = y_new[i]