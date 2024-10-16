import numpy as np

import accretingNS
import config
from geometricTask import matrix
import pathService

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

for surf in accr_col_surfs:
    # тензор косинусов между нормалью и направлением на наблюдателя размером phase x phi x theta
    # умножаем скалярно phi x theta x 3 на phase x 3 (по последнему индексу) и делаем reshape.
    rotation_psi_matrix = np.einsum('ijl,tl->tij', surf.array_normal, obs_matrix)
    # rotation_obs_matrix[k] > 0:
    # if cos_psi_range[i, j] > 0:
    #     # проверка на пересечения
    #     r = self.R_e / config.R_ns * np.sin(self.theta_range[j]) ** 2
    #
    #     origin_x = np.sin(self.theta_range[j]) * np.cos(self.phi_range[i]) * r
    #     origin_y = np.sin(self.theta_range[j]) * np.sin(self.phi_range[i]) * r
    #     origin_z = np.cos(self.theta_range[j]) * r
    #
    #     direction_x = e_obs_mu[0, 0]
    #     direction_y = e_obs_mu[0, 1]
    #     direction_z = e_obs_mu[0, 2]
    #

    # shadows
    # calc shadowed_matrix (ns + columns)
    # tau
    # calc tau_matrix

    # для разных матриц можем посчитать L и посмотреть какой вклад будет.
    # calc_L()

    # save
    # calc PF
    # calc_L_nu()
    # save

magnet_line_surfs = []
for surf in magnet_line_surfs:
    # shadows
    # calc shadowed_matrix (ns + columns)
    # tau
    # calc tau_matrix

    # УЧЕСТЬ tau в отражаемой точке!!!!!!!!!

    '''
    np.array([np.sin(theta_obs) * np.cos(phi),
              np.sin(theta_obs) * np.sin(phi),
              np.cos(theta_obs)]) - типо куда луч стреляем. нужно для нахождения tau1
    '''

    # save
    # calc nu L nu

    # calc_scatter_L_nu()

    # save
    # calc nu L nu

    # tij -> ti -> t
    # scipy.integrate(axis=-1) --- можно интегрировать по тензору.
