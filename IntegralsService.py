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

# переделать phi_0 !!!!

curr_configuration = accretingNS.AccretingPulsarConfiguration(mu, beta_mu, mc2, a_portion, phi_0)

theta_obs = 60

cur_path = pathService.PathSaver(mu, theta_obs, beta_mu, mc2, a_portion, phi_0).get_path()
print(cur_path)
e_obs = np.array([np.sin(theta_obs) * np.cos(0),
                  np.sin(theta_obs) * np.sin(0),
                  np.cos(theta_obs)])

# obs_matrix = np.zeros((config.N_theta_accretion, config.N_phi_accretion, 3))
obs_matrix = np.zeros((config.N_phase, 3))

# rotate и сохранить матрицу векторов наблюдателя
for phase_index in range(config.N_phase):
    phi_mu = config.phi_mu_0 + config.omega_ns_rad * phase_index

    A_matrix_analytic = matrix.newMatrixAnalytic(config.phi_rotate, config.betta_rotate, phi_mu, beta_mu)
    e_obs_mu = np.dot(A_matrix_analytic, e_obs)

    obs_matrix[phase_index] = e_obs_mu

# возможно хранить будем все матрицы.
pulsar_surfs = [curr_configuration.top_column.inner_surface]
for surf in pulsar_surfs:
    rotation_obs_matrix = np.einsum('ijl,tl->tij', surf.array_normal, obs_matrix)


    # shadows
    # calc shadowed_matrix (ns + columns)
    # tau
    # calc tau_matrix

    # для разных матриц можем посчитать L и посмотреть какой вклад будет.

    def calc_L():
        ...


    # save
    # calc PF
    def calc_L_nu():
        ...
    # save

magnet_line_surfs = []
for surf in magnet_line_surfs:
    # shadows
    # calc shadowed_matrix (ns + columns)
    # tau
    # calc tau_matrix

    # УЧЕСТЬ tau в отражаемой точке!!!!!!!!!

    '''
    np.array([np.sin(theta_ob) * np.cos(phi),
              np.sin(theta_ob) * np.sin(phi),
              np.cos(theta_ob)]) - типо куда луч стреляем. нужно для нахождения tau1
    '''


    def calc_scatter_L():
        ...


    # save
    # calc nu L nu

    def calc_scatter_L_nu():
        ...
    # save
    # calc nu L nu

    # tij -> ti -> t
    # scipy.integrate(axis=-1) --- можно интегрировать по тензору.
