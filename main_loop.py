import multiprocessing as mp
from itertools import repeat, cycle, product
import time
import numpy as np

import main
import config


def make_arr(theta_obs_arr, beta_mu_arr, mc2_arr, a_portion_arr, phi_0_arr):
    '''оказалось ошибкой'''
    N_big = len(theta_obs_arr) * len(beta_mu_arr) * len(mc2_arr) * len(a_portion_arr) * len(phi_0_arr)

    def func(arr, N_big):
        cur_N = N_big // len(arr)
        arr_async = np.empty(N_big)
        for i, ar in enumerate(arr):
            arr_async[i * cur_N:(i + 1) * cur_N] = ar
        return arr_async

    theta_obs_arr_async = func(theta_obs_arr, N_big)
    beta_mu_arr_async = func(beta_mu_arr, N_big)
    mc2_arr_async = func(mc2_arr, N_big)
    a_portion_arr_async = func(a_portion_arr, N_big)
    phi_0_arr_async = func(phi_0_arr, N_big)

    return theta_obs_arr_async, beta_mu_arr_async, mc2_arr_async, a_portion_arr_async, phi_0_arr_async

    # theta_obs_arr_async, beta_mu_arr_async, mc2_arr_async, a_portion_arr_async, phi_0_arr_async = make_arr(
    #     theta_obs_arr, beta_mu_arr, mc2_arr, a_portion_arr, phi_0_arr)


if __name__ == '__main__':
    mu = 0.1e31

    beta_mu_arr = [10 * i for i in range(1, 10)]
    theta_obs_arr = [10 * i for i in range(1, 10)]
    mc2_arr = [30, 60, 100]
    a_portion_arr = [0.22, 0.44, 0.66, 1]
    phi_0_arr = [0]

    theta_obs_arr = [40]
    beta_mu_arr = [60]
    # [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    mc2_arr = [10, 20, 40, 50, 60, 70, 80, 90, 110, 120, 130]
    a_portion_arr = [0.44, 0.66]
    phi_0_arr = [0]

    beta_mu_arr = [10 * i for i in range(1, 10)]
    theta_obs_arr = [10 * i for i in range(1, 10)]
    mc2_arr = [30]
    a_portion_arr = [1]
    phi_0_arr = [0]

    theta_obs_arr = [10 * i for i in range(6, 10)]
    beta_mu_arr = [10 * i for i in range(1, 10)]
    mc2_arr = [30, 60, 100]
    a_portion_arr = [0.22, 0.44, 0.66]
    phi_0_arr = [0, 40, 90]

    theta_obs_arr = [20]
    beta_mu_arr = [60]
    # [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    mc2_arr = [10, 20, 40, 50, 60, 70, 80, 90, 110, 120, 130]
    a_portion_arr = [0.44, 0.66]
    phi_0_arr = [0]

    theta_obs_arr = [20]
    beta_mu_arr = [60]
    # [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    mc2_arr = [160, 180]
    a_portion_arr = [0.44, 0.66]
    phi_0_arr = [0]

    plot_flag = True
    # ------------------------------------------------- start -------------------------------------------------------
    N_big = len(theta_obs_arr) * len(beta_mu_arr) * len(mc2_arr) * len(a_portion_arr) * len(phi_0_arr)
    if config.flag_calc_clever:

        print(f'to calculate {N_big} loops need about {100 * N_big / 3600 / config.N_cpus} hours')
        t1 = time.perf_counter()

        # print(*product([mu], theta_obs_arr, beta_mu_arr, mc2_arr, a_portion_arr, phi_0_arr, [plot_flag]))
        # for i in product([mu], theta_obs_arr, beta_mu_arr, mc2_arr, a_portion_arr, phi_0_arr, [plot_flag]):
        #     print(i)
        '''чтобы пройти по сетке берем product от массивов значений аргументов - получим декартово произведение'''
        with mp.Pool(processes=config.N_cpus) as pool:
            pool.starmap(main.calc_and_save_for_configuration,
                         (product([mu], theta_obs_arr, beta_mu_arr, mc2_arr, a_portion_arr, phi_0_arr, [plot_flag],
                                  [False])))  # False - чтобы считал на 1 процессоре каждый вызов

        #     pool.starmap(main.calc_and_save_for_configuration,
        #                  zip(repeat(mu, N_big), cycle(theta_obs_arr), cycle(beta_mu_arr), cycle(mc2_arr),
        #                      cycle(a_portion_arr), cycle(phi_0_arr), repeat(plot_flag)))

        t2 = time.perf_counter()
        print(f'{t2 - t1} seconds')

    else:
        print(f'to calculate {N_big} loops need about {20 * N_big / 3600} hours')
        for theta_obs in theta_obs_arr:
            for beta_mu in beta_mu_arr:
                for mc2 in mc2_arr:
                    for a_portion in a_portion_arr:
                        for phi_0 in phi_0_arr:
                            main.calc_and_save_for_configuration(mu, theta_obs, beta_mu, mc2, a_portion, phi_0,
                                                                 plot_flag)
