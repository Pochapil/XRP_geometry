import multiprocessing as mp
from itertools import repeat, cycle, product
import time
import numpy as np

import plot_package.romanova_phi_0 as romanova_phi_0

import main
import config

if __name__ == '__main__':
    #     romanova_tuples = [(10, 0),(20, 28), (30, 43), (40, 52), (50, 60), (60, 64), (70, 68), (80, 71), (90, 74), (100, 77)]
    #     for tuple in romanova_tuples:
    #         mc2 = tuple[0]
    #         phi_0 = tuple[1]


    mu = 0.1e31

    theta_obs_arr = [40]
    beta_mu_arr = [20]
    mc2_arr = [50]
    a_portion_arr = [0.44]
    phi_0_arr = [60]

    phi_0_arr = romanova_phi_0.phi_0_arr
    mc2_arr = romanova_phi_0.mc2_arr

    mc2_calc = mc2_arr
    phi_0_calc = phi_0_arr
    # romanova_tuples = [(10, 0), (20, 28), (30, 43), (40, 52), (50, 60), (60, 64), (70, 68), (80, 71), (90, 74),
    #                    (100, 77)]

    # mc2_calc = []
    # phi_0_calc = []
    # for tuple in romanova_tuples:
    #     mc2 = tuple[0]
    #     phi_0 = tuple[1]
    #
    #     mc2_calc.append(mc2)
    #     phi_0_calc.append(phi_0)

    plot_flag = True
    theta_obs = 60
    beta_mu = 20
    a_portion = 0.2

    # a_portion_start = 0.44
    # a_portion_stop = 0.66
    # a_portion_arr = np.linspace(a_portion_start, a_portion_stop, len(mc2_calc))
    # print(a_portion_arr)

    # a_portion_arr = [0.2]
    a_portion_start = 0.2
    a_portion_stop = 0.5
    a_portion_arr = np.linspace(a_portion_start, a_portion_stop, len(mc2_calc))

    parameter_space = zip(repeat(mu), repeat(theta_obs), repeat(beta_mu), mc2_calc, a_portion_arr, phi_0_calc,
                          repeat(plot_flag), repeat(False))

    # parameter_space = zip(repeat(mu), repeat(theta_obs), repeat(beta_mu), mc2_calc, repeat(a_portion), phi_0_calc,
    #                       repeat(plot_flag), repeat(False))

    with mp.Pool(processes=config.N_cpus) as pool:
        pool.starmap(main.calc_and_save_for_configuration, parameter_space, chunksize=1)
