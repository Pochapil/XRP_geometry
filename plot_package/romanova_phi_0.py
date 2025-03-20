import numpy as np
from matplotlib import pyplot as plt

import config

mu = 0.1e31
mc2 = 10

mc2_arr = np.array([10 * i for i in range(1, 11)])

M_accretion_rate = mc2_arr * config.L_edd / config.c ** 2
R_alfven = (mu ** 2 / (2 * M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)

# 0.43
P = 0.33 # 1.8, 2.4, 2.8
R_corot = (config.G * config.M_ns * (P / (2 * np.pi)) ** 2) ** (1 / 3)

phi_0_arr = 160 - 146 * (R_alfven / R_corot)



# phi_0_arr = [-83, -39, -17, -3, 7, 15, 21, 26, 31, 34, ]
# phi_0_arr = [-31, 4, 21, 32, 40, 46, 51, 55, 58, 61, ]
# phi_0_arr = [0, 29, 43, 52, 59, 64, 68, 72, 75, 77, ]

if __name__ == '__main__':
    phi_0_arr = 160 - 146 * (R_alfven / R_corot)

    print(R_alfven / R_corot)

    print(phi_0_arr)

    # M_accretion_rate = 5e18
    # print(M_accretion_rate*config.c**2/config.L_edd)
    B = 9e12
    mu = B * config.R_ns**3 /2
    print(mu)
    R_alfven = 1/2*(mu ** 2 / (2 * M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
    T = 9.8
    R_corot = (config.G * config.M_ns * (T / (2 * np.pi)) ** 2) ** (1 / 3)
    print(R_alfven / R_corot)

    print(f'{R_alfven=}')

    print('[', end='')
    for phi_0 in phi_0_arr:
        print(f'{round(phi_0)}', end=', ')
    print(']', end='')