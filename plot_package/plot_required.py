import matplotlib.pyplot as plt
import numpy as np
import scienceplots
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.lines as mlines
import scipy.interpolate
from scipy import stats

import pathService
import config
import newService
import save
import accretingNS

plt.style.use(['science', 'notebook', 'grid'])
# чтобы стиль был похож на теховский
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'STIXGeneral'

mu = 0.1e31

theta_obs = 40
beta_mu = 20

mc2 = 60
a_portion = 0.22
phi_0 = 0


def func(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, R='rs'):
    curr_configuration = accretingNS.AccretingPulsarConfiguration(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    R_e = curr_configuration.top_column_for_calc.R_e
    if R == 'rs':
        R = 1.1 * curr_configuration.top_column_for_calc.ksi_shock * config.R_ns
    else:
        R = 0.9 * curr_configuration.R_disk
    theta = newService.get_theta_for_R(R, R_e)

    S = newService.get_A_normal(theta, R_e, a_portion, curr_configuration.dRe_div_Re)
    v = newService.get_free_fall_velocity(theta, R_e)
    M = mc2 * config.L_edd / config.c ** 2
    rho = M / (S * v)

    # T_e = 0.1e3 / config.k_bolc  # 0.1 / 3

    T_e = 3000 * 1.6e-12 # 3kev -> erg
    T_e /= config.k_bolc # erg -> K

    h = config.h_plank_evs
    k = config.k_bolc
    nu = 3e3 / h

    # wtf ????
    k_ff = 3.7e8 * (rho / config.mass_P) ** 2 / (rho * np.sqrt(T_e) * nu ** 3) * (1 - np.exp(-h * nu / (k * T_e)))

    res = np.sqrt(3 * k_ff * (k_ff + config.k))

    # tau_t = 1.4 * 1e-18 * M * np.sqrt(R_e / config.R_ns) / (a_portion * R / config.R_ns)
    # tau_t = 5 * mc2 * np.sqrt(R / R_e) / (a_portion * R / config.R_ns)

    # tau_t = 5 * mc2 * np.sqrt(R_e / config.R_ns) / (a_portion * R / config.R_ns)
    # tau_t = 0.32 * mc2 * np.sqrt(R_e / config.R_ns) / (a_portion * R / config.R_ns)

    # print(f'old tau_t {tau_t}')

    tau_t = (config.k * M * np.sqrt(R_e))/ (4 * np.pi * a_portion * R * np.sqrt(2 * config.G * config.M_ns))

    # print(f'new tau_t {tau_t}')

    # print(k_ff)
    # print(res)
    # print(rho)
    # print(tau_t)

    return res, rho, tau_t


mu = 0.1e31

theta_obs = 40
beta_mu = 0

n = 11
mc2_arr = [10 * i for i in range(1, n)]
a_portion = 0.22
phi_0 = 0
a_portion_arr = [0.25, 0.5, 0.6, 1]

n = len(mc2_arr)
res = [0] * n
rho = [0] * n
tau_t = [0] * n

fig, ax = plt.subplot_mosaic('abc', figsize=(21, 6))

R = ['rs', 'rm']
line_style = ['-', '--']
for j, r in enumerate(R):
    for a_portion in a_portion_arr:
        for i, mc2 in enumerate(mc2_arr):
            res[i], rho[i], tau_t[i] = func(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, R=r)

        suff = '; 1.1 Rs'
        if r == 'rm':
            suff = '; 0.9 Rm'

        ax['a'].plot(mc2_arr, tau_t, label=f'{a_portion}' + suff, linestyle=line_style[j])
        x_axis_label = r'$\dot{m}$'
        y_axis_label = r'$\tau_T$'
        ax['a'].set_xlabel(x_axis_label, fontsize=24)
        ax['a'].set_ylabel(y_axis_label, fontsize=24)
        ax['a'].set_yscale('log')

        ax['b'].plot(mc2_arr, res, label=f'{a_portion}' + suff, linestyle=line_style[j])
        x_axis_label = r'$\dot{m}$'
        y_axis_label = r'$\tau_{eff}$'
        ax['b'].set_xlabel(x_axis_label, fontsize=24)
        ax['b'].set_ylabel(y_axis_label, fontsize=24)
        ax['b'].set_yscale('log')

        ax['c'].plot(mc2_arr, rho, label=f'{a_portion}' + suff, linestyle=line_style[j])
        x_axis_label = r'$\dot{m}$'
        y_axis_label = r'$\rho$'
        ax['c'].set_xlabel(x_axis_label, fontsize=24)
        ax['c'].set_ylabel(y_axis_label, fontsize=24)
        ax['c'].set_yscale('log')


# plt.legend()
ax['a'].legend()
prefix_folder = 'free-free/'
save_dir = pathService.get_dir(mu=None, theta_obs=None, beta_mu=None, mc2=None, a_portion=None, phi_0=None,
                               prefix_folder=prefix_folder)
file_name = 'tau_rho'
save.save_figure(fig, save_dir, file_name)


# plt.show()
