import matplotlib.pyplot as plt
import numpy as np
import scienceplots
import matplotlib as mpl

import pathService
import config
import newService
import save

plt.style.use(['science', 'notebook', 'grid'])
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'STIXGeneral'


def plot_L_to_mc2(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0):
    # plot_L_to_accr_rate

    # L_x_arr = np.empty(len(mc2_arr))
    # for i, mc2 in enumerate(mc2_arr):
    #     L_x_arr[i] = save.load_L_x(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    L_total_mc2 = np.empty((len(mc2_arr), config.N_phase))
    for i, mc2 in enumerate(mc2_arr):
        L_total_mc2[i] = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    # нормировка на L_max в каждом mc2
    buf = L_total_mc2 / np.max(L_total_mc2, axis=-1).ravel()[:, np.newaxis]

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=buf)
    # data_to_plot = np.empty((len(mc2_arr), config.N_phase_for_plot))
    # for i, L in enumerate(buf):
    #     data_to_plot[i] = newService.extend_arr_for_plot(L)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))

    x_axis_label = r'$\Phi$'
    # y_axis_label = r'$mean L_{iso} [erg/s]$'
    y_axis_label = r'$\dot{m}$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    im = ax['a'].pcolormesh(config.phase_for_plot, np.linspace(0, 1, len(mc2_arr)), data_to_plot)
    ax['a'].set_yticks(np.linspace(0, 1, len(mc2_arr)), mc2_arr)

    clb = plt.colorbar(im, pad=0.15, format="{x:.2}")
    clb.set_label(r'$L_{iso} \cdot max(L_{iso})^{-1}$', fontsize=24)

    folder = 'L_to_mass/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=a_portion,
                                   phi_0=phi_0, folder=folder)
    file_name = 'map_contour'
    save.save_figure(fig, save_dir, file_name)


def plot_L_to_a_portion(theta_obs, beta_mu, mc2, a_portion_arr, phi_0):
    L_total_a_portion = np.empty((len(mc2_arr), config.N_phase))
    for i, a_portion in enumerate(a_portion_arr):
        L_total_a_portion[i] = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))

    # нормировка на L_nu_avg
    buf = L_total_a_portion / np.max(L_total_a_portion, axis=-1).ravel()[:, np.newaxis]
    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=buf)

    im = ax['a'].pcolormesh(config.phase_for_plot, np.linspace(0, 1, len(a_portion_arr)), data_to_plot)
    ax['a'].set_yticks(np.linspace(0, 1, len(a_portion_arr)), np.round(a_portion_arr, 2))

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$a$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$L_{iso} \cdot max(L_{iso})^{-1}$', fontsize=24)

    folder = 'L_to_a/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=mc2, a_portion=None,
                                   phi_0=phi_0, folder=folder)
    file_name = 'map_contour'
    save.save_figure(fig, save_dir, file_name)


if __name__ == '__main__':
    mu = 0.1e31
    beta_mu = 40
    mc2 = 100
    a_portion = 0.66
    phi_0 = 0
    theta_obs = 20

    mc2_arr = [30, 60, 100]
    plot_L_to_mc2(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    a_portion_arr = [0.22, 0.44, 0.66]
    plot_L_to_a_portion(theta_obs, beta_mu, mc2, a_portion_arr, phi_0)