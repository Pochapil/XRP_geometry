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

mpl.rcParams['image.cmap'] = 'magma'


x = [['A panel', 'A panel', 'edge'],
     ['C panel', '.',       'edge']]

x = [['A', 'A', 'B', 'C'],
     ['A', 'A', 'D', 'E']]



def func(mu, beta_mu, mc2, a_portion, phi_0):
    input_phi_0 = phi_0

    L_x = save.load_L_x(mu, 10, beta_mu, mc2, a_portion, phi_0)

    theta_obs_arr = np.linspace(0, 90, 10).astype(int)
    data_array = np.empty((len(theta_obs_arr) * 2 - 1, config.N_phase))

    for i, theta_obs in enumerate(theta_obs_arr):
        L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
        data_array[i] = L_total
        # для углов больших 90 будет симметрично относительно 90 (со сдвигом на полфазы)
        if i != len(theta_obs_arr) - 1:  # wtf -1 ?? - не будем трогать 90, его не дублируем
            data_array[-i - 1] = np.roll(data_array[i], len(data_array[i]) // 2)

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_array / L_x)
    theta_obs_arr_to_plot = np.linspace(0, 180, 19).astype(int)

    mosaic = [['a', 'a', 'b', 'c', 'f', 'h'],
              ['a', 'a', 'd', 'e', 'g', '.']]

    fig, ax = plt.subplot_mosaic(mosaic, figsize=(18, 7))
    im = ax['a'].contourf(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot, levels=30)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='white', marker='*', s=300)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='white', marker='*', s=300)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='black', marker='*', s=100)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='black', marker='*', s=100)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\theta_{\rm obs} \, [^{\circ}]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    #ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    clb = plt.colorbar(im)
    clb.set_label(r'$L_{\rm iso} / L_{x}$', fontsize=24)
    #clb.ax.tick_params(labelsize=ticks_labelsize)

    # ax['a'].axhline(y=60, color='black', lw=1.5)
    ax['a'].axhline(y=60, color='white', lw=1)

    phi_0_arr = [0, 40, 60, 90, 120, 140, 180]
    axes = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
    for i, phi_0 in enumerate(phi_0_arr):
        theta_obs = 60
        L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
        data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=L_total)
        ax[axes[i]].plot(config.phase_for_plot, data_to_plot, label=r'$\phi_0 = $' + f'{phi_0}', color='black')

    prefix_folder = 'skymap_with_profile/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'{beta_mu=} {a_portion=} phi_0={input_phi_0}'
    save.save_figure(fig, save_dir, file_name)

    mosaic = 'ab'

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_array / L_x)
    theta_obs_arr_to_plot = np.linspace(0, 180, 19).astype(int)

    fig, ax = plt.subplot_mosaic(mosaic, figsize=(18, 7))
    im = ax['a'].contourf(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot, levels=30)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='white', marker='*', s=300)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='white', marker='*', s=300)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='black', marker='*', s=100)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='black', marker='*', s=100)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\theta_{\rm obs} \, [^{\circ}]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    # ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    clb = plt.colorbar(im)
    clb.set_label(r'$L_{\rm iso} / L_{x}$', fontsize=24)
    # clb.ax.tick_params(labelsize=ticks_labelsize)

    # ax['a'].axhline(y=60, color='black', lw=1.5)
    ax['a'].axhline(y=60, color='white', lw=1)

    phi_0_arr = [0, 20, 40, 60]
    axes = ['b', 'c', 'd', 'e']
    for i, phi_0 in enumerate(phi_0_arr):
        theta_obs = 60
        L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
        data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=L_total)
        ax['b'].plot(config.phase_for_plot, data_to_plot, label=r'$\phi_0 = $' + f'{phi_0}')
    ax['b'].legend()

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$L_{\rm iso}$'
    ax['b'].set_xlabel(x_axis_label, fontsize=24)
    ax['b'].set_ylabel(y_axis_label, fontsize=24)


    plt.show()



def func2(mu, beta_mu, mc2, a_portion):
    roll = True
    fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))
    phi_0_arr = [0, 40, 60, 90, 120, 140, 180]
    for i, phi_0 in enumerate(phi_0_arr):
        theta_obs = 60
        L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

        data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=L_total)
        if roll:
            # id = np.argmin(np.abs(config.phase_for_plot - phi_0 / 360))
            # id = np.argmax(data_to_plot)
            id = -np.argmax(L_total)
            data_to_plot = np.roll(data_to_plot, id-22)
        ax['a'].plot(config.phase_for_plot, data_to_plot, label=r'$\phi_0 = $' + f'{phi_0}')
    ax['a'].legend()

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$L_{\rm iso}$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    ticks_labelsize = 22
    ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    prefix_folder = 'profiles_for_sky_map/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'{beta_mu=} {a_portion=}'
    if roll:
        file_name = file_name + ' rolled_max'
    save.save_figure(fig, save_dir, file_name)





mu = 0.1e31

theta_obs = 40
beta_mu = 20

mc2 = 60
a_portion = 0.22
phi_0 = 0

a_portion_arr = [0.22, 0.44, 0.66]
for a_portion in a_portion_arr:
    func2(mu, beta_mu, mc2, a_portion)

# func(mu, beta_mu, mc2, a_portion, phi_0)
# fig, ax = plt.subplot_mosaic(x, figsize=(18, 9))
# plt.show()