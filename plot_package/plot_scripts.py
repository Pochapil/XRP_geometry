import matplotlib.pyplot as plt
import numpy as np
import scienceplots
import matplotlib as mpl

import config
import newService
import save

plt.style.use(['science', 'notebook', 'grid'])
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'STIXGeneral'


# mpl.use('Agg')  # вроде ускоряет ??


def plot_total_luminosity_of_surfaces(L_surfs, save_dir=None):
    '''рисует вклады от поверхностей + сумму. без учета рассеяния'''

    legend_labels_arr = [r'$top_{pol}$', r'$top_{eq}$', r'$bottom_{pol}$', r'$bottom_{eq}$', 'sum']
    legend_labels_arr = [r'$\rm north polar$', r'$\rm north equatorial$', r'$\rm south polar$',
                         r'$\rm south equatorial$', r'$\sum$']
    legend_labels_arr = ['north polar', 'north equatorial', 'south polar', 'south equatorial', 'sum']

    colors_arr = ['red', 'red', 'green', 'green']
    marker_arr = ['.', '*', '+', '^']
    line_style_arr = [':', '-']

    fig, ax = plt.subplot_mosaic('a', figsize=(21, 10))
    for i, L_surf in enumerate(L_surfs):
        ax['a'].plot(config.phase_for_plot, newService.extend_arr_for_plot(L_surf), label=legend_labels_arr[i],
                     color=colors_arr[i], marker=marker_arr[i], linestyle=line_style_arr[i % 2], markersize=12)

    ax['a'].plot(config.phase_for_plot, newService.extend_arr_for_plot(np.sum(L_surfs, axis=0)),
                 label=legend_labels_arr[-1], color='black')

    x_axis_label = config.symbol_phase
    y_axis_label = r'$L_{\rm{iso}}$' + ' [erg/s]'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].legend()
    # plt.show()

    if save_dir is not None:
        file_name = 'total_luminosity_of_surfaces'
        save.save_figure(fig, save_dir, file_name)

    # ax['a'].set_ylabel('')
    # ax['a'].legend_.remove()
    # ax['a'].set(xlabel=None)
    # ax['a'].tick_params(axis='both', labelsize=16)


def plot_observer_angles(observer_phi, observer_theta, save_dir=None):
    '''рисует углы наблюдателя от фазы'''

    # file_name = 'observer_angles.txt'
    # data_array = main_service.load_arr_from_txt(working_folder, file_name)
    # observer_phi = data_array[0]
    # observer_theta = data_array[1]

    observer_phi = newService.extend_arr_for_plot(observer_phi)
    observer_theta = newService.extend_arr_for_plot(observer_theta)

    fig_title = 'Observer angles'

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    ax['a'].plot(config.phase_for_plot, np.rad2deg(observer_theta), label=r'$\theta_{observer}$')
    ax['a'].plot(config.phase_for_plot, np.rad2deg(observer_phi), label=r'$\phi_{observer}$')
    ax['a'].legend()

    x_axis_label = config.symbol_phase
    y_axis_label = r'$angle ~ [deg]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    if save_dir is not None:
        file_name = 'Observer_angles'
        save.save_figure(fig, save_dir, file_name)


def plot_Teff_to_ksi(R_e, T_eff, theta_range, save_dir=None):
    '''рисует Т от кси'''
    # file_name = 'surfaces_T_eff.txt'
    # data_array = main_service.load_arr_from_txt(working_folder, file_name)
    #
    # file_name = "save_theta_range.txt"
    # theta_range = main_service.load_arr_from_txt(working_folder, file_name)

    # with open(working_folder + 'save_values.txt') as f:
    #     lines = f.readlines()
    #     R_e = float(lines[0][6:-1])

    ksi = R_e * (np.sin(theta_range)) ** 2

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    ax['a'].plot(ksi, T_eff, color='black')
    x_axis_label = r'$\xi$'
    y_axis_label = r'$T_{\rm eff} ~ [K]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    if save_dir is not None:
        file_name = 'T_eff'
        save.save_figure(fig, save_dir, file_name)


def plot_PF_to_energy(L_nu_surfs, save_dir=None):
    '''рисует PF на кажой энергии. без учета рассеянного'''
    # sum(L_nu_surfs + L_nu_scatter)

    L_nu = np.sum(L_nu_surfs, axis=0)
    # L_nu = np.sum(L_nu_surfs, axis=0) + np.sum(L_nu_scatter, axis=0)

    PF = newService.get_PF(L_nu)

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    ax['a'].plot(config.energy_arr, PF, color='black')
    x_axis_label = r'$h \nu$' + ' [keV]'
    y_axis_label = 'PF'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    ax['a'].set_ylim(0, 1)
    if save_dir is not None:
        folder = 'L_nu/'
        file_name = 'PF'
        save.save_figure(fig, save_dir / folder, file_name)


def plot_L_nu(L_nu_surfs, save_dir=None):
    '''
    отрисовка массивов графиков L_nu на каждой энергии
    попробовать распараллелить
    '''
    # sum L_nu
    L_nu_to_plot = newService.extend_arr_for_plot(L_nu_surfs)

    for i, energy in enumerate(config.energy_arr):
        fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
        # for L_nu_surf in L_nu_surfs[:, i, :]:
        L_nu_to_plot = newService.extend_arr_for_plot(np.sum(L_nu_surfs[:, i, :], axis=0))
        PF = newService.get_PF(L_nu_to_plot)
        fig_title = r'$L_{\nu}^' + r'{' + f'{energy:.2f} keV' + r'}' + r'(\Phi)$' + '  ' + f'PF = {PF:.3f}'
        ax['a'].plot(config.phase_for_plot, L_nu_to_plot, color='black',
                     label=r'$L_{\nu}^' + r'{' + f'{energy:.2f}' + r'}' + r'(\Phi)$')

        x_axis_label = r'$h \nu$' + ' [keV]'
        y_axis_label = r'$L_{\nu} \: [erg \cdot s^{-1} \cdot hz^{-1}]$'
        ax['a'].set_xlabel(x_axis_label, fontsize=24)
        ax['a'].set_ylabel(y_axis_label, fontsize=24)
        # ax['a'].legend()
        fig.suptitle(fig_title, fontsize=16)

        if save_dir is not None:
            folder = 'L_nu/'
            file_name = f'L_nu_of_energy_{energy:.2f}_KeV_of_surfaces'
            save.save_figure(fig, save_dir / folder, file_name)


def plot_L_nu_all_in_one(L_nu_surfs, save_dir=None):
    '''отрисовка графиков L_nu на каждой энергии на 1 графике'''
    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    for i, energy in enumerate(config.energy_arr):
        L_nu_to_plot = newService.extend_arr_for_plot(np.sum(L_nu_surfs[:, i, :], axis=0))
        ax['a'].plot(config.phase_for_plot, L_nu_to_plot, label=f'{energy:.2f} keV')

        x_axis_label = config.symbol_phase
        y_axis_label = r'$L_{\nu} \: [erg \cdot s^{-1} \cdot hz^{-1}]$'
        ax['a'].set_xlabel(x_axis_label, fontsize=24)
        ax['a'].set_ylabel(y_axis_label, fontsize=24)

    if save_dir is not None:
        folder = 'L_nu/'
        file_name = 'L_nu'
        save.save_figure(fig, save_dir / folder, file_name)


def plot_L_nu_on_phase(L_nu_surfs, save_dir=None, phase_index=0):
    '''рисует L_nu от энергии на заданной фазе'''
    # phase_index - индекс фазы для L_nu(nu)

    L_nu_on_phase = np.sum(L_nu_surfs[:, :, phase_index], axis=0)

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    ax['a'].plot(config.energy_arr, L_nu_on_phase, color='black')

    x_axis_label = r'$h \nu$' + ' [keV]'
    y_axis_label = r'$L_{\nu} \: [erg \cdot s^{-1} \cdot hz^{-1}]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    if save_dir is not None:
        folder = 'L_nu/'
        file_name = 'L_nu(nu)'
        save.save_figure(fig, save_dir / folder, file_name)


def plot_L_nu_avg(L_nu_surfs, save_dir=None):
    '''рисует среднее L_nu от энергии'''
    L_nu_avg_on_phase = np.mean(np.sum(L_nu_surfs, axis=0), axis=1)

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    ax['a'].plot(config.energy_arr, L_nu_avg_on_phase, color='black')

    x_axis_label = r'$h \nu$' + ' [keV]'
    y_axis_label = r'$L_{\nu} \: [erg \cdot s^{-1} \cdot hz^{-1}]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    ax['a'].set_xscale('log')
    ax['a'].set_yscale('log')
    # plt.xscale('log')

    if save_dir is not None:
        folder = 'L_nu/'
        file_name = 'L_nu(nu)_avg_log_log'
        save.save_figure(fig, save_dir / folder, file_name)


def plot_L_nu_with_bb(L_nu_surfs, T_eff, save_dir=None):
    '''рисует среднее L_nu от энергии + черное тело с максимальной Т и минимальной Т по колонке'''
    freq_arr = newService.get_frequency_from_energy(np.array(config.energy_arr))

    # file_name = "surfaces_T_eff.txt"
    # T_eff = main_service.load_arr_from_txt(config.full_file_folder, file_name)

    black_body_max = newService.plank_energy_on_frequency(freq_arr, np.max(T_eff))
    black_body_min = newService.plank_energy_on_frequency(freq_arr, np.min(T_eff))

    L_nu_avg_on_phase = np.mean(np.sum(L_nu_surfs, axis=0), axis=1)
    coeff_max = np.max(L_nu_avg_on_phase) / np.max(black_body_max)
    coeff_min = np.max(L_nu_avg_on_phase) / np.max(black_body_min)

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))

    ax['a'].plot(config.energy_arr, coeff_max * np.array(black_body_max), label=f'black body max T={np.max(T_eff):.1f}')
    ax['a'].plot(config.energy_arr, coeff_min * np.array(black_body_min), label=f'black body min T={np.min(T_eff):.1f}')
    ax['a'].plot(config.energy_arr, L_nu_avg_on_phase, label=r'$L_{\nu} \, avg$', color='black')

    x_axis_label = r'$h \nu$' + ' [keV]'
    y_axis_label = r'$L_{\nu} \: [erg \cdot s^{-1} \cdot hz^{-1}]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    ax['a'].legend()
    ax['a'].set_xscale('log')
    ax['a'].set_yscale('log')
    ax['a'].set_ylim(np.min(L_nu_avg_on_phase))

    if save_dir is not None:
        folder = 'L_nu/'
        file_name = 'L_nu(nu)_avg_and_black_body'
        save.save_figure(fig, save_dir / folder, file_name)


def plot_scatter_L(L_surfs, L_scatter, save_dir=None, log=False):
    '''рисует полный вклад каждой поверхности + рассеяния '''

    legend_labels_arr = [r'$top_{pol}$', r'$top_{eq}$', r'$bottom_{pol}$', r'$bottom_{eq}$', 'sum']
    legend_labels_arr = ['north polar', 'north equatorial', 'south polar', 'south equatorial', 'sum']

    colors_arr = ['red', 'red', 'green', 'green']
    marker_arr = ['.', '*', '+', '^']
    line_style_arr = [':', '-']

    fig, ax = plt.subplot_mosaic('a', figsize=(21, 10))
    for i, L_surf in enumerate(L_surfs):
        fillstyle = 'none'
        ax['a'].plot(config.phase_for_plot, newService.extend_arr_for_plot(L_surf), label=legend_labels_arr[i],
                     color=colors_arr[i], marker=marker_arr[i], linestyle=line_style_arr[i % 2], markersize=12,
                     fillstyle=fillstyle)

    ax['a'].plot(config.phase_for_plot, newService.extend_arr_for_plot(L_scatter[0]), label='north scatter',
                 color='purple', marker="1", markersize=12)
    ax['a'].plot(config.phase_for_plot, newService.extend_arr_for_plot(L_scatter[1]), label='south scatter',
                 color='blue', marker="2", markersize=12)

    ax['a'].plot(config.phase_for_plot,
                 newService.extend_arr_for_plot(np.sum(L_surfs, axis=0) + np.sum(L_scatter, axis=0)),
                 label=legend_labels_arr[-1], color='black')

    x_axis_label = config.symbol_phase
    y_axis_label = r'$L_{\rm{iso}}$' + ' [erg/s]'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].legend()

    if log:
        ax['a'].set_yscale('log')
        # ax['a'].set_ylim(1e12, 1e16)
        # ax['a'].set_ylim(1e35, 5e39)

    if save_dir is not None:
        folder = 'scattered_on_magnet_lines/'
        file_name = 'total_luminosity_of_surfaces_with_scatter'
        if log:
            file_name+='_log'
        save.save_figure(fig, save_dir / folder, file_name)


def plot_PF_to_energy_with_scatter(L_nu_surfs, L_nu_scatter, save_dir=None):
    '''рисуем PF от энергии с учетом рассеяния'''
    L_nu = np.sum(L_nu_surfs, axis=0) + np.sum(L_nu_scatter, axis=0)
    PF = newService.get_PF(L_nu)

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    ax['a'].plot(config.energy_arr, PF, color='black')

    x_axis_label = r'$h \nu$' + ' [keV]'
    y_axis_label = 'PF'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    ax['a'].set_ylim(0, 1)

    if save_dir is not None:
        folder = 'scattered_on_magnet_lines/' + 'L_nu/'
        file_name = 'PF'
        save.save_figure(fig, save_dir / folder, file_name)


def plot_scatter_L_nu(L_nu_surfs, L_nu_scatter, save_dir=None):
    '''
    отрисовка массивов графиков L_nu рассеяния на каждой энергии
    попробовать распараллелить
    '''
    L_nu_surfs_total = np.sum(L_nu_surfs, axis=0)
    L_nu_total = L_nu_surfs_total + np.sum(L_nu_scatter, axis=0)

    legend_labels_arr = [r'$top_{scatter}$', r'$bot_{scatter}$']
    colors_arr = ['blue', 'red']
    for i, energy in enumerate(config.energy_arr):
        fig, ax = plt.subplot_mosaic('a', figsize=(21, 10))

        for j, L_nu in enumerate(L_nu_scatter):
            ax['a'].plot(config.phase_for_plot, newService.extend_arr_for_plot(L_nu[i]), label=legend_labels_arr[j],
                         color=colors_arr[j])

        ax['a'].plot(config.phase_for_plot, newService.extend_arr_for_plot(L_nu_surfs_total[i]),
                     label='without scatter', color='green')
        ax['a'].plot(config.phase_for_plot, newService.extend_arr_for_plot(L_nu_total[i]), label='sum', color='black')

        x_axis_label = config.symbol_phase
        y_axis_label = r'$L_{\nu} \: [erg \cdot s^{-1} \cdot hz^{-1}]$'
        ax['a'].set_xlabel(x_axis_label, fontsize=24)
        ax['a'].set_ylabel(y_axis_label, fontsize=24)
        ax['a'].legend()

        PF = newService.get_PF(L_nu_total[i])
        fig_title = r'$L_{\nu}^' + r'{' + f'{energy:.2f} keV' + r'}' + r'(\Phi)$' + '  ' + f'PF = {PF:.3f}'
        fig.suptitle(fig_title, fontsize=16)

        if save_dir is not None:
            folder = 'scattered_on_magnet_lines/' + 'L_nu/'
            file_name = f'L_nu_of_energy_{energy:.2f}_KeV_of_surfaces'
            save.save_figure(fig, save_dir / folder, file_name)


if __name__ == '__main__':
    import pathService
    import save

    mu = 0.1e31
    theta_obs = 60
    beta_mu = 40
    mc2 = 100
    a_portion = 0.22
    phi_0 = 40

    L_surfs = save.load_L_surfs(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    L_scatter = save.load_L_scatter(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    cur_dir_saved = pathService.PathSaver(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    folder = 'txt/'
    # cur_path = path to save txt !!!
    cur_path = cur_dir_saved.get_path()

    plot_scatter_L(L_surfs, L_scatter, save_dir=cur_path / 'logs', log=True)
