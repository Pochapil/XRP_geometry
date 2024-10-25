import matplotlib.pyplot as plt
import numpy as np
import scienceplots
import matplotlib as mpl
import matplotlib.cm as cm

import pathService
import config
import newService
import save

plt.style.use(['science', 'notebook', 'grid'])
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'STIXGeneral'


def plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0):
    # try_sky_map(obs_i_angle_arr)
    # obs_i_angle = np.linspace(0, 180, 19)
    L_x = save.load_L_x(mu, 10, beta_mu, mc2, a_portion, phi_0)

    theta_obs_arr = np.linspace(10, 90, 9)
    data_array = np.empty((len(theta_obs_arr) * 2 - 1, config.N_phase))
    # roll -- циклическая перестановка - делаем так как симметрическая задача и для угла 180 - theta будет симметрично
    # для сдвига на полфазы -- можно расчитать только до 90 а потом для других переставить и получить для до 180
    for i, theta_obs in enumerate(theta_obs_arr):
        L_total = save.load_L_total(mu, int(theta_obs), beta_mu, mc2, a_portion, phi_0)
        data_array[i] = L_total
        if i != len(theta_obs_arr) - 1:
            data_array[-i - 1] = np.roll(data_array[i], len(data_array[i]) // 2)

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_array / L_x)
    theta_obs_arr_to_plot = np.linspace(0, 180, 17)

    fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))
    im = ax['a'].pcolormesh(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='white', marker='*', s=300)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='white', marker='*', s=300)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='black', marker='*', s=100)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='black', marker='*', s=100)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\theta_{obs} \, [^{\circ}]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    clb = plt.colorbar(im)
    clb.set_label(r'$L_{iso} \cdot L_{x}^{-1}$', fontsize=24)

    prefix_folder = 'sky_map/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion, phi_0=phi_0,
                                   prefix_folder=prefix_folder)
    file_name = 'try_map'
    save.save_figure(fig, save_dir, file_name)

    fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))
    im = ax['a'].contourf(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot, levels=30)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='white', marker='*', s=300)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='white', marker='*', s=300)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='black', marker='*', s=100)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='black', marker='*', s=100)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\theta_{obs} \, [^{\circ}]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    clb = plt.colorbar(im)
    clb.set_label(r'$L_{iso} \cdot L_{x}^{-1}$', fontsize=24)

    prefix_folder = 'sky_map/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion, phi_0=phi_0,
                                   prefix_folder=prefix_folder)
    file_name = 'try_map_contour'
    save.save_figure(fig, save_dir, file_name)


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

    prefix_folder = 'L_to_mass/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=a_portion,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
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

    prefix_folder = 'L_to_a/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=mc2, a_portion=None,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
    file_name = 'map_contour'
    save.save_figure(fig, save_dir, file_name)


def plot_masses_PF_L_nu(theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr, energy_index=8):
    # PF(L_nu) много точек, берутся линии по mc, a
    ''' рисует графики PF от nu_L_nu (nu_L_nu усреднили по фазе)

    сначала в цикле читаю PF, заношу в 2D массив
    потом в цикле по fi_0 в мапу - ключ среднее по фазе nu_Lnu значение - значение PF массив
    одновременно для L_nu запоминаю fi_0, чтобы окрасить в цвет (1 график)

    рисую 3 графика ox - nu_L_nu, oy - PF:
    1 - множество точек для разных mc2, a - папки по a - точки окрашены по fi_0
    2 - множество точек для разных mc2, a - папки по a - точки окрашены для комбинаций mc2, a
    3 - множество точек для разных mc2, a - все точки - точки окрашены для комбинаций mc2, a'''

    line_color_dict = {0: 'blue', 1: 'green', 2: 'orange', 3: 'red', 4: 'purple', 5: 'black', 6: 'yellow'}
    marker_index = 0
    line_style = ['-', '--']
    marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}
    file_name = f'L_nu_of_energy_{config.energy_arr[energy_index]:.2f}_KeV_of_surfaces'

    PF_tensor = np.empty((len(a_portion_arr), len(mc2_arr), len(phi_0_arr)))

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    for a_index, a_portion in enumerate(a_portion_arr):
        for mc2_index, mc2 in enumerate(mc2_arr):
            L_nu_data_dict = {}
            color_dict = {}
            full_dict = {}
            for phi_0_index, phi_0 in enumerate(phi_0_arr):
                PF = save.load_PF(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                PF_tensor[a_index][mc2_index][phi_0_index] = PF[energy_index]

                L_nu_total = save.load_L_nu_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

                freq_arr = newService.get_frequency_from_energy(np.array(config.energy_arr))
                nu_L_nu_total = L_nu_total * freq_arr[:, np.newaxis]

                L_nu_data_dict[np.mean(nu_L_nu_total)] = PF_tensor[a_index][mc2_index][phi_0_index]
                color_dict[np.mean(nu_L_nu_total)] = phi_0_arr[phi_0_index]
                full_dict[np.mean(nu_L_nu_total)] = (PF_tensor[a_index][mc2_index][phi_0_index], phi_0_arr[phi_0_index])

            lists = sorted(L_nu_data_dict.items())  # sorted by key, return a list of tuples
            x, y = zip(*lists)  # unpack a list of pairs into two tuples

            lists = sorted(color_dict.items())
            buffer, colors = zip(*lists)

            lists = sorted(full_dict.items(), key=lambda item: item[1])
            x_sort, y_sort = zip(*lists)
            y_sort, colors_sort = zip(*y_sort)

            lists = sorted(full_dict.items(), key=lambda item: item[1][1])
            x_sort_fi_0, y_sort_fi_0 = zip(*lists)
            y_sort_fi_0, _ = zip(*y_sort_fi_0)

            colors = (np.array(colors_sort)) / 180

            if marker_index == 0:
                ax['a'].scatter(x_sort, y_sort, s=30, facecolors='none', edgecolors=cm.jet(colors))
            else:
                ax['a'].scatter(x_sort, y_sort, marker=marker_dict[marker_index % 4], color=cm.jet(colors))

            ax['a'].plot(x_sort_fi_0, y_sort_fi_0, color='black', alpha=0.2, linestyle=line_style[mc2_index])

        marker_index = 3

        x_axis_label = r'$\nu L_{\nu}$' + ' [erg/s]'
        y_axis_label = r'$PF_{' + f'{config.energy_arr[energy_index]:.2f}' + r'}$'
        ax['a'].set_xlabel(x_axis_label, fontsize=24)
        ax['a'].set_ylabel(y_axis_label, fontsize=24)
        ax['a'].set_xscale('log')
        # ax['a'].legend()

    bounds = phi_0_arr.copy()
    cmap = mpl.cm.jet
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax['a'], orientation='vertical', label=r'$\phi_0$',
                 pad=0.01)

    prefix_folder = 'PF_to_L_nu/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)

    file_name = f'theta={theta_obs} beta_mu={beta_mu} All_PF_to_L_nu.png'
    save.save_figure(fig, save_dir, file_name)


def plot_L_to_phi_0(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr):
    # plot_L_to_new_fi_0
    # ожидаем что phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]

    L_data = np.empty(((len(phi_0_arr) - 1) * 2, config.N_phase))
    for i, phi_0 in enumerate(phi_0_arr):
        L_data[i] = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    for i in range(8):
        # симметрично для 180 - fi_0
        L_data[9 + i + 1] = L_data[9 - i - 1][::-1]

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=L_data)
    phi_0_arr_for_plot = [20 * i for i in range(18)]

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    im = ax['a'].pcolormesh(config.phase_for_plot, phi_0_arr_for_plot, data_to_plot)
    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$L_{iso}$' + ' [erg/s]', fontsize=26)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\phi_0$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'L_to_phi_0/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = 'map_contour_L_iso'
    save.save_figure(fig, save_dir, file_name)

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    im = ax['a'].contourf(config.phase_for_plot, phi_0_arr_for_plot, data_to_plot, linewidths=0.5)
    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$L_{iso}$' + ' [erg/s]', fontsize=26)

    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'L_to_phi_0/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = 'map_contour_L_iso_c'
    save.save_figure(fig, save_dir, file_name)


if __name__ == '__main__':
    mu = 0.1e31
    beta_mu = 40
    mc2 = 100
    a_portion = 0.66
    phi_0 = 0
    theta_obs = 20

    mc2_arr = [30, 60, 100]
    # plot_L_to_mc2(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    a_portion_arr = [0.22, 0.44, 0.66]
    # plot_L_to_a_portion(theta_obs, beta_mu, mc2, a_portion_arr, phi_0)

    mc2_arr = [30, 100]
    a_portion_arr = [0.22, 0.66]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    # plot_masses_PF_L_nu(theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)

    mu = 0.1e31
    beta_mu = 80
    mc2 = 30
    a_portion = 1
    phi_0 = 0
    # plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    mu = 0.1e31
    beta_mu = 40
    mc2 = 100
    a_portion = 0.66
    phi_0 = 0
    theta_obs = 20
    plot_L_to_phi_0(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr)
