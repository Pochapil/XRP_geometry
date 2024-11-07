import matplotlib.pyplot as plt
import numpy as np
import scienceplots
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as colors

import pathService
import config
import newService
import save

plt.style.use(['science', 'notebook', 'grid'])
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'STIXGeneral'

ticks_labelsize = 20
mpl.rcParams['xtick.labelsize'] = ticks_labelsize
mpl.rcParams['ytick.labelsize'] = ticks_labelsize


# plt.rcParams['axes.labelsize'] = 42
# ticks_labelsize = 18


def plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0):
    '''рисует карты неба (излучательный коэффициент в зависимости от положения наблюдателя)
    по другому - куча профилей. для разных наблюдателей. характеризует как источник излучает в пространство
    '''

    # try_sky_map(obs_i_angle_arr)
    # obs_i_angle = np.linspace(0, 180, 19)
    L_x = save.load_L_x(mu, 10, beta_mu, mc2, a_portion, phi_0)

    theta_obs_arr = np.linspace(10, 90, 9).astype(int)
    data_array = np.empty((len(theta_obs_arr) * 2 - 1, config.N_phase))
    # roll -- циклическая перестановка - делаем так как симметрическая задача и для угла 180 - theta будет симметрично
    # для сдвига на полфазы -- можно расчитать только до 90 а потом для других переставить и получить для до 180
    for i, theta_obs in enumerate(theta_obs_arr):
        L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
        data_array[i] = L_total
        # для углов больших 90 будет симметрично относительно 90 (со сдвигом на полфазы)
        if i != len(theta_obs_arr) - 1:  # wtf -1 ??
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
    y_axis_label = r'$\theta_{\rm obs} \, [^{\circ}]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    clb = plt.colorbar(im)
    clb.set_label(r'$L_{\rm iso} / L_{x}$', fontsize=24)  # r'$L_{\rm iso} \cdot L_{x}^{-1}$'

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
    y_axis_label = r'$\theta_{\rm obs} \, [^{\circ}]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    clb = plt.colorbar(im)
    clb.set_label(r'$L_{\rm iso} / L_{x}$', fontsize=24)

    prefix_folder = 'sky_map/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion, phi_0=phi_0,
                                   prefix_folder=prefix_folder)
    file_name = 'try_map_contour'
    save.save_figure(fig, save_dir, file_name)


def plot_L_to_mc2(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0):
    '''рисует 2d L от фазы и m '''
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

    fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))
    # чтобы сделать хороший масштаб и чтобы каждый профиль имел равную толщину сначала рисуем от 0 до 1 а потом
    # заменяем надписи на значения м
    im = ax['a'].pcolormesh(config.phase_for_plot, np.linspace(0, 1, len(mc2_arr)), data_to_plot)
    ax['a'].set_yticks(np.linspace(0, 1, len(mc2_arr)), mc2_arr)

    x_axis_label = r'$\Phi$'
    # y_axis_label = r'$mean L_{iso} [erg/s]$'
    y_axis_label = r'$\dot{m}$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    clb = plt.colorbar(im, pad=0.03, format="{x:.2}")  # pad=0.15
    clb.set_label(r'$\widetilde{L}_{\rm iso}$', fontsize=24)  # \cdot max(L_{\rm iso})^{-1}

    prefix_folder = 'L_to_mass/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=a_portion,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
    file_name = 'map_contour'
    save.save_figure(fig, save_dir, file_name)


def plot_L_to_a_portion(theta_obs, beta_mu, mc2, a_portion_arr, phi_0):
    '''рисует 2d L от фазы и a'''
    L_total_a_portion = np.empty((len(a_portion_arr), config.N_phase))
    for i, a_portion in enumerate(a_portion_arr):
        L_total_a_portion[i] = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    # нормировка на L_max
    buf = L_total_a_portion / np.max(L_total_a_portion, axis=-1).ravel()[:, np.newaxis]
    # buf = L_total_a_portion / np.max(np.max(L_total_a_portion, axis=-1))
    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=buf)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    # чтобы сделать хороший масштаб и чтобы каждый профиль имел равную толщину сначала рисуем от 0 до 1 а потом
    # заменяем надписи на значения а
    levels = 10 ** np.linspace(np.log10(np.min(data_to_plot)), np.log10(np.max(data_to_plot)), 10)
    im = ax['a'].pcolormesh(config.phase_for_plot, np.linspace(0, 1, len(a_portion_arr)), data_to_plot,
                            # norm=colors.LogNorm(vmin=data_to_plot.min(), vmax=data_to_plot.max())
                            norm=colors.TwoSlopeNorm(0.92, vmin=data_to_plot.min(), vmax=data_to_plot.max())
                            )
    ax['a'].set_yticks(np.linspace(0, 1, len(a_portion_arr)), np.round(a_portion_arr, 2))

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$a$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$\widetilde{L}_{\rm iso}$', fontsize=24)

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

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    for a_index, a_portion in enumerate(a_portion_arr):
        for mc2_index, mc2 in enumerate(mc2_arr):
            full_dict = {}
            for phi_0_index, phi_0 in enumerate(phi_0_arr):
                # PF_arr = save.load_PF_nu(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                # PF = PF_arr[energy_index]
                L_nu_total = save.load_L_nu_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)[energy_index]
                PF = newService.get_PF(L_nu_total)

                # freq_arr = newService.get_frequency_from_energy(np.array(config.energy_arr))
                freq = newService.get_frequency_from_energy(np.array(config.energy_arr[energy_index]))
                nu_L_nu_total = L_nu_total * freq

                full_dict[np.mean(nu_L_nu_total)] = (PF, phi_0_arr[phi_0_index])

            # lists = sorted(full_dict.items(), key=lambda item: item[1])
            # x_sort, y_sort = zip(*lists)
            # y_sort, colors_sort = zip(*y_sort)
            #
            # lists = sorted(full_dict.items(), key=lambda item: item[1][1])
            # x_sort_fi_0, y_sort_fi_0 = zip(*lists)
            # y_sort_fi_0, _ = zip(*y_sort_fi_0)

            # в словаере лежит среднее : (PF, phi_0). сортируем его по значению phi_0, чтобы посмотреть как зависит от
            # фи0 - будем рисовать по возрастанию фи0
            list_tuples = sorted(full_dict.items(), key=lambda item: item[1][1])
            x_sort_phi_0, y_sort_phi_0 = zip(*list_tuples)
            y_sort_phi_0, colors_sort_phi_0 = zip(*y_sort_phi_0)

            colors = (np.array(colors_sort_phi_0)) / np.max(colors_sort_phi_0)

            if marker_index == 0:
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=30, facecolors='none', edgecolors=cm.jet(colors))
            else:
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, marker=marker_dict[marker_index % 4], color=cm.jet(colors))

            ax['a'].plot(x_sort_phi_0, y_sort_phi_0, color='black', alpha=0.2, linestyle=line_style[mc2_index])

        marker_index = 3

    x_axis_label = r'$\nu L_{\nu}$' + r'$\rm [erg/s]$'
    y_axis_label = r'$PF_{' + f'{config.energy_arr[energy_index]:.2f}' + r'}$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].set_xscale('log')

    # ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)
    # ax['a'].xaxis.set_tick_params(labelsize=ticks_labelsize)
    # ax['a'].yaxis.set_tick_params(labelsize=ticks_labelsize)
    # ax['a'].legend()

    bounds = phi_0_arr.copy()
    cmap = mpl.cm.jet
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax['a'], orientation='vertical', pad=0.01)
    cb.set_label(label=r'$\varphi_0 ~ [^\circ]$', fontsize=24)

    prefix_folder = 'PF_to_L_nu/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)

    file_name = f'theta={theta_obs} beta_mu={beta_mu} All_PF_to_L_nu.png'
    save.save_figure(fig, save_dir, file_name)


def plot_masses_PF_L(theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr):
    marker_index = 0
    line_style = ['-', '--']
    marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    for a_index, a_portion in enumerate(a_portion_arr):
        for mc2_index, mc2 in enumerate(mc2_arr):
            full_dict = {}
            for phi_0_index, phi_0 in enumerate(phi_0_arr):
                L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                PF = newService.get_PF(L_total)
                full_dict[np.mean(L_total)] = (PF, phi_0_arr[phi_0_index])

            # в словаере лежит среднее : (PF, phi_0). сортируем его по значению phi_0
            list_tuples = sorted(full_dict.items(), key=lambda item: item[1][1])
            x_sort_phi_0, y_sort_phi_0 = zip(*list_tuples)
            y_sort_phi_0, colors_sort_phi_0 = zip(*y_sort_phi_0)

            colors = (np.array(colors_sort_phi_0)) / np.max(colors_sort_phi_0)

            if marker_index == 0:
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=30, facecolors='none', edgecolors=cm.jet(colors))
            else:
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, marker=marker_dict[marker_index % 4], color=cm.jet(colors))

            ax['a'].plot(x_sort_phi_0, y_sort_phi_0, color='black', alpha=0.2, linestyle=line_style[mc2_index])

        marker_index = 3

    x_axis_label = r'$L_{\rm iso}$' + r'$\rm [erg/s]$'
    y_axis_label = r'$PF$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].set_xscale('log')

    bounds = phi_0_arr.copy()
    cmap = mpl.cm.jet
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax['a'], orientation='vertical', pad=0.01)
    cb.set_label(label=r'$\varphi_0 ~ [^\circ]$', fontsize=24)

    prefix_folder = 'PF_to_L/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)

    file_name = f'theta={theta_obs} beta_mu={beta_mu} All_PF_to_L.png'
    save.save_figure(fig, save_dir, file_name)


def plot_L_to_phi_0(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr, flag_same_res=False):
    # plot_L_to_new_fi_0
    # ожидаем что phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]

    L_data = np.empty(((len(phi_0_arr)) * 2 - 1, config.N_phase))
    for i, phi_0 in enumerate(phi_0_arr):
        L_data[i] = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    for i in range(9):
        # симметрично для 180 - fi_0
        L_data[9 + i + 1] = L_data[9 - i - 1][::-1]

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=L_data)
    phi_0_arr_for_plot = [20 * i for i in range(19)]

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    im = ax['a'].pcolormesh(config.phase_for_plot, phi_0_arr_for_plot, data_to_plot)
    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$L_{\rm iso}$' + r'$\rm [erg/s]$', fontsize=26)

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$\varphi_0 ~ [^\circ]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    # ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    prefix_folder = 'L_to_phi_0/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = 'map_contour_L_iso'
    save.save_figure(fig, save_dir, file_name)

    # это чтобы сделать равную сетку. у нас меньше разрешение по фи0 чем по фазе поэтому можем проредить по фазе
    # мб перенести ниже и отрисовывать так только контуры
    if flag_same_res:
        step = 5
        delta = 1 if config.N_phase % step != 0 else 0
        L_data_same_res = np.empty(((len(phi_0_arr)) * 2 - 1, config.N_phase // step + delta))
        for i in range((len(phi_0_arr)) * 2 - 1):
            L_data_same_res[i] = L_data[i][::step]

        data_to_plot = np.apply_along_axis(newService.extend_arr_twice_for_plot, axis=-1, arr=L_data_same_res)
        config.phase_for_plot = np.linspace(0, 2, len(data_to_plot[0]))

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    im = ax['a'].contourf(config.phase_for_plot, phi_0_arr_for_plot, data_to_plot, levels=30)
    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$L_{\rm iso}$' + r'$\rm [erg/s]$', fontsize=26)

    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    # ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    prefix_folder = 'L_to_phi_0/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = 'map_contour_L_iso_c'
    save.save_figure(fig, save_dir, file_name)


def plot_PF_contour(mu, mc2, a_portion, phi_0):
    # PF(L_nu) много точек, берутся линии по mc, a

    theta_obs_arr = np.linspace(10, 90, 9).astype(int)
    beta_mu_arr = np.linspace(10, 90, 9).astype(int)

    final_final_array = np.empty((len(theta_obs_arr), len(beta_mu_arr)))

    for i, theta_obs in enumerate(theta_obs_arr):
        for j, beta_mu in enumerate(beta_mu_arr):
            L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            final_final_array[i][j] = newService.get_PF(L_total)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    x, y = np.meshgrid(beta_mu_arr, theta_obs_arr)
    z = final_final_array

    cs = ax['a'].contourf(x, y, z, 90)
    cbar = fig.colorbar(cs, pad=0.01)
    cbar.ax.set_ylabel('PF')

    x_axis_label = r'$\chi$'
    y_axis_label = r'$\theta_{obs}$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'PF_contour/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=mc2, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'a={a_portion}' + ' ' + f'phi_0={phi_0}'
    save.save_figure(fig, save_dir, file_name)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    cs = ax['a'].pcolormesh(x, y, z)
    cbar = fig.colorbar(cs, pad=0.01)
    cbar.ax.set_ylabel('PF')

    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'PF_contour/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=mc2, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'a={a_portion}' + ' ' + f'phi_0={phi_0} colormesh'
    save.save_figure(fig, save_dir, file_name)


def plot_PF_to_chi_theta(mu, mc2, a_portion, phi_0):
    theta_obs_arr = np.linspace(10, 90, 9).astype(int)
    beta_mu_arr = np.linspace(10, 90, 9).astype(int)

    final_final_array = np.zeros((len(theta_obs_arr), len(beta_mu_arr)))

    dict_chi_plus_theta = {}
    dict_chi_minus_theta = {}

    for i, theta_obs in enumerate(theta_obs_arr):
        for j, beta_mu in enumerate(beta_mu_arr):
            L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            final_final_array[i][j] = newService.get_PF(L_total)

            plus = theta_obs_arr[i] + beta_mu_arr[j]
            dict_chi_plus_theta.setdefault(plus, []).append(final_final_array[i][j])

            minus = beta_mu_arr[j] - theta_obs_arr[i]
            dict_chi_minus_theta.setdefault(minus, []).append(final_final_array[i][j])

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    for plus in dict_chi_plus_theta:
        arr = [plus] * len(dict_chi_plus_theta[plus])
        ax['a'].scatter(arr, dict_chi_plus_theta[plus], color='black')

    x_axis_label = r'$\chi + \theta_{obs}$'
    y_axis_label = r'$PF$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'PF_to_plus/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'a={a_portion}' + ' ' + f'mc2={mc2}' + ' ' + f'phi_0={phi_0}'
    save.save_figure(fig, save_dir, file_name)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    for minus in dict_chi_minus_theta:
        arr = [minus] * len(dict_chi_minus_theta[minus])
        ax['a'].scatter(arr, dict_chi_minus_theta[minus], color='black')

    x_axis_label = r'$\chi - \theta_{obs}$'
    y_axis_label = r'$PF$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'PF_to_minus/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'a={a_portion}' + ' ' + f'mc2={mc2}' + ' ' + f'phi_0={phi_0}'
    save.save_figure(fig, save_dir, file_name)


def plot_L_iso_to_m(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0):
    L_iso_arr = np.empty(len(mc2_arr))
    for i, mc2 in enumerate(mc2_arr):
        L_iso_arr[i] = np.mean(save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0), axis=-1)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    ax['a'].scatter(mc2_arr, L_iso_arr, s=30, facecolors='none', edgecolors='black')
    ax['a'].plot(mc2_arr, L_iso_arr, color='black', alpha=0.2, linestyle='-')

    x_axis_label = r'$\dot{m}$'
    y_axis_label = r'$L_{\rm iso}$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'L_iso_to_m/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'a={a_portion} phi_0={phi_0}'
    save.save_figure(fig, save_dir, file_name)


def plot_table(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0):
    R_e_arr = np.empty(len(mc2_arr))
    ksi_shock_arr = np.empty(len(mc2_arr))
    L_x_arr = np.empty(len(mc2_arr))

    for i, mc2 in enumerate(mc2_arr):
        R_e_arr[i] = save.load_R_e(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
        ksi_shock_arr[i] = save.load_ksi_shock(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
        L_x_arr[i] = save.load_L_x(mu, theta_obs, beta_mu, mc2, a_portion, phi_0) * 1e-38

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))

    ax['a'].scatter(mc2_arr, R_e_arr, s=30, facecolors='none', edgecolors='red', label=r'$\frac{R_e}{R_*}$')
    ax['a'].plot(mc2_arr, R_e_arr, color='black', alpha=0.2, linestyle='-')

    ax['a'].scatter(mc2_arr, ksi_shock_arr, s=30, facecolors='none', edgecolors='blue', label=r'$\xi_{s}$')
    ax['a'].plot(mc2_arr, ksi_shock_arr, color='black', alpha=0.2, linestyle='-')

    ax['a'].scatter(mc2_arr, L_x_arr, s=30, facecolors='none', edgecolors='green', label=r'$\frac{L_{\rm x}}{10^{38}}$')
    ax['a'].plot(mc2_arr, L_x_arr, color='black', alpha=0.2, linestyle='-')

    x_axis_label = r'$\dot{m}$'
    y_axis_label = r'$\xi$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)
    ax['a'].legend()

    prefix_folder = 'table/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'a={a_portion}'
    save.save_figure(fig, save_dir, file_name)


def plot_table_together(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0):
    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))

    fig1, ax1 = plt.subplot_mosaic('a', figsize=(9, 6))

    color_arr = ['red', 'green', 'blue', 'purple']

    for i, a_portion in enumerate(a_portion_arr):
        R_e_arr = np.empty(len(mc2_arr))
        ksi_shock_arr = np.empty(len(mc2_arr))
        L_x_arr = np.empty(len(mc2_arr))
        for j, mc2 in enumerate(mc2_arr):
            R_e_arr[j] = save.load_R_e(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            ksi_shock_arr[j] = save.load_ksi_shock(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            L_x_arr[j] = save.load_L_x(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

        ax['a'].scatter(mc2_arr, ksi_shock_arr, s=30, facecolors='none', edgecolors=color_arr[i],
                        label=r'$\xi_{s}$' + ' for ' + f'a={a_portion}')
        ax['a'].plot(mc2_arr, ksi_shock_arr, color='black', alpha=0.2, linestyle='-')

        ax1['a'].scatter(mc2_arr, L_x_arr, s=30, facecolors='none', edgecolors=color_arr[i],
                         label=r'$L_{\rm x}$' + ' for ' + f'a={a_portion}')
        ax1['a'].plot(mc2_arr, L_x_arr, color='black', alpha=0.2, linestyle='-')

    ax['a'].scatter(mc2_arr, R_e_arr, s=30, facecolors='none', edgecolors='black', label=r'$\frac{R_e}{R_*}$')
    ax['a'].plot(mc2_arr, R_e_arr, color='black', alpha=0.2, linestyle='-')

    x_axis_label = r'$\dot{m}$'
    y_axis_label = r'$\xi$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)
    ax['a'].legend()

    prefix_folder = 'table/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'together_a'
    save.save_figure(fig, save_dir, file_name)

    x_axis_label = r'$\dot{m}$'
    y_axis_label = r'$L_{\rm x} \rm [erg ~ s^{-1}]$'
    ax1['a'].set_xlabel(x_axis_label, fontsize=26)
    ax1['a'].set_ylabel(y_axis_label, fontsize=26)
    ax1['a'].legend()

    prefix_folder = 'table/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'together_a_L_x'
    save.save_figure(fig1, save_dir, file_name)


if __name__ == '__main__':
    mu = 0.1e31
    beta_mu = 40
    mc2 = 100
    a_portion = 0.66
    phi_0 = 0
    theta_obs = 20

    theta_obs = 20
    beta_mu = 60
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.44
    phi_0 = 0
    # plot_L_to_mc2(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    theta_obs = 20
    beta_mu = 60
    mc2 = 30
    a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    phi_0 = 0
    # plot_L_to_a_portion(theta_obs, beta_mu, mc2, a_portion_arr, phi_0)

    theta_obs = 40
    beta_mu = 40
    mc2_arr = [30, 100]
    a_portion_arr = [0.22, 0.66]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    # plot_masses_PF_L_nu(theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)
    # plot_masses_PF_L(theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)

    beta_mu = 80
    mc2 = 30
    a_portion = 1
    phi_0 = 0
    # plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    theta_obs = 40
    beta_mu = 60
    mc2 = 30
    a_portion = 0.66
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    # plot_L_to_phi_0(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr, True)

    mc2 = 30
    a_portion = 1
    phi_0 = 0
    # plot_PF_contour(mu, mc2, a_portion, phi_0)

    mc2 = 100
    a_portion = 1
    phi_0 = 0
    # plot_PF_to_chi_theta(mu, mc2, a_portion, phi_0)

    theta_obs = 40
    beta_mu = 60
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.44
    phi_0 = 0
    # plot_L_iso_to_m(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    theta_obs = 40
    beta_mu = 60
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.66
    phi_0 = 0
    # plot_table(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)
    a_portion_arr = [0.22, 0.44, 0.66, 1]
    # plot_table_together(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0)
