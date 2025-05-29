from typing import Iterable

import numpy as np
import scipy.interpolate
from scipy import stats

import matplotlib.pyplot as plt
import scienceplots
import matplotlib.colors as colors
import matplotlib.lines as mlines
import matplotlib as mpl
import matplotlib.cm as cm
import seaborn as sns

import pathService
import config
import newService
import save
from geometry import matrix

plt.style.use(['science', 'notebook', 'grid'])
# чтобы стиль был похож на теховский
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'STIXGeneral'

# стандартный размер тиков по осям
ticks_labelsize = 18
mpl.rcParams['xtick.labelsize'] = ticks_labelsize
mpl.rcParams['ytick.labelsize'] = ticks_labelsize

cmap = mpl.cm.magma
# plasma
# magma
# inferno
# 'jet'
# viridis
# twilight_shifted
mpl.rcParams['image.cmap'] = 'magma'


# plt.rcParams['axes.labelsize'] = 42
# ticks_labelsize = 18


def plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0, flag_contour=True):
    '''рисует карты неба (излучательный коэффициент в зависимости от положения наблюдателя)
    по другому - куча профилей. для разных наблюдателей. характеризует как источник излучает в пространство
    '''
    ticks_labelsize = 16
    # try_sky_map(obs_i_angle_arr)
    # obs_i_angle = np.linspace(0, 180, 19)
    L_x = save.load_L_x(mu, 10, beta_mu, mc2, a_portion, phi_0)

    theta_obs_arr = np.linspace(0, 90, 10).astype(int)
    data_array = np.empty((len(theta_obs_arr) * 2 - 1, config.N_phase))
    # roll -- циклическая перестановка - делаем так как симметрическая задача и для угла 180 - theta будет симметрично
    # для сдвига на полфазы -- можно расчитать только до 90 а потом для других переставить и получить для до 180
    for i, theta_obs in enumerate(theta_obs_arr):
        L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
        data_array[i] = L_total
        # для углов больших 90 будет симметрично относительно 90 (со сдвигом на полфазы)
        if i != len(theta_obs_arr) - 1:  # wtf -1 ?? - не будем трогать 90, его не дублируем
            data_array[-i - 1] = np.roll(data_array[i], len(data_array[i]) // 2)

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_array / L_x)
    theta_obs_arr_to_plot = np.linspace(0, 180, 19).astype(int)

    fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))
    im = ax['a'].pcolormesh(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='white', marker='*', s=300)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='white', marker='*', s=300)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='black', marker='*', s=100)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='black', marker='*', s=100)

    x_axis_label = config.symbol_phase
    y_axis_label = config.symbol_theta_obs_y
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    clb = plt.colorbar(im)
    clb.set_label(r'$L_{\rm iso} / L_{\rm x}$', fontsize=24)  # r'$L_{\rm iso} \cdot L_{\rm x}^{-1}$'
    clb.ax.tick_params(labelsize=ticks_labelsize)

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

    x_axis_label = config.symbol_phase
    y_axis_label = config.symbol_theta_obs_y
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    clb = plt.colorbar(im)
    clb.set_label(r'$L_{\rm iso} / L_{\rm x}$', fontsize=24)
    clb.ax.tick_params(labelsize=ticks_labelsize)

    prefix_folder = 'sky_map/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion, phi_0=phi_0,
                                   prefix_folder=prefix_folder)
    file_name = 'try_map_contour'
    save.save_figure(fig, save_dir, file_name)

    if flag_contour:
        def get_r(i, j):

            x = matrix.get_cartesian_from_spherical(1, np.deg2rad(theta_obs_arr_to_plot[i]),
                                                    config.phase_for_plot[j] * 2 * np.pi)
            x = np.array(x)

            y = matrix.get_cartesian_from_spherical(1, np.deg2rad(beta_mu), 0)
            y = np.array(y)

            y1 = matrix.get_cartesian_from_spherical(1, np.pi - np.deg2rad(beta_mu), np.pi)
            y1 = np.array(y1)

            z = np.arccos(x.dot(y) / np.linalg.norm(x) * np.linalg.norm(y))
            z1 = np.arccos(x.dot(y1) / np.linalg.norm(x) * np.linalg.norm(y))

            # print(x.dot(y))
            # print(np.arccos(x @ y / np.linalg.norm(x) * np.linalg.norm(y)))
            return np.min([z, z1])

        fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))
        im = ax['a'].contourf(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot, levels=30)
        ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='white', marker='*', s=300)
        ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='white', marker='*', s=300)
        ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='black', marker='*', s=100)
        ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='black', marker='*', s=100)

        x_axis_label = config.symbol_phase
        y_axis_label = config.symbol_theta_obs_y
        ax['a'].set_xlabel(x_axis_label, fontsize=24)
        ax['a'].set_ylabel(y_axis_label, fontsize=24)
        ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

        clb = plt.colorbar(im)
        clb.set_label(r'$L_{\rm iso} / L_{\rm x}$', fontsize=24)
        clb.ax.tick_params(labelsize=ticks_labelsize)

        x = config.phase_for_plot
        y = theta_obs_arr_to_plot
        X, Y = np.meshgrid(x, y)
        for i in range(data_array.shape[0]):
            for j in range(data_array.shape[1]):
                data_array[i, j] = np.rad2deg(get_r(i, j))

        Z = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_array)

        CS = ax['a'].contour(X, Y, Z, levels=10, colors='c')
        ax['a'].clabel(CS, fontsize=16)

        prefix_folder = 'sky_map/'
        save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion,
                                       phi_0=phi_0,
                                       prefix_folder=prefix_folder)
        file_name = 'try_map_contour_with'
        save.save_figure(fig, save_dir, file_name)


def plot_rvm(mu, beta_mu, mc2):
    '''рисует карты неба (излучательный коэффициент в зависимости от положения наблюдателя)
    по другому - куча профилей. для разных наблюдателей. характеризует как источник излучает в пространство
    '''
    ticks_labelsize = 16
    # try_sky_map(obs_i_angle_arr)
    # obs_i_angle = np.linspace(0, 180, 19)
    # L_x = save.load_L_x(mu, 10, beta_mu, mc2, a_portion, phi_0)

    theta_obs_arr = np.linspace(0, 90, 10).astype(int)
    data_array = np.empty((len(theta_obs_arr) * 2 - 1, config.N_phase))

    theta_obs_arr_to_plot = np.linspace(0, 180, 19).astype(int)

    first_pole = (beta_mu, 0)
    second_pole = (180 - beta_mu, np.pi)
    # second_pole = (180 - beta_mu + np.deg2rad(10), np.pi - np.deg2rad(140))
    third_pole = (beta_mu, 2 * np.pi)

    def get_r(i, j, pole):
        # return np.sqrt(2 * np.abs((np.deg2rad(theta_obs_arr_to_plot[i] - pole[0]))) ** 2 + (1 / 2 * np.abs(
        #     (config.phase_for_plot[j] * 2 * np.pi - pole[1]))) ** 2)
        x = matrix.get_cartesian_from_spherical(1, np.deg2rad(theta_obs_arr_to_plot[i]),
                                                config.phase_for_plot[j] * 2 * np.pi)
        x = np.array(x)

        y = matrix.get_cartesian_from_spherical(1, np.deg2rad(beta_mu), 0)
        y = np.array(y)

        y1 = matrix.get_cartesian_from_spherical(1, np.pi - np.deg2rad(beta_mu), np.pi)
        y1 = np.array(y1)

        z = np.arccos(x.dot(y) / np.linalg.norm(x) * np.linalg.norm(y))
        z1 = np.arccos(x.dot(y1) / np.linalg.norm(x) * np.linalg.norm(y))

        # print(x.dot(y))
        # print(np.arccos(x @ y / np.linalg.norm(x) * np.linalg.norm(y)))
        return np.min([z, z1])

    for i in range(data_array.shape[0]):
        for j in range(data_array.shape[1]):
            r1 = get_r(i, j, first_pole)
            r2 = get_r(i, j, second_pole)
            r3 = get_r(i, j, third_pole)
            data_array[i, j] = -np.min([r1, r2, r3])  # -

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_array)

    fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))
    im = ax['a'].contourf(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot, levels=30)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='white', marker='*', s=300)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='white', marker='*', s=300)
    ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='black', marker='*', s=100)
    ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='black', marker='*', s=100)

    x_axis_label = config.symbol_phase
    y_axis_label = config.symbol_theta_obs_y
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    clb = plt.colorbar(im)
    clb.set_label(r'$L_{\rm iso} / L_{\rm x}$', fontsize=24)
    clb.ax.tick_params(labelsize=ticks_labelsize)

    prefix_folder = 'sky_map_rvm/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=beta_mu, mc2=mc2, a_portion=None, phi_0=None,
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

    x_axis_label = config.symbol_phase
    # y_axis_label = r'$mean L_{iso} [erg/s]$'
    y_axis_label = config.symbol_m
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    clb = plt.colorbar(im, pad=0.03, format="{x:.2}")  # pad=0.15
    clb.set_label(r'$\widetilde{L}_{\rm iso}$', fontsize=24)  # \cdot max(L_{\rm iso})^{-1}

    prefix_folder = 'L_to_mass/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=a_portion,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
    file_name = 'map_contour'
    save.save_figure(fig, save_dir, file_name)


def plot_L_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0):
    '''рисует 2d L от фазы и a'''
    L_total_a_portion = np.empty((len(a_portion_arr), config.N_phase))
    for i, a_portion in enumerate(a_portion_arr):
        L_total_a_portion[i] = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    # нормировка на L_max
    buf = L_total_a_portion / np.max(L_total_a_portion, axis=-1).ravel()[:, np.newaxis]
    data_to_plot_norm = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=buf)
    # buf = L_total_a_portion
    # buf = L_total_a_portion / np.max(np.max(L_total_a_portion, axis=-1))
    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=np.log10(L_total_a_portion))

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    # чтобы сделать хороший масштаб и чтобы каждый профиль имел равную толщину сначала рисуем от 0 до 1 а потом
    # заменяем надписи на значения а
    levels = 10 ** np.linspace(np.log10(np.min(data_to_plot)), np.log10(np.max(data_to_plot)), 10)
    # norm - нужно сместить растянуть колорбар потому что в районе 0.9 есть перепады, но мы их не увидим из-за
    # малых значений, которые появляются после добавления а=1
    im = ax['a'].pcolormesh(config.phase_for_plot, np.linspace(0, 1, len(a_portion_arr)), data_to_plot,
                            # norm=colors.LogNorm(vmin=data_to_plot.min(), vmax=data_to_plot.max())
                            # norm=colors.TwoSlopeNorm(0.9, vmin=data_to_plot.min(), vmax=data_to_plot.max())
                            )
    ax['a'].set_yticks(np.linspace(0, 1, len(a_portion_arr)), np.round(a_portion_arr, 2))

    x_axis_label = config.symbol_phase
    y_axis_label = r'$a$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    clb = plt.colorbar(im, pad=0.01)
    # clb.set_label(r'$\widetilde{L}_{\rm iso}$', fontsize=24)
    clb.set_label(r'$\log {L}_{\rm iso}$', fontsize=24)

    prefix_folder = 'L_to_a/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=mc2, a_portion=None,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
    file_name = 'map_contour'
    save.save_figure(fig, save_dir, file_name)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    # чтобы сделать хороший масштаб и чтобы каждый профиль имел равную толщину сначала рисуем от 0 до 1 а потом
    # заменяем надписи на значения а
    # norm - нужно сместить растянуть колорбар потому что в районе 0.9 есть перепады, но мы их не увидим из-за
    # малых значений, которые появляются после добавления а=1
    im = ax['a'].pcolormesh(config.phase_for_plot, np.linspace(0, 1, len(a_portion_arr)), data_to_plot_norm,
                            # norm=colors.LogNorm(vmin=data_to_plot.min(), vmax=data_to_plot.max())
                            norm=colors.TwoSlopeNorm(0.9, vmin=data_to_plot_norm.min(), vmax=data_to_plot_norm.max())
                            )
    ax['a'].set_yticks(np.linspace(0, 1, len(a_portion_arr)), np.round(a_portion_arr, 2))

    x_axis_label = config.symbol_phase
    y_axis_label = r'$a$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$\widetilde{L}_{\rm iso}$', fontsize=24)
    # clb.set_label(r'${L}_{\rm iso}$', fontsize=24)

    prefix_folder = 'L_to_a/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=mc2, a_portion=None,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
    file_name = 'map_contour_norm'
    save.save_figure(fig, save_dir, file_name)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    ax['a'].plot(a_portion_arr, L_total_a_portion.mean(axis=1), color='black')

    x_axis_label = r'$a$'
    y_axis_label = r'${L}_{\rm iso}$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    prefix_folder = 'L_to_a/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=mc2, a_portion=None,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
    file_name = 'L_to_a'
    save.save_figure(fig, save_dir, file_name)


def plot_masses_PF_L_nu(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr, energy_index=8):
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

    global cmap

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

            # в словаере лежит среднее nu_L_nu : (PF, phi_0). сортируем его по значению phi_0,
            # чтобы посмотреть как зависит от фи0 - будем рисовать по возрастанию фи0
            list_tuples = sorted(full_dict.items(), key=lambda item: item[1][1])
            x_sort_phi_0, y_sort_phi_0 = zip(*list_tuples)
            y_sort_phi_0, colors_sort_phi_0 = zip(*y_sort_phi_0)

            colors = (np.array(colors_sort_phi_0)) / np.max(colors_sort_phi_0)

            if marker_index == 0:
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, color=cmap(colors))
                # ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=30, facecolors='none', edgecolors=cm.jet(colors))
            else:
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, marker=marker_dict[marker_index % 4],
                                color=cmap(colors))

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
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax['a'], orientation='vertical', pad=0.01)
    cb.set_label(label=config.symbol_phi_0_y, fontsize=24)

    # config.symbol_phi_0 + ' ' + config.symbol_grad

    prefix_folder = 'PF_to_L_nu/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)

    file_name = f'theta={theta_obs} beta_mu={beta_mu} All_PF_to_L_nu.png'
    save.save_figure(fig, save_dir, file_name)


def plot_masses_PF_L(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr):
    '''то же самое что plot_masses_PF_L_nu только рисуем от L_iso'''
    marker_index = 0
    line_style = ['-', '--']
    marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}

    # cmap = mpl.cm.viridis
    global cmap

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
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, color=cmap(colors))
                # ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, facecolors='none', edgecolors=cm.jet(colors))
            else:
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, marker=marker_dict[marker_index % 4],
                                color=cmap(colors))

            ax['a'].plot(x_sort_phi_0, y_sort_phi_0, color='black', alpha=0.2, linestyle=line_style[mc2_index])

        marker_index = 3

    x_axis_label = r'$L_{\rm iso}$' + r'$\rm [erg/s]$'
    y_axis_label = r'$PF$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].set_xscale('log')
    # ax['a'].set_ylim(0, 1.1)

    bounds = phi_0_arr.copy()
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax['a'], orientation='vertical', pad=0.01)
    cb.set_label(label=config.symbol_phi_0_y, fontsize=24)

    prefix_folder = 'PF_to_L/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)

    file_name = f'theta={theta_obs} beta_mu={beta_mu} All_PF_to_L.png'
    save.save_figure(fig, save_dir, file_name)


def plot_masses_PF_L_on_off(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr):
    '''то же самое что plot_masses_PF_L_nu только рисуем от L_iso'''
    marker_index = 0
    line_style = ['-', '--']
    marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}

    # cmap = mpl.cm.viridis
    # global cmap

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
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, color=cmap(colors))
                # ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, facecolors='none', edgecolors=cm.jet(colors))
            else:
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, marker=marker_dict[marker_index % 4],
                                color=cmap(colors))

            ax['a'].plot(x_sort_phi_0, y_sort_phi_0, color='black', alpha=0.2, linestyle=line_style[mc2_index])

        marker_index = 3

    x_axis_label = r'$L_{\rm iso}$' + r'$\rm [erg/s]$'
    y_axis_label = r'$PF$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].set_xscale('log')
    # ax['a'].set_ylim(0, 1.1)

    bounds = phi_0_arr.copy()
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax['a'], orientation='vertical', pad=0.01)
    cb.set_label(label=config.symbol_phi_0_y, fontsize=24)

    limits_x = ax['a'].get_xlim()
    limits_y = ax['a'].get_ylim()

    config.flag_scatter = False
    marker_index = 0

    fig1, ax1 = plt.subplot_mosaic('a', figsize=(12, 6))
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
                ax1['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, color=cmap(colors))
                # ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, facecolors='none', edgecolors=cm.jet(colors))
            else:
                ax1['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, marker=marker_dict[marker_index % 4],
                                 color=cmap(colors))

            ax1['a'].plot(x_sort_phi_0, y_sort_phi_0, color='black', alpha=0.2, linestyle=line_style[mc2_index])

        marker_index = 3

    x_axis_label = r'$L_{\rm iso}$' + r'$\rm [erg/s]$'
    y_axis_label = r'$PF$'
    ax1['a'].set_xlabel(x_axis_label, fontsize=24)
    ax1['a'].set_ylabel(y_axis_label, fontsize=24)

    ax1['a'].set_xscale('log')

    limits_x1 = ax1['a'].get_xlim()
    limits_y1 = ax1['a'].get_ylim()

    limits_xall = (min(limits_x[0], limits_x1[0]), max(limits_x[1], limits_x1[1]))
    limits_yall = (min(limits_y[0], limits_y1[0]), max(limits_y[1], limits_y1[1]))

    ax['a'].set_xlim(limits_xall)
    ax['a'].set_ylim(limits_yall)

    ax1['a'].set_xlim(limits_xall)
    ax1['a'].set_ylim(limits_yall)

    bounds = phi_0_arr.copy()
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax1['a'], orientation='vertical', pad=0.01)
    cb.set_label(label=config.symbol_phi_0_y, fontsize=24)

    prefix_folder = 'PF_to_L/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)

    file_name = f'theta={theta_obs} beta_mu={beta_mu} All_PF_to_L_on'
    save.save_figure(fig, save_dir, file_name)

    prefix_folder = 'PF_to_L/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)

    file_name = f'theta={theta_obs} beta_mu={beta_mu} All_PF_to_L_off'
    save.save_figure(fig1, save_dir, file_name)


def plot_PF_to_L_const_a_dif_m(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0_arr):
    '''то же самое что plot_masses_PF_L_nu только рисуем от L_iso'''
    marker_index = 0
    line_style = ['--']
    marker_dict = {0: '.', 1: '*', 2: '+', 3: '^', 4: "o", 5: "d", 6: "s", 7: "X"}

    global cmap
    legend_arr = [0] * len(mc2_arr)

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
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

            ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, marker=marker_dict[marker_index], color=cmap(colors))
            ax['a'].plot(x_sort_phi_0, y_sort_phi_0, color='black', alpha=0.2, linestyle=line_style[0])

        legend_arr[mc2_index] = mlines.Line2D([], [], color=cmap(colors)[2], marker=marker_dict[marker_index],
                                              markersize=15, linestyle='--', label=f'm={mc2}', alpha=0.2)
        marker_index += 1

    x_axis_label = r'$L_{\rm iso}$' + r'$\rm [erg/s]$'
    y_axis_label = r'$PF$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].set_xscale('log')

    plt.legend(handles=legend_arr)
    bounds = phi_0_arr.copy()
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax['a'], orientation='vertical', pad=0.01)
    cb.set_label(label=config.symbol_phi_0_y, fontsize=24)

    prefix_folder = 'PF_to_L_const_a/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=a_portion,
                                   phi_0=None, prefix_folder=prefix_folder)

    file_name = f'theta={theta_obs} beta_mu={beta_mu} All_PF_to_L.png'
    save.save_figure(fig, save_dir, file_name)


def plot_PF_to_L_const_a_const_phi0_dif_m(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0):
    '''то же самое что plot_masses_PF_L_nu только рисуем от L_iso'''
    marker_index = 0
    line_style = ['--']
    marker_dict = {0: '.', 1: '*', 2: '+', 3: '^', 4: "o", 5: "d", 6: "s", 7: "X"}

    # cmap = mpl.cm.plasma
    global cmap
    legend_arr = [0] * len(a_portion_arr)

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    for a_portion_index, a_portion in enumerate(a_portion_arr):
        full_dict = {}
        for mc2_index, mc2 in enumerate(mc2_arr):
            L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            PF = newService.get_PF(L_total)
            full_dict[np.mean(L_total)] = (PF, mc2)

        # в словаере лежит среднее : (PF, mc2). сортируем его по значению mc2
        list_tuples = sorted(full_dict.items(), key=lambda item: item[1][1])
        x_sort_mc2, y_sort_mc2 = zip(*list_tuples)
        y_sort_mc2, colors_sort_mc2 = zip(*y_sort_mc2)

        colors = (np.array(colors_sort_mc2)) / np.max(colors_sort_mc2)

        ax['a'].scatter(x_sort_mc2, y_sort_mc2, s=100, marker=marker_dict[a_portion_index], color=cmap(colors))
        ax['a'].plot(x_sort_mc2, y_sort_mc2, color='black', alpha=0.2, linestyle=line_style[0])

        legend_arr[a_portion_index] = mlines.Line2D([], [], color=cmap(colors)[1], marker=marker_dict[a_portion_index],
                                                    markersize=15, linestyle='--', label=f'a={a_portion}')

    x_axis_label = r'$L_{\rm iso}$' + r'$\rm [erg/s]$'
    y_axis_label = r'$PF$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].set_xscale('log')

    plt.legend(handles=legend_arr)
    bounds = mc2_arr.copy()
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax['a'], orientation='vertical', pad=0.01)
    cb.set_label(label=config.symbol_m, fontsize=24)

    prefix_folder = 'PF_to_L_const_a_diff_m/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)

    file_name = f'theta={theta_obs} beta_mu={beta_mu} All_PF_to_L.png'
    save.save_figure(fig, save_dir, file_name)


def plot_L_to_phi_0(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr, flag_same_res=False):
    '''2d L от Phi phi0'''
    # plot_L_to_new_fi_0
    # ожидаем что phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]

    L_data = np.empty(((len(phi_0_arr)) * 2 - 1, config.N_phase))
    for i, phi_0 in enumerate(phi_0_arr):
        L_data[i] = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    for i in range(9):
        # симметрично для 180 - fi_0, как будто вращаемся в другую сторону, поэтому сортировка ::-1
        L_data[9 + i + 1] = L_data[9 - i - 1][::-1]

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=L_data)
    phi_0_arr_for_plot = [20 * i for i in range(19)]

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    im = ax['a'].pcolormesh(config.phase_for_plot, phi_0_arr_for_plot, data_to_plot)
    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$L_{\rm iso}$' + r'$\rm [erg/s]$', fontsize=26)

    x_axis_label = config.symbol_phase
    y_axis_label = config.symbol_phi_0_y
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

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    im = ax['a'].contourf(np.linspace(0, 2, len(data_to_plot[0])), phi_0_arr_for_plot, data_to_plot, levels=30)
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
    '''2d PF от theta beta_mu'''
    # PF(L_nu) много точек, берутся линии по mc, a
    ticks_labelsize = 18

    theta_obs_arr = np.linspace(0, 90, 10).astype(int)
    beta_mu_arr = np.linspace(0, 80, 9).astype(int)

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

    x_axis_label = r'$\chi [^\circ]$'
    y_axis_label = config.symbol_theta_obs_y
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    cbar.ax.tick_params(labelsize=ticks_labelsize)
    ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

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

    cbar.ax.tick_params(labelsize=ticks_labelsize)
    ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    prefix_folder = 'PF_contour/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=mc2, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'a={a_portion}' + ' ' + f'phi_0={phi_0} colormesh'
    save.save_figure(fig, save_dir, file_name)


def plot_PF_to_chi_theta(mu, mc2, a_portion, phi_0):
    '''ищем закономерности в PF_theta_beta . смотрим как зависит от theta+beta и theta-beta'''
    theta_obs_arr = np.linspace(0, 90, 10).astype(int)
    beta_mu_arr = np.linspace(0, 80, 9).astype(int)

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

    x_axis_label = config.symbol_mu_ang + ' + ' + config.symbol_theta_obs_y  # r'$\chi + \theta_{\rm obs} [^\circ]$'
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

    x_axis_label = config.symbol_mu_ang + ' - ' + config.symbol_theta_obs_y  # r'$\chi - \theta_{\rm obs} [^\circ]$'
    y_axis_label = r'$PF$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'PF_to_minus/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'a={a_portion}' + ' ' + f'mc2={mc2}' + ' ' + f'phi_0={phi_0}'
    save.save_figure(fig, save_dir, file_name)


def plot_L_iso_to_m(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0):
    '''L от m для таблицы'''
    L_iso_arr = np.empty(len(mc2_arr))
    for i, mc2 in enumerate(mc2_arr):
        L_iso_arr[i] = np.mean(save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0), axis=-1)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    ax['a'].scatter(mc2_arr, L_iso_arr, s=30, facecolors='none', edgecolors='black')
    ax['a'].plot(mc2_arr, L_iso_arr, color='black', alpha=0.2, linestyle='-')

    x_axis_label = config.symbol_m
    y_axis_label = r'$L_{\rm iso}$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'L_iso_to_m/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'a={a_portion} phi_0={phi_0}'
    save.save_figure(fig, save_dir, file_name)


def plot_table(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0):
    '''R_e/R_ns от m + ksi_shock от m'''
    R_e_arr = np.empty(len(mc2_arr))
    ksi_shock_arr = np.empty(len(mc2_arr))
    L_x_arr = np.empty(len(mc2_arr))

    for i, mc2 in enumerate(mc2_arr):
        R_e_arr[i] = save.load_R_e(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
        ksi_shock_arr[i] = save.load_ksi_shock(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
        L_x_arr[i] = save.load_L_x(mu, theta_obs, beta_mu, mc2, a_portion, phi_0) * 1e-38

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))

    ax['a'].scatter(mc2_arr, R_e_arr, s=30, facecolors='none', edgecolors='red',
                    label=r'$\xi_{\rm e}$')  # \frac{R_e}{R_*}
    ax['a'].plot(mc2_arr, R_e_arr, color='black', alpha=0.2, linestyle='-')

    ax['a'].scatter(mc2_arr, ksi_shock_arr, s=30, facecolors='none', edgecolors='blue', label=r'$\xi_{\rm s}$')
    ax['a'].plot(mc2_arr, ksi_shock_arr, color='black', alpha=0.2, linestyle='-')

    ax['a'].scatter(mc2_arr, L_x_arr, s=30, facecolors='none', edgecolors='green', label=r'$\frac{L_{\rm x}}{10^{38}}$')
    ax['a'].plot(mc2_arr, L_x_arr, color='black', alpha=0.2, linestyle='-')

    x_axis_label = config.symbol_m
    y_axis_label = r'$R/R_{*}$'  # \xi
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)
    ax['a'].legend()

    prefix_folder = 'table/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'a={a_portion}'
    save.save_figure(fig, save_dir, file_name)


def plot_table_together(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0):
    '''R_e/R_ns от m + ksi_shock от m + L_x от m'''
    ticks_labelsize = 18
    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    fig1, ax1 = plt.subplot_mosaic('a', figsize=(9, 6))
    fig2, ax2 = plt.subplot_mosaic('a', figsize=(9, 6))

    color_arr = ['red', 'green', 'blue', 'purple']

    M_accretion_rate = np.array(mc2_arr) * config.L_edd / config.c ** 2
    R_alfven = (mu ** 2 / (2 * M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
    R_disk = config.ksi_param * R_alfven / config.R_ns

    for i, a_portion in enumerate(a_portion_arr):
        R_e_arr = np.empty(len(mc2_arr))
        delta_arr = np.empty(len(mc2_arr))
        ksi_shock_arr = np.empty(len(mc2_arr))
        L_x_arr = np.empty(len(mc2_arr))
        for j, mc2 in enumerate(mc2_arr):
            R_e_arr[j] = save.load_R_e(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            delta_arr[j] = newService.get_delta_distance_at_surface_NS(R_e_arr[j] * config.R_ns,
                                                                       config.dRdisk_div_Rdisk) / config.R_ns
            ksi_shock_arr[j] = save.load_ksi_shock(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            L_x_arr[j] = save.load_L_x(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

        ax['a'].scatter(mc2_arr, ksi_shock_arr, s=30, facecolors='none', edgecolors=color_arr[i],
                        label=r'$\xi_{s}$' + ' for ' + f'a={a_portion}')  # r'$\xi_{s}$' + ' for ' +
        ax['a'].plot(mc2_arr, ksi_shock_arr, color='black', alpha=0.2, linestyle='-')

        ax1['a'].scatter(mc2_arr, L_x_arr, s=30, facecolors='none', edgecolors=color_arr[i],
                         label=f'a={a_portion}')
        ax1['a'].plot(mc2_arr, L_x_arr, color='black', alpha=0.2, linestyle='-')

        ax2['a'].scatter(mc2_arr, delta_arr, s=30, facecolors='none', edgecolors=color_arr[i],
                         label=f'a={a_portion}')
        ax2['a'].plot(mc2_arr, delta_arr, color='black', alpha=0.2, linestyle='-')

    # ax['a'].scatter(mc2_arr, R_e_arr, s=30, facecolors='none', edgecolors='black', label=r'$\xi_{\rm e}$')
    # ax['a'].plot(mc2_arr, R_e_arr, color='black', alpha=0.2, linestyle='-')

    ax['a'].scatter(mc2_arr, R_disk, s=30, facecolors='none', edgecolors='black', label=r'$\xi_{\rm disk}$')
    ax['a'].plot(mc2_arr, R_disk, color='black', alpha=0.2, linestyle='-')

    ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    x_axis_label = config.symbol_m
    y_axis_label = r'$\xi$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)
    ax['a'].legend()

    prefix_folder = 'table/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'together_a'
    save.save_figure(fig, save_dir, file_name)

    x_axis_label = config.symbol_m
    y_axis_label = r'$L_{\rm x} \rm [erg ~ s^{-1}]$'
    ax1['a'].set_xlabel(x_axis_label, fontsize=26)
    ax1['a'].set_ylabel(y_axis_label, fontsize=26)
    ax1['a'].legend()

    ax1['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    prefix_folder = 'table/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'together_a_L_x'
    save.save_figure(fig1, save_dir, file_name)

    x_axis_label = config.symbol_m
    y_axis_label = r'$\delta$'
    ax2['a'].set_xlabel(x_axis_label, fontsize=26)
    ax2['a'].set_ylabel(y_axis_label, fontsize=26)
    ax2['a'].legend()

    ax2['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    prefix_folder = 'table/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'together_delta'
    save.save_figure(fig2, save_dir, file_name)


def plot_table_R_disk(mu, theta_obs, beta_mu_arr, a_portion, mc2_arr, phi_0):
    '''R_e/R_ns от m'''
    ticks_labelsize = 18
    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))

    for i, beta_mu in enumerate(beta_mu_arr):
        R_e_arr = np.empty(len(mc2_arr))
        for j, mc2 in enumerate(mc2_arr):
            R_e_arr[j] = save.load_R_e(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

        ax['a'].scatter(mc2_arr, R_e_arr, s=30, label=r'$\chi = $' + f'{beta_mu}' + r'$^\circ$')
        ax['a'].plot(mc2_arr, R_e_arr, color='black', alpha=0.2, linestyle='-')

    # ax['a'].scatter(mc2_arr, R_disk, s=30, color='black', label=r'$\frac{R_{disk}}{R_*}$')
    # ax['a'].plot(mc2_arr, R_disk, color='black', alpha=0.2, linestyle='-')

    ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)
    ax['a'].set_yscale('log')

    x_axis_label = config.symbol_m
    y_axis_label = r'$R/R_{*}$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)
    ax['a'].legend()

    prefix_folder = 'table/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'R_e'
    save.save_figure(fig, save_dir, file_name)


def plot_PF_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0):
    '''по чему усреднять? по каким параметрам'''
    L_total_a_portion = np.empty((len(a_portion_arr), config.N_phase))
    for i, a_portion in enumerate(a_portion_arr):
        L_total_a_portion[i] = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    PF = np.apply_along_axis(newService.get_PF, axis=-1, arr=L_total_a_portion)
    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    ax['a'].plot(a_portion_arr, PF, color='black')

    x_axis_label = r'$a$'
    y_axis_label = r'$PF$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    plt.show()


def plot_ksi_to_beta(mu, theta_obs, beta_mu_arr, mc2, a_portion, phi_0):
    ksi_arr = np.empty(len(beta_mu_arr))
    for i, beta_mu in enumerate(beta_mu_arr):
        ksi_arr[i] = save.load_ksi_shock(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    ax['a'].plot(beta_mu_arr, ksi_arr, color='black')

    x_axis_label = r'$\chi [^\circ]$'
    y_axis_label = r'$\xi_{\rm s}$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'ksi_to_beta/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=None, mc2=mc2, a_portion=a_portion,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
    file_name = f'ksi_to_beta'
    save.save_figure(fig, save_dir, file_name)


def plot_L_max_phase_to_m_to_a(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0):
    '''
    max phase position.
    контуры равной фазы максимума
    :return:
    '''

    # mc2_arr = [10, 30, 60, 80, 100, 130, 160]
    # a_portion_arr = [0.11, 0.22, 0.33, 0.44, 0.55, 0.66]

    phase = np.linspace(0, 1, config.N_phase)
    max_phase_idx_data = np.zeros((len(a_portion_arr), len(mc2_arr)))

    for i, a_portion in enumerate(a_portion_arr):
        for j, mc2 in enumerate(mc2_arr):
            L = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            max_idx = np.argmax(L)
            max_phase_idx_data[i, j] = phase[max_idx]

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))

    # im = ax['a'].contourf(mc2_arr, a_portion_arr, max_phase_idx_data)
    im = ax['a'].contourf(np.linspace(0, 1, len(mc2_arr)), np.linspace(0, 1, len(a_portion_arr)), max_phase_idx_data)
    ax['a'].set_xticks(np.linspace(0, 1, len(mc2_arr)), mc2_arr)
    ax['a'].set_yticks(np.linspace(0, 1, len(a_portion_arr)), a_portion_arr)

    # ax.contour(mc2_arr, a_portion_arr, max_phase_idx_data)
    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$\Phi_{max}$', fontsize=26)

    x_axis_label = config.symbol_m
    y_axis_label = r'$a$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'L_max_phase/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
    file_name = f'contour_max'
    save.save_figure(fig, save_dir, file_name)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    # im = ax['a'].pcolormesh(mc2_arr, a_portion_arr, max_phase_idx_data)
    im = ax['a'].pcolormesh(np.linspace(0, 1, len(mc2_arr)), np.linspace(0, 1, len(a_portion_arr)), max_phase_idx_data)
    ax['a'].set_xticks(np.linspace(0, 1, len(mc2_arr)), mc2_arr)
    ax['a'].set_yticks(np.linspace(0, 1, len(a_portion_arr)), a_portion_arr)
    # ax.contour(mc2_arr, a_portion_arr, max_phase_idx_data, colors='k')
    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$\Phi_{max}$', fontsize=26)

    x_axis_label = config.symbol_m
    y_axis_label = r'$a$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    prefix_folder = 'L_max_phase/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
    file_name = f'colormesh_max'
    save.save_figure(fig, save_dir, file_name)

    newpoints = 20
    xq, yq = np.linspace(min(mc2_arr), max(mc2_arr), newpoints), \
        np.linspace(min(a_portion_arr), max(a_portion_arr), newpoints)
    f = scipy.interpolate.interp2d(mc2_arr, a_portion_arr, max_phase_idx_data, kind='cubic')
    interpolate_data = f(xq, yq)

    fig, ax = plt.subplot_mosaic('a', figsize=(9, 6))
    im = ax['a'].contourf(xq, yq, interpolate_data)
    # im = ax['a'].contourf(np.linspace(0, 1, len(xq)), np.linspace(0, 1, len(yq)), interpolate_data)
    # ax['a'].set_xticks(np.linspace(0, 1, len(xq)), xq)
    # ax['a'].set_yticks(np.linspace(0, 1, len(yq)), yq)

    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$\Phi_{max}$', fontsize=26)

    x_axis_label = config.symbol_m
    y_axis_label = r'$a$'
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)
    # fig.suptitle(figure_title, fontsize=22)

    prefix_folder = 'L_max_phase/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=phi_0, prefix_folder=prefix_folder)
    file_name = f'contour_max_interpolate'
    save.save_figure(fig, save_dir, file_name)


def plot_a_restrictions():
    '''
        отрисовка пределов для а
    '''

    eps = 1e-5
    chi_arr = np.linspace(eps, np.pi / 2, 1600, endpoint=False)
    # cos >
    # ans = 1 / ((1 + delta) * np.sin(chi_arr) ** 2) - 1 / np.tan(chi_arr) ** 2

    theta_end = np.pi / 2 - chi_arr
    delta = 0.25
    delta_arr = [0.05, 0.25, 0.5, 1, 1.5]

    colors = ['black', 'blue', 'red', 'green', 'magenta', 'cyan']

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    for delta_ind, delta in enumerate(delta_arr[::-1]):

        ans = (1 / np.tan(chi_arr) ** 2 * (1 - (1 + delta) * np.sin(theta_end) ** 2) / (
                (1 + delta) * np.sin(theta_end) ** 2))

        chi_05 = np.arccos((1 / (1 + delta + eps)) ** (1 / 2))
        n_arr1 = 100
        chi_arr1 = np.linspace(eps, chi_05, n_arr1, endpoint=False)
        n_arr2 = 400
        chi_arr2 = np.linspace(chi_05 + eps, np.pi / 2, n_arr2, endpoint=False)
        chi_arr = np.concatenate((chi_arr1, chi_arr2))

        # ax['a'].axvline(x=np.rad2deg(chi_05), color='black', linestyle='--')

        theta_end = np.pi / 2 - chi_arr
        ans = (1 / np.tan(chi_arr) ** 2 * (1 - (1 + delta) * np.sin(theta_end) ** 2) / (
                (1 + delta) * np.sin(theta_end) ** 2))

        res = np.zeros(ans.shape[0])
        for i in range(ans.shape[0]):
            if ans[i] < 0:
                res[i] = 1
            elif ans[i] > 1:
                res[i] = 0
            else:
                res[i] = np.arccos(ans[i] ** (1 / 2)) / np.pi

        ax['a'].plot(np.rad2deg(chi_arr1), res[:n_arr1], label=r'$\Delta$' + f' = {delta}', color=colors[delta_ind])
        ax['a'].plot(np.rad2deg(chi_arr2), res[n_arr1:], color=colors[delta_ind])
        # ax['a'].scatter(np.rad2deg(chi_05), 0.5, color=colors[delta_ind])
        ax['a'].plot(np.rad2deg([chi_arr1[-1], chi_05 + eps]), [1, 0.5], color=colors[delta_ind], linestyle='--')

        # ax['a'].scatter(np.rad2deg(chi_arr), res, label=r'$\Delta$' + f' = {delta}')
        ax['a'].axhline(y=0.5, color='black', linestyle='--', alpha=0.1)

        x_axis_label = r'$\chi ~ [^{\circ}]$'
        y_axis_label = r'$a$'
        ax['a'].set_xlabel(x_axis_label, fontsize=24)
        ax['a'].set_ylabel(y_axis_label, fontsize=24)
        ax['a'].legend()
        ax['a'].set_ylim(0, 1.1)

    plt.show()


def plot_hist_tau(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    tensor_tau_cols = save.load_tensor(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, 'tensor_tau_cols')
    buf_tau = tensor_tau_cols.reshape(1, -1)
    buf_tau = buf_tau[buf_tau > 0]

    tensor_alpha_cols = save.load_tensor(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, 'tensor_alpha_cols')
    buf_cos = tensor_alpha_cols.reshape(1, -1)
    buf_cos = buf_cos[buf_cos > 0]

    fig1, ax1 = plt.subplot_mosaic('ab', figsize=(18, 6))
    sns.ecdfplot(np.log10(buf_tau / buf_cos), ax=ax1['a'])
    ax1['a'].set_title(r'$\log \tau/\cos \alpha_0$', fontsize=26)

    # ax['a'].hist(buf, 30, histtype='bar', rwidth=0.8)  # barstacked
    # ax['a'].set_title(r'$\tau$', fontsize=26)

    # sns.ecdfplot(buf, label='Эмпирическая CDF', ax=ax['c'])

    # kde = stats.gaussian_kde(tensor_tau_cols.reshape(1, -1))
    # x_kde = np.linspace(tensor_tau_cols.min(), tensor_tau_cols.max(), 1000)
    # ax['a'].plot(x_kde, kde(x_kde), 'r-', linewidth=2, label='KDE')

    fig, ax = plt.subplot_mosaic('abc', figsize=(18, 6))
    sns.histplot(buf_tau, kde=True, bins=100, alpha=0.5, ax=ax['a']) # stat='density',
    sns.histplot(buf_cos, kde=True, bins=100, alpha=0.5, ax=ax['b'])
    sns.histplot(np.log10(buf_tau / buf_cos), kde=True, bins=100, alpha=0.5, ax=ax['c'])

    # ax['b'].hist(buf_cos, 30, histtype='bar', rwidth=0.8)

    ax['a'].set_title(r'$\tau$', fontsize=26)
    ax['b'].set_title(r'$\cos \alpha_0$', fontsize=26)
    ax['c'].set_title(r'$\log \tau/\cos \alpha_0$', fontsize=26)

    # fig.suptitle('cols', fontsize=22)

    prefix_folder = 'tau_hist/'
    file_name = f'{theta_obs} {beta_mu} {mc2} {a_portion} {phi_0} cols'
    # file_name = str(pathService.get_args_dir(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)) + 'cols'
    save_dir = pathService.get_dir(mu=None, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    # file_name = f'contour_max_interpolate'
    save.save_figure(fig, save_dir, file_name)


    tensor_tau_cols = save.load_tensor(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, 'tensor_tau_scatter')
    buf_tau = tensor_tau_cols.reshape(1, -1)
    buf_tau = buf_tau[buf_tau > 0]

    tensor_alpha_cols = save.load_tensor(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, 'tensor_alpha_scatter')
    buf_cos = tensor_alpha_cols.reshape(1, -1)
    buf_cos = buf_cos[buf_cos > 0]

    sns.ecdfplot(np.log10(buf_tau / buf_cos), ax=ax1['b'])
    ax1['b'].set_title(r'$\log \tau/\cos \alpha_1$', fontsize=26)

    fig, ax = plt.subplot_mosaic('abc', figsize=(18, 6))
    sns.histplot(buf_tau, kde=True, bins=100, alpha=0.5, ax=ax['a'])
    sns.histplot(buf_cos, kde=True, bins=100, alpha=0.5, ax=ax['b'])
    sns.histplot(np.log10(buf_tau / buf_cos), kde=True, bins=100, alpha=0.5, ax=ax['c'])

    # ax['b'].hist(buf_cos, 30, histtype='bar', rwidth=0.8)

    ax['a'].set_title(r'$\tau$', fontsize=26)
    ax['b'].set_title(r'$\cos \alpha_1$', fontsize=26)
    ax['c'].set_title(r'$\log \tau/\cos \alpha_1$', fontsize=26)

    # fig, ax = plt.subplot_mosaic('ab', figsize=(18, 6))
    # ax['a'].hist(buf_tau, 30, histtype='bar', rwidth=0.8)  # barstacked
    # ax['a'].set_title(r'$\tau$', fontsize=26)
    #
    # ax['b'].hist(buf_cos, 30, histtype='bar', rwidth=0.8)
    # ax['b'].set_title(r'$\cos \alpha$', fontsize=26)

    # fig.suptitle('scatter', fontsize=22)

    prefix_folder = 'tau_hist/'
    file_name = f'{theta_obs} {beta_mu} {mc2} {a_portion} {phi_0} scatter'
    save_dir = pathService.get_dir(mu=None, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    save.save_figure(fig, save_dir, file_name)

    file_name = f'{theta_obs} {beta_mu} {mc2} {a_portion} {phi_0} final'
    save.save_figure(fig1, save_dir, file_name)
    # counts, bin_edges = np.histogram(buf, 20)
    # bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2.
    # err = np.random.rand(bin_centres.size) * 100
    # ax['a'].errorbar(bin_centres, counts, yerr=err, fmt='o')
    #
    # plt.show()


def plot_sky_map_NS_scatter_on_off(mu, beta_mu, mc2, a_portion, phi_0):
    '''рисует карты неба (излучательный коэффициент в зависимости от положения наблюдателя)
    по другому - куча профилей. для разных наблюдателей. характеризует как источник излучает в пространство
    '''
    ticks_labelsize = 16
    # try_sky_map(obs_i_angle_arr)
    # obs_i_angle = np.linspace(0, 180, 19)
    L_x = save.load_L_x(mu, 10, beta_mu, mc2, a_portion, phi_0)

    def get_data_for_sky_map(mu, beta_mu, mc2, a_portion, phi_0):
        theta_obs_arr = np.linspace(0, 90, 10).astype(int)
        data_array = np.empty((len(theta_obs_arr) * 2 - 1, config.N_phase))
        # roll -- циклическая перестановка - делаем так как симметрическая задача и для угла 180 - theta будет симметрично
        # для сдвига на полфазы -- можно расчитать только до 90 а потом для других переставить и получить для до 180
        for i, theta_obs in enumerate(theta_obs_arr):
            L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            data_array[i] = L_total
            # для углов больших 90 будет симметрично относительно 90 (со сдвигом на полфазы)
            if i != len(theta_obs_arr) - 1:  # wtf -1 ?? - не будем трогать 90, его не дублируем
                data_array[-i - 1] = np.roll(data_array[i], len(data_array[i]) // 2)
        return data_array

    data_all_on = get_data_for_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    config.NS_shadow_flag = False
    data_NS_off = get_data_for_sky_map(mu, beta_mu, mc2, a_portion, phi_0)
    config.NS_shadow_flag = True

    config.flag_scatter = False
    data_scatter_off = get_data_for_sky_map(mu, beta_mu, mc2, a_portion, phi_0)
    config.flag_scatter = True

    config.flag_attenuation_above_shock = False
    data_tau_off = get_data_for_sky_map(mu, beta_mu, mc2, a_portion, phi_0)
    config.flag_attenuation_above_shock = True

    config.flag_scatter = False
    config.flag_attenuation_above_shock = False
    data_scatter_off_tau_off = get_data_for_sky_map(mu, beta_mu, mc2, a_portion, phi_0)
    config.flag_scatter = True
    config.flag_attenuation_above_shock = True

    def plot_maps_same_fig(data_to_plot_arr, prefix_folder, file_name_arr):
        theta_obs_arr_to_plot = np.linspace(0, 180, 19).astype(int)
        ['all_on', 'NS_off', 'scatter_off', 'tau_off', 'scatter_off_tau_off']
        ['all_on', 'tau_off', 'scatter_off_tau_off']

        fig, ax = plt.subplot_mosaic('abc', figsize=(8 * 3, 6), gridspec_kw={'width_ratios': [1, 1, 1]})

        ax['a'].contourf(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot_arr[4], levels=30, vmin=vmin,
                         vmax=vmax)
        ax['b'].contourf(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot_arr[2], levels=30, vmin=vmin,
                         vmax=vmax)

        # from mpl_toolkits.axes_grid1 import make_axes_locatable
        # # Для графика c создаем место для colorbar заранее
        # divider = make_axes_locatable(ax['c'])
        # cax = divider.append_axes("right", size="5%", pad=0.1)

        # Рисуем график на c
        im = ax['c'].contourf(config.phase_for_plot, theta_obs_arr_to_plot,
                              data_to_plot_arr[0], levels=30, vmin=vmin, vmax=vmax)

        # Получаем позицию последнего графика (ax['c'])
        pos = ax['c'].get_position()

        # Создаем ось для colorbar с такой же высотой как у графиков
        # [left, bottom, width, height] в долях от figure
        clb_width = 0.01
        cbar_ax = fig.add_axes([pos.x1 + clb_width,  # правее последнего графика
                                pos.y0,  # нижний край как у графиков
                                clb_width,  # ширина colorbar
                                pos.height])  # высота = высоте графиков

        # Добавляем colorbar
        clb = fig.colorbar(im, cax=cbar_ax)

        # cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        # clb = fig.colorbar(im, cax=cbar_ax)
        clb.set_label(r'$L_{\rm iso} / L_{\rm x}$', fontsize=24)
        clb.ax.tick_params(labelsize=ticks_labelsize)

        # clb = plt.colorbar(im, cax=cax)
        # clb.set_label(r'$L_{\rm iso} / L_{\rm x}$', fontsize=24)
        # clb.ax.tick_params(labelsize=ticks_labelsize)

        # im = ax['c'].contourf(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot_arr[0], levels=30, vmin=vmin,
        #                       vmax=vmax)
        #
        # clb = plt.colorbar(im, )
        # clb.set_label(r'$L_{\rm iso} / L_{\rm x}$', fontsize=24)
        # clb.ax.tick_params(labelsize=ticks_labelsize)

        for ax in [ax['a'], ax['b'], ax['c']]:
            ax.scatter([0, 1, 2], [beta_mu] * 3, c='white', marker='*', s=300)
            ax.scatter([0.5, 1.5], [180 - beta_mu] * 2, c='white', marker='*', s=300)
            ax.scatter([0, 1, 2], [beta_mu] * 3, c='black', marker='*', s=100)
            ax.scatter([0.5, 1.5], [180 - beta_mu] * 2, c='black', marker='*', s=100)

            x_axis_label = config.symbol_phase
            y_axis_label = config.symbol_theta_obs_y
            ax.set_xlabel(x_axis_label, fontsize=24)
            ax.set_ylabel(y_axis_label, fontsize=24)
            ax.tick_params(axis='both', labelsize=ticks_labelsize)

        save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion,
                                       phi_0=phi_0,
                                       prefix_folder=prefix_folder)
        # plt.tight_layout(rect=[0, 0, 0.85, 1])
        # plt.tight_layout()
        save.save_figure(fig, save_dir, 'all_maps')

    def plot_maps(data_to_plot, prefix_folder, file_name, vmin=None, vmax=None, flag_remove_cbar=False):

        theta_obs_arr_to_plot = np.linspace(0, 180, 19).astype(int)

        fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))
        if 'diff' in file_name:
            data_to_plot *= 100

        im = ax['a'].contourf(config.phase_for_plot, theta_obs_arr_to_plot, data_to_plot, levels=30, vmin=vmin,
                              vmax=vmax)
        ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='white', marker='*', s=300)
        ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='white', marker='*', s=300)
        ax['a'].scatter([0, 1, 2], [beta_mu] * 3, c='black', marker='*', s=100)
        ax['a'].scatter([0.5, 1.5], [180 - beta_mu] * 2, c='black', marker='*', s=100)

        x_axis_label = config.symbol_phase
        y_axis_label = config.symbol_theta_obs_y
        ax['a'].set_xlabel(x_axis_label, fontsize=24)
        ax['a'].set_ylabel(y_axis_label, fontsize=24)
        ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

        # from mpl_toolkits.axes_grid1 import make_axes_locatable
        # divider = make_axes_locatable(ax['a'])
        # cax = divider.append_axes("right", size="5%", pad=0.1)
        # cax.set_visible(False)  # Делаем ось невидимой

        clb = plt.colorbar(im, )

        if 'diff' in file_name:
            clb.set_label(r'$\%$', fontsize=24)
        else:
            clb.set_label(r'$L_{\rm iso} / L_{\rm x}$', fontsize=24)
        clb.ax.tick_params(labelsize=ticks_labelsize)
        if flag_remove_cbar:
            # clb.remove()
            ...
            # fig.subplots_adjust(right=0.85, top=0.85)
            # cax = ax['a'].inset_axes([1.03, 0, 0.1, 1])
            # cax.clear()
            # from mpl_toolkits.axes_grid1 import make_axes_locatable
            # divider = make_axes_locatable(ax['a'])
            # cax = divider.append_axes("right", size="5%", pad=0.05)
            # cax.clear()

        save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion,
                                       phi_0=phi_0,
                                       prefix_folder=prefix_folder)

        save.save_figure(fig, save_dir, file_name)
        # plt.tight_layout(rect=[0, 0, 0.85, 1])
        # save.save_fixed_size(fig, ax['a'], save_dir, file_name, has_cbar=False)

    def plot_PF(data_to_plot_arr, prefix_folder, file_name_arr):
        theta_obs_arr_to_plot = np.linspace(0, 180, 19).astype(int)
        fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))

        line_style = ['-', ':', '--', '-.']
        marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}

        j = 0
        for data_to_plot, label_name in zip(data_to_plot_arr, file_name_arr):
            PF = []
            if label_name != 'None':
                for i in range(10):
                    PF.append(newService.get_PF(data_to_plot[i]))

                # if label_name == 'No NS eclipse':
                ax['a'].plot(theta_obs_arr_to_plot[:10], PF, label=label_name, linestyle=line_style[j])
                j += 1
                # else:
                #     ax['a'].plot(theta_obs_arr_to_plot[:10], PF, label=label_name)

        x_axis_label = config.symbol_theta_obs_y
        y_axis_label = 'PF'
        ax['a'].set_xlabel(x_axis_label, fontsize=24)
        ax['a'].set_ylabel(y_axis_label, fontsize=24)
        ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)
        ax['a'].legend()
        save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=beta_mu, mc2=mc2, a_portion=a_portion,
                                       phi_0=phi_0,
                                       prefix_folder=prefix_folder)

        save.save_figure(fig, save_dir, 'PFs')

    vmin = np.min(np.hstack((data_all_on, data_scatter_off_tau_off, data_scatter_off, data_NS_off)) / L_x)
    vmax = np.max(np.hstack((data_all_on, data_scatter_off_tau_off, data_scatter_off, data_NS_off)) / L_x)

    PF_data_to_plot = []
    eps = 1e-8

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_all_on / L_x)
    plot_maps(data_to_plot, 'sky_map_difference/', 'all_on', vmin=vmin, vmax=vmax)
    # plot_maps(np.log10(data_to_plot + eps), 'sky_map_difference/', 'all_on_log', vmin=np.log10(vmin + eps),
    #           vmax=np.log10(vmax))
    plot_maps(np.log10(data_to_plot + eps), 'sky_map_difference/', 'all_on_log', vmin=np.log10(vmin + eps),
              vmax=np.log10(vmax))
    PF_data_to_plot.append(data_to_plot)

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_NS_off / L_x)
    plot_maps(data_to_plot, 'sky_map_difference/', 'NS_off', vmin=vmin, vmax=vmax)
    plot_maps(np.log10(data_to_plot + eps), 'sky_map_difference/', 'NS_off_log', vmin=np.log10(vmin + eps),
              vmax=np.log10(vmax))
    PF_data_to_plot.append(data_to_plot)
    # plot_PF(data_to_plot, 'sky_map_difference/', 'NS_off')

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1,
                                       arr=np.abs(data_all_on - data_NS_off) / data_all_on)
    plot_maps(data_to_plot, 'sky_map_difference/', 'NS_off_diff')

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_scatter_off / L_x)
    plot_maps(data_to_plot, 'sky_map_difference/', 'scatter_off', vmin=vmin, vmax=vmax)
    plot_maps(np.log10(data_to_plot + eps), 'sky_map_difference/', 'scatter_off_log', vmin=np.log10(vmin + eps),
              vmax=np.log10(vmax))
    PF_data_to_plot.append(data_to_plot)
    # plot_PF(data_to_plot, 'sky_map_difference/', 'scatter_off')

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1,
                                       arr=np.abs(data_all_on - data_scatter_off) / L_x)
    plot_maps(data_to_plot, 'sky_map_difference/', 'scatter_off_diff')

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_tau_off / L_x)
    plot_maps(data_to_plot, 'sky_map_difference/', 'tau_off', vmin=vmin, vmax=vmax)
    plot_maps(np.log10(data_to_plot + eps), 'sky_map_difference/', 'tau_off_log', vmin=np.log10(vmin + eps),
              vmax=np.log10(vmax))
    PF_data_to_plot.append(data_to_plot)
    # plot_PF(data_to_plot, 'sky_map_difference/', 'tau_off')

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1,
                                       arr=np.abs(data_all_on - data_tau_off) / L_x)
    plot_maps(data_to_plot, 'sky_map_difference/', 'tau_off_diff')

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_scatter_off_tau_off / L_x)
    plot_maps(data_to_plot, 'sky_map_difference/', 'scatter_off_tau_off', vmin=vmin, vmax=vmax)
    plot_maps(np.log10(data_to_plot + eps), 'sky_map_difference/', 'scatter_off_tau_off_log', vmin=np.log10(vmin + eps),
              vmax=np.log10(vmax))
    PF_data_to_plot.append(data_to_plot)
    # plot_PF(data_to_plot, 'sky_map_difference/', 'scatter_off_tau_off')

    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1,
                                       arr=np.abs(data_all_on - data_scatter_off_tau_off) / L_x)
    plot_maps(data_to_plot, 'sky_map_difference/', 'scatter_off_tau_off_diff')

    plot_maps_same_fig(PF_data_to_plot, 'sky_map_difference/',
                       ['all_on', 'NS_off', 'scatter_off', 'tau_off', 'scatter_off_tau_off'])

    plot_PF(PF_data_to_plot, 'sky_map_difference/',
            ['Total', 'No NS eclipse', 'No FF reflection', 'None', 'No FF effect'])
    # ['all_on', 'NS_off', 'scatter_off', 'tau_off', 'scatter_off_tau_off']
    # ['all_on', 'No NS eclipse', 'scatter_off', 'tau_off', 'scatter_off_tau_off']


def plot_PF_obs(mu, beta_mu, mc2_arr, a_portion_arr, phi_0):
    '''PF to obs'''
    marker_index = 0
    line_style = ['-', '--']
    marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}

    # cmap = mpl.cm.viridis
    global cmap

    theta_obs_arr = [10 * i for i in range(0, 10)]

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    for a_index, a_portion in enumerate(a_portion_arr):
        for mc2_index, mc2 in enumerate(mc2_arr):
            full_dict = {}
            for theta_obs_index, theta_obs in enumerate(theta_obs_arr):
                L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                PF = newService.get_PF(L_total)
                full_dict[np.mean(L_total)] = (PF, theta_obs_arr[theta_obs_index])

            # в словаере лежит среднее : (PF, phi_0). сортируем его по значению phi_0
            list_tuples = sorted(full_dict.items(), key=lambda item: item[1][1])
            x_sort_phi_0, y_sort_phi_0 = zip(*list_tuples)
            y_sort_phi_0, colors_sort_phi_0 = zip(*y_sort_phi_0)

            colors = (np.array(colors_sort_phi_0)) / np.max(colors_sort_phi_0)

            if marker_index == 0:
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, color=cmap(colors))
                # ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, facecolors='none', edgecolors=cm.jet(colors))
            else:
                ax['a'].scatter(x_sort_phi_0, y_sort_phi_0, s=100, marker=marker_dict[marker_index % 4],
                                color=cmap(colors))

            ax['a'].plot(x_sort_phi_0, y_sort_phi_0, color='black', alpha=0.2, linestyle=line_style[mc2_index])

        marker_index = 3

    x_axis_label = r'$L_{\rm iso}$' + r'$\rm [erg/s]$'
    y_axis_label = r'$PF$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].set_xscale('log')
    # ax['a'].set_ylim(0, 1.1)

    bounds = theta_obs_arr.copy()
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax['a'], orientation='vertical', pad=0.01)
    cb.set_label(label=config.symbol_theta_obs_y, fontsize=24)

    prefix_folder = 'PF_to_obs/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=phi_0, prefix_folder=prefix_folder)

    file_name = f'beta_mu={beta_mu} phi_0={phi_0} All_PF_to_obs.png'
    save.save_figure(fig, save_dir, file_name)


def plot_romanova(mu, theta_obs, beta_mu, mc2_arr, a_portion):
    # romanova_tuples = [(10, 0), (20, 28), (30, 43), (40, 52), (50, 60), (60, 64), (70, 68), (80, 71), (90, 74),
    #                    (100, 77)]
    # for i, tuple in enumerate(romanova_tuples):
    #     mc2 = tuple[0]
    #     phi_0 = tuple[1]
    # mc2_arr = np.array([10 * i for i in range(1, len(romanova_tuples)+1)])

    import romanova_phi_0
    phi_0_arr = romanova_phi_0.phi_0_arr
    mc2_arr = romanova_phi_0.mc2_arr

    data_array = np.empty((len(phi_0_arr), config.N_phase))

    if not isinstance(a_portion, np.ndarray):
        a_portion = np.full(10, a_portion)
        suffix = ''
    else:
        suffix = f'a_range {a_portion[0]} - {a_portion[-1]}'

    for i in range(len(phi_0_arr)):
        phi_0 = phi_0_arr[i]
        mc2 = mc2_arr[i]
        L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion[i], phi_0)
        data_array[i] = L_total

    avg = 1
    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_array)
    data_to_plot = data_to_plot / avg

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    im = ax['a'].contourf(config.phase_for_plot, mc2_arr, data_to_plot, levels=30)
    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$\langle L_{\rm iso} \rangle$', fontsize=26)

    x_axis_label = config.symbol_phase
    y_axis_label = config.symbol_m
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    # ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    prefix_folder = 'romanova/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=a_portion[0],
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'P={romanova_phi_0.P} contour_L_iso' + suffix
    save.save_figure(fig, save_dir, file_name)

    avg = np.mean(data_array, axis=1)[:, np.newaxis]
    data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=data_array)
    data_to_plot = data_to_plot / avg

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    im = ax['a'].contourf(config.phase_for_plot, mc2_arr, data_to_plot, levels=30)
    clb = plt.colorbar(im, pad=0.01)
    clb.set_label(r'$\langle L_{\rm iso} \rangle$', fontsize=26)

    x_axis_label = config.symbol_phase
    y_axis_label = config.symbol_m
    ax['a'].set_xlabel(x_axis_label, fontsize=26)
    ax['a'].set_ylabel(y_axis_label, fontsize=26)

    # ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)

    prefix_folder = 'romanova/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=a_portion[0],
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'P={romanova_phi_0.P} contour_L_iso_avg' + suffix
    save.save_figure(fig, save_dir, file_name)


def plot_show_phase_shift(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr):
    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))

    line_style = ['-', '--']
    for i, phi_0 in enumerate(phi_0_arr):
        L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
        data_to_plot = np.apply_along_axis(newService.extend_arr_for_plot, axis=-1, arr=L_total)
        ax['a'].plot(config.phase_for_plot, data_to_plot, label=r'$\phi_0 = $' + f'{phi_0}', linestyle=line_style[i],
                     color='black')
    ax['a'].legend()
    prefix_folder = 'phase_shift_phi/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=theta_obs, beta_mu=beta_mu, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'{mc2=} {a_portion=} phi_0={phi_0_arr}'
    save.save_figure(fig, save_dir, file_name)


def plot_Teff_to_ksi_diff_a(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0):
    '''рисует Т от кси'''

    color_arr = ['red', 'green', 'blue', 'purple']
    line_style = ['-', '--']
    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    ksi_max = [0]
    for j, mc2 in enumerate(mc2_arr):
        for i, a_portion in enumerate(a_portion_arr):
            # R_e_arr = np.empty(len(mc2_arr))
            R_e = save.load_R_e(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            # ksi_shock = save.load_ksi_shock(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            T_eff = save.load_T_eff(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            theta_range = save.load_theta_range(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
            ksi = R_e * (np.sin(theta_range)) ** 2

            # ax['a'].scatter(ksi, T_eff, s=30, facecolors='none', edgecolors=color_arr[i], label=f'a={a_portion}')
            ax['a'].plot(ksi, T_eff, color=color_arr[i], alpha=0.8, linestyle=line_style[j],
                         label=f'a={a_portion}' + config.symbol_m + f'={mc2}')

            y1 = np.log(T_eff)
            x1 = np.log(ksi)
            n = 40
            a = stats.linregress(x1[:n], y1[:n])
            print(a)
            if ksi[-1] > ksi_max[-1]:
                ksi_max = ksi
            # ax['a'].plot(ksi, ksi ** (a[0]) * max(T_eff), color='black', alpha=0.8, linestyle='-', label=f'regr')
        # ax['a'].plot(ksi, T_eff, color='black')
    ax['a'].plot(ksi_max, ksi_max ** (-5 / 8) * max(T_eff), color='black', alpha=0.8, linestyle='-', label=f'-5/8')
    x_axis_label = r'$\xi$'
    y_axis_label = r'$T_{\rm eff} ~ [K]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].legend()

    ax['a'].set_xscale('log')
    ax['a'].set_yscale('log')

    prefix_folder = 'table/'
    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                   phi_0=None, prefix_folder=prefix_folder)
    file_name = f'T_eff'
    save.save_figure(fig, save_dir, file_name)

    # ax['a'].set_xlim(1, 8)
    file_name = f'T_eff_lim'
    save.save_figure(fig, save_dir, file_name)


def plot_coeff_gamma():
    import scipy.special as special
    gamma = np.linspace(0, 1, 100)
    ksi_shock_arr = [2, 3, 5, 8, 10]
    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    for ksi_shock in ksi_shock_arr:
        beta = 1 - gamma * np.exp(gamma) * (special.expn(1, gamma) - special.expn(1, gamma * ksi_shock))
        coef = -3 / 8 + gamma / 4 - 1 / (4 * beta)
        ax['a'].plot(gamma, coef, label=r'$\xi_s=$' + f'{ksi_shock}', linestyle='-')
        coef = -3 / 8 + gamma / 4
        # ax['a'].plot(gamma, coef, label=r'$\xi_s=$' + f'{ksi_shock}', linestyle='--')
        coef = -3 / 8 - 1 / (4 * beta)
        # ax['a'].plot(gamma, coef, label=r'$\xi_s=$' + f'{ksi_shock}', linestyle='-.')

    x_axis_label = r'$\gamma$'
    y_axis_label = r'$p$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    ax['a'].legend()
    plt.show()


def plot_PF_to_theta(mu, beta_mu_arr, mc2_arr, a_portion, phi_0):
    theta_obs_arr = np.linspace(0, 90, 10).astype(int)

    fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))
    for beta_mu in beta_mu_arr:
        for mc2 in mc2_arr:
            data = np.zeros((len(theta_obs_arr)))
            for i, theta_obs in enumerate(theta_obs_arr):
                L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                data[i] = newService.get_PF(L_total)
            label_name = config.symbol_mu_ang + f'={beta_mu}' + ' ' + config.symbol_m + f'={mc2}'
            ax['a'].plot(theta_obs_arr, data, label=label_name)

    x_axis_label = config.symbol_theta_obs_y
    y_axis_label = 'PF'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)
    ax['a'].tick_params(axis='both', labelsize=ticks_labelsize)
    ax['a'].legend()

    ax['a'].set_ylim(0, 1)

    prefix_folder = 'PF_to_theta'

    save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=a_portion,
                                   phi_0=phi_0,
                                   prefix_folder=prefix_folder)

    save.save_figure(fig, save_dir, 'PFs')


if __name__ == '__main__':

    # plot_a_restrictions()
    mu = 0.1e31
    theta_obs = 10
    beta_mu = 20
    mc2 = 100
    a_portion = 0.66
    phi_0 = 0

    theta_obs = 60
    beta_mu = 20
    mc2 = 60
    a_portion = 0.8
    phi_0 = 0
    # plot_hist_tau(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)  # active

    beta_mu = 20
    mc2 = 60
    a_portion = 0.75
    phi_0 = 0
    # plot_sky_map_NS_scatter_on_off(mu, beta_mu, mc2, a_portion, phi_0)

    theta_obs = 20
    beta_mu = 20
    mc2 = 30
    a_portion = 0.22
    phi_0_arr = [0, 180]
    # plot_show_phase_shift(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr)

    theta_obs = 60
    beta_mu = 20
    a_portion = 0.22

    a_portion_start = 0.44
    a_portion_stop = 0.66
    a_portion_arr = np.linspace(a_portion_start, a_portion_stop, 10)

    # plot_romanova(mu, theta_obs, beta_mu, None, a_portion)

    theta_obs = 60
    beta_mu = 20
    mc2_arr = [30, 100]
    a_portion_arr = [0.25, 0.75]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    # plot_masses_PF_L_on_off(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)

    beta_mu = 70
    mc2_arr = [30, 100]
    a_portion_arr = [0.22, 0.66]
    phi_0 = 0
    # plot_PF_obs(mu, beta_mu, mc2_arr, a_portion_arr, phi_0)

    mu = 0.1e31
    beta_mu = 40
    mc2 = 100
    a_portion = 0.66
    phi_0 = 0
    theta_obs = 20

    # ---------------------- $PF(Liso)$ for $a=const$
    theta_obs = 40
    beta_mu = 20
    mc2_arr = [20, 40, 60, 80, 100, 120]
    a_portion_arr = [0.22, 0.44, 0.66]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    for a_portion in a_portion_arr:
        # plot_PF_to_L_const_a_dif_m(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0_arr)
        ...
    #
    theta_obs = 40
    beta_mu = 20
    mc2_arr = [20, 40, 60, 80, 100, 120]
    a_portion_arr = [0.22, 0.44, 0.66]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    # plot_PF_to_L_const_a_const_phi0_dif_m(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0)

    theta_obs = 20
    beta_mu_arr = [10 * i for i in range(0, 9)]
    mc2 = 30
    a_portion = 0.66
    phi_0 = 0
    # plot_ksi_to_beta(mu, theta_obs, beta_mu_arr, mc2, a_portion, phi_0)

    theta_obs = 60
    beta_mu = 20
    mc2_arr = [30, 100]
    a_portion_arr = [0.22, 0.66]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    # plot_masses_PF_L(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)

    # -------------------
    theta_obs = 40
    beta_mu = 20
    mc2_arr = [20, 40, 60, 80, 100, 120]
    a_portion_arr = [0.22, 0.44, 0.66, 1]
    phi_0 = 0
    # plot_L_max_phase_to_m_to_a(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0)
    theta_obs = 40

    beta_mu = 20
    mc2 = 100
    a_portion = 1
    phi_0 = 0
    # plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0, True)

    theta_obs = 20
    beta_mu = 40
    mc2 = 60
    a_portion = 0.66
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    # plot_L_to_phi_0(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr, True)

    theta_obs = 20
    beta_mu = 20
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.66
    phi_0 = 20
    # (mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    theta_obs = 60
    beta_mu = 20
    mc2 = 60
    # a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    a_portion_arr = np.linspace(0.1, 1, 19)
    phi_0 = 0
    # plot_L_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)

    theta_obs = 60
    beta_mu = 20
    mc2_arr = [30, 100]
    a_portion_arr = [0.22, 0.66]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]

    theta_obs = 40
    beta_mu_arr = [10 * i for i in range(0, 9)]
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.44
    phi_0 = 0

    # plot_table_R_disk(mu, theta_obs, beta_mu_arr, a_portion, mc2_arr, phi_0)

    # plot_masses_PF_L_nu(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)

    theta_obs = 20
    beta_mu = 40
    mc2 = 100
    a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    phi_0 = 0
    # plot_PF_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)  # wtf ???

    # theta_obs_arr = [20, 60, 80]
    # beta_mu = 20
    # mc2 = 30
    # a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    # phi_0 = 0
    # for theta_obs in theta_obs_arr:
    # plot_L_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)

    # ----- old --------
    # beta_mu = 60
    # mc2 = 30
    # a_portion = 1
    # phi_0 = 0
    # # plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    # theta_obs = 20
    # beta_mu = 60
    # mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    # a_portion = 0.44
    # phi_0 = 0
    # # plot_L_to_mc2(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    # theta_obs = 20
    # beta_mu = 40
    # mc2 = 100
    # a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    # phi_0 = 0
    # # plot_L_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)
    # # plot_PF_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)

    # theta_obs = 40
    # beta_mu = 60
    # mc2 = 30
    # a_portion = 0.22
    # phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    # # plot_L_to_phi_0(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr, True)

    # theta_obs = 20
    # beta_mu = 60
    # mc2_arr = [30, 100]
    # a_portion_arr = [0.22, 0.66]
    # phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    # # plot_masses_PF_L_nu(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)
    # # plot_masses_PF_L(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)

    # mc2 = 30
    # a_portion = 1
    # phi_0 = 0
    # # plot_PF_contour(mu, mc2, a_portion, phi_0)
    #
    # mc2 = 30
    # a_portion = 1
    # phi_0 = 0
    # # plot_PF_to_chi_theta(mu, mc2, a_portion, phi_0)

    # theta_obs = 40
    # beta_mu = 60
    # mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    # a_portion = 0.44
    # phi_0 = 0
    # # plot_L_iso_to_m(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)
    #

    theta_obs = 20
    beta_mu = 0
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.66
    phi_0 = 0
    # plot_table(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    a_portion_arr = [0.25, 0.5, 0.75, 1]
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    # plot_table_together(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0)
    # plot_Teff_to_ksi_diff_a(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0)

    # plot_coeff_gamma()
    plot_rvm(mu, 20, 60)

    beta_mu_arr = [40, 80]
    a_portion = 1
    mc2_arr = [30, 60, 100]
    phi_0 = 0
    # plot_PF_to_theta(mu, beta_mu_arr, mc2_arr, a_portion, phi_0)
