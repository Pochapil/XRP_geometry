import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

import config

# plt.style.use(['science', 'notebook', 'grid'])
# чтобы стиль был похож на теховский
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'STIXGeneral'

# стандартный размер тиков по осям
ticks_labelsize = 18
mpl.rcParams['xtick.labelsize'] = ticks_labelsize
mpl.rcParams['ytick.labelsize'] = ticks_labelsize


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


if __name__ == '__main__':

    # phi_range = np.linspace(-np.pi * a_portion, np.pi * a_portion, config.N_phi_accretion) + np.deg2rad(phi_0)
    #
    # disk_min_theta_angle = np.min(np.pi / 2 - np.arctan(np.tan(np.deg2rad(beta_mu)) * np.cos(phi_range)))
    # theta_end = disk_min_theta_angle
    #
    # ans = (1 / np.tan(chi) ** 2 * (1 - (1 + delta) * np.sin(theta_end) ** 2) / (
    #         (1 + delta) * np.sin(theta_end) ** 2))
    #
    # cos_arr = np.cos(phi_range) ** 2
    #
    # if (cos_arr >= ans).all():
    #     print(f'ok {a_portion}')
    # else:
    #     print(f'NO! {a_portion}')


    beta_mu = 30
    a_portion = 0.4
    phi_0 = 60

    delta = 0.25

    chi = np.deg2rad(beta_mu)


    def func(chi_arr, phi_0):
        # chi_arr = [10 * i for i in range(1, 9)]

        res = []
        eps = 1e-3
        for chi in chi_arr:
            left = eps
            right = 1
            while right - left > eps:
                middle = (left + right) / 2
                a_portion = middle

                phi_range = np.linspace(-np.pi * a_portion, np.pi * a_portion, config.N_phi_accretion) + np.deg2rad(phi_0)

                disk_min_theta_angle = np.min(np.pi / 2 - np.arctan(np.tan(chi) * np.cos(phi_range)))
                theta_end = disk_min_theta_angle

                ans = (1 / np.tan(chi) ** 2 * (1 - (1 + delta) * np.sin(theta_end) ** 2) / (
                        (1 + delta) * np.sin(theta_end) ** 2))

                cos_arr = np.cos(phi_range) ** 2

                if (cos_arr >= ans).all():
                    left = middle
                    # print(f'ok {a_portion}')
                else:
                    right = middle
                    # print(f'NO! {a_portion}')
            res.append(right)
        return res


    eps = 1e-5
    chi_arr = np.linspace(eps, np.pi / 2, 200, endpoint=False)
    res = func(chi_arr, 0)
    print(res)

    marker_dict = {0: '.', 1: '*', 2: '+', 3: '^'}
    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    phi_0_arr = [0, 20, 40, 60, 80, 90, 100, 120, 140, 160, 180]

    line_style = ['-', '--']
    for phi_0 in phi_0_arr:
        line_style = '-'
        if phi_0 in [0, 90, 180]:
            line_style = '--'
        res = func(chi_arr, phi_0)
        ax['a'].plot(np.rad2deg(chi_arr), res, linestyle=line_style,  label=f'{phi_0=}')

    x_axis_label = r'$\chi ~ [^{\circ}]$'
    y_axis_label = r'$a$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    ax['a'].grid()
    plt.legend(fontsize=18)
    plt.show()