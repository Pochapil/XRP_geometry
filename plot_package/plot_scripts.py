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


def plot_total_luminosity_of_surfaces(L_surfs, save_dir=None):
    legend_labels_arr = [r'$top_{pol}$', r'$top_{eq}$', r'$bottom_{pol}$', r'$bottom_{eq}$', 'sum']

    colors_arr = ['red', 'red', 'green', 'green']
    marker_arr = ['.', '*', '+', '^']
    line_style_arr = [':', '-']

    phi_for_plot = np.linspace(0, 2, config.N_phase_for_plot)

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    for i, L_surf in enumerate(L_surfs):
        ax['a'].plot(phi_for_plot, newService.extend_arr_for_plot(L_surf), label=legend_labels_arr[i],
                     color=colors_arr[i], marker=marker_arr[i], linestyle=line_style_arr[i % 2], markersize=12)

    ax['a'].plot(phi_for_plot, newService.extend_arr_for_plot(np.sum(L_surfs, axis=0)), label=legend_labels_arr[-1],
                 color='black')

    x_axis_label = r'$\Phi$'
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
    # file_name = 'observer_angles.txt'
    # data_array = main_service.load_arr_from_txt(working_folder, file_name)
    # observer_phi = data_array[0]
    # observer_theta = data_array[1]

    observer_phi = newService.extend_arr_for_plot(observer_phi)
    observer_theta = newService.extend_arr_for_plot(observer_theta)

    fig_title = 'Observer angles'

    phi_for_plot = np.linspace(0, 2, config.N_phase_for_plot)
    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    ax['a'].plot(phi_for_plot, np.rad2deg(observer_theta), label=r'$\theta_{observer}$')
    ax['a'].plot(phi_for_plot, np.rad2deg(observer_phi), label=r'$\phi_{observer}$')
    ax['a'].legend()

    x_axis_label = r'$\Phi$'
    y_axis_label = r'$angle ~ [deg]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    if save_dir is not None:
        file_name = 'Observer_angles'
        save.save_figure(fig, save_dir, file_name)


def plot_Teff_to_ksi(R_e, T_eff, theta_range, save_dir=None):
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
    y_axis_label = r'$T_{eff} \: [K]$'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    if save_dir is not None:
        file_name = 'T_eff'
        save.save_figure(fig, save_dir, file_name)


def f():
    folder = 'L_nu/'

    phi_for_plot = np.linspace(0, 2, config.N_phase_for_plot)

    x_axis_label = r'$h \nu$' + ' [keV]'
    y_axis_label = r'$L_{\nu} \: [erg \cdot s^{-1} \cdot hz^{-1}]$'


def plot_PF_to_energy(L_nu_surfs, L_nu_scatter, save_dir=None):
    # sum(L_nu_surfs + L_nu_scatter)

    folder = 'L_nu/'

    energy_arr = config.energy_arr
    L_nu = np.sum(L_nu_surfs, axis=0) + np.sum(L_nu_scatter, axis=0)
    PF = newService.get_PF(L_nu)

    fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
    ax['a'].plot(energy_arr, PF, color='black')
    x_axis_label = r'$h \nu$' + ' [keV]'
    y_axis_label = 'PF'
    ax['a'].set_xlabel(x_axis_label, fontsize=24)
    ax['a'].set_ylabel(y_axis_label, fontsize=24)

    if save_dir is not None:
        file_name = 'PF'
        save.save_figure(fig, save_dir / folder, file_name)


def plot_L_nu(L_nu_surfs, save_dir):
    # sum L_nu
    folder = 'L_nu/'

    L_nu_to_plot = newService.extend_arr_for_plot(L_nu_surfs)

    phi_for_plot = np.linspace(0, 2, config.N_phase_for_plot)

    x_axis_label = r'$h \nu$' + ' [keV]'
    y_axis_label = r'$L_{\nu} \: [erg \cdot s^{-1} \cdot hz^{-1}]$'

    for i, energy in enumerate(config.energy_arr):
        fig, ax = plt.subplot_mosaic('a', figsize=(12, 6))
        # for L_nu_surf in L_nu_surfs[:, i, :]:
        L_nu_to_plot = newService.extend_arr_for_plot(np.sum(L_nu_surfs[:, i, :], axis=0))
        PF = newService.get_PF(L_nu_to_plot)
        fig_title = r'$L_{\nu}^' + r'{' + f'{energy:.2f} keV' + r'}' + r'(\Phi)$' + '  ' + f'PF = {PF:.3f}'
        ax['a'].plot(phi_for_plot, L_nu_to_plot, color='black',
                     label=r'$L_{\nu}^' + r'{' + f'{energy:.2f}' + r'}' + r'(\Phi)$')
        ax['a'].set_xlabel(x_axis_label, fontsize=24)
        ax['a'].set_ylabel(y_axis_label, fontsize=24)
        # ax['a'].legend()
        fig.suptitle(fig_title, fontsize=16)

        if save_dir is not None:
            file_name = f'L_nu_of_energy_{energy:.2f}_KeV_of_surfaces'
            save.save_figure(fig, save_dir / folder, file_name)


if __name__ == '__main__':
    L1 = [[0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
           0.00000000e+00, 0.00000000e+00, 7.74675923e+35, 5.18107356e+37,
           6.07466173e+37, 6.60328489e+37, 6.88307236e+37, 6.92754072e+37,
           6.74831813e+37, 6.36379964e+37, 5.80172783e+37, 5.09673155e+37,
           4.27782056e+37, 3.38483264e+37, 2.44768525e+37, 1.60911792e+37,
           8.97480927e+36, 4.11467726e+36, 2.33266697e+36, 2.33266697e+36,
           4.11467726e+36, 8.97480927e+36, 1.60911792e+37, 2.44768525e+37,
           3.38483264e+37, 4.27782056e+37, 5.09673155e+37, 5.80172783e+37,
           6.36379964e+37, 6.74831813e+37, 6.92754072e+37, 6.88307236e+37,
           6.60328489e+37, 6.07466173e+37, 5.18107356e+37, 7.74675923e+35,
           0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
           0.00000000e+00],
          [2.92594683e+38, 2.94732594e+38, 3.01881578e+38, 3.10697842e+38,
           3.07191100e+38, 2.94506643e+38, 2.75274600e+38, 2.58320916e+38,
           2.76913805e+38, 2.98650499e+38, 3.20326055e+38, 3.40945845e+38,
           3.60553250e+38, 3.77968799e+38, 3.93196516e+38, 4.06667819e+38,
           4.18484405e+38, 4.28017196e+38, 4.36540182e+38, 4.42701025e+38,
           4.46857739e+38, 4.49579977e+38, 4.47871846e+38, 4.47871846e+38,
           4.49579977e+38, 4.46857739e+38, 4.42701025e+38, 4.36540182e+38,
           4.28017196e+38, 4.18484405e+38, 4.06667819e+38, 3.93196516e+38,
           3.77968799e+38, 3.60553250e+38, 3.40945845e+38, 3.20326055e+38,
           2.98650499e+38, 2.76913805e+38, 2.58320916e+38, 2.75274600e+38,
           2.94506643e+38, 3.07191100e+38, 3.10697842e+38, 3.01881578e+38,
           2.94732594e+38],
          [1.41206174e+38, 1.43003828e+38, 1.50506024e+38, 1.68606097e+38,
           2.04719317e+38, 2.71544862e+38, 4.32422422e+38, 4.54576215e+38,
           4.72142996e+38, 4.86508073e+38, 5.03539513e+38, 5.12118880e+38,
           5.15559606e+38, 5.16109324e+38, 5.14762860e+38, 5.11607258e+38,
           5.06573173e+38, 5.00286918e+38, 4.93617197e+38, 4.87289966e+38,
           4.82118014e+38, 4.78837437e+38, 4.77928546e+38, 4.77928546e+38,
           4.78837437e+38, 4.82118014e+38, 4.87289966e+38, 4.93617197e+38,
           5.00286918e+38, 5.06573173e+38, 5.11607258e+38, 5.14762860e+38,
           5.16109324e+38, 5.15559606e+38, 5.12118880e+38, 5.03539513e+38,
           4.86508073e+38, 4.72142996e+38, 4.54576215e+38, 4.32422422e+38,
           2.71544862e+38, 2.04719317e+38, 1.68606097e+38, 1.50506024e+38,
           1.43003828e+38],
          [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
           0.00000000e+00, 0.00000000e+00, 2.72409853e+36, 1.60094353e+35,
           0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
           0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
           0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
           3.73605687e+34, 6.51295841e+34, 1.18419567e+34, 1.18419567e+34,
           6.51295841e+34, 3.73605687e+34, 0.00000000e+00, 0.00000000e+00,
           0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
           0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
           0.00000000e+00, 0.00000000e+00, 1.60094353e+35, 2.72409853e+36,
           0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
           0.00000000e+00]
          ]

    buf = L1.copy()
    L1[0] = buf[1]
    L1[1] = buf[0]
    L1[2] = buf[3]
    L1[3] = buf[2]

    plot_total_luminosity_of_surfaces(L1)
