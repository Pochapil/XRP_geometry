import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

import geometricTask.matrix as matrix
import config
import main_service

plt.style.use(['science', 'notebook', 'grid'])  # для красивых графиков

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'


# mpl.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

def add_to_ax(ax, i_angle, betta_mu, mc2, a_portion, fi_0, count):
    label = f'i={i_angle} chi={betta_mu} mc2={mc2} a={a_portion}'
    full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2, a_portion, fi_0)

    file_name = "total_luminosity_of_surfaces.txt"
    L_arr = main_service.load_arr_from_txt(full_file_folder, file_name)[4]

    buf = main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/',
                                         'scattered_energy_top.txt')
    if not np.isnan(buf).any():
        L_arr += buf

    buf = main_service.load_arr_from_txt(full_file_folder + 'scattered_on_magnet_lines/',
                                         'scattered_energy_bot.txt')
    if not np.isnan(buf).any():
        L_arr += buf

    y_data = L_arr.copy()
    mean = np.mean(y_data)
    std = np.std(y_data)
    # y_data = (y_data - mean) / std

    max_val = np.max(y_data)
    # y_data = y_data / mean

    mean_label = np.mean(y_data) / 1e38
    # y_data *= mean_label
    y_data /= 1e38
    max_idx = np.argmax(y_data)
    y_data = np.roll(y_data, -max_idx)

    y_data_plot = main_service.extend_arr_for_phase(y_data)

    PF = (np.max(y_data) - np.min(y_data)) / (np.max(y_data) + np.min(y_data))

    label = label + f' mean/1e38 = {mean_label:.2f} PF = {PF:.2f}'

    count = count // 7 + 1
    alpha = 1/count

    ax.plot(x_data_plot, y_data_plot, alpha=alpha, label=label)
    # ax.plot(x_data_plot, fit_all_plot, label='modelled')


x_data = np.linspace(0, 2 * np.pi, 45)
x_data_plot = np.linspace(0, 2, 90)

# figlegend = plt.figure()
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
# ax.plot(x_data_plot, y_data_plot + noise, label='real')

i_angle_arr = [20, 40]
betta_mu_arr = [20, 70]
mc2_arr = [30, 100]
a_portion_arr = [0.22, 0.44]
fi_0_arr = [0]

count = 1

for i_angle in i_angle_arr:
    for betta_mu in betta_mu_arr:
        for mc2 in mc2_arr:
            for a_portion in a_portion_arr:
                for fi_0 in fi_0_arr:
                    fi_0 += config.fi_0_dict[a_portion]
                    fi_0 = fi_0 % 360

                    add_to_ax(ax, i_angle, betta_mu, mc2, a_portion, fi_0, count)
                    count += 1


# i_angle, betta_mu, mc2, a_portion, fi_0 = 20, 20, 100, 0.22, 320
#
#
# i_angle, betta_mu, mc2, a_portion, fi_0 = 20, 20, 100, 0.22, 40
# add_to_ax(ax, i_angle, betta_mu, mc2, a_portion, fi_0)
#
# i_angle, betta_mu, mc2, a_portion, fi_0 = 20, 20, 30, 0.22, 320
# add_to_ax(ax, i_angle, betta_mu, mc2, a_portion, fi_0)
#
# i_angle, betta_mu, mc2, a_portion, fi_0 = 40, 70, 100, 0.66, 240
# add_to_ax(ax, i_angle, betta_mu, mc2, a_portion, fi_0)
#
# i_angle, betta_mu, mc2, a_portion, fi_0 = 20, 20, 100, 0.44, 280
# add_to_ax(ax, i_angle, betta_mu, mc2, a_portion, fi_0)
#
# i_angle, betta_mu, mc2, a_portion, fi_0 = 20, 70, 100, 0.44, 280
# add_to_ax(ax, i_angle, betta_mu, mc2, a_portion, fi_0)
#
# i_angle, betta_mu, mc2, a_portion, fi_0 = 40, 70, 100, 0.44, 280
# add_to_ax(ax, i_angle, betta_mu, mc2, a_portion, fi_0)


x_axis_label = r'$\Phi$'
y_axis_label = r'$relative \, L_{iso}$'

ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)

# plt.legend()
plt.show()
# noise = np.random.normal(loc=0.0, scale=0.2, size=90)
# figlegend.legend(ax.get_legend_handles_labels())
# plt.show()
