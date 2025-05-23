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


obs_i_angle_arr_to_plot = np.linspace(10, 170, 3)


def fit_func_all(x, a_1, b_1, a_2, b_2, a_3, b_3, c):
    return a_1 * np.sin(x + b_1) + a_2 * np.sin(2 * x + b_2) + a_3 * np.sin(3 * x + b_3) + c

i_angle, betta_mu, mc2, a_portion, fi_0 = 20, 40, 100, 0.66, 240

full_file_folder = config.get_folder_with_args(i_angle, betta_mu, mc2, a_portion, fi_0)

file_name = "total_luminosity_of_surfaces.txt"
L_arr = main_service.load_arr_from_txt(full_file_folder, file_name)[4]

y_data = L_arr
mean = np.mean(y_data)
std = np.std(y_data)
# y_data = (y_data - mean) / std

max_val = np.max(y_data)
y_data = y_data / mean

x_data = np.linspace(0, 2 * np.pi, 45)

# popt, pcov = curve_fit(fit_func_all, x_data, y_data)
# fit_all = fit_func_all(x_data, *popt)
# fit_params = popt

a_1 = 1.3671454589735008
b_1 = 1.5595338205691776
a_2 = 2.63740593267352524
b_2 = 1.001138446292847
a_3 = 1.14638485845640548
b_3 = 2.3590588942813622
c = -0.03633449368283235

# y_data = fit_func_all(x_data, *fit_params)
fit_all = np.roll(y_data[::1], 30)

fit_all_plot = main_service.extend_arr_for_phase(fit_all)
y_data_plot = main_service.extend_arr_for_phase(y_data)
x_data_plot = np.linspace(0, 2, 90)

# noise = np.random.normal(loc=0.0, scale=0.2, size=90)

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)
# ax.plot(x_data_plot, y_data_plot + noise, label='real')
ax.plot(x_data_plot, y_data_plot, color='black', label='real')
# ax.plot(x_data_plot, fit_all_plot, label='modelled')

x_axis_label = r'$\Phi$'
y_axis_label = r'$relative \, L_{iso}$'

ax.set_xlabel(x_axis_label, fontsize=24)
ax.set_ylabel(y_axis_label, fontsize=24)

# plt.legend()
plt.show()