# from tkinter import *
import tkinter as Tk
import matplotlib
import matplotlib.pyplot as plt

import pathService
import save

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np

import config

mu = 0.1e31

beta_mu_deg = 60
beta_mu = np.deg2rad(beta_mu_deg)
mc2 = 10
a_portion = 1
phi_0 = 0

colors = {'magnet': 'blue', 'disc': 'brown', 'alfv': 'black'}

phi_range = np.linspace(np.deg2rad(phi_0), 2 * np.pi, 200)
theta_range = np.zeros_like(phi_range)

for i in range(len(phi_range)):
    theta_range[i] = np.pi / 2 - np.arctan(np.tan(beta_mu) * np.cos(phi_range[i]))
    # theta_range[i] = np.pi / 2 - betta_mu * np.cos(phi_range[i])

mc2 = 30
M_accretion_rate = mc2 * config.L_edd / config.c ** 2

R_alfven = (mu ** 2 / (2 * M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)
R_alfven = R_alfven / config.R_ns

# R_e = R_alfven * config.ksi_param
R_disk = config.ksi_param * R_alfven

disk_min_theta_angle = np.min(np.pi / 2 - np.arctan(np.tan((beta_mu)) * np.cos(phi_range)))
idx = np.argmin(np.pi / 2 - np.arctan(np.tan((beta_mu)) * np.cos(phi_range)))
disk_min_phi_angle = phi_range[idx]
# print(f'{self.disk_min_theta_angle =}')

# радиус на котором лежат внутренние дипольные линии
R_e = R_disk / ((np.sin(disk_min_theta_angle)) ** 2)
# R_e = R_alfven * config.ksi_param

fig, ax = plt.subplot_mosaic('a', figsize=(8, 6))

# внутренний радиус колонки
r = R_e * np.sin(theta_range) ** 2
r1 = r
# r1 = r * np.sin(theta_range)
x = r1 * np.cos(phi_range)
y = r1 * np.sin(phi_range)
ax['a'].plot(x, y, colors['magnet'], label='inner R')

# внешний радиус колонки
dRe_div_Re = 0.25
r = (1 + dRe_div_Re) * R_e * np.sin(theta_range) ** 2
r1 = r
# r1 = r * np.sin(theta_range)
x = r1 * np.cos(phi_range)
y = r1 * np.sin(phi_range)
ax['a'].plot(x, y, colors['magnet'], label='outer R')

# круг R_e
x = R_disk * np.cos(phi_range)
y = R_disk * np.sin(phi_range)
ax['a'].plot(x, y, colors['disc'], label='R_m circle')

x = 1.25 * R_disk * np.cos(phi_range)
y = 1.25 * R_disk * np.sin(phi_range)
ax['a'].plot(x, y, colors['disc'], label='1.25 R_m circle')

# попытка менять R_e - 69 стр - там ksi
'''  The process also requires selecting a
value for the normalization constant ξ appearing in Equation (22), such that the minimum
value of the Alfv´en radius equals the maximum value of R2, disk, corresponding to α = 0'''

# H = 2 * mu / config.R_ns ** 3
# B_disk = H / 2 * (3 * np.cos(theta_range) ** 2 + 1) ** (1 / 2)
r = R_disk * ((3 * np.cos(theta_range) ** 2 + 1)) ** (2 / 7)
ksi = R_disk / np.max(r)
r = r * ksi

# c = 1
# r = ((3 * np.cos(theta_range) ** 2 + 1) / c) ** (2 / 7)
# ksi = max(R_disk * np.sin((theta_range)) ** 2) / min(r)
# r *= ksi
# r *= R_disk

# ksi = 2 * max(np.sin((theta_range)) ** 2)
# print(ksi)
# r = ksi * R_e * 1 / 2 * (1 + 3 * np.cos(theta_range) ** 2) ** (1 / 2) / (R_e * np.sin(theta_range) ** 3)
rotate_angle = 0
rotate_angle = np.deg2rad(rotate_angle)
x = r * np.cos(phi_range + rotate_angle)
y = r * np.sin(phi_range + rotate_angle)
# eq22 https://arxiv.org/pdf/1612.02411
ax['a'].plot(x, y, colors['alfv'], label='fake R_alfv')

# вектор магнитного поля
origin = np.array([0, 0])
vector = np.array([R_disk * 2.25, 0])
# plt.quiver(origin, vector, angles='xy', scale_units='xy',color='black', scale=10)
ax['a'].arrow(0, 0, *vector, head_width=0.05, head_length=1, color='black')

ax['a'].legend()
ax['a'].set_xlabel('r/R_ns', fontsize=24)
ax['a'].set_ylabel('r/R_ns', fontsize=24)
ax['a'].axis('equal')

# plt.title(r'$\beta_{\mu}$' + f'={beta_mu_deg}')
# plt.axis('equal')


# x_val, y_val = X[i, j], Y[i, j]
#         if x_val >= 0 and y_val >= 0 and x_val <= 6 and y_val <= 6:  # Ограничение по границам сетки
#             if x_val**2 + y_val**2 >= 25:  # Точки выше окружности (внешняя область)
#                 plt.fill_between([x_val, x_val + 1], [y_val, y_val], [y_val + 1, y_val + 1],
#                                 color='lightgray', alpha=0.3)






prefix_folder = '2d/'
save_dir = pathService.get_dir(mu=mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                               phi_0=None, prefix_folder=prefix_folder)
file_name = f'{beta_mu_deg}'
save.save_figure(fig, save_dir, file_name)

# plt.axis('scaled')
# ax.set_ylim(-R_disk*1.25, 1.25*R_disk)
# ax.set_xlim(-R_disk*1.25, 1.25*R_disk)
# plt.show()

crit_betta_mu = np.arctan((1 / 12) ** (1 / 2))
crit_betta_mu = np.rad2deg(crit_betta_mu)
betta_mu_arr = np.array([crit_betta_mu + 5 * i for i in range(18)])

phi_arr = np.arccos(1 / (12 * np.tan(np.deg2rad(betta_mu_arr)) ** 2))
phi_arr = np.rad2deg(phi_arr)
a_arr = 2 * phi_arr / 360

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
ax.plot(betta_mu_arr, a_arr)

ax.set_xlabel(r'$\beta_{\mu}$', fontsize=24)
ax.set_ylabel('a', fontsize=24)

# plt.show()
