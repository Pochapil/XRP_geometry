from numpy import pi
import numpy as np

# Parameters

# глобальные постоянные
M_Sun = 1.9891e33  # масса молнца [г]
G = 6.67e-8  # гравитационная постоянная [см3·с−2·г−1]
c = 2.99792458e10  # скорость света [см/с]
sigm_Stf_Bolc = 5.67e-5  # постоянная Стефана Больцмана в сгс
a_rad_const = 7.5657e-15  # радиационная константа p=aT**4 [эрг см-3 К-4] постоянная излучения - Плотность энергии и давление равновесного излучения
sigma_T = 6.652e-25  # сечение томсона [см-2]
mass_P = 1.67e-24  # масса протона [г]
h_plank_ergs = 6.62607015e-27  # постоянная Планка в [эрг * с]
h_plank_evs = 4.135667669e-15  # постоянная Планка в [эв * с]
k_bolc = 1.380649e-16  # постоянная Больцмана [эрг/К]

# параметры НЗ
M_ns = 1.4 * M_Sun  # масса нз [г]
R_ns = 1e6  # радиус нз [см]
# H = 2 * 10 ** 13  # магнитное поле стр 19 над формулой 37
# mu магнитный момент [Гаусс * см3] H = 2 * mu / R_ns ** 3
# p_spin = 3.62  # период вращения, [с]

# параметры аккреционного потока
dRe_div_Re = 0.25  # взял просто число
# M_accretion_rate = 10 ** 38 * R_ns / G / MSun  # темп аккреции
ksi_rad = 3 / 2
ksi_param = 0.5  # между 1 и 2 формулой в статье - размер магнитосферы
k = 0.35  # opacity непрозрачность [см**2 / г]

new_magnet_lines_flag = True
tau_flag = True
NS_shadow_flag = True

FLAG_PHI_0_OLD = False
FLAG_R_E_OLD = False
outer_R_e_ksi_flag = False  # if true = обрезаем по ksi; False - по тета (тета iner == тета outer)
ASYNC_FLAG = True
old_path_flag = False

L_nu_flag = False

tau_cutoff = 0
opacity_above_shock = 0  # непрозрачность вещества над ударной волной: 0 - полностью прозрачное, 1 - непрозрачное
# L_ed = M_ns / MSun * 10 ** 38
L_edd = 4 * pi * G * M_ns * c / k

# таблица 1

# a - в азимутальном направлении поток занимает фиксированную долю a полного круга 2πR sinθ

# верхний предел по phi
# нижний предел по phi - phi_0 !!!!
# нижний предел по phi

# цикл для поворотов, сколько точек на графике интегралов - для фазы от 0 до 2 - с перекрытием чтобы форму макс

N_cpus = 10

N_phase = 45  # цикл для поворотов, сколько точек для фазы от 0 до 1 (полного поворота)
omega_ns_deg = 360 / 45  # скорость вращения НЗ - будет меняться только угол phi_mu!
omega_ns_rad = np.deg2rad(omega_ns_deg)
N_phase_for_plot = 2 * N_phase  # сколько точек на графике интегралов - для фазы от 0 до 2 - с перекрытием чтобы форму макс

# количество шагов
N_phi_accretion = 100
N_theta_accretion = 100
N_wavelength_range = 10
N_frequency_range = 100

phase_for_plot = np.linspace(0, 2, N_phase_for_plot)

N_energy = 20
energy_min = 1  # [КэВ]
energy_max = 40  # [КэВ]

# лог сетка по энергии: e_n+1 = e_n * const
energy_step = (energy_max / energy_min) ** (1 / (N_energy - 1))
energy_arr = list(energy_min * energy_step ** i for i in range(N_energy - 1))
# чтобы убрать погрешности и закрыть массив точным числом
energy_arr.append(energy_max)
# energy_arr = np.array(energy_arr)

# obs_i_angle_deg угол между нормалью к двойной системе и наблюдателем

# угол между осью вращения системы и собственным вращением НЗ (берем ось z сонаправленно с осью собств. вращения omega)
betta_rotate = np.deg2rad(0)
phi_rotate = np.deg2rad(0)
# угол между собственным вращением НЗ и магнитной осью

phi_mu_0 = np.deg2rad(0)

# для pretty графиков - индекс по энергии и сколько subfigures
N_column_plot = 5
energy_indexes = [0, 12, 15, 17, 19]
# меньше на 1 т.к. в диапазоне энергий => меньше на 1 чем точек разбиения
energy_indexes_luminosity = [0, 12, 14, 16, 18]

fi_0_dict = {0.11: 340, 0.165: 330, 0.22: 320, 0.275: 310, 0.33: 300, 0.385: 290, 0.44: 280, 0.5: 270, 0.55: 260,
             0.605: 250, 0.66: 240, 0.715: 230, 0.77: 220, 0.825: 210, 0.25: 320, 0.65: 240, 1: 0}

# val_ksi = 100
# val_n = 1e15
# # val = 2 ** (1 / 2) * (G * M_ns) ** (3 / 2) / (val_ksi * R_ns ** (7 / 2) * M_accretion_rate)
# val = G * M_ns * M_accretion_rate / ((val_ksi * R_ns) ** 3 * val_n)
# print(val)
# print(H/10**11)
