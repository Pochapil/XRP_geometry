import numpy as np

from geometry import matrix
import integralsService

# phi x theta x 3
#


phi_arr = np.array([0, np.deg2rad(30), np.deg2rad(60), np.deg2rad(90), np.deg2rad(45), np.deg2rad(30)])
theta_arr = np.array([0, np.deg2rad(30), np.deg2rad(60), np.deg2rad(90)])

z = matrix.newE_fi_n(phi_arr)
print(z)

z = matrix.newE_n_n(phi_arr, theta_arr)
print(z)

array_normal = np.zeros((phi_arr.shape[0], theta_arr.shape[0]), dtype=object)
count = 1
for i in range(phi_arr.shape[0]):
    for j in range(theta_arr.shape[0]):
        # array_normal.append(coefficient * matrix.newE_n(phi_range[i], theta_range[j]))
        # array_normal[i, j] = 1 * matrix.newE_r(phi_arr[i], theta_arr[j])
        array_normal[i, j] = 1 * matrix.newE_n(phi_arr[i], theta_arr[j])
        print(array_normal[i, j] - z[i, j, :])
        print(count)
        count += 1

print(np.pi - 1.382)
print(np.pi + 1.382)

# x = np.

z1 = matrix.get_xyz_coord_angles(1, phi_arr, theta_arr)

print(z1)
print(z1.shape)

# L = np.array([[1, 2, 4, 12], [4, 2, 9, 13]])
# print(L.shape)
# z = integralsService.get_PF(L)


phi_range = np.linspace(1, 10, 10)
print('hello')
print(phi_range.shape[0])
phi_range_plus = phi_range + 2 * np.pi
phi_range_minus = phi_range - 2 * np.pi
phi_range = np.hstack([phi_range, phi_range_plus, phi_range_minus])
print(phi_range)

print(np.linspace(0.1,1,19))



a = np.ones((2,3,4))

print(a)

a[0, 2, 1] = 10
a[0, 1, 2] = 10
a[1, 1, 2] = 10
a[1, 1, 3] = 10

print(a[a > 1].shape)
# print(array_normal)
# print(array_normal


import numpy as np

array = np.zeros((5, 5))
k = -1

rows, cols = np.indices(array.shape)
row_values = np.diag(rows, k=k)
col_values = np.diag(cols, k=k)
array[row_values, col_values] = 1

print(array)


import config
delta_R = 1
C1 = 1.4e-27 * 1.2 / (config.mass_P)**2
print(f'{C1=}')

C2 = (4 * np.pi)**2 * 2 * config.G * config.M_ns / (config.k**2 * np.pi * delta_R)
print(f'{C2=}')

T = 3000 * 1.6e-12/config.k_bolc # эрг


print(f'{C1 * C2}')




mc2 = 60
M_accretion_rate = mc2 * config.L_edd / config.c ** 2
mu = 0.1e31
R_m = (mu ** 2 / (2 * M_accretion_rate * (2 * config.G * config.M_ns) ** (1 / 2))) ** (2 / 7)

C2 = 1 / (np.pi * delta_R)
C3 = (2  * (2 * config.G * config.M_ns) ** (1 / 2)) ** (6 / 7) / ((1/2) ** 3)

print(C1 * C2 * C3)

M_accretion_rate = 1e18
mu = 0.1e31
T = 3000 * 1.6e-12/config.k_bolc # эрг
print('ans:')
print(C1 * C2 * C3 * M_accretion_rate ** (20/7) / mu ** (12/7) * T**(1/2))
print(T)

print(2 * mu / config.R_ns ** 3)


print('final')
C1 = 1.4e-27 * 1.2
C2 = C1 / (config.mass_P)**2
C3 = C2 * (np.pi / (4 * np.pi * np.sqrt(2 * config.G * config.M_ns))**2)
print(C3)
print(C3 * M_accretion_rate**2 * T**(1/2))


print('total')
C1 = 1.4e-27 * 1.2
C2 = C1 * (np.pi / (config.k * config.mass_P)**2)
print(C2)
tau=10
R=1e8
print(C2 * tau**2 * T**(1/2)*R)

tau=80
R=0.27e8
print(C2 * tau**2 * T**(1/2)*R)