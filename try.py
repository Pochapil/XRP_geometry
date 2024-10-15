import numpy as np

from geometricTask import matrix

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
# print(array_normal)
# print(array_normal