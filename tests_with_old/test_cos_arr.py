'''normals'''
# array_normal = np.zeros((phi_arr.shape[0], theta_arr.shape[0]), dtype=object)
# count = 1
# for i in range(phi_arr.shape[0]):
#     for j in range(theta_arr.shape[0]):
#         # array_normal.append(coefficient * matrix.newE_n(phi_range[i], theta_range[j]))
#         # array_normal[i, j] = 1 * matrix.newE_r(phi_arr[i], theta_arr[j])
#         array_normal[i, j] = 1 * matrix.newE_n(phi_arr[i], theta_arr[j])
#         print(array_normal[i, j] - z[i, j, :])
#         print(count)
#         count += 1
# # print(array_normal)
# # print(array_normal


'''cos'''
# cos_psi_range = np.empty([config.N_phi_accretion, config.N_theta_accretion])
# for i in range(config.N_phi_accretion):
#     for j in range(config.N_theta_accretion):
#         # умножать на N_theta
#         # cos_psi_range[i, j] = np.dot(e_obs_mu, self.array_normal[i * config.N_theta_accretion + j])
#         cos_psi_range[i, j] = np.dot(obs_matrix[6], surf.array_normal[i][j])
#
# print(rotation_obs_matrix - cos_psi_range)

'''in cycle! all empty'''
# for k in range(config.N_phase):
#     cos_psi_range = np.empty([config.N_phi_accretion, config.N_theta_accretion])
#     for i in range(config.N_phi_accretion):
#         for j in range(config.N_theta_accretion):
#             # умножать на N_theta
#             # cos_psi_range[i, j] = np.dot(e_obs_mu, self.array_normal[i * config.N_theta_accretion + j])
#             cos_psi_range[i, j] = np.dot(obs_matrix[k], surf.array_normal[i][j])
#
#     print(rotation_obs_matrix[k][(rotation_obs_matrix[k] - cos_psi_range) > 0])