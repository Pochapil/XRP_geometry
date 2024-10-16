'''normals'''
# array_normal = np.zeros((phi_arr.shape[0], theta_arr.shape[0]), dtype=object)
# count = 1
# for phi_index in range(phi_arr.shape[0]):
#     for theta_index in range(theta_arr.shape[0]):
#         # array_normal.append(coefficient * matrix.newE_n(phi_range[phi_index], theta_range[theta_index]))
#         # array_normal[phi_index, theta_index] = 1 * matrix.newE_r(phi_arr[phi_index], theta_arr[theta_index])
#         array_normal[phi_index, theta_index] = 1 * matrix.newE_n(phi_arr[phi_index], theta_arr[theta_index])
#         print(array_normal[phi_index, theta_index] - z[phi_index, theta_index, :])
#         print(count)
#         count += 1
# # print(array_normal)
# # print(array_normal


'''cos'''
# cos_psi_range = np.empty([config.N_phi_accretion, config.N_theta_accretion])
# for phi_index in range(config.N_phi_accretion):
#     for theta_index in range(config.N_theta_accretion):
#         # умножать на N_theta
#         # cos_psi_range[phi_index, theta_index] = np.dot(e_obs_mu, self.array_normal[phi_index * config.N_theta_accretion + theta_index])
#         cos_psi_range[phi_index, theta_index] = np.dot(obs_matrix[6], surf.array_normal[phi_index][theta_index])
#
# print(rotation_obs_matrix - cos_psi_range)

'''in cycle! all empty'''
# for phase_index in range(config.N_phase):
#     cos_psi_range = np.empty([config.N_phi_accretion, config.N_theta_accretion])
#     for phi_index in range(config.N_phi_accretion):
#         for theta_index in range(config.N_theta_accretion):
#             # умножать на N_theta
#             # cos_psi_range[phi_index, theta_index] = np.dot(e_obs_mu, self.array_normal[phi_index * config.N_theta_accretion + theta_index])
#             cos_psi_range[phi_index, theta_index] = np.dot(obs_matrix[phase_index], surf.array_normal[phi_index][theta_index])
#
#     print(rotation_obs_matrix[phase_index][(rotation_obs_matrix[phase_index] - cos_psi_range) > 0])
