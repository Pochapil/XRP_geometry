import config
import save
import pandas as pd
import numpy as np

# ------------------------------------------------------------------
cols_model_args = ['theta_obs', 'beta_mu', 'mc2', 'a_portion', 'phi_0']
cols_L = [f'L_{i}' for i in range(config.N_phase)]
cols_vals = ['L_x', 'PF', 'R_e', 'ksi_shock', 'beta', 'gamma']
cols = cols_model_args + cols_vals + cols_L

df = pd.DataFrame(columns=cols)

# ------------------------------------------------------------------
mu = 0.1e31

theta_obs_arr = [10 * i for i in range(0, 10)]
beta_mu_arr = [10, 20]
# mc2_arr = [10, 20, 40, 50, 70, 80, 90, 110, 120, 130]
mc2_arr = [30, 60, 100]
a_portion_arr = np.linspace(0.1, 1, 9, endpoint=False)
phi_0_arr = [20 * i for i in range(0, 10)]
phi_0_arr = [0, 40, 90]  # 140, 180

L_total_arr = []
PF_arr = []
ksi_shock_arr = []
R_e_arr = []
L_x_arr = []
beta_arr = []
gamma_arr = []

model_args_arr = []

pd_index = 0
for theta_obs in theta_obs_arr:
    for beta_mu in beta_mu_arr:
        for mc2 in mc2_arr:
            for a_portion in a_portion_arr:
                for phi_0 in phi_0_arr:
                    flag = save.check_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                    if flag:
                        try:
                            # print(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            # print(flag)
                            # print()
                            L_total_arr.append(save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0))

                            PF_arr.append(save.load_PF(mu, theta_obs, beta_mu, mc2, a_portion, phi_0))
                            ksi_shock_arr.append(save.load_ksi_shock(mu, theta_obs, beta_mu, mc2, a_portion, phi_0))
                            R_e_arr.append(save.load_R_e(mu, theta_obs, beta_mu, mc2, a_portion, phi_0))
                            L_x_arr.append(save.load_L_x(mu, theta_obs, beta_mu, mc2, a_portion, phi_0))

                            beta_arr.append(save.load_beta(mu, theta_obs, beta_mu, mc2, a_portion, phi_0))
                            gamma_arr.append(save.load_gamma(mu, theta_obs, beta_mu, mc2, a_portion, phi_0))

                            model_args_arr.append((theta_obs, beta_mu, mc2, a_portion, phi_0))

                            # L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            # PF = save.load_PF(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            # ksi_shock = save.load_ksi_shock(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            # R_e = save.load_R_e(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            # L_x = save.load_L_x(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            # model_args = (theta_obs, beta_mu, mc2, a_portion, phi_0)

                            # df.loc[pd_index, cols_L] = L_total
                            # df.loc[pd_index, cols_model_args] = model_args
                            # df.loc[pd_index, cols_vals] = (L_x, PF, R_e, ksi_shock)  # 'L_x ', 'PF', 'R_e', 'ksi_shock'

                            # L_total = L_total / np.max(L_total)
                            # print(L_total)
                            # print(flag)
                            # print()
                            pd_index += 1
                        except FileNotFoundError:
                            pass
                            # print('ERROR')
                            # print(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            # print('NOT EXIST')

df1 = pd.DataFrame(model_args_arr, columns=cols_model_args)
# df1 = pd.DataFrame(dict(zip(cols_model_args, model_args_arr)))
df2 = pd.DataFrame(dict(zip(cols_vals, (L_x_arr, PF_arr, R_e_arr, ksi_shock_arr, beta_arr, gamma_arr))))  #
# df2 = pd.DataFrame((L_x_arr, PF_arr, R_e_arr, ksi_shock_arr, beta_arr, gamma_arr), columns=cols_vals) #
df3 = pd.DataFrame(L_total_arr, columns=cols_L)
# df3 = pd.DataFrame(dict(zip(cols_L, L_total_arr))) #

df = pd.concat([df1, df2, df3], axis=1)

# dict_args = dict(zip(cols_model_args, model_args_arr))
# dict_vals = dict(zip(cols_vals, (L_x_arr, PF_arr, R_e_arr, ksi_shock_arr, beta_arr, gamma_arr)))
# dict_L = dict(zip(cols_L, L_total_arr))
#
# full_dict = {}
# full_dict.update(dict_args)
# full_dict.update(dict_vals)
# full_dict.update(dict_L)

# merged_dict = dict_args | dict_vals | dict_L

# df = pd.DataFrame(full_dict)

# df = pd.DataFrame({
#     cols_model_args: model_args_arr,
#     cols_vals: (L_x_arr, PF_arr, R_e_arr, ksi_shock_arr, beta_arr, gamma_arr),
#     cols_L: L_total_arr,
# })

df.to_parquet('data.parquet')
