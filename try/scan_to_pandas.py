import save
import pandas as pd
import numpy as np

# ------------------------------------------------------------------
cols_model_args = ['theta_obs', 'beta_mu', 'mc2', 'a_portion', 'phi_0']
cols_L = [f'L_{i}' for i in range(45)]
cols_vals = ['L_x', 'PF', 'R_e', 'ksi_shock']
cols = cols_model_args + cols_vals + cols_L

df = pd.DataFrame(columns=cols)

# ------------------------------------------------------------------
mu = 0.1e31

theta_obs_arr = [10 * i for i in range(0, 10)]
beta_mu_arr = [10 * i for i in range(1, 9)]
# mc2_arr = [10, 20, 40, 50, 70, 80, 90, 110, 120, 130]
mc2_arr = [30, 60, 100]
a_portion_arr = np.linspace(0.1, 1, 19)
phi_0_arr = [20 * i for i in range(0, 10)]

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
                            L_total = save.load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            PF = save.load_PF(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            ksi_shock = save.load_ksi_shock(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            R_e = save.load_R_e(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            L_x = save.load_L_x(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                            model_args = (theta_obs, beta_mu, mc2, a_portion, phi_0)

                            df.loc[pd_index, cols_L] = L_total
                            df.loc[pd_index, cols_model_args] = model_args
                            df.loc[pd_index, cols_vals] = (L_x, PF, R_e, ksi_shock) # 'L_x ', 'PF', 'R_e', 'ksi_shock'

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

df.to_parquet('data.parquet')
