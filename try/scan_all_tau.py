import config
import save
import pandas as pd
import numpy as np

# ------------------------------
cols_model_args = ['theta_obs', 'beta_mu', 'mc2', 'a_portion', 'phi_0']

mu = 0.1e31

lin_cos = True

theta_obs_arr = [10 * i for i in range(0, 10)]
if lin_cos:
    theta_obs_arr = [0, 60, 90, 25.84, 36.87, 45.57, 53.13, 66.42, 72.54, 78.46, 84.26]
theta_obs_arr.sort()

beta_mu_arr = [10, 20]
mc2_arr = [30, 60, 100]  # 130
a_portion_arr = np.linspace(0.1, 1, 10)  # 9, endpoint=False
phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]  # 90 # [20 * i for i in range(0, 10)]

tensor_tau_cols_scatter_arr = []
tensor_alpha_cols_scatter_arr = []

idx = 0

for theta_obs in theta_obs_arr:
    tensor_tau_cols_arr = []
    tensor_alpha_cols_arr = []
    model_args_arr = []
    for beta_mu in beta_mu_arr:
        for mc2 in mc2_arr:
            for a_portion in a_portion_arr:
                for phi_0 in phi_0_arr:
                    flag = save.check_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
                    if flag:
                        print(f'{theta_obs=} {beta_mu=} {mc2=} {a_portion=} {phi_0=}')
                        try:

                            model_args_arr.append((theta_obs, beta_mu, mc2, a_portion, phi_0))

                            tensor_tau_cols = save.load_tensor(mu, theta_obs, beta_mu, mc2, a_portion, phi_0,
                                                               'tensor_tau_cols')
                            tensor_tau_cols = tensor_tau_cols[tensor_tau_cols > 0]

                            # print(tensor_tau_cols.shape)
                            cur_str = ','.join(f'{val:.2f}' for val in tensor_tau_cols)

                            # tensor_tau_cols_arr.append(list(np.round(tensor_tau_cols, 2)))
                            tensor_tau_cols_arr.append(cur_str)

                            tensor_alpha_cols = save.load_tensor(mu, theta_obs, beta_mu, mc2, a_portion, phi_0,
                                                                 'tensor_alpha_cols')
                            tensor_alpha_cols = tensor_alpha_cols[tensor_alpha_cols > 0]
                            cur_str = ','.join(f'{val:.2f}' for val in tensor_alpha_cols)
                            tensor_alpha_cols_arr.append(cur_str)

                            # tensor_alpha_cols_arr.append(list(tensor_tau_cols))

                            #
                            #
                            # tensor_tau_cols_scatter = save.load_tensor(mu, theta_obs, beta_mu, mc2, a_portion, phi_0,
                            #                                            'tensor_tau_scatter')
                            # tensor_tau_cols_scatter = tensor_tau_cols_scatter[tensor_tau_cols_scatter > 0]
                            # tensor_tau_cols_scatter_arr.append(tensor_tau_cols_scatter)
                            #
                            # tensor_alpha_cols_scatter = save.load_tensor(mu, theta_obs, beta_mu, mc2, a_portion, phi_0,
                            #                                              'tensor_alpha_scatter')
                            # tensor_alpha_cols_scatter = tensor_alpha_cols_scatter[tensor_alpha_cols_scatter > 0]
                            # tensor_alpha_cols_scatter_arr.append(tensor_alpha_cols_scatter)


                        except FileNotFoundError:
                            pass
                            # print('wtf')

    # Создаём DataFrame, дополняя короткие строки NaN
    # df_tau_cols = pd.DataFrame([pd.Series(row) for row in tensor_tau_cols_arr])
    # df_alpha_cols = pd.DataFrame([pd.Series(row) for row in tensor_alpha_cols_arr])

    # Альтернативный вариант (через словарь)
    # df_tau_cols = pd.DataFrame({i: pd.Series(row) for i, row in enumerate(tensor_tau_cols_arr)}).T
    # df_alpha_cols = pd.DataFrame({i: pd.Series(row) for i, row in enumerate(tensor_alpha_cols_arr)}).T

    # df_tau_cols = pd.DataFrame(tensor_tau_cols_arr, columns=["tau"])
    # df_alpha_cols = pd.DataFrame(tensor_alpha_cols_arr, columns=["alpha"])

    # df = pd.concat([df, df_tau_cols, df_alpha_cols], axis=1)

    # df_tau_cols.to_parquet(f'df_tau_cols_{idx}.parquet')
    # df_alpha_cols.to_parquet(f'df_alpha_cols_{idx}.parquet')
    # df_args.to_parquet(f'df_args_{idx}.parquet')
    #
    # df = pd.concat([df_args, df_tau_cols, df_alpha_cols], axis=1)
    # df.to_parquet(f'df_all_{idx}.parquet')

    df_args = pd.DataFrame(model_args_arr, columns=cols_model_args)
    df_tau_cols = pd.DataFrame({"tau": tensor_tau_cols_arr})
    df_alpha_cols = pd.DataFrame({"alpha": tensor_alpha_cols_arr})

    df_tau_cols.to_parquet(f'df_tau_cols_{idx}.parquet')
    df_alpha_cols.to_parquet(f'df_alpha_cols_{idx}.parquet')
    df_args.to_parquet(f'df_args_{idx}.parquet')

    df = pd.concat([df_args, df_tau_cols, df_alpha_cols], axis=1)
    df.to_parquet(f'df_all_{idx}.parquet')

    idx += 1
