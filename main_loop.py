import main

if __name__ == '__main__':
    mu = 0.1e31

    beta_mu_arr = [10 * i for i in range(1, 10)]
    theta_obs_arr = [10 * i for i in range(1, 10)]
    mc2_arr = [30, 60, 100]
    a_portion_arr = [0.22, 0.44, 0.66, 1]
    phi_0_arr = [0]

    theta_obs_arr = [40]
    beta_mu_arr = [60]
    # [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    mc2_arr = [10, 20, 40, 50, 60, 70, 80, 90, 110, 120, 130]
    a_portion_arr = [0.44, 0.66]
    phi_0_arr = [0]

    flag = False
    # ------------------------------------------------- start -------------------------------------------------------
    N_big = len(theta_obs_arr) * len(beta_mu_arr) * len(mc2_arr) * len(a_portion_arr) * len(phi_0_arr)
    print(f'to calculate {N_big} loops need about {30 * N_big / 3600} hours')
    for theta_obs in theta_obs_arr:
        for beta_mu in beta_mu_arr:
            for mc2 in mc2_arr:
                for a_portion in a_portion_arr:
                    for phi_0 in phi_0_arr:
                        main.calc_and_save_for_configuration(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, flag)
