import plot_many_characteristics

mu = 0.1e31


def plot_sky_map_as_in_art():
    beta_mu = 20
    mc2 = 60
    a_portion = 0.22
    phi_0 = 0
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.44
    phi_0 = 0
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.66
    phi_0 = 0
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.22
    phi_0 = 20
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.44
    phi_0 = 20
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.66
    phi_0 = 20
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.22
    phi_0 = 40
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.44
    phi_0 = 40
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.66
    phi_0 = 40
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.22
    phi_0 = 60
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.44
    phi_0 = 60
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)

    beta_mu = 20
    mc2 = 60
    a_portion = 0.66
    phi_0 = 60
    plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)


def plot_L_to_mc2_as_in_art():
    theta_obs = 20
    beta_mu = 60
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.44
    phi_0 = 0
    plot_many_characteristics.plot_L_to_mc2(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    theta_obs = 40
    beta_mu = 60
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.44
    phi_0 = 0
    plot_many_characteristics.plot_L_to_mc2(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    theta_obs = 20
    beta_mu = 60
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.66
    phi_0 = 0
    plot_many_characteristics.plot_L_to_mc2(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    theta_obs = 40
    beta_mu = 60
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.66
    phi_0 = 0
    plot_many_characteristics.plot_L_to_mc2(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)


def plot_L_to_a_portion_as_in_art():
    theta_obs = 20
    beta_mu = 40
    mc2 = 100
    a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    phi_0 = 0
    plot_many_characteristics.plot_L_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)

    theta_obs = 20
    beta_mu = 60
    mc2 = 30
    a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    phi_0 = 0
    plot_many_characteristics.plot_L_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)

    theta_obs = 40
    beta_mu = 20
    mc2 = 30
    a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    phi_0 = 0
    plot_many_characteristics.plot_L_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)

    theta_obs = 40
    beta_mu = 20
    mc2 = 100
    a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    phi_0 = 0
    plot_many_characteristics.plot_L_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)

    theta_obs = 60
    beta_mu = 20
    mc2 = 30
    a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    phi_0 = 0
    plot_many_characteristics.plot_L_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)

    theta_obs = 60
    beta_mu = 20
    mc2 = 100
    a_portion_arr = [0.165, 0.22, 0.275, 0.33, 0.385, 0.44, 0.5, 0.55, 0.605, 0.66, 0.715, 0.77, 0.825, 1]
    phi_0 = 0
    plot_many_characteristics.plot_L_to_a_portion(mu, theta_obs, beta_mu, mc2, a_portion_arr, phi_0)


def plot_L_to_phi_0_as_in_art():
    theta_obs = 40
    beta_mu = 60
    mc2 = 30
    a_portion = 0.22
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    plot_many_characteristics.plot_L_to_phi_0(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr, True)

    # plt.clf()
    # plt.cla()
    # plt.close()

    theta_obs = 40
    beta_mu = 60
    mc2 = 30
    a_portion = 0.66
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    plot_many_characteristics.plot_L_to_phi_0(mu, theta_obs, beta_mu, mc2, a_portion, phi_0_arr, True)


def plot_masses_PF_L_nu_as_in_art():
    theta_obs = 20
    beta_mu = 40
    mc2_arr = [30, 100]
    a_portion_arr = [0.22, 0.66]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    plot_many_characteristics.plot_masses_PF_L_nu(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)
    plot_many_characteristics.plot_masses_PF_L(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)

    theta_obs = 20
    beta_mu = 60
    mc2_arr = [30, 100]
    a_portion_arr = [0.22, 0.66]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    plot_many_characteristics.plot_masses_PF_L_nu(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)
    plot_many_characteristics.plot_masses_PF_L(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)

    theta_obs = 40
    beta_mu = 40
    mc2_arr = [30, 100]
    a_portion_arr = [0.22, 0.66]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    plot_many_characteristics.plot_masses_PF_L_nu(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)
    plot_many_characteristics.plot_masses_PF_L(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)

    theta_obs = 40
    beta_mu = 60
    mc2_arr = [30, 100]
    a_portion_arr = [0.22, 0.66]
    phi_0_arr = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]
    plot_many_characteristics.plot_masses_PF_L_nu(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)
    plot_many_characteristics.plot_masses_PF_L(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0_arr)


def plot_PF_contour_as_in_art():
    mc2 = 30
    a_portion = 1
    phi_0 = 0
    plot_many_characteristics.plot_PF_contour(mu, mc2, a_portion, phi_0)
    plot_many_characteristics.plot_PF_to_chi_theta(mu, mc2, a_portion, phi_0)

    mc2 = 100
    a_portion = 1
    phi_0 = 0
    plot_many_characteristics.plot_PF_contour(mu, mc2, a_portion, phi_0)
    plot_many_characteristics.plot_PF_to_chi_theta(mu, mc2, a_portion, phi_0)


def plot_table_as_in_art():
    theta_obs = 40
    beta_mu = 20
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.44
    phi_0 = 0
    plot_many_characteristics.plot_L_iso_to_m(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    theta_obs = 40
    beta_mu = 20
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    a_portion = 0.66
    phi_0 = 0
    plot_many_characteristics.plot_table(mu, theta_obs, beta_mu, mc2_arr, a_portion, phi_0)

    theta_obs = 40
    beta_mu = 20
    mc2_arr = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130]
    phi_0 = 0
    a_portion_arr = [0.22, 0.44, 0.66, 1]
    plot_many_characteristics.plot_table_together(mu, theta_obs, beta_mu, mc2_arr, a_portion_arr, phi_0)

if __name__ == "__main__":
    beta_mu = 20
    mc2 = 30
    a_portion_arr = [0.22, 0.66]
    phi_0_arr = [20, 40, 60]
    for a_portion in a_portion_arr:
        for phi_0 in phi_0_arr:
            pass
            # plot_many_characteristics.plot_sky_map(mu, beta_mu, mc2, a_portion, phi_0)


    # plot_sky_map_as_in_art()
    # plot_sky_map_as_in_art()
    # plot_L_to_phi_0_as_in_art()
    # plot_L_to_mc2_as_in_art()
    # plot_L_to_a_portion_as_in_art()
    # plot_masses_PF_L_nu_as_in_art()
    # plot_PF_contour_as_in_art()
    # plot_table_as_in_art()