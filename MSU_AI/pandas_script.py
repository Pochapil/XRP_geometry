import pandas as pd
import numpy as np
from pathlib import Path

import config
import main_service

i_angle = [10 * i for i in range(0, 10)]
betta_mu = [10 * i for i in range(0, 10)]
a_portion = [0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 1]

mc2 = np.linspace(10, 130, 13)
mc2 = list(map(int, mc2))

fi_0 = [0, 20, 40, 60, 80, 90, 100, 120, 140, 160, 180]

old_fi_0 = fi_0.copy()

df = pd.DataFrame()

data = []

for i_angle_index in range(len(i_angle)):
    for betta_mu_index in range(len(betta_mu)):
        for i in range(len(mc2)):
            for j in range(len(a_portion)):
                for k in range(len(fi_0)):
                    # fi_0[k] = 360 - 180 * a_portion[j]

                    # fi_0[k] = (config.fi_0_dict[a_portion[j]] + 40) % 360
                    fi_0[k] = (config.fi_0_dict[a_portion[j]] + old_fi_0[k]) % 360

                    # fi_0[k] = (config.fi_0_dict[a_portion[j]] + 20 * (k)) % 360

                    config.set_e_obs(i_angle[i_angle_index], 0)
                    config.set_betta_mu(betta_mu[betta_mu_index])

                    config.M_rate_c2_Led = mc2[i]
                    config.a_portion = a_portion[j]
                    config.phi_accretion_begin_deg = fi_0[k]

                    config.update()

                    s = Path(config.full_file_folder)
                    print(s)
                    if s.exists():
                        if Path(config.full_file_folder + "total_luminosity_of_surfaces.txt").exists() and Path(config.full_file_folder + 'scattered_on_magnet_lines/' + 'scattered_energy_bot.txt').exists():
                            file_name = "total_luminosity_of_surfaces.txt"
                            L_arr = main_service.load_arr_from_txt(config.full_file_folder, file_name)[4]

                            buf = main_service.load_arr_from_txt(config.full_file_folder + 'scattered_on_magnet_lines/',
                                                                 'scattered_energy_top.txt')
                            if not np.isnan(buf).any():
                                L_arr += buf

                            buf = main_service.load_arr_from_txt(config.full_file_folder + 'scattered_on_magnet_lines/',
                                                                 'scattered_energy_bot.txt')
                            if not np.isnan(buf).any():
                                L_arr += buf

                            L_avg = np.mean(L_arr)

                            PF = (np.max(L_arr) - np.min(L_arr)) / (np.max(L_arr) + np.min(L_arr))
                            # y_data = np.roll(y_data, -max_idx)
                            data.append(
                                (i_angle[i_angle_index], betta_mu[betta_mu_index], mc2[i], a_portion[j], fi_0[k],
                                 L_arr, PF))
df = pd.DataFrame(data)
main_service.create_file_path('pandas/')
df.to_csv('pandas/my_data')
