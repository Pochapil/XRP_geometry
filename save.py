import pathlib
import numpy as np
import matplotlib.pyplot as plt

import newService
import pathService
import config


def create_file_path_from_str(file_path_str):
    # parents - создает пропущенные папки если их нет, exist_ok=True - не выдает ошибку если существует папка
    pathlib.Path(file_path_str).mkdir(parents=True, exist_ok=True)

    # p = pathlib.Path("temp/")
    # p.mkdir(parents=True, exist_ok=True)
    # fn = "test.txt"  # I don't know what is your fn
    # filepath = p / fn


def create_file_path(file_path):
    # parents - создает пропущенные папки если их нет, exist_ok=True - не выдает ошибку если существует папка
    file_path.mkdir(parents=True, exist_ok=True)


def save_arr_as_txt(arr, file_folder, file_name):
    full_file_name = file_folder / pathlib.Path(file_name)
    create_file_path(file_folder)
    np.savetxt(full_file_name, arr)  # encoding='utf-8'


def save_arr_as_txt_from_str(arr, file_folder, file_name):
    full_file_name = file_folder + file_name
    create_file_path(file_folder)
    np.savetxt(full_file_name, arr)  # encoding='utf-8'


def load_arr_from_txt(file_folder, file_name):
    return np.loadtxt(file_folder / file_name)


def load_arr_from_txt_from_str(file_folder, file_name):
    return np.loadtxt(file_folder + file_name)


def save_figure(fig, file_path, file_name):
    create_file_path(file_path)
    full_file_name = file_path / (file_name + '.png')
    fig.savefig(full_file_name, dpi=fig.dpi)
    # fig.savefig(full_file_name, dpi=200)
    # plt.clf()
    plt.close(fig)


def save_figure_txt(fig, file_path, file_name):
    create_file_path(file_path)
    full_file_name = file_path + file_name
    fig.savefig(full_file_name, dpi=fig.dpi)
    # fig.savefig(full_file_name, dpi=200)
    plt.close(fig)


def get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    cur_dir_saved = pathService.PathSaver(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    return cur_dir_saved.get_path() / 'txt/'


def load_L_x(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    with open(data_folder / 'save_values.txt') as f:
        lines = f.readlines()
    return float(lines[3][12:20]) * 10 ** float(lines[3][27:29])


def load_R_e(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    with open(data_folder / 'save_values.txt') as f:
        lines = f.readlines()
    return float(lines[0][6:-1])


def load_ksi_shock(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    with open(data_folder / 'save_values.txt') as f:
        lines = f.readlines()
    return float(lines[1][12:-1])


def load_L_surfs(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    folder = 'surfs/'
    file_name = "L_surfs.txt"
    L_surfs = load_arr_from_txt(data_folder / folder, file_name)
    L_surfs[np.isnan(L_surfs)] = 0
    return L_surfs


def load_L_scatter(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    folder = 'scatter/'
    file_name = "L_scatter.txt"
    L_scatter = load_arr_from_txt(data_folder / folder, file_name)
    L_scatter[np.isnan(L_scatter)] = 0
    return L_scatter


def load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    '''returns array len config.N_phase'''
    # sum of surfs
    L_surfs = load_L_surfs(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    # if not np.isnan(buf).any():
    L_scatter = load_L_scatter(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    return np.sum(L_surfs, axis=0) + np.sum(L_scatter, axis=0)


def load_L_nu_surfs(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, energy):
    '''для конкретной энергии'''
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    folder = 'surfs/'
    file_name = f'L_nu_surfs_of_energy_{energy:.2f}_KeV.txt'
    L_nu_surfs = load_arr_from_txt(data_folder / folder, file_name)
    L_nu_surfs[np.isnan(L_nu_surfs)] = 0
    return L_nu_surfs


def load_L_nu_scatter(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, energy):
    '''для конкретной энергии'''
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    folder = 'scatter/'
    file_name = f'L_nu_scatter_of_energy_{energy:.2f}_KeV.txt'
    L_nu_scatter = load_arr_from_txt(data_folder / folder, file_name)
    L_nu_scatter[np.isnan(L_nu_scatter)] = 0
    return L_nu_scatter


def load_L_nu_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    '''для массива энергий (лежат в конфиге). возвращает'''
    L_nu_total = np.empty((config.N_energy, config.N_phase))
    for i, energy in enumerate(config.energy_arr):
        L_nu_surfs = load_L_nu_surfs(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, energy)
        L_nu_scatter = load_L_nu_scatter(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, energy)
        L_nu_total[i] = np.sum(L_nu_surfs, axis=0) + np.sum(L_nu_scatter, axis=0)
    return L_nu_total


def load_PF_nu(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    L_nu_total = load_L_nu_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    return newService.get_PF(L_nu_total)


def load_PF(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    L_total = load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    return newService.get_PF(L_total)


def load_tensor(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, file_name):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    file_name = file_name + '.npy'
    return np.load(data_folder / file_name)


def load_T_eff(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    file_name = 'surfaces_T_eff.txt'
    T_eff = load_arr_from_txt(data_folder, file_name)
    return T_eff


def load_theta_range(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    file_name = 'save_theta_range.txt'
    theta_range = load_arr_from_txt(data_folder, file_name)
    return theta_range

