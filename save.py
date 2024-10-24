import pathlib
import numpy as np
import matplotlib.pyplot as plt

import pathService


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


def load_L_surfs(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    folder = 'surfs/'
    file_name = "L_surfs.txt"
    return load_arr_from_txt(data_folder / folder, file_name)


def load_L_scatter(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    folder = 'scatter/'
    file_name = "L_scatter.txt"
    return load_arr_from_txt(data_folder / folder, file_name)


def load_L_total(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    # sum of surfs
    L_surfs = load_L_surfs(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    # if not np.isnan(buf).any():
    L_scatter = load_L_scatter(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    return np.sum(L_surfs, axis=0) + np.sum(L_scatter, axis=0)
