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


# Функция для сохранения с фиксированными границами
def save_fixed_size(fig, ax, file_path, file_name, has_cbar=False):
    # попытка от deepseek
    create_file_path(file_path)
    full_file_name = file_path / (file_name + '.png')

    dpi = 300
    # ax.set_aspect('equal')
    plt.tight_layout()

    # Жестко фиксируем границы
    if has_cbar:
        # Для графиков с colorbar
        fig.savefig(full_file_name, dpi=dpi, bbox_inches='tight', pad_inches=0.1)
    else:
        # Для графиков без colorbar добавляем искусственные отступы
        # Добавляем прозрачный прямоугольник справа для выравнивания
        from PIL import Image

        fig.savefig(full_file_name, dpi=dpi, bbox_inches='tight',
                    pad_inches=0.1, bbox_extra_artists=[plt.Rectangle((0, 0), 1, 1, fill=False)])

        content_width = 8  # Ширина контента в дюймах
        content_height = 6  # Высота контента в дюймах
        cbar_width = 0.5

        img = Image.open(full_file_name)
        new_img = Image.new('RGBA', (int((content_width + cbar_width) * dpi) - 20, int(content_height * dpi)),
                            (0, 0, 0, 0))
        new_img.paste(img, (0, 0))
        new_img.save(full_file_name)

    plt.close(fig)


def save_figs_same_w(fig, ax, file_path, file_name, has_cbar=False):
    # Общие параметры
    content_width = 8  # Ширина контента в дюймах
    content_height = 6  # Высота контента в дюймах
    cbar_width = 0.5  # Ширина colorbar в дюймах
    dpi = 300

    def save_plot(ax, filename, has_cbar=False):
        # Жестко фиксируем размеры контента
        ax.set_position([0.1, 0.1, 0.8, 0.8])

        if has_cbar:
            # Для графиков с colorbar
            fig = ax.figure
            fig.set_size_inches(content_width + cbar_width, content_height)
            cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.05)
            fig.savefig(filename, dpi=dpi, bbox_inches=Bbox([[0, 0], [content_width + cbar_width, content_height]]))
        else:
            # Для графиков без colorbar
            fig = ax.figure
            fig.set_size_inches(content_width, content_height)
            fig.savefig(filename, dpi=dpi, bbox_inches=Bbox([[0, 0], [content_width, content_height]]))
            # Добавляем прозрачный прямоугольник справа для выравнивания
            from PIL import Image
            img = Image.open(filename)
            new_img = Image.new('RGBA', (int((content_width + cbar_width) * dpi), int(content_height * dpi)),
                                (0, 0, 0, 0))
            new_img.paste(img, (0, 0))
            new_img.save(filename)


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


def load_L_x_calc(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    with open(data_folder / 'save_values.txt') as f:
        lines = f.readlines()
    return float(lines[5][41:49]) * 10 ** float(lines[5][56:58])


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


def load_beta(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    with open(data_folder / 'save_values.txt') as f:
        lines = f.readlines()
    return float(lines[2][7:-1])


def load_gamma(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    with open(data_folder / 'save_values.txt') as f:
        lines = f.readlines()
    return float(lines[13][8:-1])


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


def check_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    data_folder = get_data_folder(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    if data_folder.exists():
        return True
    else:
        return False
