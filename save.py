import pathlib
import numpy as np
import matplotlib.pyplot as plt


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
    np.savetxt(full_file_name, arr)


def load_arr_from_txt(file_folder, file_name):
    return np.loadtxt(file_folder + file_name)


def save_figure(fig, file_path, file_name):
    create_file_path(file_path)
    full_file_name = file_path + file_name
    fig.savefig(full_file_name, dpi=fig.dpi)
    # fig.savefig(full_file_name, dpi=200)
    plt.close(fig)
