import pathlib
import config

print_flag = False


def get_prefix_dir():
    PROJECT_DIR = pathlib.Path(__file__).resolve().parent
    prefix_dir = PROJECT_DIR / 'data'
    if not config.NS_shadow_flag:
        prefix_dir = prefix_dir / 'NS_shadow_off'
    if not config.flag_scatter:
        prefix_dir = prefix_dir / 'scatter_off'
    if not config.flag_attenuation_above_shock:
        prefix_dir = prefix_dir / 'tau_off'

    return prefix_dir


def get_args_dir(mu=None, theta_obs=None, beta_mu=None, mc2=None, a_portion=None, phi_0=None):
    if mu is not None:
        buf = mu
        count = 1
        while buf > 1:
            count += 1
            buf //= 10
    else:
        count = 31

    args_dir = pathlib.Path(f'mu=0.1e{count}')
    if theta_obs is not None:
        args_dir = args_dir / f'theta_obs={theta_obs}'
    if beta_mu is not None:
        args_dir = args_dir / f'beta_mu={beta_mu}'
    if mc2 is not None:
        args_dir = args_dir / f'mc2={mc2}'
    if a_portion is not None:
        args_dir = args_dir / f'a={a_portion:.2f}'
    if phi_0 is not None:
        args_dir = args_dir / f'phi_0={phi_0}'
    return args_dir


def get_dir(mu, theta_obs=None, beta_mu=None, mc2=None, a_portion=None, phi_0=None, prefix_folder=None,
            suffix_folder=None):
    prefix_dir = get_prefix_dir()
    if prefix_folder is not None:
        prefix_dir = prefix_dir / prefix_folder
    args_dir = get_args_dir(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    if suffix_folder is not None:
        args_dir = args_dir / suffix_folder
    return prefix_dir / args_dir


class PathSaver:

    def __init__(self, mu, theta_obs, beta_mu, mc2, a_portion, phi_0, NS_shadow_flag=None, flag_scatter=None, flag_attenuation_above_shock=None):
        '''
            поменять название betta_mu на beta, i на theta_obs, fi_0 на phi_0 !!!!
            'NS_shadow_off/'
        '''
        # возможно убрать 1 parent!
        self.PROJECT_DIR = pathlib.Path(__file__).resolve().parent

        self.old_prefix_dir = self.PROJECT_DIR / 'new_magnet_lines' / 'new_condition' / 'new_data'

        buf = mu
        count = 1
        while buf > 1:
            count += 1
            buf //= 10

        self.old_prefix_dir = self.old_prefix_dir / f'mu=0.1e{count}' / 'tau' / 'figs' / 'loop'
        fi_0_old = -10000
        # fi_0_old = (phi_0 + config.fi_0_dict[a_portion]) % 360
        # поменять название betta_mu на beta. i на theta!
        self.old_save_dir = self.old_prefix_dir / f'i={theta_obs} betta_mu={beta_mu}'
        # поменять fi_0 на phi_0
        self.old_save_dir = self.old_save_dir / f'mc2={mc2}' / f'a={a_portion:.2f} fi_0={fi_0_old}'

        # ------------------------ new way ---------------------------

        self.prefix_dir = self.PROJECT_DIR / 'data'

        if NS_shadow_flag is None:
            NS_shadow_flag = config.NS_shadow_flag
        if flag_scatter is None:
            flag_scatter = config.flag_scatter
        if flag_attenuation_above_shock is None:
            flag_attenuation_above_shock = config.flag_attenuation_above_shock

        # if NS_shadow_flag and flag_scatter:
        #     pass
        # else:
        #     if not NS_shadow_flag and not flag_scatter:
        #         self.prefix_dir = self.prefix_dir / 'NS_shadow_off_scatter_off'
        #     elif not NS_shadow_flag:
        #         self.prefix_dir = self.prefix_dir / 'NS_shadow_off'
        #     else:
        #         self.prefix_dir = self.prefix_dir / 'scatter_off'

        if not flag_scatter:
            self.prefix_dir = self.prefix_dir / 'scatter_off'
        if not NS_shadow_flag:
            self.prefix_dir = self.prefix_dir / 'NS_shadow_off'
        if not flag_attenuation_above_shock:
            self.prefix_dir = self.prefix_dir / 'tau_off'

        # if config.tau_flag:
        #     self.prefix_dir = self.prefix_dir / 'tau'
        # else:
        #     self.prefix_dir = self.prefix_dir / 'transperent'
        self.save_dir = self.prefix_dir / f'mu=0.1e{count}'
        self.save_dir = self.save_dir / f'theta_obs={theta_obs}' / f'beta_mu={beta_mu}'
        # поменять fi_0 на phi_0
        self.save_dir = self.save_dir / f'mc2={mc2}' / f'a={a_portion:.2f}' / f'phi_0={phi_0}'

        if print_flag:
            print(f'PROJECT_DIR: {self.PROJECT_DIR}')
            print(f'prefix_dir: {self.prefix_dir}')
            print(f'save_dir: {self.save_dir}')
        # print('D:\MyProgs\Python\Diplom\\new_magnet_lines\\new_condition\\new_data\mu=0.1e31\\tau\\figs\loop\i=60 betta_mu=40\mc2=100\\a=0.44 fi_0=280')
        # filepath = p / fn

        # print(PROJECT_DIR)
        # print(PROJECT_DIR.resolve())

        # print(pathlib.Path(__file__).resolve().parent.parent.absolute())
        # mkdir()
        # print([x for x in PROJECT_DIR.iterdir() if x.is_dir()])

        # print(p.is_dir())
        # print(p.exists())

        # print(p)
        # print(self.PROJECT_DIR.joinpath('new_log'))
        # p = pathlib.Path('.').resolve()
        # print(p)
        # print([x for x in p.iterdir() if x.is_dir()])

        # with open(path, mode='wt') as config:
        #     config.write('# config goes here')

        # os.chdir
        # from pathlib import Path
        # from os import chdir
        #
        # parent = Path('..')
        # chdir(parent)

    def get_path(self):
        return self.save_dir

    def get_save_args_dir(self):
        return self.save_dir

    def get_old_path(self):
        return self.old_save_dir


if __name__ == '__main__':
    mu = 0.1e31
    beta_mu = 40
    mc2 = 100
    a_portion = 0.44
    phi_0 = 0
    theta_obs = 60

    path1 = PathSaver(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    print(get_args_dir(mu, theta_obs, beta_mu, mc2, a_portion))

    print(get_dir(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, 'aggr/'))
