import numpy as np
import math
import matplotlib.pyplot as plt
from mayavi import mlab
import time

from traits.api import HasTraits, Range, Button, Instance, on_trait_change
from traitsui.api import View, Item, HGroup, VGroup
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene

import config
import newService
import geometry
from geometry import matrix
from mayavi_plots import mayavi_geometry


def lim_axes(ax, lim_value):
    ax.set_xlim([-lim_value, lim_value])
    ax.set_ylim([-lim_value, lim_value])
    ax.set_zlim([-lim_value, lim_value])
    ax.set_aspect("equal")


def hide_axes_and_background(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    ax.set_axis_off()


def plot_NS_mayavi(phi_range_column, theta_range_column):
    '''рисуем в масштабах магнитосферы R_e - так проще'''
    # рисуем звезду
    theta_range = np.linspace(0, np.pi, num=config.N_theta_accretion, endpoint=True)
    phi_range = np.linspace(0, 2 * np.pi, num=config.N_phi_accretion, endpoint=True)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    mlab.mesh(x, y, z, color=(0, 0, 1))


def get_data_for_NS(theta_range_column):
    theta_range = np.linspace(0, np.pi, num=config.N_theta_accretion, endpoint=True)
    phi_range = np.linspace(0, 2 * np.pi, num=config.N_phi_accretion, endpoint=True)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    return x, y, z


def get_data_for_accretion_disc_side_surface(accretion_disc_con_val):
    # cylinder
    radius = 4
    z = np.linspace(-accretion_disc_con_val * radius, accretion_disc_con_val * radius, config.N_theta_accretion)
    phi = np.linspace(0, 2 * np.pi, config.N_phi_accretion)
    phi_grid, z = np.meshgrid(phi, z)
    x = radius * np.cos(phi_grid)
    y = radius * np.sin(phi_grid)

    return x, y, z


def get_data_for_accretion_disc(accretion_disc_con_val):
    r, phi = np.mgrid[1:4:100j, 0:2 * np.pi:100j]
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    z = accretion_disc_con_val * (x ** 2 + y ** 2) ** (1 / 2)
    # z, z1 = np.mgrid[-0.003:0.003:100j, -0.003:0.003:100j]
    return x, y, z


def plot_accr_columns_mayavi(phi_range_column, theta_range_column):
    # рисуем силовые линии

    # верх
    theta_range = theta_range_column
    phi_range = phi_range_column

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    mlab.mesh(x, y, z, color=(1, 0, 0))
    # ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)

    # низ
    mlab.mesh(-x, -y, -z, color=(0, 1, 0))


def get_data_for_accretion_columns(phi_range_column, theta_range_column):
    # phi_range = np.linspace(0, 2 * np.pi, num=config.N_phi_accretion, endpoint=True)

    phi_range = phi_range_column
    theta_range = theta_range_column

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    return x, y, z


def get_data_for_accretion_columns_outer(phi_range_column, theta_range_column):
    # phi_range = phi_range_column
    # theta_range = theta_range_column
    #
    # r, p = np.meshgrid((1 + config.dRe_div_Re) * np.sin(theta_range) ** 2, phi_range)
    # r1 = r * np.sin(theta_range)
    # x = r1 * np.cos(p)
    # y = r1 * np.sin(p)
    # z = r * np.cos(theta_range)
    #
    # return x, y, z

    x, y, z = get_data_for_accretion_columns(phi_range_column, theta_range_column)
    return x * (1 + config.dRe_div_Re), y * (1 + config.dRe_div_Re), z * (1 + config.dRe_div_Re)


def get_data_for_accretion_columns_hat(phi_range_column, theta_range_column):
    phi_range = phi_range_column

    r1 = np.sin(theta_range_column[-1]) ** 3
    r2 = (np.sin(theta_range_column[-1]) / (1 + config.dRe_div_Re) ** (1 / 2)) ** 3 * (1 + config.dRe_div_Re)
    r = np.linspace(r1, r2, 100)

    r, phi = np.meshgrid(r, phi_range)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    z = np.full_like(y, np.sin(theta_range_column[-1]) ** 2 * np.cos(theta_range_column[-1]))
    # z, z1 = np.mgrid[-0.003:0.003:100j, -0.003:0.003:100j]
    return x, y, z


def get_data_for_magnet_lines_with_mask(phi_range_magnet_lines, theta_range_magnet_lines, mask_magnet_lines):
    '''смог реализовать только с помощью маски
    -1 индекс так как начало магнитных линий = конец колонки'''
    phi_range = phi_range_magnet_lines
    theta_range = theta_range_magnet_lines

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    return x, y, z, mask_magnet_lines


def get_data_for_magnet_lines_outer_with_mask(phi_range_magnet_lines, theta_range_magnet_lines, mask_magnet_lines):
    x, y, z, mask_magnet_lines = get_data_for_magnet_lines_with_mask(phi_range_magnet_lines, theta_range_magnet_lines,
                                                                     mask_magnet_lines)
    return x * (1 + config.dRe_div_Re), y * (1 + config.dRe_div_Re), z * (1 + config.dRe_div_Re), mask_magnet_lines


def plot_magnet_lines(ax, phi_range_column):
    step_theta_accretion = (np.pi / 2) / (config.N_theta_accretion - 1)
    theta_range = np.array([step_theta_accretion * j for j in range(config.N_theta_accretion)])

    # верх
    phi_range = phi_range_column

    r, p = np.meshgrid(np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    ax.plot_wireframe(x, y, z, rstride=4, cstride=4, color="blue", alpha=0.06)
    # ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)

    # низ
    ax.plot_wireframe(-x, -y, -z, rstride=4, cstride=4, color="blue", alpha=0.06)


import accretingNS


def plot_main(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, flag_do_not_draw):
    curr_configuration = accretingNS.AccretingPulsarConfiguration(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)

    class Visualization(HasTraits):
        slider_phi_0 = Range(0, 360, phi_0)
        slider_theta_obs = Range(0, 90, theta_obs)
        slider_beta_mu = Range(0, 90, beta_mu)
        slider_phase = Range(0., 2., 0.)
        slider_distance = Range(0.1, 10, 1)

        button_magnet_line = Button('draw_magnet_lines')
        button_accr_disc = Button('accr_disc_omega_mu')
        button_cut_magnet_lines = Button('cut_magnet_lines')
        button_hide_accr_disc = Button('hide_accr_disc')
        button_check_data = Button('check_data')
        button_animate = Button('animate')

        scene = Instance(MlabSceneModel, ())

        color_accretion_column_top = (1, 0, 0)
        color_accretion_column_bot = (0, 1, 0)

        color_accretion_column_top_outer = (1, 1, 0)
        color_accretion_column_bot_outer = (0, 1, 1)

        color_magnet_lines_top = (0, 0, 1)
        color_magnet_lines_top_outer = (0, 0.5, 1)
        color_magnet_lines_bot = (0, 0, 1)
        color_magnet_lines_bot_outer = (0, 0, 1)

        mu_vector_tube_radius = 0.005
        omega_vector_tube_radius = 0.005

        mu_vector_color = (1, 0, 0)
        omega_vector_color = (0, 0, 0)

        def draw_accretion_disc(self):
            self.flag_accretion_disc_hide = False
            disc_color = (160, 82, 45)
            disc_color = tuple(rgb / 255 for rgb in disc_color)

            accretion_disc_con_val = 0.01
            x, y, z = get_data_for_accretion_disc(accretion_disc_con_val)

            # сначала отрисовали в СК beta_mu
            self.accretion_disc_top = mlab.mesh(x, y, z, color=disc_color)
            self.accretion_disc_bot = mlab.mesh(x, y, -z, color=disc_color)

            self.flag_accretion_disc_omega_mu = True  # True = omega

            x, y, z = get_data_for_accretion_disc_side_surface(accretion_disc_con_val)

            # боковая поверхность диска (как боковая поверхность цилиндра)
            self.accretion_disc_side_surface = mlab.mesh(x, y, z, color=disc_color)

            # поворачиваем диск на -betta чтобы было перпендикулярно omega -> переходим в СК omega
            # deg ??
            self.accretion_disc_rotate_angle = beta_mu
            self.rotate_accretion_disc(-self.accretion_disc_rotate_angle)

        def __init__(self, curr_configuration: accretingNS.AccretingPulsarConfiguration):
            # Do not forget to call the parent's __init__
            HasTraits.__init__(self)
            self.scene.background = (1, 1, 1)
            # mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0., 0., 0.))
            # NS
            x, y, z = get_data_for_NS(curr_configuration.top_column.inner_surface.theta_range)
            self.NS = self.scene.mlab.mesh(x, y, z, color=(0, 0, 0))

            # рисуем колонки
            x, y, z = get_data_for_accretion_columns(curr_configuration.top_column.inner_surface.phi_range,
                                                     curr_configuration.top_column.inner_surface.theta_range)
            # верх
            self.accretion_column_top = self.scene.mlab.mesh(x, y, z, color=self.color_accretion_column_top)
            # низ
            self.accretion_column_bot = self.scene.mlab.mesh(-x, -y, -z, color=self.color_accretion_column_bot)

            # внешние
            x, y, z = get_data_for_accretion_columns_outer(curr_configuration.top_column.outer_surface.phi_range,
                                                           curr_configuration.top_column.outer_surface.theta_range)
            # верх
            self.accretion_column_top_outer = self.scene.mlab.mesh(x, y, z, color=self.color_accretion_column_top_outer)
            # низ
            self.accretion_column_bot_outer = self.scene.mlab.mesh(-x, -y, -z,
                                                                   color=self.color_accretion_column_bot_outer)

            # попытка отрисовать боковые
            # x, y, z = get_data_for_accretion_columns_hat(theta_range_column, phi_range_column,
            #                                              config.phi_accretion_begin_deg)

            # self.accretion_column_top_outer_hat = self.scene.mlab.mesh(x, y, z, color=self.color_accretion_column_top)

            self.flag_draw_magnet_lines = True
            self.flag_cut_magnet_lines = False
            # x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, config.phi_accretion_begin_deg)
            x, y, z, mask = get_data_for_magnet_lines_with_mask(curr_configuration.top_magnet_lines.phi_range,
                                                                curr_configuration.top_magnet_lines.theta_range,
                                                                curr_configuration.top_magnet_lines.mask_array)
            opacity_for_magnet_line = 0.1

            if not flag_do_not_draw:
                # .visible = False - чтобы сделать невидимым
                self.magnet_lines_top = self.scene.mlab.mesh(x, y, z, color=(0, 0, 1), opacity=opacity_for_magnet_line,
                                                             representation='wireframe', mask=mask)

                # self.magnet_lines_top = self.scene.mlab.surf(x, y, z, color=(1, 0, 0), warp_scale=0.3,
                #                                              representation='wireframe', line_width=0.5)
                self.magnet_lines_bot = self.scene.mlab.mesh(-x, -y, -z, color=self.color_magnet_lines_bot,
                                                             opacity=opacity_for_magnet_line,
                                                             representation='wireframe', mask=mask)

            x, y, z, mask = get_data_for_magnet_lines_outer_with_mask(curr_configuration.top_magnet_lines.phi_range,
                                                                      curr_configuration.top_magnet_lines.theta_range,
                                                                      curr_configuration.top_magnet_lines.mask_array)

            if not flag_do_not_draw:
                self.magnet_lines_top_outer = self.scene.mlab.mesh(x, y, z,
                                                                   color=self.color_magnet_lines_top_outer,
                                                                   opacity=opacity_for_magnet_line,
                                                                   representation='wireframe', mask=mask)

                self.magnet_lines_bot_outer = self.scene.mlab.mesh(-x, -y, -z,
                                                                   color=self.color_magnet_lines_bot_outer,
                                                                   opacity=opacity_for_magnet_line,
                                                                   representation='wireframe', mask=mask)

            # mu_vector
            # mlab.plot3d([0, 0], [0, 0], [0, 1.0], color=self.mu_vector_color, tube_radius=self.mu_vector_tube_radius,
            #             tube_sides=6)

            self.mu_vector = mlab.quiver3d(0, 0, 1, mode='2ddash', scale_factor=1, color=self.mu_vector_color)
            self.mu_vector_1 = mlab.quiver3d(0, 0, -1, mode='2ddash', scale_factor=1, color=self.mu_vector_color)

            omega_vector = [np.sin(np.deg2rad(-beta_mu)) * np.cos(0),
                            np.sin(np.deg2rad(-beta_mu)) * np.sin(0),
                            np.cos(np.deg2rad(-beta_mu))]

            omega_vector = matrix.get_cartesian_from_spherical(1, np.deg2rad(-beta_mu), 0)
            # omega_vector
            # self.omega_vector = mlab.plot3d([0, omega_vector[0]], [0, omega_vector[1]], [0, omega_vector[2]],
            #                                 color=self.omega_vector_color, tube_radius=self.omega_vector_tube_radius,
            #                                 tube_sides=4)

            self.omega_vector = mlab.quiver3d(omega_vector[0], omega_vector[1], omega_vector[2], mode='2ddash',
                                              scale_factor=1, color=self.omega_vector_color)
            self.omega_vector_1 = mlab.quiver3d(-omega_vector[0], -omega_vector[1], -omega_vector[2], mode='2ddash',
                                                scale_factor=1, color=self.omega_vector_color)

            # рисуем аккреционный диск
            self.draw_accretion_disc()

            # self.check_data()

        def rotate_accretion_disc(self, rotate_angle):
            self.accretion_disc_top.actor.actor.rotate_y(rotate_angle)
            self.accretion_disc_bot.actor.actor.rotate_y(rotate_angle)
            self.accretion_disc_side_surface.actor.actor.rotate_y(rotate_angle)

        def update_accretion_disc_rotate_angle(self):
            if self.flag_accretion_disc_omega_mu:
                self.rotate_accretion_disc(self.accretion_disc_rotate_angle)
                self.accretion_disc_rotate_angle = beta_mu
                self.rotate_accretion_disc(-self.accretion_disc_rotate_angle)
            else:
                self.accretion_disc_rotate_angle = beta_mu

        '''def update_magnet_lines(self):
            opacity_for_magnet_line = 0.1
            x, y, z, mask = get_data_for_magnet_lines_with_mask(theta_range_column, phi_range_column,
                                                                self.slider_phi_0)


            self.magnet_lines_top.mlab_source.trait_set(x=[0], y=[0], z=[0])
            self.magnet_lines_bot.mlab_source.trait_set(x=[0], y=[0], z=[0])
            # new
            self.magnet_lines_top = self.scene.mlab.mesh(x, y, z, color=self.color_magnet_lines_top,
                                                         opacity=opacity_for_magnet_line,
                                                         representation='wireframe', mask=mask)
            self.magnet_lines_bot = self.scene.mlab.mesh(-x, -y, -z, color=self.color_magnet_lines_bot,
                                                         opacity=opacity_for_magnet_line,
                                                         representation='wireframe', mask=mask)

            x, y, z, mask = get_data_for_magnet_lines_outer_with_mask(theta_range_column, phi_range_column,
                                                                      self.slider_phi_0)

            self.magnet_lines_top_outer.mlab_source.trait_set(x=[0], y=[0], z=[0])
            self.magnet_lines_bot_outer.mlab_source.trait_set(x=[0], y=[0], z=[0])
            # new

            self.magnet_lines_top_outer = self.scene.mlab.mesh(x, y, z,
                                                               color=self.color_magnet_lines_top_outer,
                                                               opacity=opacity_for_magnet_line,
                                                               representation='wireframe', mask=mask)

            self.magnet_lines_bot_outer = self.scene.mlab.mesh(-x, -y, -z,
                                                               color=self.color_magnet_lines_bot_outer,
                                                               opacity=opacity_for_magnet_line,
                                                               representation='wireframe', mask=mask)'''

        def view_phase(self, phase=0):

            e_obs = matrix.get_cartesian_from_spherical(1, np.deg2rad(theta_obs), 0)

            A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, np.deg2rad(phase), np.deg2rad(beta_mu))
            e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

            azimuth, elevation = matrix.vec_to_angles(e_obs_mu)
            roll_angle = mayavi_geometry.calculate_roll_angle(theta_obs, beta_mu, e_obs_mu, phase)
            # print(roll_angle)

            # ax.view_init(90 - elevation / config.grad_to_rad, azimuth / config.grad_to_rad, roll=roll_angle)
            distance = self.slider_distance
            self.scene.mlab.view(azimuth=np.rad2deg(azimuth), elevation=np.rad2deg(elevation),
                                 distance=distance, focalpoint=[0, 0, 0])
            # roll angle считаем для плоскости камеры - поэтому roll там
            # только здесь нашел про камеру https://docs.enthought.com/mayavi/mayavi/mlab_figures_decorations.html
            camera = self.scene.camera
            camera.roll(roll_angle)

        '''def check_data(self):
            # columns = {0: top_column, 1: bot_column}
            def get_cos_magnet_lines():
                # magnet_lines
                top_column_cos_array, bot_column_cos_array = [], []
                magnet_lines_cos_file_folder = 'data/magnet_cos/' + config.file_folder_angle_args + config.file_folder_accretion_args
                file_name_for_magnet_lines_cos_of_columns = {0: 'top_column', 1: 'bot_column'}
                cos_array_dict = {0: top_column_cos_array, 1: bot_column_cos_array}
                for key, column_name in file_name_for_magnet_lines_cos_of_columns.items():
                    full_magnet_line_cos_file_folder = magnet_lines_cos_file_folder + \
                                                       file_name_for_magnet_lines_cos_of_columns[key] + '/'
                    for cos_index in range(config.t_max):
                        file_name = 'save_magnet_lines_cos_' + column_name + ('_%d_phase' % cos_index) + '.txt'
                        cos_range = main_service.load_arr_from_txt(full_magnet_line_cos_file_folder, file_name)
                        cos_array_dict[key].append(cos_range)
                return cos_array_dict

            def get_cos_surfaces():
                # accret_columns
                top_column_outer_surface_cos_array, top_column_inner_surface_cos_array = [], []
                bot_column_outer_surface_cos_array, bot_column_inner_surface_cos_array = [], []

                cos_array_dict = {0: top_column_outer_surface_cos_array, 1: top_column_inner_surface_cos_array,
                                  2: bot_column_outer_surface_cos_array, 3: bot_column_inner_surface_cos_array}

                file_name_for_cos_of_surfaces = {0: 'top_outer', 1: 'top_inner', 2: 'bot_outer', 3: 'bot_inner'}
                # cos_file_folder = 'data/cos/' + config.file_folder_angle_args + config.file_folder_accretion_args
                cos_file_folder = config.PROJECT_DIR + 'data/cos/' + config.file_folder_angle_args + config.file_folder_accretion_args
                for key, surface_name in file_name_for_cos_of_surfaces.items():
                    full_cos_file_folder = cos_file_folder + file_name_for_cos_of_surfaces[key] + '/'
                    for cos_index in range(config.t_max):
                        file_name = 'save_cos_' + surface_name + ('_%d_phase' % cos_index) + '.txt'

                        cos_range = main_service.load_arr_from_txt(full_cos_file_folder, file_name)
                        cos_array_dict[key].append(cos_range)
                return cos_array_dict

            # self.cos_array_dict_magnet_lines = get_cos_magnet_lines()
            self.cos_array_dict_surfaces = get_cos_surfaces()'''

        # def try_check_data(self):
        #     '''получить фазу, по нужной фазе достать косинусы, закрасить те участки которые не видны (или наоборот)'''
        #     phase = (360 * self.slider_phase) % 360
        #     index = math.floor(phase // config.omega_ns)
        #
        #     x, y, z = get_data_for_accretion_columns(theta_range_column, phi_range_column,
        #                                              config.phi_accretion_begin_deg)
        #     mask = np.zeros_like(x).astype(bool)
        #     for i in range(config.N_phi_accretion):
        #         for j in range(config.N_theta_accretion):
        #
        #             # if self.cos_array_dict_magnet_lines[0][index][i][j] > 0:
        #             if self.cos_array_dict_surfaces[0][index][i][j] < 0.15:
        #                 mask[i][j] = True
        #
        #                 # self.scene.mlab.points3d(
        #                 #     self.magnet_lines_top.mlab_source.x[i][j],
        #                 #     self.magnet_lines_top.mlab_source.y[i][j],
        #                 #     self.magnet_lines_top.mlab_source.z[i][j],
        #                 #     color=(1, 0, 1),
        #                 #     scale_factor=0.01
        #                 # )
        #     self.accretion_column_top.mlab_source.trait_set(x=[0], y=[0], z=[0])
        #     self.accretion_column_top = self.scene.mlab.mesh(x, y, z, color=self.color_accretion_column_top, mask=mask)
        #     # smth.mlab_source.x[i][j] = 0

        @on_trait_change('slider_theta_obs')
        def func_change_slider_i_angle_slider(self):
            theta_obs = self.slider_theta_obs

            curr_configuration.theta_obs = theta_obs

            self.update_accretion_disc_rotate_angle()

            phase = 360 * self.slider_phase
            self.view_phase(phase)

        @on_trait_change('slider_beta_mu')
        def func_change_slider_betta_mu(self):
            beta_mu = self.slider_beta_mu

            curr_configuration.beta_mu = beta_mu

            omega_vector = matrix.get_cartesian_from_spherical(1, np.deg2rad(-beta_mu), 0)

            self.omega_vector.mlab_source.vectors = np.reshape([omega_vector[0], omega_vector[1], omega_vector[2]],
                                                               (1, 3))
            self.omega_vector_1.mlab_source.vectors = np.reshape([-omega_vector[0], -omega_vector[1], -omega_vector[2]],
                                                                 (1, 3))

            self.update_accretion_disc_rotate_angle()

            if self.flag_draw_magnet_lines:
                # x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, self.slider_phi_0,
                #                                     self.flag_cut_magnet_lines)

                # обновляем т.к. зависит от beta_mu
                self.update_magnet_lines()

            phase = 360 * self.slider_phase
            self.view_phase(phase)

        @on_trait_change('slider_phi_0')
        def func_change_slider_fi_0(self):
            global curr_configuration

            curr_configuration = accretingNS.AccretingPulsarConfiguration(curr_configuration.mu,
                                                                          curr_configuration.theta_obs,
                                                                          curr_configuration.beta_mu,
                                                                          curr_configuration.mc2,
                                                                          curr_configuration.a_portion,
                                                                          self.slider_phi_0)

            x, y, z = get_data_for_accretion_columns(curr_configuration.top_column.inner_surface.phi_range,
                                                     curr_configuration.top_column.inner_surface.theta_range)

            self.accretion_column_top.mlab_source.trait_set(x=x, y=y, z=z, color=self.color_accretion_column_top)
            self.accretion_column_bot.mlab_source.trait_set(x=-x, y=-y, z=-z, color=self.color_accretion_column_bot)

            x, y, z = get_data_for_accretion_columns_outer(curr_configuration.top_column.inner_surface.phi_range,
                                                           curr_configuration.top_column.inner_surface.theta_range)
            # верх
            self.accretion_column_top_outer.mlab_source.trait_set(x=x, y=y, z=z,
                                                                  color=self.color_accretion_column_top_outer)
            # низ
            self.accretion_column_bot_outer.mlab_source.trait_set(x=-x, y=-y, z=-z,
                                                                  color=self.color_accretion_column_bot_outer)

            if self.flag_draw_magnet_lines:
                # x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, self.slider_phi_0,
                #                                     self.flag_cut_magnet_lines)

                self.update_magnet_lines()
            # self.magnet_lines_top.mlab_source.trait_set(x=x, y=y, z=z, mask=mask)
            # self.magnet_lines_bot.mlab_source.trait_set(x=-x, y=-y, z=-z, mask=mask)

            phase = 360 * self.slider_phase
            self.view_phase(phase)

        @ on_trait_change('slider_phase, slider_distance')
        def func_change_slider_phase_slider_distance(self):
            phase = 360 * self.slider_phase
            self.view_phase(phase)


        # @on_trait_change('button_magnet_line')
        # def func_change_button_magnet_line(self):
        #     self.flag_draw_magnet_lines = not self.flag_draw_magnet_lines
        #     if self.flag_draw_magnet_lines:
        #         x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, self.slider_phi_0,
        #                                             self.flag_cut_magnet_lines)
        #         self.magnet_lines_top.mlab_source.trait_set(x=x, y=y, z=z, color=self.color_magnet_lines_top)
        #         self.magnet_lines_bot.mlab_source.trait_set(x=-x, y=-y, z=-z, color=self.color_magnet_lines_bot)
        #     else:
        #         self.magnet_lines_top.mlab_source.trait_set(x=[0], y=[0], z=[0], color=self.color_magnet_lines_top)
        #         self.magnet_lines_bot.mlab_source.trait_set(x=[0], y=[0], z=[0], color=self.color_magnet_lines_bot)


        @on_trait_change('button_cut_magnet_lines')
        def func_change_button_cut_magnet_lines(self):
            self.flag_cut_magnet_lines = not self.flag_cut_magnet_lines
            # if self.flag_draw_magnet_lines:
            #     x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, self.slider_phi_0,
            #                                         self.flag_cut_magnet_lines)
            #     self.magnet_lines_top.mlab_source.trait_set(x=x, y=y, z=z, color=(0, 0, 1))
            #     self.magnet_lines_bot.mlab_source.trait_set(x=-x, y=-y, z=-z, color=(0, 0, 1))
            # else:
            #     self.magnet_lines_top.mlab_source.trait_set(x=[0], y=[0], z=[0], color=(0, 0, 1))
            #     self.magnet_lines_bot.mlab_source.trait_set(x=[0], y=[0], z=[0], color=(0, 0, 1))

            theta = np.pi / 2 - config.betta_mu

            r, phi = np.mgrid[np.sin(theta) ** 2:4:100j, 0:2 * np.pi:100j]
            x = r * np.cos(phi)
            y = r * np.sin(phi)
            z, z1 = np.mgrid[-0.003:0.003:100j, -0.003:0.003:100j]

            # сначала отрисовали в beta_mu
            self.accretion_disc_top.mlab_source.trait_set(x=x, y=y)


        @on_trait_change('button_accr_disc')
        def func_change_button_accr_disc(self):
            if self.flag_accretion_disc_omega_mu:
                self.rotate_accretion_disc(self.accretion_disc_rotate_angle)
            else:
                self.rotate_accretion_disc(-self.accretion_disc_rotate_angle)
            self.flag_accretion_disc_omega_mu = not self.flag_accretion_disc_omega_mu


        @on_trait_change('button_hide_accr_disc')
        def func_change_button_hide_accr_disc(self):
            if self.flag_accretion_disc_hide:
                # self.draw_accretion_disc()
                self.flag_accretion_disc_hide = not self.flag_accretion_disc_hide
                self.accretion_disc_top.visible = True
                self.accretion_disc_bot.visible = True
                self.accretion_disc_side_surface.visible = True
            else:
                self.flag_accretion_disc_hide = not self.flag_accretion_disc_hide
                self.accretion_disc_top.visible = False
                self.accretion_disc_bot.visible = False
                self.accretion_disc_side_surface.visible = False


        @on_trait_change('button_check_data')
        def func_change_button_check_data(self):
            self.try_check_data()


        @on_trait_change('button_animate')
        def anim(self):
            self.scene.reset_zoom()
            # save_folder = config.PROJECT_DIR_ORIGIN + 'mayavi_figs/'
            N = 100
            for i in range(N):
                self.slider_distance = 3.89764
                self.slider_phase = 1 * i / (N - 1)
                self.scene.save_png('mayavi_figs/' + 'anim%02d.png' % i)


        # @mlab.animate
        # def anim(self):
        #     for i in range(10):
        #         self.slider_phase = 2 * i / 10
        #         yield

        # the layout of the dialog created
        view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                         height=250, width=300, show_label=False),
                    VGroup(
                        HGroup(
                            '_', 'slider_theta_obs', 'slider_beta_mu'
                        ),
                        HGroup(
                            '_', 'slider_phi_0', 'slider_phase'
                        ),
                        HGroup('button_magnet_line', 'button_accr_disc', 'button_cut_magnet_lines',
                               'button_hide_accr_disc', 'slider_distance', 'button_animate')
                    )
                    )

    visualization = Visualization(curr_configuration)
    visualization.configure_traits()
    # visualization.anim()


if __name__ == "__main__":
    mu = 0.1e31
    beta_mu = 20
    mc2 = 100
    a_portion = 0.44
    phi_0 = 0

    theta_obs = 80

    plot_main(mu, theta_obs, beta_mu, mc2, a_portion, phi_0, False)