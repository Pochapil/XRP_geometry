import numpy as np
import math
import matplotlib.pyplot as plt
from mayavi import mlab
import time

from traits.api import HasTraits, Range, Button, Instance, on_trait_change, Float
from traits.trait_types import false
from traitsui.api import View, Item, HGroup, VGroup
from tvtk.pyface.scene_editor import SceneEditor
from mayavi.tools.mlab_scene_model import MlabSceneModel
from mayavi.core.ui.mayavi_scene import MayaviScene

import config
import newService
from geometry import matrix
from mayavi_plots import mayavi_geometry
import accretingNS


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
    '''
    рисуем в масштабах магнитосферы R_e - так проще
    поэтому радиус звезды = np.sin(theta_range_column[0]) ** 2 то есть начало колонки
    '''

    theta_range = np.linspace(0, np.pi, num=config.N_theta_accretion, endpoint=True)
    phi_range = np.linspace(0, 2 * np.pi, num=config.N_phi_accretion, endpoint=True)

    u, v = np.meshgrid(phi_range, theta_range)
    r1 = np.sin(theta_range_column[0]) ** 2
    # print(r1)
    x = r1 * np.sin(v) * np.cos(u)
    y = r1 * np.sin(v) * np.sin(u)
    z = r1 * np.cos(v)

    return x, y, z


def get_data_for_accretion_disc_side_surface(accretion_disc_con_val):
    # cylinder - боковая поверхность у диска
    # accretion_disc_con_val = насколько меняется радиус ??
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
    # передаем углы колонок, получаем сетку xyz в 3d
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


def get_data_for_accretion_columns_with_mask(accretion_surface: accretingNS.Surface):
    '''теперь тоже с маской - чтобы сделать обрезы'''
    # передаем углы колонок, получаем сетку xyz в 3d
    # phi_range = np.linspace(0, 2 * np.pi, num=config.N_phi_accretion, endpoint=True)

    phi_range = accretion_surface.phi_range
    theta_range = accretion_surface.theta_range

    coef = 1
    if accretion_surface.surface_type == accretingNS.surface_surf_types['outer']:
        coef = (1 + config.dRe_div_Re)

    r, p = np.meshgrid(coef * np.sin(theta_range) ** 2, phi_range)
    r1 = r * np.sin(theta_range)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range)

    return x, y, z, accretion_surface.mask_array


def get_data_for_accretion_columns_hat(phi_range_column, theta_range_column):
    # переписать так как радиусы поменялись в сравнении со старой версией
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
    '''смог реализовать только с помощью маски'''

    r, p = np.meshgrid(np.sin(theta_range_magnet_lines) ** 2, phi_range_magnet_lines)
    r1 = r * np.sin(theta_range_magnet_lines)
    x = r1 * np.cos(p)
    y = r1 * np.sin(p)
    z = r * np.cos(theta_range_magnet_lines)

    return x, y, z, mask_magnet_lines


def get_data_for_magnet_lines_outer_with_mask(phi_range_magnet_lines, theta_range_magnet_lines, mask_magnet_lines):
    # просто удлинить в (1 + config.dRe_div_Re)
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


class Visualization(HasTraits):
    '''внутри класса буду пересоздавать конфигурацию и хранить в ней главные параметры'''

    slider_theta_obs = Range(0, 90, 0)
    slider_beta_mu = Range(0, 90, 0)
    slider_mc2 = Range(1, 300, 0)
    slider_a_portion = Range(0., 1., 0)
    slider_phi_0 = Range(0, 360, 0)

    R_e_val = Float(0)
    ksi_val = Float(0)
    R_disk_val = Float(0)

    slider_distance_init_val = 100

    slider_phase = Range(0., 2., 0.)
    slider_distance = Range(1, 1e5, slider_distance_init_val)

    button_magnet_line = Button('draw_magnet_lines')
    button_accr_disc = Button('accr_disc_omega_mu')
    button_cut_magnet_lines = Button('cut_magnet_lines')
    button_hide_accr_disc = Button('hide_accr_disc')
    button_check_data = Button('check_data')
    button_animate = Button('animate')

    scene = Instance(MlabSceneModel, ())  # wtf ?? - так было в гайде

    color_accretion_column_top = (1, 0, 0)
    color_accretion_column_bot = (0, 1, 0)

    color_accretion_column_top_outer = (1, 1, 0)
    color_accretion_column_bot_outer = (0, 1, 1)
    color_accretion_column_for_calc = (0, 0, 0)

    color_magnet_lines_top = (0, 0, 1)
    color_magnet_lines_top_outer = (0, 0.5, 1)
    color_magnet_lines_bot = (0, 0, 1)
    color_magnet_lines_bot_outer = (0, 0, 1)

    mu_vector_tube_radius = 0.005
    omega_vector_tube_radius = 0.005

    mu_vector_color = (1, 0, 0)
    omega_vector_color = (0, 0, 0)

    opacity_for_magnet_line = 0.1

    def __init__(self, curr_configuration: accretingNS.AccretingPulsarConfiguration):

        self.slider_theta_obs = curr_configuration.theta_obs
        self.slider_beta_mu = curr_configuration.beta_mu
        self.slider_mc2 = curr_configuration.mc2
        self.slider_a_portion = curr_configuration.a_portion
        self.slider_phi_0 = curr_configuration.phi_0

        self.R_e_val = curr_configuration.top_column_for_plot.R_e / config.R_ns
        self.ksi_val = curr_configuration.top_column_for_plot.ksi_shock
        self.R_disk_val = curr_configuration.R_disk / config.R_ns

        self.R_e_for_distance = self.R_e_val

        self.curr_configuration = curr_configuration

        self.slider_distance = self.ksi_val * 4
        # Do not forget to call the parent's __init__
        HasTraits.__init__(self)

        self.scene.background = (1, 1, 1)
        # self.scene.mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0., 0., 0.))

        # ---------------------------------------------------NS---------------------------------------------------
        self.draw_NS()
        # -----------------------------------------------рисуем колонки-----------------------------------------------
        self.draw_accretion_columns()

        # попытка отрисовать боковые
        # x, y, z = get_data_for_accretion_columns_hat(theta_range_column, phi_range_column,
        #                                              config.phi_accretion_begin_deg)

        # self.accretion_column_top_outer_hat = self.scene.mlab.mesh(x, y, z, color=self.color_accretion_column_top)

        # --------------------------------------------рисуем магнитные линии-------------------------------------------
        self.flag_draw_magnet_lines = True  # для чего ? вроде убрать меш сетку
        self.flag_cut_magnet_lines = False  # для чего ?
        # x, y, z = get_data_for_magnet_lines(theta_range_column, phi_range_column, config.phi_accretion_begin_deg)
        self.draw_magnet_lines()

        # --------------------------------------------- рисуем вектора ------------------------------------------
        self.draw_mu_omega()
        # ------------------------------------------ рисуем аккреционный диск ------------------------------------------
        self.draw_accretion_disc()

        # -----
        self.plot_proj_on_plane()

        # self.plot_plane()

    def plot_plane(self):
        normal = np.array([0, 0, 1])

        # a plane is a*x+b*y+c*z+d=0
        # [a,b,c] is the normal. Thus, we have to calculate
        # d and we're set

        # create x,y
        plane_max = self.R_disk_val * 1.5
        xx, yy = np.meshgrid(np.linspace(-plane_max, plane_max, 100),
                             np.linspace(-plane_max, plane_max, 100))

        # calculate corresponding z
        z = (-normal[0] * xx - normal[1] * yy) * 1. / normal[2]

        # contour_surf
        # mlab.mesh
        # mlab.surf
        # color=(0, 0, 1), transparent=True,
        self.scene.mlab.mesh(xx, yy, z, mode='2dsquare', transparent=True, opacity=0.4)

    def plot_proj_on_plane(self):
        # r_min = ksi
        # self.curr_configuration.top_column.inner_surface.theta_range

        R_disc = self.curr_configuration.R_disk / config.R_ns

        # r_min = self.curr_configuration.top_column.R_e / config.R_ns * np.sin(
        #     np.pi / 2 - np.deg2rad(self.curr_configuration.beta_mu)) ** 2

        # * np.cos(np.deg2rad(self.curr_configuration.beta_mu))
        phi = np.linspace(0, 2 * np.pi, 100)
        # - так как проекция на плоскость то х надо уменьшить
        x = R_disc * np.cos(phi) #  * np.cos(np.deg2rad(self.curr_configuration.beta_mu))
        y = R_disc * np.sin(phi)
        z = np.zeros_like(phi)

        tube_radius = 0.05
        color = (0, 0, 0)

        self.scene.mlab.plot3d(x, y, z, color=color, tube_radius=tube_radius)
        self.scene.mlab.quiver3d(2 * x[0], y[0], z[0], mode='2ddash', scale_factor=1, color=color)  # X

        flag = True
        if flag:
            phi = np.linspace(self.curr_configuration.top_column_for_plot.inner_surface.phi_range[0],
                              self.curr_configuration.top_column_for_plot.inner_surface.phi_range[-1], 100)

            R_column = self.curr_configuration.top_column_for_plot.inner_surface.surf_R_e / config.R_ns
            theta = self.curr_configuration.top_column_for_plot.inner_surface.theta_range[-1]
            r = R_column * np.sin(theta) ** 3  # * np.cos(theta)
            x = r * np.cos(phi)  # * np.cos(np.deg2rad(self.curr_configuration.beta_mu))
            y = r * np.sin(phi)
            z = np.zeros_like(phi)

            self.scene.mlab.plot3d(x, y, z, color=(0, 0, 0), tube_radius=tube_radius)
            self.scene.mlab.plot3d(-x, -y, z, color=(0, 0, 0), tube_radius=tube_radius)

            phi_center = phi[50]
            x = r * np.cos(phi_center)  # * np.cos(np.deg2rad(self.curr_configuration.beta_mu))
            y = r * np.sin(phi_center)
            z = 0
            scale = 1
            self.scene.mlab.quiver3d(scale * x, scale * y, scale * z, mode='2ddash', scale_factor=1, color=(0, 0, 0))
            self.scene.mlab.quiver3d(-scale * x, -scale * y, -scale * z, mode='2ddash', scale_factor=1, color=(0, 0, 0))


            x = r * np.cos(phi[0])
            y = r * np.sin(phi[0])
            z = 0

            self.scene.mlab.quiver3d(x, y, z, mode='2ddash', scale_factor=1, color=(0, 0, 0))
            self.scene.mlab.quiver3d(-x, -y, z, mode='2ddash', scale_factor=1, color=(0, 0, 0))

            x = r * np.cos(phi[-1])
            y = r * np.sin(phi[-1])
            z = 0

            self.scene.mlab.quiver3d(x, y, z, mode='2ddash', scale_factor=1, color=(0, 0, 0))
            self.scene.mlab.quiver3d(-x, -y, z, mode='2ddash', scale_factor=1, color=(0, 0, 0))

            # plot the surface
            # plt3d = plt.figure().gca(projection='3d')
            # plt3d.plot_surface(xx, yy, z)
            # plt.show()

            # self.check_data()

    def draw_mu_omega(self):
        radi = self.ksi_val * 3

        self.mu_vector = mlab.quiver3d(0, 0, radi, mode='2ddash', scale_factor=1, color=self.mu_vector_color)
        # self.mu_vector_1 = mlab.quiver3d(0, 0, -radi, mode='2ddash', scale_factor=1, color=self.mu_vector_color)

        omega_vector = matrix.get_cartesian_from_spherical(radi, np.deg2rad(-self.curr_configuration.beta_mu), 0)
        # omega_vector
        # self.omega_vector = mlab.plot3d([0, omega_vector[0]], [0, omega_vector[1]], [0, omega_vector[2]],
        #                                 color=self.omega_vector_color, tube_radius=self.omega_vector_tube_radius,
        #                                 tube_sides=4)
        # mu_vector
        # mlab.plot3d([0, 0], [0, 0], [0, 1.0], color=self.mu_vector_color, tube_radius=self.mu_vector_tube_radius,
        #             tube_sides=6)

        self.omega_vector = mlab.quiver3d(omega_vector[0], omega_vector[1], omega_vector[2], mode='2ddash',
                                          scale_factor=1, color=self.omega_vector_color)
        # self.omega_vector_1 = mlab.quiver3d(-omega_vector[0], -omega_vector[1], -omega_vector[2], mode='2ddash',
        #                                     scale_factor=1, color=self.omega_vector_color)

    def update_mu_omega(self):
        radi = self.ksi_val * 2

        omega_vector = matrix.get_cartesian_from_spherical(radi, np.deg2rad(-self.curr_configuration.beta_mu), 0)
        self.omega_vector.mlab_source.vectors = np.reshape([omega_vector[0], omega_vector[1], omega_vector[2]], (1, 3))
        # self.omega_vector_1.mlab_source.vectors = np.reshape([-omega_vector[0], -omega_vector[1], -omega_vector[2]],
        #                                                      (1, 3))

        self.mu_vector.mlab_source.vectors = np.reshape([0, 0, radi], (1, 3))
        # self.mu_vector_1.mlab_source.vectors = np.reshape([0, 0, -radi], (1, 3))

    def get_data_for_accretion_columns_with_mask(self, accretion_surface: accretingNS.Surface):
        '''теперь тоже с маской - чтобы сделать обрезы'''
        # передаем углы колонок, получаем сетку xyz в 3d
        # phi_range = np.linspace(0, 2 * np.pi, num=config.N_phi_accretion, endpoint=True)

        phi_range = accretion_surface.phi_range
        theta_range = accretion_surface.theta_range

        coef = 1
        if accretion_surface.surface_type == accretingNS.surface_surf_types['outer']:
            coef = (1 + self.curr_configuration.dRe_div_Re)

        radi = self.curr_configuration.top_column_for_plot.R_e / config.R_ns

        r, p = np.meshgrid(radi * coef * np.sin(theta_range) ** 2, phi_range)
        r1 = r * np.sin(theta_range)
        x = r1 * np.cos(p)
        y = r1 * np.sin(p)
        z = r * np.cos(theta_range)

        return x, y, z, accretion_surface.mask_array

    def get_data_for_accretion_columns_with_mask_for_calc(self, accretion_surface: accretingNS.Surface):
        '''теперь тоже с маской - чтобы сделать обрезы'''
        # передаем углы колонок, получаем сетку xyz в 3d
        # phi_range = np.linspace(0, 2 * np.pi, num=config.N_phi_accretion, endpoint=True)

        phi_range = accretion_surface.phi_range
        theta_range = accretion_surface.theta_range

        coef = 1

        radi = accretion_surface.surf_R_e / config.R_ns

        r, p = np.meshgrid(radi * coef * np.sin(theta_range) ** 2, phi_range)
        r1 = r * np.sin(theta_range)
        x = r1 * np.cos(p)
        y = r1 * np.sin(p)
        z = r * np.cos(theta_range)

        return x, y, z, accretion_surface.mask_array

    def get_data_for_magnet_lines_with_mask(self, magnet_surface: accretingNS.MagnetLine):
        '''смог реализовать только с помощью маски'''

        radi = self.R_e_val
        r, p = np.meshgrid(radi * np.sin(magnet_surface.theta_range) ** 2, magnet_surface.phi_range)
        r1 = r * np.sin(magnet_surface.theta_range)
        x = r1 * np.cos(p)
        y = r1 * np.sin(p)
        z = r * np.cos(magnet_surface.theta_range)

        return x, y, z, magnet_surface.mask_array

    def get_data_for_magnet_lines_with_mask_for_calc(self, magnet_surface: accretingNS.MagnetLine):
        '''смог реализовать только с помощью маски'''

        radi = self.curr_configuration.top_column_for_calc.R_e / config.R_ns
        r, p = np.meshgrid(radi * np.sin(magnet_surface.theta_range) ** 2, magnet_surface.phi_range)
        r1 = r * np.sin(magnet_surface.theta_range)
        x = r1 * np.cos(p)
        y = r1 * np.sin(p)
        z = r * np.cos(magnet_surface.theta_range)

        return x, y, z, magnet_surface.mask_array

    def get_data_for_magnet_lines_outer_with_mask(self, magnet_surface: accretingNS.MagnetLine):
        # просто удлинить в (1 + config.dRe_div_Re)
        x, y, z, mask_magnet_lines = self.get_data_for_magnet_lines_with_mask(magnet_surface)
        coef = (1 + self.curr_configuration.dRe_div_Re)
        return x * coef, y * coef, z * coef, mask_magnet_lines

    def get_data_for_NS(self):
        '''
        рисуем в масштабах магнитосферы R_e - так проще
        поэтому радиус звезды = np.sin(theta_range_column[0]) ** 2 то есть начало колонки
        '''

        theta_range = np.linspace(0, np.pi, num=config.N_theta_accretion, endpoint=True)
        phi_range = np.linspace(0, 2 * np.pi, num=config.N_phi_accretion, endpoint=True)

        u, v = np.meshgrid(phi_range, theta_range)
        r1 = 1
        # np.sin(np.pi / 2 - np.deg2rad(self.curr_configuration.beta_mu)) ** 2
        # print(r1)
        x = r1 * np.sin(v) * np.cos(u)
        y = r1 * np.sin(v) * np.sin(u)
        z = r1 * np.cos(v)

        return x, y, z

    def draw_NS(self):
        # x, y, z = self.get_data_for_NS(self.curr_configuration.top_column.inner_surface.theta_range)
        x, y, z = self.get_data_for_NS()
        self.NS = self.scene.mlab.mesh(x, y, z, color=(0, 0, 0))

    def update_NS(self):
        self.NS.mlab_source.trait_set(x=[0], y=[0], z=[0])
        self.NS = 0
        self.draw_NS()

    def get_data_for_accretion_disc(self, accretion_disc_con_val):

        # r_min = ksi
        # self.curr_configuration.top_column.inner_surface.theta_range
        r_min = self.curr_configuration.R_disk / config.R_ns

        # r_min = self.curr_configuration.top_column.R_e / config.R_ns * np.sin(
        #     np.pi / 2 - np.deg2rad(self.curr_configuration.beta_mu)) ** 2

        r_max = 1.25 * r_min
        r, phi = np.mgrid[r_min:r_max:100j, 0:2 * np.pi:100j]
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        z = accretion_disc_con_val * (x ** 2 + y ** 2) ** (1 / 2)
        # z, z1 = np.mgrid[-0.003:0.003:100j, -0.003:0.003:100j]
        return x, y, z

    def get_data_for_accretion_disc_side_surface(self, accretion_disc_con_val):
        # cylinder - боковая поверхность у диска
        # accretion_disc_con_val = насколько меняется радиус ??
        r_min = self.curr_configuration.R_disk / config.R_ns
        # self.curr_configuration.top_column.R_e / config.R_ns * np.sin(
        # np.pi / 2 - np.deg2rad(self.curr_configuration.beta_mu)) ** 2
        r_max = 1.25 * r_min
        radius = r_max
        z = np.linspace(-accretion_disc_con_val * radius, accretion_disc_con_val * radius, config.N_theta_accretion)
        phi = np.linspace(0, 2 * np.pi, config.N_phi_accretion)
        phi_grid, z = np.meshgrid(phi, z)
        x = radius * np.cos(phi_grid)
        y = radius * np.sin(phi_grid)

        return x, y, z

    def draw_accretion_disc(self):
        self.flag_accretion_disc_hide = False
        disc_color = (160, 82, 45)
        disc_color = tuple(rgb / 255 for rgb in disc_color)

        accretion_disc_con_val = 0.01
        x, y, z = self.get_data_for_accretion_disc(accretion_disc_con_val)

        # сначала отрисовали в СК beta_mu
        disc_opacity = 1
        self.accretion_disc_top = mlab.mesh(x, y, z, color=disc_color, opacity=disc_opacity)
        self.accretion_disc_bot = mlab.mesh(x, y, -z, color=disc_color, opacity=disc_opacity)

        # флаг позволяет сделать плоскость диска не омега а mu и обратно
        self.flag_accretion_disc_omega_mu = True  # True = omega

        # боковая поверхность диска (как боковая поверхность цилиндра)
        x, y, z = self.get_data_for_accretion_disc_side_surface(accretion_disc_con_val)
        self.accretion_disc_side_surface = mlab.mesh(x, y, z, color=disc_color, opacity=disc_opacity)

        # поворачиваем диск на -beta чтобы было перпендикулярно omega -> переходим в СК omega
        # deg ?? вроде да в градусах
        self.accretion_disc_rotate_angle = self.curr_configuration.beta_mu
        self.rotate_accretion_disc(-self.accretion_disc_rotate_angle)

    def update_accretion_disc(self):
        self.accretion_disc_top.mlab_source.trait_set(x=[0], y=[0], z=[0])
        self.accretion_disc_bot.mlab_source.trait_set(x=[0], y=[0], z=[0])
        self.accretion_disc_side_surface.mlab_source.trait_set(x=[0], y=[0], z=[0])
        if not self.flag_accretion_disc_hide:
            self.draw_accretion_disc()

    def rotate_accretion_disc(self, rotate_angle):
        self.accretion_disc_top.actor.actor.rotate_y(rotate_angle)
        self.accretion_disc_bot.actor.actor.rotate_y(rotate_angle)
        self.accretion_disc_side_surface.actor.actor.rotate_y(rotate_angle)

    def update_accretion_disc_rotate_angle(self):
        if self.flag_accretion_disc_omega_mu:
            self.rotate_accretion_disc(self.accretion_disc_rotate_angle)
            self.accretion_disc_rotate_angle = self.curr_configuration.beta_mu
            self.rotate_accretion_disc(-self.accretion_disc_rotate_angle)
        else:
            self.accretion_disc_rotate_angle = self.curr_configuration.beta_mu

    def draw_magnet_lines(self):
        # new
        # x, y, z, mask = get_data_for_magnet_lines_with_mask(self.curr_configuration.top_magnet_lines.phi_range,
        #                                                     self.curr_configuration.top_magnet_lines.theta_range,
        #                                                     self.curr_configuration.top_magnet_lines.mask_array)
        # self.get_data_for_magnet_lines_with_mask_for_calc(self.curr_configuration.top_magnet_lines_for_calc)
        x, y, z, mask = self.get_data_for_magnet_lines_with_mask(
            self.curr_configuration.top_magnet_lines_for_plot)  # top_magnet_lines_for_plot
        # если маска для магнитных линий вся из True почему то программа падает, для этого проверяем
        if not (mask == True).all():
            # flag_do_not_draw = True
            self.magnet_lines_top = self.scene.mlab.mesh(x, y, z, color=self.color_magnet_lines_top,
                                                         opacity=self.opacity_for_magnet_line,
                                                         representation='wireframe', mask=mask)
            self.magnet_lines_bot = self.scene.mlab.mesh(-x, -y, -z, color=self.color_magnet_lines_bot,
                                                         opacity=self.opacity_for_magnet_line,
                                                         representation='wireframe', mask=mask)
        else:
            self.magnet_lines_top = self.scene.mlab.mesh([[0]], [[0]], [[0]])
            self.magnet_lines_bot = self.scene.mlab.mesh([[0]], [[0]], [[0]])

        # self.mlab_source.trait_set(x=-x, y=-y, z=-z, color=self.color_accretion_column_bot_outer)

        # new
        x, y, z, mask = self.get_data_for_magnet_lines_outer_with_mask(
            self.curr_configuration.top_magnet_lines_for_plot_outer)
        # curr_configuration.top_magnet_lines_for_plot_outer

        if not (mask == True).all():
            self.magnet_lines_top_outer = self.scene.mlab.mesh(x, y, z,
                                                               color=self.color_magnet_lines_top_outer,
                                                               opacity=self.opacity_for_magnet_line,
                                                               representation='wireframe', mask=mask)
            self.magnet_lines_bot_outer = self.scene.mlab.mesh(-x, -y, -z,
                                                               color=self.color_magnet_lines_bot_outer,
                                                               opacity=self.opacity_for_magnet_line,
                                                               representation='wireframe', mask=mask)
        else:
            self.magnet_lines_top_outer = self.scene.mlab.mesh([[0]], [[0]], [[0]])
            self.magnet_lines_bot_outer = self.scene.mlab.mesh([[0]], [[0]], [[0]])

    def update_magnet_lines(self):
        # убираем старые линии
        self.magnet_lines_top.mlab_source.trait_set(x=[0], y=[0], z=[0])
        self.magnet_lines_bot.mlab_source.trait_set(x=[0], y=[0], z=[0])
        # убираем старые линии
        self.magnet_lines_top_outer.mlab_source.trait_set(x=[0], y=[0], z=[0])
        self.magnet_lines_bot_outer.mlab_source.trait_set(x=[0], y=[0], z=[0])
        self.draw_magnet_lines()

    def draw_accretion_columns(self):
        # внутренние
        x, y, z, mask = self.get_data_for_accretion_columns_with_mask(
            self.curr_configuration.top_column_for_plot.inner_surface)
        if not (mask == True).all():
            self.accretion_column_top = self.scene.mlab.mesh(x, y, z, color=self.color_accretion_column_top, mask=mask)
            self.accretion_column_bot = self.scene.mlab.mesh(-x, -y, -z, color=self.color_accretion_column_bot,
                                                             mask=mask)
        # внешние
        x, y, z, mask = self.get_data_for_accretion_columns_with_mask(
            self.curr_configuration.top_column_for_plot.outer_surface)
        if not (mask == True).all():
            # ---------верх
            self.accretion_column_top_outer = self.scene.mlab.mesh(x, y, z, color=self.color_accretion_column_top_outer,
                                                                   mask=mask)
            # ---------низ
            self.accretion_column_bot_outer = self.scene.mlab.mesh(-x, -y, -z,
                                                                   color=self.color_accretion_column_bot_outer,
                                                                   mask=mask)
        # middle
        x, y, z, mask = self.get_data_for_accretion_columns_with_mask_for_calc(
            self.curr_configuration.top_column_for_calc.outer_surface)
        if not (mask == True).all():
            # ---------верх
            self.accretion_column_top_center = self.scene.mlab.mesh(x, y, z, color=self.color_accretion_column_for_calc,
                                                                    mask=mask)
            # ---------низ
            self.accretion_column_bot_center = self.scene.mlab.mesh(-x, -y, -z,
                                                                    color=self.color_accretion_column_for_calc,
                                                                    mask=mask)

    def update_accretion_columns(self):
        # рисуем колонки
        # x, y, z = get_data_for_accretion_columns(self.curr_configuration.top_column.inner_surface.phi_range,
        #                                          self.curr_configuration.top_column.inner_surface.theta_range)
        # self.accretion_column_top.mlab_source.trait_set(x=x, y=y, z=z, color=self.color_accretion_column_top, mask=mask)
        # self.accretion_column_bot.mlab_source.trait_set(x=-x, y=-y, z=-z, color=self.color_accretion_column_bot,
        #                                                 mask=mask)

        self.accretion_column_top.mlab_source.trait_set(x=[0], y=[0], z=[0])
        self.accretion_column_bot.mlab_source.trait_set(x=[0], y=[0], z=[0])

        self.accretion_column_top_outer.mlab_source.trait_set(x=[0], y=[0], z=[0])
        self.accretion_column_bot_outer.mlab_source.trait_set(x=[0], y=[0], z=[0])

        self.accretion_column_top_center.mlab_source.trait_set(x=[0], y=[0], z=[0])
        self.accretion_column_bot_center.mlab_source.trait_set(x=[0], y=[0], z=[0])

        self.draw_accretion_columns()

    def view_phase(self, phase=0):
        # поворот на фазу; устанавливаем конфигурацию на эту фазу
        # начальный наблюдатель
        e_obs = matrix.get_cartesian_from_spherical(1, np.deg2rad(self.curr_configuration.theta_obs), 0)
        # поворот на фазу с помощью матрицы поворотов
        A_matrix_analytic = matrix.newMatrixAnalytic(0, config.betta_rotate, np.deg2rad(phase),
                                                     np.deg2rad(self.curr_configuration.beta_mu))
        e_obs_mu = np.dot(A_matrix_analytic, e_obs)  # переход в магнитную СК

        azimuth, elevation = matrix.vec_to_angles(e_obs_mu)
        roll_angle = mayavi_geometry.calculate_roll_angle(self.curr_configuration.theta_obs,
                                                          self.curr_configuration.beta_mu, e_obs_mu, phase)
        # print(roll_angle)

        # ax.view_init(90 - elevation / config.grad_to_rad, azimuth / config.grad_to_rad, roll=roll_angle)
        distance = self.slider_distance
        self.scene.mlab.view(azimuth=np.rad2deg(azimuth), elevation=np.rad2deg(elevation),
                             distance=distance, focalpoint=[0, 0, 0])
        # roll angle считаем для плоскости камеры - поэтому roll там
        # только здесь нашел про камеру https://docs.enthought.com/mayavi/mayavi/mlab_figures_decorations.html
        camera = self.scene.camera
        camera.roll(roll_angle)

    def update_curr_configuration(self, mu=None, theta_obs=None, beta_mu=None, mc2=None, a_portion=None,
                                  phi_0=None):
        if mu is not None:
            self.curr_configuration.mu = mu
        if theta_obs is not None:
            self.curr_configuration.theta_obs = theta_obs
        if beta_mu is not None:
            self.curr_configuration.beta_mu = beta_mu
        if mc2 is not None:
            self.curr_configuration.mc2 = mc2
        if a_portion is not None:
            self.curr_configuration.a_portion = a_portion
        if phi_0 is not None:
            self.curr_configuration.phi_0 = phi_0

        self.update_accretion_disc_rotate_angle()

        if theta_obs is None:
            self.curr_configuration = accretingNS.AccretingPulsarConfiguration(self.curr_configuration.mu,
                                                                               self.curr_configuration.theta_obs,
                                                                               self.curr_configuration.beta_mu,
                                                                               self.curr_configuration.mc2,
                                                                               self.curr_configuration.a_portion,
                                                                               self.curr_configuration.phi_0)
            if self.flag_draw_magnet_lines:
                self.update_magnet_lines()
            if mc2 is not None:
                self.update_accretion_disc()
            # upate
            # self.update_accretion_columns()
            self.update_all()
            self.update_mu_omega()
        phase = 360 * self.slider_phase
        self.view_phase(phase)

    def update_ksi_R_e(self):
        self.R_e_val = self.curr_configuration.top_column_for_plot.R_e / config.R_ns
        self.ksi_val = self.curr_configuration.top_column_for_plot.ksi_shock
        self.slider_distance = max(self.ksi_val * 4, self.slider_distance)
        self.R_disk_val = self.curr_configuration.R_disk / config.R_ns

    def update_all(self):
        self.update_accretion_columns()
        self.update_ksi_R_e()
        self.update_NS()
        self.update_accretion_disc()
        self.update_magnet_lines()

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
    def func_change_slider_theta_obs(self):
        self.update_curr_configuration(theta_obs=self.slider_theta_obs)
        # self.update_accretion_disc_rotate_angle()

    @on_trait_change('slider_beta_mu')
    def func_change_slider_betta_mu(self):
        self.update_curr_configuration(beta_mu=self.slider_beta_mu)

    @on_trait_change('slider_mc2')
    def func_change_slider_mc2(self):
        self.update_curr_configuration(mc2=self.slider_mc2)

    @on_trait_change('slider_a_portion')
    def func_change_slider_a_portion(self):
        self.update_curr_configuration(a_portion=self.slider_a_portion)

    @on_trait_change('slider_phi_0')
    def func_change_slider_fi_0(self):
        self.update_curr_configuration(phi_0=self.slider_phi_0)

    @on_trait_change('slider_phase, slider_distance')
    def func_change_slider_phase_slider_distance(self):
        phase = 360 * self.slider_phase
        self.view_phase(phase)

    @on_trait_change('button_magnet_line')
    def func_change_button_magnet_line(self):
        self.flag_draw_magnet_lines = not self.flag_draw_magnet_lines
        self.magnet_lines_top.visible = not self.magnet_lines_top.visible
        self.magnet_lines_top_outer.visible = not self.magnet_lines_top_outer.visible
        self.magnet_lines_bot.visible = not self.magnet_lines_bot.visible
        self.magnet_lines_bot_outer.visible = not self.magnet_lines_bot_outer.visible

    @on_trait_change('button_accr_disc')
    def func_change_button_accr_disc(self):
        if self.flag_accretion_disc_omega_mu:
            self.rotate_accretion_disc(self.accretion_disc_rotate_angle)
        else:
            self.rotate_accretion_disc(-self.accretion_disc_rotate_angle)
        self.flag_accretion_disc_omega_mu = not self.flag_accretion_disc_omega_mu

    @on_trait_change('button_hide_accr_disc')
    def func_change_button_hide_accr_disc(self):
        self.flag_accretion_disc_hide = not self.flag_accretion_disc_hide
        self.accretion_disc_top.visible = not self.accretion_disc_top.visible
        self.accretion_disc_bot.visible = not self.accretion_disc_bot.visible
        self.accretion_disc_side_surface.visible = not self.accretion_disc_side_surface.visible
        if not self.flag_accretion_disc_hide:
            self.update_accretion_disc()

    @on_trait_change('button_check_data')
    def func_change_button_check_data(self):
        self.try_check_data()

    @on_trait_change('button_animate')
    def anim(self):
        self.scene.reset_zoom()
        N = 100
        for i in range(N):
            self.slider_distance = 3.89764  # wtf???
            self.slider_phase = 1 * i / (N - 1)
            self.scene.save_png('mayavi_figs/' + f'anim{i:02d}.png')
            # mlab.figure(size = (1024,768), bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5))
            # mlab.save('pure_beauty.png')

    '''сделать видео:
    
    from mayavi import mlab
    from mayavi.core.lut_manager import LUTManager 
    from tvtk.util import ctf
    from moviepy import editor
    from tqdm import tqdm
    import os
    
    # Set misc variables for composing a video from the .pngs
    input_ext = ".png"
    output_name = "test.mp4"
    fps = 60
    codec = "libx264"
    # Get a list of image paths
    imgs = sorted(os.listdir(log_dir))
    images = [img for img in imgs if img.endswith(input_ext)]
    image_paths = [os.path.join(log_dir, img) for img in images]
    # Ganarate a video and save it on local
    video_name = os.path.join(log_dir, output_name)
    clip = editor.ImageSequenceClip(image_paths, fps=fps)
    clip.write_videofile(video_name, codec=codec, fps=fps)
    '''

    # @mlab.animate
    # def anim(self):
    #     for i in range(10):
    #         self.slider_phase = 2 * i / 10
    #         yield

    # the layout of the dialog created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene), height=250, width=300, show_label=False),
                VGroup(
                    HGroup(
                        '_', 'slider_theta_obs', 'slider_beta_mu', 'slider_mc2', 'slider_a_portion', 'slider_phi_0'
                    ),
                    HGroup('_', 'slider_phase', 'slider_distance'),
                    HGroup('button_magnet_line', 'button_hide_accr_disc', 'button_accr_disc', 'button_animate'),
                    HGroup('R_e_val', 'ksi_val', 'R_disk_val')
                )
                )


def plot_main(mu, theta_obs, beta_mu, mc2, a_portion, phi_0):
    curr_configuration = accretingNS.AccretingPulsarConfiguration(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
    visualization = Visualization(curr_configuration)
    visualization.configure_traits()
    # visualization.anim()


if __name__ == "__main__":
    mu = 0.1e31

    theta_obs = 10
    beta_mu = 60
    mc2 = 100
    a_portion = 0.22
    phi_0 = 0

    theta_obs = 10
    beta_mu = 20
    mc2 = 60
    a_portion = 0.66
    phi_0 = 0

    theta_obs = 40
    beta_mu = 30
    mc2 = 60
    a_portion = 0.2
    phi_0 = 60

    # theta_obs = 80
    # beta_mu = 20
    # mc2 = 100
    # a_portion = 0.82
    # phi_0 = 80

    plot_main(mu, theta_obs, beta_mu, mc2, a_portion, phi_0)
