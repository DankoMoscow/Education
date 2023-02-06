import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import os
from raschet import *


num_steps = 100  # количество шагов
l = np.empty(num_steps + 2, dtype=np.int16)

proc_time = 100*3600 # 100 часов сушки
c_bound = 0.
c_init = 1.

density_ips = 785.1  # кг/м3
density_co2 = 468  # кг/м3
class scd_apparatus():
    def __init__(self, width, length, height, diff_coef, key, key_sch, number_samples):
        self.width = width
        self.length = length
        self.height = height
        self.diff_coef = diff_coef
        self.key = key
        self.key_sch = key_sch
        self.number_samples = number_samples
        self.c_init_list = np.zeros(num_steps + 2)
        self.matrix_of_c = np.zeros((n_t, len(self.c_init_list)))
        self.list_of_mass = np.zeros(n_t)
        self.c_app = np.zeros(n_t)
        self.density_ips = 785.1 #кг/м3
        self.density_co2 = 468 #кг/м3

        for i in range(num_steps + 2):
            self.c_init_list[i] = c_init
            if i == num_steps + 1:
                self.c_init_list[i] = 0

    def __str__(self):
        print(
            f'width: {self.width}, length: {self.length}, diff_coef: {self.diff_coef}, number_samples: {self.number_samples}')

    def fick_conc(self, diff_coefficient, density, c, c_bound, dr, dt, r):
        global sverka_method
        sverka_method = 0
        stab_cond = dt / dr ** 2  # условие устойчивости
        alfa_i = np.zeros(num_steps + 2)
        betta_i = np.zeros(num_steps + 2)
        alfa_i[0] = 1  # прогоночный коэффициент альфа на нулевом шаге
        betta_i[0] = 0  # прогоночный коэффициент бетта на нулевом шаге
        c_temp = []
        c_temp = np.copy(c)
        y_ips = np.linspace(1, 0, num_steps + 2)

        if self.key == 'one_dim':
            if self.key_sch == 'explicit':
                if stab_cond > 1 / (2 * self.diff_coef):
                    sverka_method = 5
                    pass
                else:
                    c_temp[1:-1] = c[1:-1] + self.diff_coef * (dt / dr ** 2) * (c[2:] - 2 * c[1:-1] + c[0: -2])
                    c_temp[-1] = c_bound
                    c_temp[0] = c_temp[1]
                    return c_temp

            elif self.key_sch == 'implicit':
                a = -self.diff_coef * dt / (dr) ** 2  # коэффициент 'a'
                b = 1 + 2 * self.diff_coef * dt / (dr) ** 2  # коэффицент b
                c_koef = -self.diff_coef * dt / (dr) ** 2  # коэффициент c

                for i in range(1, len(l)):
                    alfa_i[i] = (-a) / (b + c_koef * alfa_i[i - 1])
                    betta_i[i] = (c[i] - c_koef * betta_i[i - 1]) / (b + c_koef * alfa_i[i - 1])
                alfa_i[-1] = 0
                betta_i[-1] = 0
                c_temp[1:-1] = c[2:] * alfa_i[1:-1] + betta_i[1:-1]

                c_temp[-1] = c_bound
                c_temp[0] = c_temp[1]
                return c_temp

        elif self.key == 'cyl':
            if self.key_sch == 'explicit':
                if 2 * self.diff_coef * stab_cond <= 1:
                    sverka_method = 5
                    pass
                else:
                    diff_coef = []
                    density_list = []
                    # for i in range(len(r)):
                    #     # TODO доделать уравнение с плотностью
                    #     # TODO пересмотреть коэффициенты
                    #     x_ips = (y_ips[i] * M_ips) / ((1 - y_ips[i]) * M_co2 + y_ips[i] * M_ips)
                    #     x_ips_list.append(x_ips)
                    diff_coefficient[1:-1] = (D_coef_co2_ips ** x_ips_list[1:-1] * D_coef_ips_co2 ** (1 - x_ips_list[1:-1]))
                    for i in range(len(r)):
                        density[i] = (y_ips[i] / density_ips + (1 - y_ips[i]) / density_co2) ** -1
                        density_list.append(density)
                    print('x_ips_list', x_ips_list)
                    print('diff_coefficient', diff_coefficient)
                    #print('density', density)

                    density[0] = density[1]
                    density[-1] = density_co2


                    diff_coefficient[-1] = D_coef_co2_ips
                    diff_coefficient[0] = diff_coefficient[1]
                    for i in range(len(r)):
                        c_temp[1:-1] = c[1:-1] + porosity * dt / (tau_izv * dr * r[i]) * \
                                       ((diff_coefficient[2:] - diff_coefficient[1:-1]) / 2 * (
                                                   density[2:] - density[1:-1]) / 2 * (r[2:] + r[i]) / 2 * (
                                                    y_ips[2:] - y_ips[1:-1]) / dr - (
                                                    diff_coefficient[0:-2] - diff_coefficient[1:-1]) / 2 *
                                        (density[0:-2] - density[1:-1]) / 2 * (r[0:-2] + r[i]) / 2 * (
                                                    y_ips[1:-1] - y_ips[0:-2]) / dr)
                        if i==1:
                            print('c_temp', c_temp)
                            print('c_temp[1:-1]', c_temp[1:-1])
                    c_temp[-1] = c_bound
                    c_temp[0] = c_temp[1]
                    return c_temp

            elif self.key_sch == 'implicit':  # должна быть абсолютно устойчива это с ЛКР
                for j in range(1, len(r)):
                    a = -self.diff_coef * dt / (dr) ** 2
                    b = 1 + 2 * dt * self.diff_coef / (dr) ** 2 - dt * self.diff_coef / (dr * r[j])
                    c_koef = - dt * self.diff_coef / (dr) ** 2 + dt * self.diff_coef / (dr * r[j])
                for i in range(1, len(l)):
                    alfa_i[i] = (-a) / (b + c_koef * alfa_i[i - 1])
                    betta_i[i] = (c[i] - c_koef * betta_i[i - 1]) / (b + c_koef * alfa_i[i - 1])
                alfa_i[-1] = 0
                betta_i[-1] = 0
                c_temp[1:-1] = c[2:] * alfa_i[1:-1] + betta_i[1:-1]
                c_temp[-1] = c_bound
                c_temp[0] = c_temp[1]
                return c_temp

        elif self.key == 'sphere':
            if self.key_sch == 'explicit':
                if 2 * self.diff_coef * stab_cond <= 1:  # diff_coef*(dt/dr + 2 * stab_cond) > 1:
                    sverka_method = 5
                    pass
                else:
                    for i in range(len(r)):
                        c_temp[1:-1] = c[1:-1] + self.diff_coef * dt * (
                                    (c[2:] - 2 * c[1:-1] + c[0:-2]) / dr ** 2 +  (c[2:] - c[0:-2]) / (
                                        r[i] * dr))  # явная разностная схема с ЦКР работает
                    c_temp[-1] = c_bound
                    c_temp[0] = c_temp[1]
                    return c_temp

            elif self.key_sch == 'implicit':
                for j in range(1, len(r)):
                    a = -self.diff_coef * dt / (dr) ** 2 - (1 / r[j]) * self.diff_coef * dt / dr
                    b = 1 + 2 * dt * self.diff_coef / (dr) ** 2
                    c_koef = -dt * self.diff_coef / (dr) ** 2 + (1 / r[j]) * self.diff_coef * dt / dr

                for i in range(1, len(l)):
                    alfa_i[i] = (-a) / (b + c_koef * alfa_i[i - 1])
                    betta_i[i] = (c[i] - c_koef * betta_i[i - 1]) / (b + c_koef * alfa_i[i - 1])

                alfa_i[-1] = 0
                betta_i[-1] = 0
                c_temp[1:-1] = c[2:] * alfa_i[1:-1] + betta_i[1:-1]
                c_temp[-1] = c_bound
                c_temp[0] = c_temp[1]
                return c_temp

    def fick_mass(self, c, length, width):
        m = 0.
        for i in range(1, len(c)):
            if self.key == 'sphere':
                m += c[i] * 4 / 3 * np.pi * ((i * dr) ** 3 - ((i - 1) * dr) ** 3)
            elif self.key == 'cyl':
                m += c[i] * length * np.pi * ((i * dr) ** 2 - ((i - 1) * dr) ** 2)
            elif self.key == 'one_dim':
                m += c[i] * ((2 * i * dr) - ((i - 1) * 2 * dr)) * self.length * self.width

        return m

    def time_iteration(self,diffusion_list_init, density_init_list, c_init_list, volume, flowrate, n_t, dt, dr, key, key_sch):
        global method_value
        method_value = 0  # костыль для определения и не вылетания объёма аппарата
        c_app = np.zeros(n_t)
        residence_time = volume / flowrate

        c_matrix = np.zeros((n_t, len(c_init_list)))
        density_matrix = np.zeros((n_t, len(density_init_list)))
        diffusion_matrix = np.zeros((n_t, len(diffusion_list_init)))
        mass_list = np.zeros(n_t)
        c_matrix[0] = c_init_list
        density_matrix[0] = density_init_list
        diffusion_matrix[0] = diffusion_list_init
        mass_list[0] = self.fick_mass(c_matrix[0], self.length, self.width)
        c_app[0] = 0.
        for i in range(1, n_t):
            c_bound = c_app[i - 1]
            c_matrix[i] = self.fick_conc(diffusion_matrix[i-1], density_matrix[i-1], c_matrix[i - 1], c_bound, dr, dt, r)

            if volume < self.number_samples * 4 / 3 * np.pi * (R) ** 3 and self.key == 'sphere':
                method_value = 5
                pass

            elif volume < self.number_samples * self.length * np.pi * R ** 2 and self.key == 'cyl':
                method_value = 5
                pass

            elif volume < self.number_samples * self.length * self.width * R and self.key == 'one_dim':
                method_value = 5
                pass

            else:
                if method_value != 5:
                    mass_list[i] = self.fick_mass(c_matrix[i], self.length, self.width)
                    delta_mass = - self.number_samples * (mass_list[i] - mass_list[i - 1])
                    c_app[i] = self.ideal_mixing(c_app[i - 1], 0, residence_time, dt, volume, delta_mass)

        return c_matrix, mass_list, c_app

    def ideal_mixing(self, c, c_inlet, residence_time, dt, volume, delta_mass):
        c_mixing = c + dt / residence_time * (c_inlet - c) + dt * delta_mass / volume
        return c_mixing


def main(width, length, height, volume, flowrate, dt, diff_coef, number_samples, value, key_sch, working, working_scheme):
    global n_t, R, dr, c_r, r
    R = height / 2  # meters
    dr = R / num_steps  # шаг по радиусу meters
    #n_t = int(proc_time / dt) + 1  # количество шагов с учетом нулевого шага
    n_t = 1
    c_init_list = np.zeros(num_steps + 2)
    density_init_list = np.zeros(num_steps+2)
    diffusion_list_init = np.zeros(num_steps+2)


    c_r = np.zeros(num_steps + 2)
    r = np.linspace(0, R, num_steps + 2)
    # создаю равномерный список изменения TODO передаю его в fick_conc
    for i in range(num_steps + 2):
        c_init_list[i] = c_init
        density_init_list[i] = density_ips
        diffusion_list_init[i] = D_coef_ips_co2
        if i == num_steps + 1:
            c_init_list[i] = 0
            density_init_list[i] = density_co2
            diffusion_list_init[i] = D_coef_co2_ips

    object1 = scd_apparatus(width, length, height, diff_coef, value, key_sch, number_samples)
    object1.__str__()

    print('n_t:', n_t, 'proc_time:', proc_time, 'variable of item',value)
    time = np.linspace(0, proc_time, n_t)
    value = ['one_dim', 'cyl', 'sphere']
    key_sch = ['explicit', 'implicit']
    for i in value:
        for j in key_sch:
            matrix_of_c, list_of_mass, c_app = object1.time_iteration(diffusion_list_init, density_init_list, c_init_list, volume, flowrate, n_t, dt, dr, key = i, key_sch = j)

    #print(sverka_method)
    return matrix_of_c, list_of_mass, c_app, time, i, r, method_value, sverka_method
