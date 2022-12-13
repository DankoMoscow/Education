import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import os

file_path = os.path.dirname(os.path.abspath(__file__))
img_path = os.path.join(file_path, 'Images')

num_steps = 100  # количество шагов
l = np.empty(num_steps + 2, dtype=np.int16)
height = 0.02  # высота образца meters
R = height / 2  # meters
dr = R / num_steps  # шаг по радиусу meters
proc_time = 100000
# dt = 10  # шаг по времени seconds
# n_t = int(proc_time / dt) + 1  # количество шагов с учетом нулевого шага

c_bound = 0.
c_init = 1.

c_r = np.zeros(num_steps + 2)
r = np.linspace(0, R, num_steps + 2)
c_init_list = np.zeros(num_steps + 2)


class scd_apparatus():
    def __init__(self, width, length, diff_coef, key, key_sch, number_samples):
        self.width = width
        self.length = length
        self.diff_coef = diff_coef
        self.key = key
        self.key_sch = key_sch
        self.number_samples = number_samples

        self.matrix_of_c = np.zeros((n_t, len(c_init_list)))
        self.list_of_mass = np.zeros(n_t)
        self.c_app = np.zeros(n_t)

        self.c_init_list = np.zeros(num_steps + 2)
        for i in range(num_steps + 2):
            self.c_init_list[i] = c_init
            if i == num_steps + 1:
                self.c_init_list[i] = 0

    def __str__(self):
        print(
            f'width: {self.width}, length: {self.length}, diff_coef: {self.diff_coef}, number_samples: {self.number_samples}')

    def fick_conc(self, c, c_bound, dr, dt, r):
        global sverka_method
        sverka_method = 0
        stab_cond = dt / dr ** 2  # условие устойчивости
        alfa_i = np.zeros(num_steps + 2)
        betta_i = np.zeros(num_steps + 2)
        alfa_i[0] = 1  # прогоночный коэффициент альфа на нулевом шаге
        betta_i[0] = 0  # прогоночный коэффициент бетта на нулевом шаге
        c_temp = []
        c_temp = np.copy(c)

        if self.key == 'one_dim':
            if self.key_sch == 'explicit':
                if stab_cond > 1 / (2 * self.diff_coef):
                    sverka_method = 5
                    #print('Не выполняется условие устойчивости')
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
                if 2 * self.diff_coef * stab_cond <= 1:  # (2*diff_coef*dt-dt)/ dr**2 <= 1:     #diff_coef*(dt/dr + 2 * stab_cond) > 1:страница 84 методички
                    sverka_method = 5
                    #print('Не выполняется условие устойчивости')
                    pass
                else:
                    for i in range(len(r)):
                        c_temp[1:-1] = c[1:-1] + self.diff_coef * dt * (
                                    (c[2:] - 2 * c[1:-1] + c[0:-2]) / dr ** 2 + 1 / r[i] * (c[1:-1] - c[0:-2]) / dr)
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
                    #print('Не выполняется условие устойчивости')
                    pass
                else:
                    for i in range(len(r)):
                        c_temp[1:-1] = c[1:-1] + self.diff_coef * dt * (
                                    (c[2:] - 2 * c[1:-1] + c[0:-2]) / dr ** 2 + (c[2:] - c[0:-2]) / (
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
            else:
                print("Key input error!")
                return None
        return m

    def time_iteration(self, c_init_list, volume, flowrate, n_t, dt, dr, key, key_sch):
        global method_value
        method_value = 0  # костыль для определения и не вылетания объёма аппарата
        c_app = np.zeros(n_t)
        # volume = 0.00025 #кубические метры из диссертации
        # flow_rate = 0.0000017 #кубических метров в секунду из диссертации
        residence_time = volume / flowrate

        c_matrix = np.zeros((n_t, len(c_init_list)))
        mass_list = np.zeros(n_t)
        c_matrix[0] = c_init_list

        mass_list[0] = self.fick_mass(c_matrix[0], self.length, self.width)
        c_app[0] = 0.
        for i in range(1, n_t):
            c_bound = c_app[i - 1]
            c_matrix[i] = self.fick_conc(c_matrix[i - 1], c_bound, dr, dt, r)

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


"""    def plot_conc(self, r_list, time, c_list):
        time_ratio = 100 #с какой частотой писать легенду для графика
        plt.figure()
        c_list = c_list.T
        conc_fig = plt.plot(r_list, c_list[:, ::time_ratio])
        plt.legend(iter(conc_fig), time[::time_ratio], loc = 1, fontsize = 8)
        plt.title('Concentration change profile of alcohol', fontsize=18)
        plt.xlabel('Radius, m')
        plt.grid(True)
        plt.ylabel('Concentration of alcohol')
        return

    def plot_mass(self, r_list, mass_list, legend_key):
        plt.figure(2)
        plt.plot(r_list, mass_list, label=legend_key)
        plt.legend()
        plt.title('Mass change profile of alcohol', fontsize=18)
        plt.xlabel('Time, second')
        plt.grid(True)
        plt.ylabel('Mass of alcohol, kg')
        return

    def plot_3D(self, r_list, time, c_list, name):
        fig = plt.figure(figsize=(7, 4))
        xgrid, ygrid = np.meshgrid(r_list, time)
        zgrid = c_list
        ax_3d = Axes3D(fig)
        fig.add_axes()
        ax_3d.plot_surface(xgrid, ygrid, zgrid, cmap=cm.jet)
        plt.title(('3D concentration ' + name), fontsize=12)
        plt.xlabel('Radius, m')
        plt.grid(True)
        plt.ylabel('Time, second')
        return

    def ideal_mixing_plot(self, time, c_mixing):
        plt.figure()
        plt.plot(time, c_mixing)
        plt.title('Concentration change profile of alcohol for ideal mixing', fontsize=16)
        plt.xlabel('Time, second')
        plt.grid(True)
        plt.ylabel('Concentration of alcohol')
        return    """


def main(width, length, volume, flowrate, dt, diff_coef, number_samples, value, key_sch, working, working_scheme):
    global n_t

    n_t = int(proc_time / dt) + 1  # количество шагов с учетом нулевого шага
    c_init_list = np.zeros(num_steps + 2)
    for i in range(num_steps + 2):
        c_init_list[i] = c_init
        if i == num_steps + 1:
            c_init_list[i] = 0

    # img_path_key = os.path.join(img_path, str(i))

    object1 = scd_apparatus(width, length, diff_coef, value, key_sch, number_samples)
    object1.__str__()

    print('n_t:', n_t, 'proc_time:', proc_time, 'variable of item',value)
    time = np.linspace(0, proc_time, n_t)
    r_list = np.linspace(0, R, num_steps + 2)
    value = ['one_dim', 'cyl', 'sphere']
    key_sch = ['explicit', 'implicit']
    for i in value:
        for j in key_sch:
            matrix_of_c, list_of_mass, c_app = object1.time_iteration(c_init_list, volume, flowrate, n_t, dt, dr, key = i, key_sch = j)

    """object1.plot_mass(time, list_of_mass, i)
    plt.savefig("plot_mass.png")
    object1.plot_conc(r_list, time-1, matrix_of_c)
    plt.savefig("plot_conc.png")
    object1.ideal_mixing_plot(time, c_app)
    plt.savefig("plot_mixing.png")
    object1.plot_3D(r_list, time, matrix_of_c, str(i))
    plt.savefig("plot_3D.png")"""
    print(sverka_method)
    return matrix_of_c, list_of_mass, c_app, time, i, r_list, method_value, sverka_method
