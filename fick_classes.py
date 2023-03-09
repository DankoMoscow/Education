import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import os
import math
from scipy import optimize
from scipy.optimize import minimize

target_material_porosity = 120 #кг/м3
porosity = 0.95
tau_izv = 11 #5.5 #извилистость пор

Vcr_ips = 220 #см3/моль - критический молярный объём ips
Vcr_co2 = 94.07
Pcr_ips = 47.62 #критическаое давлние ипс бар
Pcr_co2 = 73.74

T = 333.15 #температура в кельвинах
P = 120 #давление в барах

PSI_ips = 0.665 #ацентрический фактор ипс
PSI_co2 = 0.225

Par_ips = 164.4 #парахора ипс
Par_co2 = 44.8

ips_volume_norm = 0.285*Vcr_ips**1.048 #нормальный молярынй объём ипс
co2_volume_norm = 0.285*Vcr_co2**1.048

M_ips = 60.09 #г на см3
M_co2 = 44.01

Tcr_ips = 508.3 #критическая температура ипс
Tcr_co2 = 304.12

a1 = 8.07652E+15
a2 = 8.91648E+13
a3 = -2.89429E+13
a4 = 89749140000
a5 = 230022.4
dynamic_viscosity_ips = (a1 + a2 * P / 10) / (a3 + a4 * T + a5 * T ** 3 + P / 10) / 1000000

D_coef_co2_ips = 8.93*10**(-12)*(co2_volume_norm/ips_volume_norm**2)**(1/6)*(Par_ips/Par_co2)**0.6*T/(dynamic_viscosity_ips*1000)    #коэфф диффузии co2 in ips

density_co2_coef=14258.88822-84.97903074*T+11.30377033*P+0.017974145*T*P+0.135119422*T**2-0.071358164*P**2-4.73474E-05*T**3+0.000110024*P**3
koef = Tcr_co2*Vcr_co2/(1000*M_co2)
A_coef = 14.882+5.908*koef+2.0821*koef**2
Vmol_co2 = M_co2/density_co2_coef * 1000 # Молярный объем диоксида углерода при норм температуре кипения, см3/моль
Vv_co2=Vmol_co2/Vcr_co2

D_coef_ips_co2 = A_coef * 10**(-9)*(T/M_ips)**0.5*math.exp(-0.3887/(Vv_co2-0.23))  #коэффициент диффузии ипс в CO2

y_ips = np.linspace(1, 0, 102)

density_ips = 785.1  # кг/м3
density_co2 = 468  # кг/м3


file_path = os.path.dirname(os.path.abspath(__file__))
img_path = os.path.join(file_path, 'Images')

num_steps = 100  # количество шагов
l = np.empty(num_steps + 2, dtype=np.int16)

proc_time = 15*360
c_bound = 0.
c_init = 10.


class scd_apparatus():
    def __init__(self,volume, flowrate, width, length, height, diff_coef, key, key_sch, number_samples):
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
        self.volume = volume
        self.flowrate =flowrate

        if self.key == 'one_dim':
            self.V_gels = self.height * self.width * self.length * self.number_samples

        elif self.key == 'cyl':
            self.V_gels = self.number_samples * self.length * np.pi * R ** 2

        elif self.key == 'sphere':
            self.V_gels = self.number_samples * 4 / 3 * np.pi * (R) ** 3

        self.V_apparat_free = self.volume - self.V_gels #объём свободного пространства аппарата
        self.m_ips_gramms = self.V_gels * target_material_porosity * density_ips * 1000 #масса ипс

        # for i in range(num_steps + 2):
        #     self.c_init_list[i] = c_init
        #     if i == num_steps + 1:
        #         self.c_init_list[i] = 0

    def __str__(self):
        print(
            f'width: {self.width}, length: {self.length}, diff_coef: {self.diff_coef}, number_samples: {self.number_samples}')
        print('Объём гелей', self.V_gels)
    """
    эта функция отвечает за создание коэффициента диффузии в зависимости от массовой доли
    """
    def diffusion(self, x):
        diff_coef = (D_coef_co2_ips ** x) * (D_coef_ips_co2 ** (1 - x))
        return diff_coef

    def density(selfself, ro):
        density = (ro / density_ips + (1 - ro) / density_co2) ** -1
        return density

    def golden_method(self, x):
        numerator = abs(
            (self.V_apparat_free + self.V_gels * target_material_porosity) *
            ((x / density_ips + (1 - x) / density_co2) ** -1) * x
            - (self.m_ips_gramms / 1000))
        return numerator / (self.m_ips_gramms / 1000)

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
                if dt> dr/(2 *self.diff_coef*porosity/tau_izv):
                #if 2 * self.diff_coef * stab_cond <= 1:  # (2*diff_coef*dt-dt)/ dr**2 <= 1:     #diff_coef*(dt/dr + 2 * stab_cond) > 1:страница 84 методички
                    sverka_method = 5
                    pass
                else:
                    for i in range(1, len(r) - 1):
                        # y_start = optimize.golden(f, maxiter=10000)  # массовая доля, кг/кгсм
                        #
                        # V_ips = density_co2 * y_start / (density_ips * (1 - y_start) + density_co2 * y_start)

                        c_temp[1:-1] = c[1:-1] + self.diff_coef * (dt / dr ** 2) * (c[2:] - 2 * c[1:-1] + c[0: -2])
                        # c_temp[1:-1] = c[1:-1] + porosity * dt / (tau_izv * dr * r[i]) * ((diff_coef_now * density_now * r[i] * (y_ips[i] - y_ips[i-1])/ dr) - \
                        #                                                                    (diff_coef_past * density_past * r[i-1] * (y_ips[i] - y_ips[i-1])/ dr))
                        # c_temp[1:-1] = c[1:-1] + porosity * dt / (tau_izv * dr * r[i]) * \
                        #            ((diff_coef_fut + diff_coef_now) / 2  * (density_fut + density_now) /2 * (r[i+1] + r[i]) / 2 * (
                        #                     y_ips[i+1] - y_ips[i]) / dr - (
                        #                     diff_coef_past + diff_coef_now) / 2 * (density_past + density_now) /2 * (r[i-1] + r[i]) / 2 * (
                        #                     y_ips[i] - y_ips[i-1]) / dr)

                    c_temp[-1] = c_bound
                    c_temp[0] = c_temp[1]

                    return c_temp

            elif self.key_sch == 'implicit':  # должна быть абсолютно устойчива это с ЛКР
                k_coef = dt * porosity / (tau_izv * dr)
                # TODO надо как-то подумать, чтобы диапазон был до len(r)
                for i in range(1, len(r) - 1):
                    x_ips_past = (y_ips[i - 1] * M_ips) / ((1 - y_ips[i - 1]) * M_co2 + y_ips[i - 1] * M_ips)
                    x_ips_now = (y_ips[i] * M_ips) / ((1 - y_ips[i]) * M_co2 + y_ips[i] * M_ips)
                    x_ips_fut = (y_ips[i + 1] * M_ips) / ((1 - y_ips[i + 1]) * M_co2 + y_ips[i + 1] * M_ips)

                    density_past = self.density(y_ips[i - 1])
                    density_now = self.density(y_ips[i])
                    density_fut = self.density(y_ips[i + 1])

                    diff_coef_past = self.diffusion(x_ips_past)
                    diff_coef_now = self.diffusion(x_ips_now)
                    diff_coef_fut = self.diffusion(x_ips_fut)

                    a_coef = - k_coef / r[i] * ((diff_coef_fut + diff_coef_now) / 2 * (r[i+1] + r[i]) / 2 * (density_fut +  density_now) / 2) / dr
                    b_coef = 1 + k_coef / r[i] * ((diff_coef_fut + diff_coef_now) / 2 * (r[i+1] + r[i]) / 2 * (density_fut +  density_now) / 2) / dr  +  ((diff_coef_past + diff_coef_now) / 2 * (r[i-1] + r[i]) / 2 * (density_past +  density_now) / 2) / dr
                    c_koef = -((diff_coef_past + diff_coef_now) / 2 * (r[i-1] + r[i]) / 2 * (density_past +  density_now) / 2) / dr

                # for j in range(1, len(r)):
                #     a = -self.diff_coef * dt / (dr) ** 2
                #     b = 1 + 2 * dt * self.diff_coef / (dr) ** 2 - dt * self.diff_coef / (dr * r[j])
                #     c_koef = - dt * self.diff_coef / (dr) ** 2 + dt * self.diff_coef / (dr * r[j])

                for i in range(1, len(l) - 1):
                    alfa_i[i] = (-a_coef) / (b_coef + c_koef * alfa_i[i - 1])
                    betta_i[i] = (c[i] - c_koef * betta_i[i - 1]) / (b_coef + c_koef * alfa_i[i - 1])

                # TODO тут проверить
                alfa_i[-1] = 0
                betta_i[-1] = 0
                #c_temp[1:-1] = c[2:] * alfa_i[1:-1] + betta_i[1:-1]
                c_temp = c[2:] * alfa_i[1:-1] + betta_i[1:-1]
                c_temp.append(c_bound)
                #c_temp[-1] = c_bound
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

    def time_iteration(self, c_init_list, volume, flowrate, n_t, dt, dr, key, key_sch):
        global method_value
        method_value = 0  # костыль для определения и не вылетания объёма аппарата
        c_app = np.zeros(n_t)
        residence_time = volume / flowrate

        c_matrix = np.zeros((n_t, len(c_init_list)))
        mass_list = np.zeros(n_t)
        c_matrix[0] = c_init_list

        mass_list[0] = self.fick_mass(c_matrix[0], self.length, self.width)
        c_app[0] = 0.
        for i in range(1, n_t):
            c_bound = c_app[i - 1]
            c_matrix[i] = self.fick_conc(c_matrix[i - 1], c_bound, dr, dt, r)

            if volume < self.V_gels:
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
    n_t = int(proc_time / dt) + 1  # количество шагов с учетом нулевого шага
    c_init_list = np.zeros(num_steps + 2)

    c_r = np.zeros(num_steps + 2)
    r = np.linspace(0, R, num_steps + 2)

    for i in range(num_steps + 2):
        c_init_list[i] = c_init
        if i == num_steps + 1:
            c_init_list[i] = 0


    object1 = scd_apparatus(volume, flowrate, width, length, height, diff_coef, value, key_sch, number_samples)
    object1.__str__()

    print('n_t:', n_t, 'proc_time:', proc_time, 'variable of item',value)
    time = np.linspace(0, proc_time, n_t)
    value = ['one_dim', 'cyl', 'sphere']
    key_sch = ['explicit', 'implicit']
    for i in value:
        for j in key_sch:
            matrix_of_c, list_of_mass, c_app = object1.time_iteration(c_init_list, volume, flowrate, n_t, dt, dr, key = i, key_sch = j)

    print(sverka_method)
    return matrix_of_c, list_of_mass, c_app, time, i, r, method_value, sverka_method
