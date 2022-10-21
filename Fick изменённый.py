import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


num_steps = 100 # количество шагов
l = np.empty(num_steps+2, dtype=np.int16)
height = 0.02   # высота образца meters
R = height / 2   # meters
dr = R / num_steps  # шаг по радиусу meters
proc_time = 100000
dt = 50  # шаг по времени seconds
n_t = int(proc_time / dt) + 1 # количество шагов с учетом нулевого шага
diff_coef = 1e-8
width = 0.1 #metrs
length = 0.5 #metrs


c_bound = 0.
c_init = 1.

c_r = np.zeros(num_steps + 2)
r = np.linspace(0, R, num_steps + 2)


#key = input('Which equation is considered: one_dim, cyl or  sphere? ' )
#key_sch = input('what scheme will you use: explicit or implicit: ' )

def fick_conc(c, key, key_sch, c_bound, dr, dt, r):
    stab_cond = dt/dr**2   #условие устойчивости
    alfa_i = np.zeros(num_steps + 2)
    betta_i = np.zeros(num_steps + 2)
    alfa_i[0] = 1 #прогоночный коэффициент альфа на нулевом шаге
    betta_i[0] = 0 #прогоночный коэффициент бетта на нулевом шаге
    c_temp = []
    c_temp = np.copy(c)

    if key == 'one_dim':
        if key_sch == 'explicit':
            if stab_cond > 1 / (2*diff_coef):
                print('Не выполняется условие устойчивости')
                raise SystemExit
            c_temp[1:-1] = c[1:-1] + diff_coef * (dt / dr**2) * (c[2:] - 2 * c[1:-1] + c[0: -2])
            c_temp[-1] = c_bound
            c_temp[0] = c_temp[1]
            return c_temp
        
        elif key_sch == 'implicit':            
            a = -diff_coef * dt / (dr)**2 #коэффициент 'a'
            b = 1 + 2 * diff_coef * dt / (dr)**2 #коэффицент b
            c_koef = -diff_coef * dt / (dr)**2 #коэффициент c
            
            for i in range(1, len(l)):
                alfa_i[i] = (-a)/(b + c_koef * alfa_i[i-1])
                betta_i[i] = (c[i] - c_koef * betta_i[i-1])/(b + c_koef * alfa_i[i-1])
            alfa_i[-1] = 0
            betta_i[-1] = 0
            c_temp[1:-1] = c[2:] * alfa_i[1:-1] + betta_i[1:-1]
            
            c_temp[-1] = c_bound
            c_temp[0] = c_temp[1]
            return c_temp
        
    elif key == 'cyl':
        if key_sch == 'explicit':
            if 2*diff_coef* stab_cond <= 1: #(2*diff_coef*dt-dt)/ dr**2 <= 1:     #diff_coef*(dt/dr + 2 * stab_cond) > 1:страница 84 методички
                print('Не выполняется условие устойчивости')
                raise SystemExit
            for i in range(len(r)):
                c_temp[1:-1] = c[1:-1] + diff_coef * dt  *( (c[2:] - 2*c[1:-1] + c[0:-2]) / dr**2 +  1/ r[i] * (c[1:-1]-c[0:-2])/dr  )        
            c_temp[-1] = c_bound
            c_temp[0] = c_temp[1]
            return c_temp
        
     
        elif key_sch == 'implicit' : #должна быть абсолютно устойчива это с ЛКР
            for j in range(1, len(r)):
                a = -diff_coef * dt / (dr)**2
                b = 1 + 2 * dt * diff_coef / (dr)**2 - dt * diff_coef / (dr * r[j])
                c_koef = - dt * diff_coef / (dr)**2 + dt * diff_coef / (dr * r[j])
            for i in range(1, len(l)):
                alfa_i[i] = (-a)/(b+c_koef * alfa_i[i-1])
                betta_i[i] = (c[i] - c_koef * betta_i[i-1])/(b + c_koef * alfa_i[i-1])
            alfa_i[-1] = 0
            betta_i[-1] = 0
            c_temp[1:-1] = c[2:]*alfa_i[1:-1] + betta_i[1:-1]
            c_temp[-1] = c_bound
            c_temp[0] = c_temp[1]
            return c_temp    
            
    elif key == 'sphere':
        if key_sch == 'explicit':
            if 2*diff_coef * stab_cond <= 1: #diff_coef*(dt/dr + 2 * stab_cond) > 1:
                print('Не выполняется условие устойчивости')
                raise SystemExit
            for i in range(len(r)):
                c_temp[1:-1] = c[1:-1] + diff_coef * dt  *( (c[2:] - 2* c[1:-1] + c[0:-2] ) / dr**2 + (c[2:] - c[0:-2])/(r[i] * dr)) #явная разностная схема с ЦКР работает
            c_temp[-1] = c_bound
            c_temp[0] = c_temp[1]
            return c_temp
        
        elif key_sch == 'implicit':
            for j in range(1, len(r)):
                a = -diff_coef*dt/(dr)**2 - (1/r[j]) * diff_coef * dt/dr
                b = 1 + 2 * dt * diff_coef/(dr)**2
                c_koef = -dt * diff_coef/(dr)**2 + (1/r[j]) * diff_coef*dt/dr
            
            for i in range(1, len(l)):
                alfa_i[i] = (-a)/(b+c_koef * alfa_i[i-1])
                betta_i[i] = (c[i] - c_koef * betta_i[i-1])/(b + c_koef * alfa_i[i-1])
                
            alfa_i[-1] = 0
            betta_i[-1] = 0
            c_temp[1:-1] = c[2:]*alfa_i[1:-1] + betta_i[1:-1]
            c_temp[-1] = c_bound
            c_temp[0] = c_temp[1]
            return c_temp


def fick_mass(c, length, width, key='sphere'):
    m = 0.
    for i in range(1, len(c)):
        if key == 'sphere':
            m += c[i] * 4 / 3 * np.pi * ((i * dr)**3 - ((i-1) * dr)**3)
        elif key == 'cyl':
            m += c[i] * length * np.pi *((i * dr)**2 - ((i-1) * dr)**2)
        elif key == 'one_dim':
            m += c[i] * ((2 * i * dr) - ((i-1) * 2 * dr)) * length * width
        else:
            print("Key input error!")
            return None
    return m


def way_doing(c_r, c_init):
    # initialization of list of concentrations
    for i in range(num_steps + 2):
        c_r[i] = c_init
        if i == num_steps + 1:
            c_r[i] = 0
            
    mass_list = np.zeros(n_t)
    time = np.linspace(0, proc_time, n_t)
    mass_list[0] = fick_mass(c_r, length, width)
    for i in range(1, n_t):
        c_r = fick_conc(c_r, key, key_sch, c_bound, dr, dt, r)
        mass_list[i] = fick_mass(c_r, length, width)
        plt.figure(1)
        if i % 100 == 0:
            
            plt.plot(r, c_r)
            plt.title(f'Concentration change profile for {key} {key_sch} ', fontsize=18)
            plt.xlabel('Radius, m')
            plt.grid(True)
            plt.ylabel('Concentration of alcohol ')
    
    plt.figure(2)
    plt.plot(time, mass_list)
    plt.title(f'Mass change profile for {key} {key_sch}', fontsize=18)
    plt.xlabel('Time, second')
    plt.grid(True)
    plt.ylabel('Mass of alcohol, kg')


def time_iteration(c_init_list, n_t, dt, dr, key='sphere', key_sch='implicit'):
    c_app = np.zeros(n_t)
    #volume = 0.00025 #кубические метры из диссертации
    #flow_rate = 0.0000017 #кубических метров в секунду из диссертации
    volume = 1.  # cubic meter
    flowrate = 0.0001  # cubic meter per second
    residence_time = volume / flowrate
    c_matrix = np.zeros((n_t, len(c_init_list)))
    mass_list = np.zeros(n_t)
    c_matrix[0] = c_init_list
    mass_list[0] = fick_mass(c_matrix[0], length, width)
    c_app[0] = 0.
    for i in range(1, n_t):
        c_bound = c_app[i - 1]
        c_matrix[i] = fick_conc(c_matrix[i-1], key, key_sch, c_bound, dr, dt, r)
        # TODO добавить количество образцов в расчет массы
        N = 1000
        
        if volume < N * 4 / 3 * np.pi * (R)**3 and key == 'sphere' :
            print('Объём аппарата превышен ')
            raise SystemExit()
            
        elif volume < N * length * np.pi *R**2 and key == 'cyl' :
            print('Объём аппарата превышен')
            raise SystemExit()
            
        elif volume < N * length * width *R and key == 'one_dim' :
            print('Объём аппарата превышен')
            raise SystemExit()
        
        else:   
            mass_list[i] = fick_mass(c_matrix[i], length, width)
            delta_mass = - N * (mass_list[i] - mass_list[i - 1])
            c_app[i] = ideal_mixing(c_app[i - 1], 0, residence_time, volume, delta_mass)


        # TODO добавить расчет идеального смешения с учетом прибыли массы из высушиваемых частиц

    return c_matrix, mass_list, c_app


def plot_conc(r_list, time, c_list):
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


def plot_mass(r_list, mass_list, legend_key):
    plt.figure(2)
    plt.plot(r_list, mass_list, label=legend_key)
    plt.legend()
    plt.title('Mass change profile of alcohol', fontsize=18)
    plt.xlabel('Time, second')
    plt.grid(True)
    plt.ylabel('Mass of alcohol, kg')
    return


def plot_3D(r_list, time, c_list, name):
    fig = plt.figure(figsize=(7, 4))
    xgrid, ygrid = np.meshgrid(r_list, time)
    zgrid = c_list
    ax_3d = Axes3D(fig)
    ax_3d.plot_surface(xgrid, ygrid, zgrid, cmap=cm.jet)

    plt.title(('3D concentration ' + name), fontsize=12)
    plt.xlabel('Radius, m')
    plt.grid(True)
    plt.ylabel('Time, second')
    #plt.zlabel('Concentration of CO2')
    return


def ideal_mixing(c, c_inlet, residence_time, volume, delta_mass):

    c_mixing = c + dt * 1 / residence_time * (c_inlet - c) + dt * delta_mass / volume
    return c_mixing


def ideal_mixing_plot(time, c_mixing):
    plt.figure()
    plt.plot(time, c_mixing)
    plt.title('Concentration change profile of alcohol for ideal mixing', fontsize=18)
    plt.xlabel('Time, second')
    plt.grid(True)
    plt.ylabel('Concentration of alcohol')
    return    


# главный цикл
def main():
    #fick_conc(c_r, key, key_sch, c_bound, dr, dt, r)
    #way_doing(c_r, c_init)
    
    c_init_list = np.zeros(num_steps + 2)
    for i in range(num_steps + 2):
        c_init_list[i] = c_init
        if i == num_steps + 1:
            c_init_list[i] = 0

    #matrix_of_c, list_of_mass = time_iteration(c_init_list, n_t, dt, dr)
    print('n_t:', n_t, 'proc_time:', proc_time)
    time = np.linspace(0, proc_time, n_t)
    #time = np.linspace(0, proc_time, num_steps + 2)
    r_list = np.linspace(0, R, num_steps + 2)

    #прохожусь по списку фигур
    key_list = ['sphere', 'cyl', 'one_dim']

    for i in key_list:
        matrix_of_c, list_of_mass, c_app = time_iteration(c_init_list, n_t, dt, dr, key=i)
        plot_mass(time, list_of_mass, i)
        
        plot_conc(r_list, time-1, matrix_of_c)
        ideal_mixing_plot(time, c_app)
        plot_3D(r_list, time, matrix_of_c, i)
        
    plt.show()

    
main()