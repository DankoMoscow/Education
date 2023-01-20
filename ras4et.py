import numpy as np
from raschet import *
density_ips = 785.1  # кг/м3
density_co2 = 468  # кг/м3

density = []
density_new = []
y = np.linspace(1, 0, 102)
x_ips = []
x_ips_list = []
print(y)
for i in range(len(y)):
    density = (y[i]/density_ips + (1-y[i])/density_co2)**-1
    density_new.append(density)

    x_ips = (y[i] * M_ips) / ((1 - y[i]) * M_co2 + y[i] * M_ips)
    x_ips_list.append(x_ips)

diff_coef_list = []
for k in range(len(x_ips_list)):
    diff_coef = (diff_coef_co2 ** x_ips_list[k] * diff_coef_ips ** (1 - x_ips_list[k]))
    diff_coef_list.append(diff_coef)
    print(x_ips_list[k])

print(diff_coef_list)
