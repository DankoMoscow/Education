import signal
import panel as pn
import holoviews as hv
import panel.command
from bokeh.plotting import figure
from holoviews import opts
import numpy as np
import fick_classes
import os, sys
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from selenium import webdriver

#driver = webdriver.Chrome()

pid = os.getpid()
pn.extension()

float_width = pn.widgets.FloatSlider(name='Width of samples', start=0.1, end=0.5, step=0.05, value=0.1)
float_length = pn.widgets.FloatSlider(name='Length of samples', start=0.1, end=1, step=0.1, value=0.5)
float_diff_coef = pn.widgets.FloatSlider(name='Diffusion coef', start=0.000000001, end=0.00000001, step=0.000000001, value=0.00000001)
int_number_samples = pn.widgets.IntSlider(name='Number of samples', start=0, end=1500, step=100, value=1000)

group_of_actions = pn.widgets.RadioButtonGroup(
    name='Выбор необходимого графика', options=['График изменения массы', 'График изменения концентрации', 'График 3D'], button_type='success', orientation = 'vertical')

group_of_key = pn.widgets.RadioButtonGroup(
    name='Выбор необходимого вида образца', options=['one_dim', 'cyl', 'sphere'], button_type='success', orientation = 'vertical')

button = pn.widgets.Button(name='Нажмите для запуска расчёта', button_type='primary')
button_exit = pn.widgets.Button(name='Нажмите для выхода', button_type='primary')

#виджеты для отображения
main_column = pn.Column(float_width, float_length, float_diff_coef, int_number_samples, group_of_actions, group_of_key, button, button_exit)

iter = 3
plots = []
abfs = []

plot = hv.Curve([0]).opts(width=600)
#общий виджет
abf = pn.Row(pn.Spacer(width=100), main_column, pn.Spacer(width = 200), plot, sizing_mode='stretch_width')

def onClose(event):
    print(pid)
    """
    driver.close()
    driver.get()"""

    #os.system("taskkill /F /IM Crome.url")
    # убивает процесс
    panel.command.die('By')

    #os.kill(pid, signal.SIGINT)

def run(event):
    global abf, float_width, float_length, float_diff_coef, int_number_samples

    matrix_of_c, list_of_mass, c_app, time, i, r_list = fick_classes.main(float_width.value, float_length.value, float_diff_coef.value, int_number_samples.value, group_of_key.value)
    #plot = hv.Curve(list(matrix_of_c), sizing_mode='stretch_width').opts(width=600)

    if group_of_actions.value == 'График изменения массы':
        plot = hv.Curve(list(list_of_mass)).opts( title="График изменения массы",  width=300)

    elif group_of_actions.value == 'График изменения концентрации':
        plot = hv.Curve(list(matrix_of_c)).opts( title = 'График изменения концентрации',      width=300)

    elif group_of_actions.value == 'График 3D':
        X = r_list
        Y = time
        X, Y = np.meshgrid(X, Y)
        Z = matrix_of_c

        #r = [X, Y, Z]
        #plot = hv.Surface(r, bounds=(-5, -5, 5, 5)).opts(colorbar=True, width=500, height=500)

    abf[0] = main_column
    abf[1] = plot
    abf


button.on_click(run)

button_exit.on_click(onClose)
abf.show()