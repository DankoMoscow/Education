import signal
import panel as pn
import holoviews as hv
import panel.command
import pandas as pd
from bokeh.plotting import figure
from holoviews import opts
import numpy as np
import fick_classes
import os, sys
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

pid = os.getpid()
pn.extension()

float_width = pn.widgets.FloatSlider(name='Width of samples', start=0.1, end=0.5, step=0.05, value=0.1)
float_length = pn.widgets.FloatSlider(name='Length of samples', start=0.1, end=1, step=0.1, value=0.5)
float_diff_coef = pn.widgets.FloatSlider(name='Diffusion coef', start=0.000000001, end=0.00000001, step=0.000000001, value=0.00000001)
int_number_samples = pn.widgets.IntSlider(name='Number of samples', start=0, end=1500, step=100, value=1000)
text_place = pn.widgets.TextInput(name='Поле для вывода ошибок', value=' ')

group_of_actions = pn.widgets.RadioButtonGroup(
    name='Выбор необходимого графика', options=['График изменения массы', 'График изменения концентрации', 'График 3D'], button_type='success', orientation = 'vertical')

group_of_key = pn.widgets.RadioButtonGroup(
    name='Выбор необходимого вида образца', options=['one_dim', 'cyl', 'sphere'], button_type='success', orientation = 'vertical')

group_of_ways = pn.widgets.RadioButtonGroup(
    name = 'Выбор разностной схемы', options = ['implicit', 'explicit'], button_type = 'success', orientation = 'vertical')

button = pn.widgets.Button(name='Нажмите для запуска расчёта', button_type='primary')
button_exit = pn.widgets.Button(name='Нажмите для выхода', button_type='primary')

#виджеты для отображения
main_column = pn.Column(float_width, float_length, float_diff_coef, int_number_samples, group_of_actions, group_of_key,group_of_ways, text_place,  button, button_exit)
plot = hv.Curve([0]).opts(width=300)

#общий виджет
abf = pn.Row(pn.Spacer(width=100), main_column, pn.Spacer(width = 200), plot, sizing_mode='stretch_width')

def onClose(event):
    # убивает процесс
    text_place.value = 'До новых встреч'
    panel.command.die('By')


def work_process():
    if podskazka == 5:
        text_place.value = 'Ошибка. Объём аппарата превышен'
    else:
        text_place.value = 'Всё работает правильно'


def run(event):
    global abf, float_width, float_length, float_diff_coef, int_number_samples, podskazka
    matrix_of_c, list_of_mass, c_app, time, i, r_list, podskazka = fick_classes.main(float_width.value, float_length.value, float_diff_coef.value, int_number_samples.value, group_of_key.value, group_of_ways.value, text_place.value)

    if group_of_actions.value == 'График изменения массы':
        plot = hv.Curve(list(list_of_mass),('x', 'r_list'), ('y', 'mass')).opts( title="График изменения массы", width=300).\
            redim(x = hv.Dimension('radius', range=(0, 2001))).\
            redim(y =hv.Dimension('334', range=(list_of_mass[-1], list_of_mass[0])))
        work_process()

    elif group_of_actions.value == 'График изменения концентрации':
        cols = []
        list_of_curves = []
        list_of_c = []
        time_ratio = 100
        matrix_of_c = matrix_of_c.T
        for l in matrix_of_c[:, ::time_ratio]:
            list_of_c.append(l)
        for col in range(len(list_of_c[0])):
            cols.append(str(col))
        df_x = pd.DataFrame(list_of_c, columns=cols)
        df_y = pd.DataFrame(r_list, columns=['r_list'])
        df = pd.concat([df_x, df_y], axis=1)
        for col in df_x.columns:
            list_of_curves.append(hv.Curve(df[['r_list', col]]))
        plot = hv.Overlay(list_of_curves).opts(height=300, width=600, title="График изменения концентрации", xlabel='r_list',
            ylabel='concentration', legend_position='right')
        work_process()

    elif group_of_actions.value == 'График 3D':
        X = r_list
        Y = time
        X, Y = np.meshgrid(X, Y)
        Z = matrix_of_c

        #r = [X, Y, Z]
        #plot = hv.Surface(r, bounds=(-5, -5, 5, 5)).opts(colorbar=True, width=500, height=500)
        work_process()

    abf[0] = main_column
    abf[1] = plot
    abf


button.on_click(run)
button_exit.on_click(onClose)
abf.show()