import panel as pn
import holoviews as hv
import panel.command
import pandas as pd
from bokeh.plotting import figure
from holoviews import opts
import numpy as np
import fick_classes
import os, sys
import plotly.graph_objs as go
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from bokeh.models.formatters import PrintfTickFormatter
import time
import plotly.express as px
from datetime import datetime
# from holoviews import dim, opts
# from holoviews.plotting.bokeh.styles import (line_properties, fill_properties, text_properties)

pid = os.getpid()
hv.extension('plotly')

float_width = pn.widgets.FloatSlider(name='Ширина образцов', start=0.01, end=0.1, step=0.01, value=0.01,
                                     format=PrintfTickFormatter(format='%.2f м'))
float_length = pn.widgets.FloatSlider(name='Длина образцов', start=0.1, end=1, step=0.1, value=0.5,
                                      format=PrintfTickFormatter(format='%.2f м'))

float_height = pn.widgets.FloatSlider(name='Высота образца', start=0.01, end=0.2, step=0.01, value=0.05,
                                      format=PrintfTickFormatter(format='%.2f м'))

float_volume = pn.widgets.FloatSlider(name='Объём аппарата', start=0.1, end=1, step=0.1, value=0.3,
                                      format=PrintfTickFormatter(format='%.2f кубических метров'))
float_flowrate = pn.widgets.FloatSlider(name='Объёмный расход', start=0.0001, end=0.01, step=0.001, value=0.001,
                                        format=PrintfTickFormatter(format='%.4f кубических метров в секунду'))

float_dt = pn.widgets.FloatSlider(name='Шаг по времени', start=10, end=200, step=10, value=50)

float_diff_coef = pn.widgets.FloatSlider(name='Коэффициент диффузии', start=0.000000001, end=0.00000001, step=0.000000001,
                         value=0.00000001, format=PrintfTickFormatter(format='%.0e м2/сек'))

int_number_samples = pn.widgets.IntSlider(name='Количество образцов', start=100, end=1000, step=100, value=200,
                                          format=PrintfTickFormatter(format='%.1f штук'))

group_of_key = pn.widgets.RadioButtonGroup(
    name='Выбор необходимого вида образца', options=['one_dim', 'cyl', 'sphere'], button_type='success',
    orientation='vertical')

group_of_ways = pn.widgets.RadioButtonGroup(
    name='Выбор разностной схемы', options=['implicit', 'explicit'], button_type='success', orientation='vertical')

groups = pn.Row(group_of_key, group_of_ways)

static_cond = pn.widgets.StaticText(name='Условие устойчивости', value=' ')
static_text = pn.widgets.StaticText(name='Поле для вывода ошибок', value=' ')
static_time = pn.widgets.StaticText(name='Время расчёта', value=' ')
static_time_process = pn.widgets.StaticText(name='Расчётное время проведения сушки в часах', value=' ')

button = pn.widgets.Button(name='Нажмите для запуска расчёта', button_type='primary')
button_exit = pn.widgets.Button(name='Нажмите для выхода', button_type='primary')

# виджеты для отображения
main_column = pn.Column('# Расчёт процесса сверхкритической сушки', float_width, float_length, float_height,
                        float_diff_coef, float_volume, float_flowrate, float_dt, \
                        int_number_samples, groups, static_text, static_cond, static_time,
                        static_time_process, button, button_exit)

plot = hv.Curve([0]).opts(width=300)
# общий виджет
main_window = pn.Row(pn.Spacer(width=100), main_column, pn.Spacer(width=50), plot, pn.Spacer(width=50), pn.Spacer(width=50), sizing_mode='stretch_width')


def onClose(event):
    # убивает процесс
    static_text.value = 'До новых встреч'
    panel.command.die('By')


def get_condition():
    if cond_scheme == 5:
        static_cond.value = 'Условие устойчивости не выполняется '
    else:
        static_cond.value = 'Условие выполянется'


def get_time_drying(n):
    static_time_process.value = n * float_dt.value / 600


def work_process():
    if podskazka == 5:
        static_text.value = 'Ошибка. Объём аппарата превышен'
    else:
        static_text.value = 'Всё работает правильно'


def get_time():
    static_time.value = delta_time


def show_time_process(list):
    condition_stop = list[0] * 0.05
    for n, value in enumerate(list):
        if value <= condition_stop:
            break
    return n

def run(event):
    start_time = datetime.now()
    global main_window, float_width,float_dt, float_length, float_diff_coef, int_number_samples, podskazka, delta_time, cond_scheme
    matrix_of_c, list_of_mass, c_app, time, i, r_list, podskazka, cond_scheme = fick_classes.main(float_width.value, float_length.value, float_height.value,
        float_volume.value, float_flowrate.value, float_dt.value, float_diff_coef.value, int_number_samples.value, group_of_key.value, group_of_ways.value, static_text.value,static_cond.value)

    n = show_time_process(list_of_mass)
    get_time_drying(n)

    plot_mass = hv.Curve(list(list_of_mass), ('x', 'r_list'), ('y', 'mass')).opts( title="График изменения массы",   height=500,
                                                                                   bgcolor='black'). \
        redim(x=hv.Dimension('Количество шагов', range=(0, fick_classes.n_t))). \
        redim(y=hv.Dimension('Масса спирта в образцах, кг', range=(list_of_mass[-1], list_of_mass[0])))

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

    plot_conc = hv.Overlay(list_of_curves).opts(height=500, width=500, margin=dict(l=100, r=50, b=65, t=90),
      title="График изменения концентрации",  bgcolor= 'black', xlabel='Радиус, м', ylabel='Концентрация, кг/м3')

    template = "plotly_dark"
    plot_3d = go.Figure(data=[go.Surface(x=r_list, y=time, z=matrix_of_c)])
    plot_3d.update_layout(title="График отображения 3D концентрации",  template=template,  scene=dict(xaxis_title='Радиус, м',yaxis_title='Время, с', zaxis_title='Концентрация спирта, кг/метр3'),
                          width=500, height=500, margin=dict(l=100, r=50, b=65, t=90),  font=dict(family="Franklin Gothic", size=10))

    get_condition()
    work_process()
    main_window[0] = main_column
    main_window[1] = plot_mass
    main_window[2] = pn.Spacer(width = 10)
    main_window[3] = plot_conc
    main_window[4] = pn.Spacer(width = 10)
    main_window[5] = plot_3d

    end_time = datetime.now()
    delta_time = end_time - start_time
    get_time()

def main():
    button.on_click(run)
    button_exit.on_click(onClose)
    main_window.show()


main()