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
from plotly.subplots import make_subplots
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from bokeh.models.formatters import PrintfTickFormatter

pid = os.getpid()
hv.extension('plotly')

float_width = pn.widgets.FloatSlider(name='Width of samples', start=0.01, end=0.1, step=0.01, value=0.01, format=PrintfTickFormatter(format='%.2f м'))
float_length = pn.widgets.FloatSlider(name='Length of samples', start=0.1, end=1, step=0.1, value=0.5, format=PrintfTickFormatter(format='%.2f м'))

float_volume = pn.widgets.FloatSlider(name='Объём аппарата', start=0.1, end=1, step=0.1, value=0.3, format=PrintfTickFormatter(format='%.2f кубических метров'))
float_flowrate = pn.widgets.FloatSlider(name='Объёмный расход', start=0.0001, end=0.01, step=0.001, value=0.001, format=PrintfTickFormatter(format='%.4f кубических метров в секунду'))

float_diff_coef = pn.widgets.FloatSlider(name='Diffusion coef', start=0.000000001, end=0.00000001, step=0.000000001, value=0.00000001, format=PrintfTickFormatter(format='%.9f квадратных метров в секунду'))
int_number_samples = pn.widgets.IntSlider(name='Number of samples', start=100, end=1000, step=100, value=200, format=PrintfTickFormatter(format='%.1f штук'))
text_place = pn.widgets.TextInput(name='Поле для вывода ошибок', value=' ')


group_of_actions = pn.widgets.RadioButtonGroup(
    name='Выбор необходимого графика', options=['График изменения массы', 'График изменения концентрации', 'График 3D'], button_type='success', orientation = 'vertical')

group_of_key = pn.widgets.RadioButtonGroup(
    name='Выбор необходимого вида образца', options=['one_dim', 'cyl', 'sphere'], button_type='success', orientation = 'vertical')

group_of_ways = pn.widgets.RadioButtonGroup(
    name = 'Выбор разностной схемы', options = ['implicit', 'explicit'], button_type = 'success', orientation = 'vertical')

static_text = pn.widgets.StaticText(name='Поле для вывода ошибок', value=' ')

button = pn.widgets.Button(name='Нажмите для запуска расчёта', button_type='primary')
button_exit = pn.widgets.Button(name='Нажмите для выхода', button_type='primary')

#виджеты для отображения
main_column = pn.Column('# Модель расчёта массы и концентрации спирта в образцах', float_width, float_length, float_diff_coef, float_volume, float_flowrate,\
                        int_number_samples, group_of_actions, group_of_key,group_of_ways, static_text, button, button_exit)
plot_list = []
plot = hv.Curve([0]).opts(width=300)

#общий виджет
main_window = pn.Row(pn.Spacer(width=100), main_column, pn.Spacer(width = 200), plot, sizing_mode='stretch_width')


spisok_key = ['one_dim', 'cyl', 'sphere']

def onClose(event):
    # убивает процесс
    text_place.value = 'До новых встреч'
    static_text.value = 'До новых встреч'
    panel.command.die('By')

def pole_name():
    if podskazka == 5:
        static_text.value = 'Ошибка. Объём аппарата превышен'
    else:
        static_text.value = 'Всё работает правильно'


def work_process():
    if podskazka == 5:
        text_place.value = 'Ошибка. Объём аппарата превышен'
    else:
        text_place.value = 'Всё работает правильно'


def run(event):
    global main_window, float_width, float_length, float_diff_coef, int_number_samples, podskazka
    list_of_plots = make_subplots(rows=3, cols=1, shared_yaxes=True, specs=[[{'type': 'surface'}],
                                                                            [{'type': 'surface'}],
                                                                            [{'type': 'surface'}]])
    for index, key in enumerate(spisok_key):
        matrix_of_c, list_of_mass, c_app, time, i, r_list, podskazka = fick_classes.main(float_width.value, float_length.value, float_volume.value, float_flowrate.value, float_diff_coef.value, int_number_samples.value, key, group_of_ways.value, text_place.value)
        if group_of_actions.value == 'График изменения массы':

            plot = hv.Curve(list(list_of_mass), ('x', 'r_list'), ('y', 'mass')).opts(title="График изменения массы", height=600,  width=600).\
                redim(x=hv.Dimension('Количество шагов', range=(0, 2001))).\
                redim(y=hv.Dimension('Масса спирта в образцах, кг', range=(list_of_mass[-1], list_of_mass[0])))
            pole_name()

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
            plot = hv.Overlay(list_of_curves).opts(height=600, width=600, title="График изменения концентрации", xlabel='Радиус, м',
                ylabel='Концентрация, кг/кубический метр', legend_position='right')

            pole_name()

        elif group_of_actions.value == 'График 3D':
            x = r_list
            y = time
            z = matrix_of_c

            # plot = go.Figure(data=[go.Surface(x = x, y = y, z=z)])
            # plot.update_layout(title="График отображения 3D концентрации" + "_" + key,
            #                    scene=dict(xaxis_title='Радиус, м',
            #                        yaxis_title='Время, с',
            #                        zaxis_title='Концентрация спирта, кг/кубический метр'),
            #                    width=500, height =500, margin=dict(l=65, r=50, b=65, t=90))

            list_of_plots.add_trace(go.Surface(x = x, y = y, z=z), row=index + 1, col=1)
            list_of_plots.update_layout(title_text='3D subplots with different colorscales',
                                        height=800, width=800)
            pole_name()

        #plot
    main_window[1] = main_column
    #main_window[3] = plot
    main_window[3] = list_of_plots

def main():
    button.on_click(run)
    button_exit.on_click(onClose)
    main_window.show()
main()