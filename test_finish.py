import panel as pn
import holoviews as hv
from bokeh.plotting import figure

import fick_classes
pn.extension()

float_width = pn.widgets.FloatSlider(name='Width of samples', start=0.1, end=0.5, step=0.05, value=0.1)
float_length = pn.widgets.FloatSlider(name='Length of samples', start=0.1, end=1, step=0.1, value=0.5)
float_diff_coef = pn.widgets.FloatSlider(name='Diffusion coef', start=0.000000001, end=0.00000001, step=0.000000001, value=0.00000001)
int_number_samples = pn.widgets.IntSlider(name='Number of samples', start=0, end=1500, step=100, value=1000)

group_of_actions = pn.widgets.RadioButtonGroup(
    name='Выбор необходимого графика', options=['График изменения массы', 'График изменения концентрации', 'График 3D'], button_type='success', orientation = 'vertical')

group_of_key = pn.widgets.RadioButtonGroup(
    name='Выбор необходимого вида образца', options=['Плоскопараллельное тело', 'Цилиндр', 'Сфера'], button_type='success', orientation = 'vertical')


button = pn.widgets.Button(name='Нажмите для запуска расчёта', button_type='primary')

#виджеты для отображения
main_column = pn.Column(float_width, float_length, float_diff_coef, int_number_samples, group_of_actions, group_of_key, button)

iter = 3
plots = []
abfs = []
for j in range(iter):
    plot = figure(width=300, height=300)
    plots.append(plot)

    #общий виджет
    abf = pn.Row(pn.Spacer(width=100), main_column, pn.Spacer(width = 200), plot, sizing_mode='stretch_width')
    abfs.append(abf)



def run(event):
    global abf, float_width, float_length, float_diff_coef, int_number_samples

    matrix_of_c, list_of_mass, c_app, time, i, r_list = fick_classes.main(float_width.value, float_length.value, float_diff_coef.value, int_number_samples.value)
    #plot = hv.Curve(list(matrix_of_c), sizing_mode='stretch_width').opts(width=600)

    plot.line(r_list, list_of_mass)
    abf


button.on_click(run)
abf.show()