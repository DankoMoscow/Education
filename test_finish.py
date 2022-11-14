import panel as pn
import holoviews as hv
from bokeh.plotting import figure

import fick_classes
pn.extension()

float_width = pn.widgets.FloatSlider(name='Width of samples', start=0.1, end=0.5, step=0.05, value=0.1)
float_length = pn.widgets.FloatSlider(name='Length of samples', start=0.1, end=1, step=0.1, value=0.5)
float_diff_coef = pn.widgets.FloatSlider(name='Diffusion coef', start=0.000000001, end=0.00000001, step=0.000000001, value=0.00000001)
int_number_samples = pn.widgets.IntSlider(name='Number of samples', start=0, end=1500, step=100, value=1000)

button = pn.widgets.Button(name='Click me', button_type='primary')

#виджеты для отображения
main_colomn = pn.Column(float_width, float_length, float_diff_coef, int_number_samples, button)
plot = figure(width=300, height=300)

#общий виджет
abf = pn.Row(pn.Spacer(width=200), main_colomn, plot, sizing_mode='stretch_width')

def run(event):
    global abf, float_width, float_length, float_diff_coef, int_number_samples

    matrix_of_c, list_of_mass, c_app, time, i, r_list = fick_classes.main(float_width.value, float_length.value, float_diff_coef.value, int_number_samples.value)
    #plot = hv.Curve(list(list_of_mass), sizing_mode='stretch_width').opts(width=600)
    plot = figure(r_list, time-1, matrix_of_c)
    abf[1] = plot


button.on_click(run)
abf.show()