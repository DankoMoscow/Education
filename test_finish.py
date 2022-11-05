import panel as pn

import fick_classes
pn.extension()

float_width = pn.widgets.FloatSlider(name='Width of samples', start=0.1, end=0.5, step=0.05, value=0.1)
float_length = pn.widgets.FloatSlider(name='Length of samples', start=0.1, end=1, step=0.1, value=0.5)
float_diff_coef = pn.widgets.FloatSlider(name='Diffusion coef', start=0.000000001, end=0.00000001, step=0.000000001, value=0.00000001)
int_number_samples = pn.widgets.IntSlider(name='Number of samples', start=0, end=1500, step=100, value=1000)

button = pn.widgets.Button(name='Click me', button_type='primary')

abf = pn.Column(float_width, float_length, float_diff_coef, int_number_samples, button)
abf.show()
if True:
    button.param.watch(fick_classes.main(float_width.value, float_length.value, float_diff_coef.value, int_number_samples.value), 'clicks')