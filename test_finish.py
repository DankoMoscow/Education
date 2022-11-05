import panel as pn
from panel.interact import interact
import fick_classes
from fick_classes import scd_apparatus, main
pn.extension()

float_width = pn.widgets.FloatSlider(name = 'Width of samples', start=0.1, end=0.5, step=0.05, value=0.1)
float_length = pn.widgets.FloatSlider(name = 'Length of samples', start=0.1, end=1, step=0.1, value=0.5)
float_diff_coef = pn.widgets.FloatSlider(name='Diffusion coef', start = 0.000000001, end = 0.00000001, step=0.000000001, value=0.00000001)
int_number_samples = pn.widgets.IntSlider(name='Number of samples', start=0, end=1500, step=100, value=1000)

button = pn.widgets.Button(name='Click me', button_type='primary')
text = pn.widgets.TextInput(value='Ready')

def b(event):
    text.value = 'Clicked {0} times'.format(button.clicks)
button.on_click(b)

abf = pn.Column(float_width,float_length, float_diff_coef, int_number_samples, button,text)
abf.show()
if True:
    pn.interact(main, width=float_width, length =float_length, diff_coef = float_diff_coef, number_samples = int_number_samples)
    button.clicks