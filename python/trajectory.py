import sys, os
import pandas as pd
import numpy as np
import bokehplot as bkp
from bokeh.layouts import gridplot

##########################################################
######## Constants #######################################
##########################################################
NFIGS=2
FIGDIMS = 900, 300
STEPS = FIGDIMS[0]/100, FIGDIMS[1]/100, 
MODE=sys.argv[1]
if os.environ.get('CONDA_PREFIX'):
    HTML_NAME = os.path.join(os.environ['CONDA_PREFIX'], os.path.basename(__file__)[:-2] + 'html')
else:
    HTML_NAME = os.path.join(os.environ['HOME'], os.path.basename(__file__)[:-2] + 'html')

##########################################################
######## Data Processing #################################
##########################################################
df = pd.read_csv('data/track_' + MODE + '.csv')

##########################################################
######## Plotting ########################################
##########################################################

b = bkp.BokehPlot( HTML_NAME, nfigs=NFIGS, nwidgets=0,
                   fig_width=FIGDIMS[0], fig_height=FIGDIMS[1])

additional_kwargs = {'text_font_size': '9pt', 'text_font_style': 'italic',
                     'x_units': 'screen', 'y_units': 'screen',
                     'angle': np.pi/2}
dlabel = dict(x=FIGDIMS[0]/2+5*STEPS[0],
              y=FIGDIMS[1]-38*STEPS[1], **additional_kwargs)

figkw = {'x.axis_label': 'Z [cm]',
         'y.fixed_location': 0.,
         't.text': 'Trajectory'}
b.graph(idx=0, data=[df.z,df.x], style='circle',
        color='orange', size=1, line=False,
        fig_kwargs=figkw)
xaxis_label = 'X [cm]'
b.label(xaxis_label, idx=0, **dlabel)

figkw.update(x_range=b.get_figure(idx=0).x_range)
b.graph(idx=1, data=[df.z,df.y], style='circle',
        color='orange', size=1, line=False,
        fig_kwargs=figkw)
xaxis_label = 'Y [cm]'
b.label(xaxis_label, idx=1, **dlabel)


##########################################################
######## Saving ##########################################
##########################################################
b.save_frame(ncols=1, show=True)
#b.save_figs(path='.', mode='png')

