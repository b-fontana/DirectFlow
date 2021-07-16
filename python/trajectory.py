import sys, os
import pandas as pd
import numpy as np
import bokehplot as bkp

##########################################################
######## Constants #######################################
##########################################################
NFIGS=1
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
xaxis_label = 'X [cm]'
additional_kwargs = {'text_font_size': '9pt', 'text_font_style': 'italic',
                     'x_units': 'screen', 'y_units': 'screen',
                     'angle': np.pi/2}
b.label(xaxis_label, idx=0, x=FIGDIMS[0]/2+3.5*STEPS[0], y=FIGDIMS[1]-28*STEPS[1], **additional_kwargs)

figkw = {'y.axis_label': 'X [cm]',
         'x.axis_label': 'Z [cm]',
         'y.fixed_location': 0.,
         't.text': 'Trajectory'}
b.graph(idx=0, data=[df.z,df.x], style='circle',
        color='orange', size=1, line=False,
        fig_kwargs=figkw)


##########################################################
######## Saving ##########################################
##########################################################
b.save_frame(show=False)
#b.save_figs(path='.', mode='png')

