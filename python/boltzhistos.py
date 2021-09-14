import sys, os
import glob
import pandas as pd
import numpy as np
import bokehplot as bkp
from bokeh.models import Range1d
import argparser
import argparse

from latex import LatexLabel
from bokeh.layouts import gridplot
    
def run(FLAGS):
    ##########################################################
    ######## Constants #######################################
    ##########################################################
    NFIGS=1
    FIGDIMS = 600, 600
    STEPS = FIGDIMS[0]/100, FIGDIMS[1]/100, 

    BASE = os.environ['PWD']
    HTML_NAME = os.path.join(BASE, os.path.basename(__file__)[:-3] + '_' + \
                             str(FLAGS.x) + 'X_' + str(FLAGS.x) + 'Y_' + \
                             str(FLAGS.energy) + 'En.html')

    ##########################################################
    ######## File reading ####################################
    ##########################################################
    search_str = os.path.join(BASE, 'data/' + FLAGS.dataset + '.csv')

    l = glob.glob(search_str)
    if len(l) != 1:
        print('Pattern: ', search_str)
        print('NHits: ', len(l))
        raise RuntimeError('Only one data file can match the pattern. Check the input arguments and the existing data files.')
    else:
        print('Data: ', l[0])
    
    df = pd.read_csv( l[0] )
    ##########################################################
    ######## Plotting ########################################
    ##########################################################

    b = bkp.BokehPlot( HTML_NAME, nfigs=NFIGS, nwidgets=0,
                       fig_width=FIGDIMS[0], fig_height=FIGDIMS[1])

    additional_kwargs = {'text_font_size': '9pt', 'text_font_style': 'italic',
                         'x_units': 'screen', 'y_units': 'screen',
                         'angle': np.pi/2}
    dlabel = dict(x=FIGDIMS[0]/2+3*STEPS[0],
                  y=FIGDIMS[1]-32*STEPS[1], **additional_kwargs)

    figkw = {'y.axis_label': 'Counts',
             'x_range': Range1d(0,1.5)}
    b.histogram(idx=0, data=np.histogram(df.data, bins=5000),
                color='blue', fig_kwargs=figkw)
    
    
    ##########################################################
    ######## Saving ##########################################
    ##########################################################
    b.save_frame(show=True)
    #b.save_figs(path='.', mode='png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    FLAGS, _ = argparser.add_args(parser, 'gen')
    run(FLAGS)
