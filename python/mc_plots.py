import os
import glob
import pandas as pd
import numpy as np
import bokehplot as bkp
import argparser
import argparse
import re

from latex import LatexLabel
from bokeh.layouts import gridplot
from bokeh.embed import components
from bokeh.models import Range1d
    
def run(FLAGS):
    ##########################################################
    ######## Constants #######################################
    ##########################################################
    NFIGS = 2
    FIGDIMS = 800, 600
    STEPS = FIGDIMS[0]/100, FIGDIMS[1]/100, 

    BASE = os.environ['PWD']
    HTML_BASE_NAME = os.path.join(BASE, os.path.basename(__file__)[:-3] )

    ##########################################################
    ######## File reading ####################################
    ##########################################################
    l1 = ['data/collision_prob.csv']

    df = pd.read_csv( l1[0] )
    g = df.groupby('Theta').mean()
    print(df.head())
    print(g.head())
    
    ##########################################################
    ######## Plotting ########################################
    ##########################################################    
    b = bkp.BokehPlot( HTML_BASE_NAME + '.html',
                       nfigs=NFIGS, nwidgets=0,
                       fig_width=FIGDIMS[0], fig_height=FIGDIMS[1] )

    additional_kwargs = {'text_font_size': '9pt', 'text_font_style': 'italic',
                         'x_units': 'screen', 'y_units': 'screen',
                         'angle': np.pi/2}
    dlabel = dict(x=3*STEPS[0],
                  y=FIGDIMS[1]/2+10*STEPS[1], **additional_kwargs)

    phistr   = '\u03D5'
    thetastr = '\u03B8'

    #limit = 0.04
    
    figkw = {'x.axis_label': thetastr + ' [rad]',
             'y.axis_label': 'Distance being close'}
             #'y_range': Range1d(0,limit)}

    ### FIGURE 0 ###
    b.graph(idx=0, data=[g.index,g.distProxy],
            color='red', fig_kwargs=figkw)

    figkw.update({'x.axis_label': 'Distance being close',
                  'y.axis_label': thetastr + ' [rad]'})
    #figkw.pop('y_range')
    b.histogram(idx=1,
                data=np.histogram2d(df.distProxy, df.Theta, bins=50),
                #style='quad%Viridis',
                fig_kwargs=figkw)


    ##########################################################
    ######## Saving ##########################################
    ##########################################################
    #b.save_frame(nrows=5, ncols=4, show=True)
    b.save_frame(layout=[[0,1]], show=True)
    #b.save_figs(path='.', mode='png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    FLAGS, _ = argparser.add_args(parser, 'sim')
    run(FLAGS)
