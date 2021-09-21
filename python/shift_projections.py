import sys, os
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
    l1 = ['data/histo_euler_0p000000X_0p800000Y_1380p00000En_1p000EnScale_1p000WScale_0p000YShift_10p000000SZ_0p938000MInt_1NP.csv',
         'data/histo_euler_0p000000X_0p800000Y_1380p00000En_1p000EnScale_0p000WScale_-0p10YShift_10p000000SZ_0p938000MInt_1NP.csv',
         'data/histo_euler_0p000000X_0p800000Y_1380p00000En_1p000EnScale_0p000WScale_0p000YShift_10p000000SZ_0p938000MInt_1NP.csv',
         'data/histo_euler_0p000000X_0p800000Y_1380p00000En_1p000EnScale_0p000WScale_0p100YShift_10p000000SZ_0p938000MInt_1NP.csv',
         ]

    l2 = []
    for item in FLAGS.fermi_shifts:
        l2.append('data/histo_euler_0p000000X_0p800000Y_1380p00000En_1p000EnScale_1p000WScale_0p000YShift_' + '{:<05}'.format(str(item)).replace('.', 'p') + 'FShift_10p000000SZ_0p938000MInt_1NP.csv')

    running_var = []
    for f in l2:
        running_var.append( float(re.findall(r'.+_(.*)FShift.+.csv', f)[0].replace('p', '.')) )

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

    figkw = {'x.axis_label': 'X [cm]',
             'y.axis_label': 'PDF',
             'l.click_policy': 'hide'}

    dlatex = dict(x=FIGDIMS[0]/2-1*STEPS[0],
                  y=FIGDIMS[1]/2+36*STEPS[1],
                  x_units="screen",
                  y_units="screen",
                  text_font_size='10pt')

    phistr = '\u03D5'
    etastr = '\u03B7'
    
    ### FIGURE 0 ###
    limit = 2.5
    colors = ['blue', 'red', 'orange', 'pink']
    leglabels = ['all', 'Y=0.7', 'Y=0.8', 'Y=0.9']
    for idx,f in enumerate(l1):
        df = pd.read_csv( f )

        b.histogram(idx=0, data=np.histogram(df.XHit, bins=100,
                                             density=True, range=[-limit, limit]),
                    legend_label=leglabels[idx],
                    color=colors[idx], style='step', fig_kwargs=figkw)

    ### FIGURE 1 ###
    for idx,f in enumerate(l2):
        df = pd.read_csv( f )
        b.histogram(idx=1, data=np.histogram(df.XHit, bins=100,
                                             density=True, range=[-limit, limit]),
                    legend_label='P={}GeV'.format(running_var[idx],5),
                    color=colors[idx], style='step', fig_kwargs=figkw)


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
