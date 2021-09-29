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
    NFIGS = 1
    FIGDIMS = 800, 600
    STEPS = FIGDIMS[0]/100, FIGDIMS[1]/100, 

    BASE = os.environ['PWD']
    HTML_BASE_NAME = os.path.join(BASE, os.path.basename(__file__)[:-3] )

    ##########################################################
    ######## File reading ####################################
    ##########################################################
    l1 = ['data/area_of_intersection.csv']

    df = pd.read_csv( l1[0] )
    g = df.groupby('Theta').mean()
    # print(df.head())
    # print(g.head())
    
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
             'y.axis_label': 'Area [barn]',
             'l.click_policy': 'hide',
             'y_axis_type': 'log'}
             #'y_range': Range1d(0,limit)}

    ### FIGURE 0 ###
    nvals = 1
    colors = ('blue', 'orange', 'green', 'dark orange', 'brown')
    legends = ('bW=6.68fm','bW=2*leadR','bW=3*leadR','bW=4*leadR','bW=5*leadR')
    myvars = tuple('Area' + str(x) for x in range(nvals))
    ylims = [0., 1e10]
    for idx,ivar in enumerate(myvars):
        #total_area = np.trapz(g[ivar])
        b.graph(idx=0, data=[g.index,g[ivar]],
                color=colors[idx], fig_kwargs=figkw,
                legend_label=legends[idx])

        if np.max(g[ivar])<ylims[1]:
            ylims[1] = np.max(g[ivar])
        if np.min(g[ivar])>ylims[0]:
            ylims[0] = np.min(g[ivar])
            

    beamShift = 0.8 # [cm]
    beamWidth = 0.1 # [cm]
    deflection_dist = 5000 # [cm]
    alpha_min = np.arctan( (beamShift - beamWidth) / deflection_dist )
    alpha_max = np.arctan( (beamShift + beamWidth) / deflection_dist )
    theta_min = 2*alpha_min
    theta_max = 2*alpha_max
    b.line(idx=0, x=[theta_min,theta_min], y=[ylims[0], ylims[1]],
                              color='red', width=3.)
    b.line(idx=0, x=[theta_max,theta_max], y=[ylims[0], ylims[1]],
                              color='red', width=3.)


    ### FIGURE 1 ###
    # figkw.pop('l.click_policy')
    # b.histogram(idx=1,
    #             data=np.histogram2d(df.Theta, df.Area0, bins=50),
    #             style='quad%Viridis',
    #             continuum_color="grey",
    #             continuum_value=0.,
    #             fig_kwargs=figkw)

    ##########################################################
    ######## Saving ##########################################
    ##########################################################
    b.save_frame(nrows=1, show=True)
    #b.save_frame(layout=[[0,1]], show=True)
    #b.save_figs(path='.', mode='png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    FLAGS, _ = argparser.add_args(parser, 'sim')
    run(FLAGS)
