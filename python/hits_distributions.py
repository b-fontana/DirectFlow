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
    NFRAMES=1
    FIGDIMS = 400, 400
    STEPS = FIGDIMS[0]/100, FIGDIMS[1]/100, 

    BASE = os.environ['PWD']
    HTML_BASE_NAME = os.path.join(BASE, os.path.basename(__file__)[:-3] + '_' + \
                                  str(FLAGS.x) + 'X_' + str(FLAGS.x) + 'Y_' + \
                                  str(FLAGS.energy) + 'En_' + \
                                  str(FLAGS.step_size) + 'SZ_' + \
                                  str(FLAGS.npartons) + 'NP')

    ##########################################################
    ######## File reading ####################################
    ##########################################################
    posstr = str(FLAGS.x).replace('.', 'p') + "*X*" + \
        str(FLAGS.y).replace('.', 'p') + "*Y*" + \
        str(FLAGS.energy).replace('.', 'p') + "*En*" + \
        str(FLAGS.step_size).replace('.', 'p') + "*SZ*" + \
        str(FLAGS.npartons)+"NP"
    search_str = os.path.join(BASE, 'data/histo_' + FLAGS.mode + '_' + posstr + '*.csv')

    l = glob.glob(search_str)
    NFIGS = 4*len(l) + 2
    print('Pattern: ', search_str)
    print('Number of files: ', NFIGS)
    running_var = []
    for f in l:
        running_var.append( float(re.findall(r'.+_(.*)WScale.+.csv', f)[0].replace('p', '.')) )

    l = [ x for _,x in sorted(zip(running_var,l), key=lambda pair: pair[0]) ] #overwrite
    running_var = sorted(running_var)
    for f in l:
        print('-- ', f)

    ##########################################################
    ######## Plotting ########################################
    ##########################################################    
    b = bkp.BokehPlot( [HTML_BASE_NAME + 'Cat' + str(x) + '.html' for x in range(NFRAMES)],
                       nframes=NFRAMES, nfigs=NFIGS, nwidgets=0,
                       fig_width=FIGDIMS[0], fig_height=FIGDIMS[1] )

    additional_kwargs = {'text_font_size': '9pt', 'text_font_style': 'italic',
                         'x_units': 'screen', 'y_units': 'screen',
                         'angle': np.pi/2}
    dlabel = dict(x=3*STEPS[0],
                  y=FIGDIMS[1]/2+10*STEPS[1], **additional_kwargs)

    figkw = {'x.axis_label': 'X [cm]',
             'y.axis_label': 'Y [cm]'}

    dlatex = dict(x=FIGDIMS[0]/2-1*STEPS[0],
                  y=FIGDIMS[1]/2+36*STEPS[1],
                  x_units="screen",
                  y_units="screen",
                  text_font_size='10pt')

    phistr = '\u03D5'
    etastr = '\u03B7'

    def add_latex(idx):
        b.get_figure(idx).add_layout(
            LatexLabel(
                text="W_{{scale}} = {}".format(running_var[(idx-ashift)%len(l)],5),
                render_mode="css",
                background_fill_alpha=0,
                **dlatex
            ) )

    ashift = 2
    for idx,f in enumerate(l):
        df = pd.read_csv( f )

        add_latex(ashift+idx)
        add_latex(ashift+idx+len(l))
        add_latex(ashift+idx+2*len(l))
        add_latex(ashift+idx+3*len(l))
        
        limit = 4
        climbefore, climafter = 2500, 500

        if 'y_range' in figkw:
            figkw.pop('y_range')
        figkw.update({'title': 'Before the boost'})
        b.histogram(idx=ashift+idx,
                    data=np.histogram2d(df.XHitNoBoost, df.YHitNoBoost, bins=100,
                                        range=[[-limit,limit],[-limit,limit]]),
                    #style='quad%Viridis',
                    fig_kwargs=figkw)

        figkw.update({'title': 'After the boost'})
        b.histogram(idx=ashift+idx+len(l),
                    data=np.histogram2d(df.XHit, df.YHit, bins=100,
                                        range=[[-limit,limit],[-limit,limit]]),
                    #style='quad%Viridis',
                    fig_kwargs=figkw)
        
        figkw.update({'title': 'Before the boost',
                      'y.axis_label': 'Counts',
                      'y_range': Range1d(0,climbefore)})
        b.histogram(idx=ashift+idx+2*len(l),
                    data=np.histogram(df.XHitNoBoost, bins=100,
                                      range=[-limit,limit]),
                    color='red',
                    fig_kwargs=figkw)

        figkw.update({'title': 'After the boost',
                      'y_range': Range1d(0,climafter)})
        b.histogram(idx=ashift+idx+3*len(l),
                    data=np.histogram(df.XHit, bins=100,
                                      range=[-limit,limit]),
                    color='red',
                    fig_kwargs=figkw)

        
    figkw.update({'x.axis_label': 'Pz [GeV]',
                  'y.axis_label': 'Counts',
                  'title': 'Before the boost'})
    if 'y_range' in figkw:
        figkw.pop('y_range')
    b.histogram(idx=0,
                data=np.histogram(df.FermiPzBeforeBoost, bins=100),
                fig_kwargs=figkw)
    figkw.update({'title': 'After the boost'})
    b.histogram(idx=1,
                data=np.histogram(df.FermiPzAfterBoost, bins=100),
                fig_kwargs=figkw)

    ##########################################################
    ######## Saving ##########################################
    ##########################################################
    #b.save_frame(nrows=5, ncols=4, show=True)
    b.save_frame(iframe=0, layout=[[0,1],[2,3,4,5],[10,11,12,13],[6,7,8,9],[14,15,16,17]],
                 show=True)
    #b.save_figs(path='.', mode='png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    FLAGS, _ = argparser.add_args(parser, 'sim')
    run(FLAGS)
