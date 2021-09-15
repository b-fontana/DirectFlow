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
    NFIGS = 2*len(l)
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

    dlatex = dict(x=3*STEPS[0],
                  y=FIGDIMS[1]/2+42*STEPS[1],
                  x_units="screen",
                  y_units="screen",
                  text_font_size='10pt')

    phistr = '\u03D5'
    etastr = '\u03B7'

    def add_latex(idx):
        b.get_figure(idx).add_layout(
            LatexLabel(
                text="W_{{scale}} = {}".format(running_var[idx],5),
                render_mode="css",
                background_fill_alpha=0,
                **dlatex
            ) )

    for idx,f in enumerate(l):
        df = pd.read_csv( f )

        add_latex(idx)

        #figkw.update({'x_range': Range1d(-1,1)})
        b.histogram(idx=idx,
                    data=np.histogram2d(df.XHit, df.YHit, bins=100),
                    #style='quad%Viridis',
                    fig_kwargs=figkw)

        b.histogram(idx=idx+len(l),
                    data=np.histogram2d(df.XHit, df.YHit, bins=100, range=[[-2,2],[-2,2]]),
                    #style='quad%Viridis',
                    fig_kwargs=figkw)


    ##########################################################
    ######## Saving ##########################################
    ##########################################################
    b.save_frame(nrows=2, ncols=4, show=True)
    #b.save_frame(iframe=0, layout=[[5,3,4],[0,1,2]], show=True)
    #b.save_figs(path='.', mode='png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    FLAGS, _ = argparser.add_args(parser, 'sim')
    run(FLAGS)
