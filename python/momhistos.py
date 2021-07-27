import sys, os
import glob
import pandas as pd
import numpy as np
import bokehplot as bkp
import argparser
import argparse

from latex import LatexLabel
from bokeh.layouts import gridplot
    
def run(FLAGS):
    ##########################################################
    ######## Constants #######################################
    ##########################################################
    NFIGS=3
    FIGDIMS = 600, 600
    STEPS = FIGDIMS[0]/100, FIGDIMS[1]/100, 

    BASE = os.environ['PWD']
    HTML_NAME = os.path.join(BASE, os.path.basename(__file__)[:-3] + '_' + \
                             str(FLAGS.x) + 'X_' + str(FLAGS.x) + 'Y_' + \
                             str(FLAGS.energy) + 'En.html')

    ##########################################################
    ######## File reading ####################################
    ##########################################################
    posstr = str(FLAGS.x).replace('.', 'p') + "*" + \
        str(FLAGS.y).replace('.', 'p') + "*" + \
        str(FLAGS.energy).replace('.', 'p')
    search_str = os.path.join(BASE, 'data/histo_' + FLAGS.mode + '_' + posstr + '*.csv')

    l = glob.glob(search_str)
    if len(l) != 1:
        print('Pattern: ', search_str)
        print('NHits: ', len(l))
        raise RuntimeError('Only one data file can match the pattern. Check the input arguments and the existing data files.')
    else:
        print('Data: ', l[0])
    
    df = pd.read_csv( l[0] )
    df['sumMomXAbs'] = np.abs(df.sumMomX)
    df['sumMomYAbs'] = np.abs(df.sumMomY)
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

    figkw = {'y.axis_label': 'Counts'}
    figkw.update({'x.axis_label': 'X momentum sum [GeV]'})
    b.histogram(idx=0, data=np.histogram(df.sumMomXAbs, bins=100), style='v%0.8%red',
            color='orange', fig_kwargs=figkw)
    figkw.update({'x.axis_label': 'Y momentum sum [GeV]'})
    b.histogram(idx=1, data=np.histogram(df.sumMomYAbs, bins=100), style='v%0.8%green',
            color='orange', fig_kwargs=figkw)
    figkw.update({'x.axis_label': 'Z momentum sum [GeV]'})
    b.histogram(idx=2, data=np.histogram(df.sumMomZ, bins=100), style='o%0.8%purple',
            color='orange', fig_kwargs=figkw)

    
    dlatex = dict(x=FIGDIMS[0]/2+9.5*STEPS[0],
                  y=FIGDIMS[1]-10*STEPS[1],
                  x_units="screen",
                  y_units="screen",
                  text_font_size='10pt')
    latex1 = LatexLabel(
        text="x_{{init}} = {}, y_{{init}} = {}".format(round(FLAGS.x,5), round(FLAGS.y,5)),
        render_mode="css",
        background_fill_alpha=0,
        **dlatex
    )
    for i in range(NFIGS):
        b.get_figure(i).add_layout(latex1)
    
    ##########################################################
    ######## Saving ##########################################
    ##########################################################
    b.save_frame(show=True)
    #b.save_figs(path='.', mode='png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    FLAGS, _ = argparser.add_args(parser)
    run(FLAGS)
