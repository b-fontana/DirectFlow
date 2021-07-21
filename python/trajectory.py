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
    NFIGS=2
    FIGDIMS = 900, 300
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
    search_str = os.path.join(BASE, 'data/track_' + FLAGS.mode + '_' + posstr + '*.csv')

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

    figkw = {'x.axis_label': 'Z [cm]',
             'y.fixed_location': 0.}
    b.graph(idx=0, data=[df.z1,df.x1], style='circle',
            color='orange', size=1, line=False,
            fig_kwargs=figkw)
    b.graph(idx=0, data=[df.z2,df.x2], style='circle',
            color='blue', size=1, line=False,
            fig_kwargs=figkw)
    xaxis_label = 'X [cm]'
    b.label(xaxis_label, idx=0, **dlabel)

    figkw.update(x_range=b.get_figure(idx=0).x_range)
    b.graph(idx=1, data=[df.z1,df.y1], style='circle',
            color='orange', size=1, line=False,
            fig_kwargs=figkw)
    b.graph(idx=1, data=[df.z2,df.y2], style='circle',
            color='blue', size=1, line=False,
            fig_kwargs=figkw)
    xaxis_label = 'Y [cm]'
    b.label(xaxis_label, idx=1, **dlabel)

    dlatex = dict(x=FIGDIMS[0]-16*STEPS[0],
                  y=FIGDIMS[1]-25*STEPS[1],
                  x_units="screen",
                  y_units="screen",
                  text_font_size='10pt')
    latex1 = LatexLabel(
        text="x_{{init}} = {}, y_{{init}} = {}".format(round(FLAGS.x,5), round(FLAGS.y,5)),
        render_mode="css",
        background_fill_alpha=0,
        **dlatex
    )
    b.get_figure(0).add_layout(latex1)

    # latex = LatexLabel(
    #     text="f = \\sum_{n=1}^\\infty\\frac{-e^{i\\pi}}{2^n}!",
    #     x_units="screen",
    #     y_units="screen",
    #     render_mode="css",
    #     text_font_size="21px",
    #     background_fill_alpha=0,
    #     **dlatex
    # )
    # b.get_figure(0).add_layout(latex)
    
    ##########################################################
    ######## Saving ##########################################
    ##########################################################
    b.save_frame(ncols=1, show=True)
    #b.save_figs(path='.', mode='png')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    FLAGS, _ = argparser.add_args(parser)
    run(FLAGS)
