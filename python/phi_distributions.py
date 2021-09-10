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
                                  str(FLAGS.step_size) + 'SZ_')

    ##########################################################
    ######## File reading ####################################
    ##########################################################
    posstr = str(FLAGS.x).replace('.', 'p') + "*" + \
        str(FLAGS.y).replace('.', 'p') + "*" + \
        str(FLAGS.energy).replace('.', 'p') + "*" + \
        str(FLAGS.step_size).replace('.', 'p')
    search_str = os.path.join(BASE, 'data/histo_' + FLAGS.mode + '_' + posstr + '*.csv')

    l = glob.glob(search_str)
    NFIGS=len(l)
    print('Pattern: ', search_str)
    print('Number of files: ', NFIGS)
    mass_interaction = []
    for f in l:
        mass_interaction.append( float(re.findall(r'.+_(.*)MInt.csv', f)[0].replace('p', '.')) )

    l = [ x for _,x in sorted(zip(mass_interaction,l), key=lambda pair: pair[0]) ] #overwrite
    mass_interaction = sorted(mass_interaction)
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

    figkw = {'y.axis_label': 'Counts'}

    dlatex = dict(x=3*STEPS[0],
                  y=FIGDIMS[1]/2+36*STEPS[1],
                  x_units="screen",
                  y_units="screen",
                  text_font_size='10pt')

    phistr = '\u03D5'

    for idx,f in enumerate(l):
        df = pd.read_csv( f )

        b.get_figure(idx).add_layout(
            LatexLabel(
                text="m = {}GeV".format(mass_interaction[idx],5),
                render_mode="css",
                background_fill_alpha=0,
                **dlatex
            ) )
        
        figkw.update({'x.axis_label': phistr + ' [rad]',
                      'x_range': Range1d(0,2*np.pi)})
        b.histogram(idx=idx,
                    data=np.histogram(df.Phi, bins=100),
                    color='purple', fig_kwargs=figkw)

    ##########################################################
    ######## Saving ##########################################
    ##########################################################
    b.save_frame(nrows=3, ncols=4, show=True)
    #b.save_frame(iframe=0, layout=[[5,3,4],[0,1,2]], show=True)
    #b.save_figs(path='.', mode='png')

    script, div = components( b.get_figure(8) )
    with open('histos.html', 'a') as f:
        f.write(script)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    FLAGS, _ = argparser.add_args(parser, 'sim')
    run(FLAGS)
