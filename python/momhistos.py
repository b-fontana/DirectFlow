import sys, os
import glob
import pandas as pd
import numpy as np
import bokehplot as bkp
import argparser
import argparse

from latex import LatexLabel
from bokeh.layouts import gridplot
from bokeh.embed import components
    
def run(FLAGS):
    ##########################################################
    ######## Constants #######################################
    ##########################################################
    NFIGS=9
    FIGDIMS1 = 600, 600
    FIGDIMS2 = 800, 700
    STEPS = FIGDIMS1[0]/100, FIGDIMS1[1]/100, 

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
    # df['sumMomXAbs'] = np.abs(df.sumMomX)
    # df['sumMomYAbs'] = np.abs(df.sumMomY)
    df['sumMomXAbs'] = df.sumMomX
    df['sumMomYAbs'] = df.sumMomY
    ##########################################################
    ######## Plotting ########################################
    ##########################################################

    figwidths = [FIGDIMS1[0] for _ in range(NFIGS-1)]
    figwidths.extend([FIGDIMS2[0]])
    figheights = [FIGDIMS1[0] for _ in range(NFIGS-1)]
    figheights.extend([FIGDIMS2[1]])
    
    b = bkp.BokehPlot( HTML_NAME, nfigs=NFIGS, nwidgets=0,
                       fig_width=figwidths,
                       fig_height=figheights )

    additional_kwargs = {'text_font_size': '9pt', 'text_font_style': 'italic',
                         'x_units': 'screen', 'y_units': 'screen',
                         'angle': np.pi/2}
    dlabel = dict(x=FIGDIMS1[0]/2+3*STEPS[0],
                  y=FIGDIMS1[1]-32*STEPS[1], **additional_kwargs)

    figkw = {'y.axis_label': 'Counts'}
    figkw.update({'x.axis_label': 'X momentum sum [GeV]'})
    b.histogram(idx=0, data=np.histogram(df.sumMomXAbs, bins=100),
            color='orange', fig_kwargs=figkw)
    figkw.update({'x.axis_label': 'Y momentum sum [GeV]'})
    b.histogram(idx=1, data=np.histogram(df.sumMomYAbs, bins=100),
            color='orange', fig_kwargs=figkw)
    figkw.update({'x.axis_label': 'Z momentum sum [GeV]'})
    b.histogram(idx=2, data=np.histogram(df.sumMomZ, bins=100),
            color='orange', fig_kwargs=figkw)

    
    dlatex = dict(x=FIGDIMS1[0]/2+3.*STEPS[0],
                  y=FIGDIMS1[1]-10*STEPS[1],
                  x_units="screen",
                  y_units="screen",
                  text_font_size='10pt')
    latex1 = LatexLabel(
        text="x_{{init}} = {}cm, y_{{init}} = {}cm".format(round(FLAGS.x,5),round(FLAGS.y,5)),
        render_mode="css",
        background_fill_alpha=0,
        **dlatex
    )
    for i in range(NFIGS-1):
        b.get_figure(i).add_layout(latex1)

    psistr = '\u03A8'
    phistr = '\u03D5'
    pistr = '\u03C0'
    figkw.update({'x.axis_label': psistr + '=' + psistr + 'A +' + psistr + 'B +' + pistr +  ' [rad]'})
    b.histogram(idx=3, data=np.histogram(df.Psi, bins=100),
                color='purple', fig_kwargs=figkw)

    figkw.update({'x.axis_label': phistr + ' [rad]'})
    b.histogram(idx=4, data=np.histogram(df.Phi, bins=100),
                color='purple', fig_kwargs=figkw)

    figkw.update({'x.axis_label': 'Cos(' + phistr + '+' + psistr + ') [rad]'})
    b.histogram(idx=5, data=np.histogram(np.cos(df.Psi+df.Phi), bins=100),
                color='purple', fig_kwargs=figkw)


    #psiA and psiB
    figkw.update({'x.axis_label': psistr + 'A [rad]'})
    b.histogram(idx=6, data=np.histogram(df.PsiA, bins=100),
                color='red', fig_kwargs=figkw)
    figkw.update({'x.axis_label': psistr + 'B [rad]'})
    b.histogram(idx=7, data=np.histogram(df.PsiB, bins=100),
                color='red', fig_kwargs=figkw)

    figkw.update({'x.axis_label': psistr + 'A + ' + pistr + ' [rad]',
                  'y.axis_label': psistr + 'B [rad]'})
    b.histogram(idx=8, data=np.histogram2d(df.PsiA+np.pi, df.PsiB, bins=50),
                style='quad%Viridis',
                fig_kwargs=figkw)

    ##########################################################
    ######## Saving ##########################################
    ##########################################################
    #b.save_frame(nrows=3, ncols=3, show=True)
    b.save_frame(layout=[[0,1,2],[8,6,7],[3,4,5]], show=True)
    #b.save_figs(path='.', mode='png')

    script, div = components( b.get_figure(8) )
    with open('histos.html', 'a') as f:
        f.write(script)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    FLAGS, _ = argparser.add_args(parser)
    run(FLAGS)
