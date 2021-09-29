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
from bokeh.models import Range1d
    
def run(FLAGS):
    ##########################################################
    ######## Constants #######################################
    ##########################################################
    NFRAMES=4
    NFIGS=6
    FIGDIMS1 = 600, 600
    FIGDIMS2 = 800, 700
    STEPS = FIGDIMS1[0]/100, FIGDIMS1[1]/100, 

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

    figwidths = [FIGDIMS1[0] for _ in range(NFRAMES)]
    figheights = [FIGDIMS1[0] for _ in range(NFRAMES)]
    
    b = bkp.BokehPlot( [HTML_BASE_NAME + 'Cat' + str(x) + '.html' for x in range(NFRAMES)],
                       nframes=NFRAMES,
                       nfigs=[NFIGS if i!=NFRAMES-1 else 3 for i in range(NFRAMES)],
                       nwidgets=[0 for _ in range(NFRAMES)],
                       fig_width=figwidths,
                       fig_height=figheights )

    additional_kwargs = {'text_font_size': '9pt', 'text_font_style': 'italic',
                         'x_units': 'screen', 'y_units': 'screen',
                         'angle': np.pi/2}
    dlabel = dict(x=FIGDIMS1[0]/2+3*STEPS[0],
                  y=FIGDIMS1[1]-32*STEPS[1], **additional_kwargs)

    figkw = {'y.axis_label': 'Counts'}
    figkw.update({'x.axis_label': 'X momentum sum [GeV]'})
    b.histogram(idx=0, iframe=3,
                data=np.histogram(df.sumMomXAbs, bins=100),
                color='orange', fig_kwargs=figkw)
    figkw.update({'x.axis_label': 'Y momentum sum [GeV]'})
    b.histogram(idx=1, iframe=3,
                data=np.histogram(df.sumMomYAbs, bins=100),
                color='orange', fig_kwargs=figkw)
    figkw.update({'x.axis_label': 'Z momentum sum [GeV]'})
    b.histogram(idx=2, iframe=3,
                data=np.histogram(df.sumMomZ, bins=100),
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
    for i in range( b.get_nfigs(iframe=3) ):
        b.get_figure(idx=i, iframe=3).add_layout(latex1)

    psistr = '\u03A8'
    phistr = '\u03D5'
    pistr = '\u03C0'

    selections = [ (df.cat1>=0) & (df.cat1<=2), df.cat1==1, df.cat1==2]

    for iframe,sel in enumerate(selections):       
        figkw.update({'x.axis_label': psistr + ': Angle/2 between ' + psistr + 'A and ' + psistr + 'B [rad]',
                      'x_range': Range1d(0,np.pi)})
        b.histogram(idx=0, iframe=iframe,
                    data=np.histogram(df.Psi[sel], bins=100),
                    color='purple', fig_kwargs=figkw)

        figkw.update({'x.axis_label': phistr + ' [rad]',
                      'x_range': Range1d(1,6)})
        b.histogram(idx=1, iframe=iframe,
                    data=np.histogram(df.Phi[sel], bins=100),
                    color='purple', fig_kwargs=figkw)
        
        figkw.update({'x.axis_label': 'Cos(' + phistr + '-' + psistr + ') [rad]',
                      'x_range': Range1d(-1,1)})
        b.histogram(idx=2, iframe=iframe,
                    data=np.histogram(df.Cos[sel], bins=100),
                    color='purple', fig_kwargs=figkw)
        figkw.pop('x_range')
        
        #psiA and psiB
        figkw.update({'x.axis_label': psistr + 'A [rad]'})
        b.histogram(idx=3, iframe=iframe,
                    data=np.histogram(df.PsiA[sel], bins=100),
                    color='red', fig_kwargs=figkw)
        figkw.update({'x.axis_label': psistr + 'B [rad]'})
        b.histogram(idx=4, iframe=iframe,
                    data=np.histogram(df.PsiB[sel], bins=100),
                    color='red', fig_kwargs=figkw)

        figkw.update({'x.axis_label': psistr + 'A [rad]',
                      'y.axis_label': psistr + 'B [rad]'})
        b.histogram(idx=5, iframe=iframe,
                    data=np.histogram2d(df.PsiA[sel], df.PsiB[sel], bins=50),
                    style='quad%Viridis',
                    fig_kwargs=figkw)

    ##########################################################
    ######## Saving ##########################################
    ##########################################################
    b.save_frame(iframe=0, layout=[[5,3,4],[0,1,2]], show=True)
    b.save_frame(iframe=1, layout=[[5,3,4],[0,1,2]], show=False)
    b.save_frame(iframe=2, layout=[[5,3,4],[0,1,2]], show=False)
    b.save_frame(iframe=3, nrows=1, ncols=3, show=True)
    #b.save_figs(path='.', mode='png')

    script, div = components( b.get_figure(8) )
    with open('histos.html', 'a') as f:
        f.write(script)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    FLAGS, _ = argparser.add_args(parser, 'sim')
    run(FLAGS)
