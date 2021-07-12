import os, sys
import pandas as pd
import numpy as np 
from scipy import optimize
from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.layouts import row
from bokeh.models import ColumnDataSource as CDS, Range1d, Label
output_file( os.path.basename(__file__)[:-2] + 'html' )

#define constants
DEBUG = True
MOMTRUTH = 1500
PROTONMASS = 0.938
XMIN, XMAX = 1500, 2700
YMIN, YMAX = 1000, 2200
assert(XMAX-XMIN == YMAX-YMIN)
CENTER_EST = 2100, 1600

# get the data and shift it to get positive values
df = pd.read_csv('data/track_' + sys.argv[1] + '.csv')
df.x = df.x + 2000
df.z = df.z + 9000
dffilt = df # selection can be done here
df['xfilt'] = dffilt.x
df['zfilt'] = dffilt.z

#fit the circle
def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((df['zfilt']-xc)**2 + (df['xfilt']-yc)**2)

def f_2(c):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

center, ier = optimize.leastsq(f_2, CENTER_EST)

xc, yc   = center
Ri       = calc_R(*center)
R        = Ri.mean()
residues = sum((Ri - R)**2)
if DEBUG:
    print('Center: ', center)
    print('Radius: ', R)
    print('Residues: ', residues)

source = CDS( df )

p1 = figure()
p1.title = sys.argv[1]
p1.circle(x=xc, y=yc, color='red', fill_color=None, alpha=0.9, radius=R, legend_label='fit', radius_units='data')
#p1.circle(x='z', y='x', color='blue', source=source, size=1, legend_label='full data')
p1.circle(x='zfilt', y='xfilt', color='green', source=source, size=1, legend_label='fitted data')

p2 = figure()
p2.circle(x='index', y='energy', legend_label='energies', size=2, source=source)

Rstring = Label(x=xc-200, y=yc,
                 text='Radius: '+str(round(R,2))+' cm',
                 text_color='red')

mom = 0.2997924580 * R/100. * 3.529 * 350. # divide radius for meters
MomEstStr = Label(x=xc-350, y=yc-70,
                  text='Momentum estimate: '+ str(round(mom,2))+' GeV',
                  text_color='black')
MomTruthStr = Label(x=xc-350, y=yc-140,
                    text='Momentum truth: '+ str(round(np.sqrt(MOMTRUTH*MOMTRUTH+PROTONMASS*PROTONMASS),2))+' GeV',
                    text_color='black')

p1.add_layout(Rstring)
p1.add_layout(MomEstStr)
p1.add_layout(MomTruthStr)

#plot circle center and radius
p1.xaxis.axis_label = 'Z [cm]'
p1.yaxis.axis_label = 'X [cm]'
p1.x_range = Range1d(XMIN,XMAX)
p1.y_range = Range1d(YMIN,YMAX)

show(row(p1,p2))


