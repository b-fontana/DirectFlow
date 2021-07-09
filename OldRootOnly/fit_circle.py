import os
import pandas as pd
import numpy as np 
from scipy import optimize
from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource as CDS, Range1d, Label
output_file( os.path.basename(__file__)[:-2] + 'html' )

#define constants
DEBUG = True
XMIN, XMAX = 1000, 3200
YMIN, YMAX = 200, 2400
assert(XMAX-XMIN == YMAX-YMIN)
CENTER_EST = 2100, 1700

# get the data and shift it to get positive values
df = pd.read_csv('track_pos.csv')
df.x = df.x + 2000
df.z = df.z + 9000
dffilt = df
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
    print('Resiudes: ', residues)

source = CDS( df )

fig = figure()
fig.circle(x=xc, y=yc, color='red', fill_color=None, alpha=0.9, radius=R, legend_label='fit', radius_units='data')
fig.circle(x='z', y='x', color='blue', source=source, size=1, legend_label='full data')
#fig.triangle(x='zfilt', y='xfilt', color='green', source=source, legend_label='fitted data')

Rstring = Label(x=xc-200, y=yc,
                 text='Radius: '+str(round(R,2))+' cm',
                 text_color='red')

mom = 0.2997924580 * R/100. * 3.529 * 350. # divide radius for meters
MomentumStr = Label(x=xc-350, y=yc-100,
                    text='Momentum estimate: '+str(round(mom,2))+' GeV',
                 text_color='black')

fig.add_layout(Rstring)
fig.add_layout(MomentumStr)

#plot circle center and radius
fig.xaxis.axis_label = 'Z [cm]'
fig.yaxis.axis_label = 'X [cm]'
#fig.x_range = Range1d(XMIN,XMAX)
#fig.y_range = Range1d(YMIN,YMAX)

show(fig)


