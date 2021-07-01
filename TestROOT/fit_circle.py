import os
import pandas as pd
import numpy as np 
from scipy import optimize
from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource as CDS, Range1d, Label
output_file( os.path.basename(__file__)[:-2] + 'html' )

#define constants
XMIN, XMAX = -1000, 50
YMIN, YMAX = -7300, -6250

# get the data 
df = pd.read_csv('track_pos.csv')
df = df[ (df.x > -900) & (df.x < -50) ]

#fit the circle
def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((df['x']-xc)**2 + (df['z']-yc)**2)

def f_2(c):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

CENTER_EST = (XMAX+XMIN)/2, (YMAX+YMIN)/2

center, ier = optimize.leastsq(f_2, CENTER_EST)

xc, yc   = center
Ri       = calc_R(*center)
R        = Ri.mean()
residues = sum((Ri - R)**2)
print(center)
print(Ri)
print(R)
print(residues)

source = CDS( df )

fig = figure()
fig.circle(x=xc, y=yc, color='red', fill_color=None, alpha=0.9, size=R, legend_label='fit')
fig.circle(x='x', y='z', color='blue', source=source, legend_label='data')
Rstring = Label(x=CENTER_EST[0]-50, y=CENTER_EST[1],
                text='Radius: '+str(round(R,2)),
                text_color='red')
fig.add_layout(Rstring)

#plot circle center and radius
fig.xaxis.axis_label = 'x'
fig.yaxis.axis_label = 'z'
fig.x_range = Range1d(XMIN,XMAX)
fig.y_range = Range1d(YMIN,YMAX)

show(fig)


