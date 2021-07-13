import os, sys
import pandas as pd
import numpy as np

from scipy import optimize
from scipy.interpolate import CubicSpline, interp1d

from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource as CDS, Range1d, Label
output_file( os.path.basename(__file__)[:-2] + 'html' )

#define constants
DEBUG = True
MOMTRUTH = 1
BSCALE = 1
BSTRENGTH = 0.5 # T
PROTONMASS = 0.938

XSHIFT, ZSHIFT = 2000, 9000
XMIN, XMAX = 650, 2050
YMIN, YMAX = 1300, 2700
assert(XMAX-XMIN == YMAX-YMIN)
# XMIN, XMAX = 1900, 2500
# YMIN, YMAX = 1700, 2050

# get the data and shift it to get positive values
mode = sys.argv[1]
df = pd.read_csv('data/track_' + mode + '.csv')
df.x = df.x + XSHIFT
df.z = df.z + ZSHIFT
spline = interp1d(df.x, df.z)

def get_uniformly_spaced_xarray(x, y, threshold):
    """Extracts uniformly spaced x points from non-uniformly distributed arrays."""
    assert(len(x) == len(y))
    def distance(xold,yold,xnew,ynew):
        a = xnew-xold
        b = ynew-yold
        return np.sqrt(a*a+b*b)

    new_x = []
    xold, yold = x[0], y[0]
    
    for ix,iy in zip(x,y):
        d = distance(xold,yold,ix,iy)
        if d > threshold:
            new_x.append(ix)
            xold, yold = ix, iy
        else:
            continue
        
    return new_x

NSELPOINTS = 7000
dffilt = df # selection can be done here
df['xfilt'] = dffilt.iloc[:NSELPOINTS].x
df['zfilt'] = dffilt.iloc[:NSELPOINTS].z

def fit(x, y):
    """performs circular fit"""
    center_est_ = 1400, XSHIFT
    
    def calc_R_(xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((x-xc)**2 + (y-yc)**2)

    def calc_dist_(c):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R_(*c)
        return Ri - Ri.mean()

    center_, _ = optimize.leastsq(calc_dist_, center_est_, ftol=1e-10, xtol=1e-10)
    Ri_       = calc_R_(*center_)
    R_        = Ri_.mean()
    residues_ = sum((Ri_ - R_)**2)
    if DEBUG:
        print('Center: ', center_)
        print('Radius: ', R_)
        print('Residues: ', residues_)
        
    return center_, R_

#########################################################################
### SOURCES #############################################################
#########################################################################
filt_data = CDS( dict(x=df.xfilt, z=df.zfilt) )

xuniform = get_uniformly_spaced_xarray(df.xfilt, df.zfilt, 10)
interp_data = CDS( dict(x=xuniform, z=spline(xuniform)) )

#########################################################################
##########################################################################
p0 = figure()
p0.circle(x='x', y='z', size=.5, source=filt_data,
          color='blue', legend_label='data')
p0.circle(x='x', y='z', size=2, source=interp_data,
          color='orange', legend_label='spline')
p0.title = 'Uniformly-distributed circular 1D spline'
    
p1 = figure()
p1.title = mode[0].upper() + mode[1:] + ' method'
(xc, yc), R = fit(x=xuniform, y=spline(xuniform))

thetas = np.linspace(0, 2*np.pi, len(xuniform))
xfit = (R*np.cos(thetas)) + xc
yfit = (R*np.sin(thetas)) + yc
interp_data.data['xfit'] = xfit
interp_data.data['yfit'] = yfit

p1.circle(x='xfit', y='yfit', color='red', size=2, 
          source=interp_data,
          legend_label='fit')
p1.circle(x='x', y='z', color='green', size=1,
          source=interp_data,
          legend_label='interpolated data (spline)')

Rstring = Label(x=xc-200, y=yc,
                 text='Radius: '+str(round(R,4))+' cm',
                 text_color='red')

mom = 0.2997924580 * R/100. * BSTRENGTH * BSCALE # divide radius for meters
MomEstStr = Label(x=xc-350, y=yc-70,
                  text='Momentum estimate: '+ str(round(mom,4))+' GeV',
                  text_color='black')
MomTruthStr = Label(x=xc-350, y=yc-140,
                    text='Momentum truth: '+ str(round(MOMTRUTH,4))+' GeV',
                    text_color='black')

p1.add_layout(Rstring)
p1.add_layout(MomEstStr)
p1.add_layout(MomTruthStr)

p3 = figure()
source = CDS( df )
p3.circle(x='index', y='energy', legend_label='energies', size=.5, source=source)
##########################################################################

#plot circle center and radius
# for p in [p0,p1]:
#     p.xaxis.axis_label = 'X [cm]'
#     p.yaxis.axis_label = 'Z [cm]'
#     p.x_range = Range1d(XMIN,XMAX)
#     p.y_range = Range1d(YMIN,YMAX)

show( gridplot( [[p0,p1,p3]] ) )
