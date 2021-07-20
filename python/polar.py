# Plot polar data using Bokeh
# Author: Jonathan J. Helmus

import numpy as np
from bokeh.plotting import figure, output_file, show
from bokeh.palettes import brewer
output_file("sample.html")

# Data in polar coordinates
azimuths = np.linspace(0, 2 * np.pi, 200)
ranges = np.array([2, 4, 6, 8])

# convert to x and y meshes
azimuths_m, ranges_m = np.meshgrid(azimuths, ranges)
xx = ranges_m * np.cos(azimuths_m)
yy = ranges_m * np.sin(azimuths_m)

# array of color values
C = np.random.randint(0, 10, xx.shape)
colors = brewer["Spectral"][10]


def pcolor(C, xx, yy):
    """ Create a pseudocolor like-plot. Note the last points in C ignored. """
    xx_count, yy_count = xx.shape

    f = figure()
    for j in range(yy_count-1):
        for i in range(xx_count-1):

            ys = (yy[i, j], yy[i, j+1], yy[i+1, j+1], yy[i+1, j])
            xs = (xx[i, j], xx[i, j+1], xx[i+1, j+1], xx[i+1, j])
            fill_color = colors[C[i, j]]
            f.patch(ys, xs, fill_color=fill_color)
    return f

# create plot

#hold()
f = pcolor(C, xx, yy)
show(f)
