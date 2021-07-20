#!/usr/bin/env python3
from std_lib import *

import os
import numpy as np
import cv2

from scipy import interpolate
import matplotlib.pyplot as plt

CUR_DIR = os.getcwd()
CIRCLE_FILE = "circle.txt"
CURVE_FILE  = "curve.txt"
SQUARE_FILE = "square.txt"

#test
CIRCLE_NAME = "circle"
CURVE_NAME  = "curve"
SQUARE_NAME = "square"

SYS_TOKEN_CNT = 2   # x, y

total_pt_cnt = 0        # total no. of points      
x_arr = np.array([])    # x position set
y_arr = np.array([])    # y position set

def convert_coord_to_array(file_path):
    global total_pt_cnt
    global x_arr
    global y_arr

    if file_path == "":
        return FALSE

    with open(file_path) as f:
        content = f.readlines()

    content = [x.strip() for x in content] 

    total_pt_cnt = len(content)

    if (total_pt_cnt <= 0):
        return FALSE

    ##
    x_arr = np.empty((0, total_pt_cnt))
    y_arr = np.empty((0, total_pt_cnt))

    #compare the first and last x 
    # if ((content[0][0]) > (content[-1])):
        # is_reverse = TRUE

    for x in content:
        token_cnt = get_token_cnt(x, ',') 

        if (token_cnt != SYS_TOKEN_CNT):
            return FALSE

        for idx in range(token_cnt):
            token_string = get_token_string(x, ',', idx)
            token_string = token_string.strip()
            if (not token_string.isdigit()): 
                return FALSE

            # save x, y set
            if (idx == 0):
                x_arr = np.append(x_arr, int(token_string))
            else:
                y_arr = np.append(y_arr, int(token_string))

    return TRUE

def linear_interpolation(fig, axs):
    xnew = np.linspace(x_arr.min(), x_arr.max(), len(x_arr))
    f = interpolate.interp1d(xnew , y_arr)

    axs.plot(xnew, f(xnew))
    axs.set_title('linear')

def cubic_interpolation(fig, axs):
    xnew = np.linspace(x_arr.min(), x_arr.max(), len(x_arr))
    f = interpolate.interp1d(xnew , y_arr, kind='cubic')

    axs.plot(xnew, f(xnew))
    axs.set_title('cubic')

def cubic_spline_interpolation(fig, axs):
    xnew = np.linspace(x_arr.min(), x_arr.max(), len(x_arr))
    tck = interpolate.splrep(x_arr, y_arr, s=0) #always fail (ValueError: Error on input data)
    ynew = interpolate.splev(xnew, tck, der=0)

    axs.plot(xnew, ynew)
    axs.set_title('cubic spline')

def parametric_spline_interpolation(fig, axs):
    xnew = np.linspace(x_arr.min(), x_arr.max(), len(x_arr))
    tck, u = interpolate.splprep([x_arr, y_arr], s=0)
    out = interpolate.splev(xnew, tck)

    axs.plot(out[0], out[1])
    axs.set_title('parametric spline')

def univariate_spline_interpolated(fig, axs):   
    s = interpolate.InterpolatedUnivariateSpline(x_arr, y_arr)# ValueError: x must be strictly increasing
    xnew = np.linspace(x_arr.min(), x_arr.max(), len(x_arr))
    ynew = s(xnew)

    axs.plot(xnew, ynew)
    axs.set_title('univariate spline')

def rbf(fig, axs):
    xnew = np.linspace(x_arr.min(), x_arr.max(), len(x_arr))
    rbf = interpolate.Rbf(x_arr, y_arr) # numpy.linalg.linalg.LinAlgError: Matrix is singular.
    fi = rbf(xnew)

    axs.plot(xnew, fi)
    axs.set_title('rbf')

def interpolation():
    fig, axs = plt.subplots(nrows=4)
    axs[0].plot(x_arr, y_arr, 'r-')
    axs[0].set_title('org')

    cubic_interpolation(fig, axs[1])
    # cubic_spline_interpolation(fig, axs[2]) 
    parametric_spline_interpolation(fig, axs[2])
    # univariate_spline_interpolated(fig, axs[3])
    # rbf(fig, axs[3])        
    linear_interpolation(fig, axs[3])

    plt.show()

#------- main -------
if __name__ == "__main__":
    interpolation()
    # np.seterr(divide='ignore', invalid='ignore')
