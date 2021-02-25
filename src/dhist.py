# dhist is part of ljpy for Lennard Jones simulations.                      #
# Copyright (C) 2021 Thomas Allen Knotts IV - All Rights Reserved          	#
#																		   	#
# This program is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation, either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                          	#
# This program is distributed in the hope that it will be useful,   	 	#
# but WITHOUT ANY WARRANTY; without even the implied warranty of           	#
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            	#
# GNU General Public License for more details.                             	#
#                                                                          	#
# You should have received a copy of the GNU General Public License        	#
# along with this program.  If not, see <http://www.gnu.org/licenses/>.   	#

# ========================================================================= #
# dhist.py                                                      	        #
#                                                                          	#
# Thomas A. Knotts IV                                                      	#
# Brigham Young University                                                 	#
# Department of Chemical Engineering                                       	#
# Provo, UT  84606                                                         	#
# Email: thomas.knotts@byu.edu                                             	#
# ========================================================================= #
# Version 1.0 - February 2021                                              	#
# ========================================================================= #

"""
This module is a library to dynamically create histograms. It is based off
the gsl_histogram library.
"""
import numpy as np

class hist(object):   
    def __init__(self, xmin, xmax, N):
        # xmin is the minimum value of x accumulated in the histogram
        # xmax is the maximum value of x accumulated in the histogram
        # bin_width is the amount of x covered by each bin
        # bin holds the accumulated value (the y value) of each bin
        # range[i] is the left x value of each bin
        # mrange[i] is the middle x value of each bin. Used for graphing.
        # bin[i] corresponds to range[i] <= x < range[i+1], so the histogram 
        # includes xmin but does not include xmax.
        self.xmin=xmin
        self.xmax=xmax
        self.N=N
        self.bin_width=(xmax-xmin)/N
        self.range=np.zeros(N)
        self.mrange=np.zeros(N)        
        for i in range(N):
            self.range[i]=xmin+self.bin_width*i
            self.mrange[i]=xmin+self.bin_width*(i+0.5)
        self.bin=np.zeros(N)
        
        return
    
    def accumulate(self, x, weight):
        print("x=",x)
        if x>=self.xmin and x<self.xmax:
            index=int((x-self.xmin)/self.bin_width)
            print("index =",index)
            self.bin[index]+=weight 
        return
        
    def increment(self, x):     
        self.accumulate(x,1.0)
        return

    def write(self,fp):
        for i in range(self.N): 
            fp.write("{:>10}\t{:13.6f}\t{:13.6f}\n".format(i+1, \
                     self.mrange[i], self.bin[i]))

    def print(self):
        for i in range(self.N): 
            print("{:>10}\t{:13.6f}\t{:13.6f}\n".format(i+1, \
                  self.mrange[i], self.bin[i]))

def clone(src):
    xmin=src.xmin
    xmax=src.xmax
    N=len(src.range)
    h=hist(xmin,xmax,N)
    h.range=src.range
    h.mrange=src.mrange
    h.bin=src.bin
    return(h)

def equalbins(h1, h2):
    if(h1.N != h2.N): return(0)
    for i in range(h1.N):
        if(h1.range[i] != h2.range[i]): return(0)
    return(1)

def div(h1, h2):
    if not equalbins(h1,h2):
        print("The histograms have differen binning.\n")
        return(None)
    h=clone(h1)
    for i in range(h.N):
        if h2.bin[i] != 0.0:
            h.bin[i] = h1.bin[i] / h2.bin[i]
        else:
            h.bin[i] = 0.0
    return(h)