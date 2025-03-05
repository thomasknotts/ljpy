# rdf is part of ljpy for Lennard Jones simulations.                        #
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
# rdf.py                                                                 	#
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
This module is part of ljpy. It contains functionality to calculate the 
radial distrubtion function from the simulation.
"""
# import relevant libraries
import numpy as np
from numba import njit

#@njit
def rdf_accumulate(sim, atom, h):
    # Variables
    hL=sim.length*0.5   # half the box length
    N=np.int64(sim.N)
    
    # Loop around the atoms to calculate the distances
    for i in range(N-1):
        for j in range(i+1, N):
            # Calculate the distance between sites i and j
            dx=atom[i].x-atom[j].x
            dy=atom[i].y-atom[j].y
            dz=atom[i].z-atom[j].z
            
            # Minimum image convention
            if np.abs(dx)>hL:
                if dx < 0.0: dx=dx+sim.length
                else: dx=dx-sim.length
            if np.abs(dy)>hL:
                if dy < 0.0: dy=dy+sim.length
                else: dy=dy-sim.length
            if np.abs(dz)>hL:
                if dz < 0.0: dz=dz+sim.length
                else: dz=dz-sim.length
            
            dr=np.sqrt(dx*dx+dy*dy+dz*dz)
            
            # Increment the histogram for the calculated
            # value of dr
            h.increment(dr)

def rdf_finalize(sim, h, Ncalls):
    # Variables
    sphere=4.0/3.0*np.pi*sim.rho    # constant to calculate volume for a shell
    
    # Loop over all entries in the rdf histogram
    # to calculate the rdf for each entry.
    for i in range(h.N):
        r1=h.range[i]                       # small radius for shell
        r2=h.range[i]+h.bin_width           # larger radius for shell
        nideal=sphere*(r2*r2*r2 - r1*r1*r1) # number of particle in shell of sim.rho
        h.bin[i]=h.bin[i]/Ncalls/nideal/sim.N*2.0
    
            
            

