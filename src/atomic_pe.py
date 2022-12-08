# atomic_pe is part of ljpy for Lennard Jones simulations.                  #
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
# atomic_pe.py                                                    	        #
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
This module is part of ljpy. It has a function that calculates the 
potential energy of each atom with each of its neighbors.  It is used in MC
to calculate the difference in energies for use in the metropolis criterion.    
"""
# Import relevant libraries
import numpy as np
import numba as nb

# This functions take a simulation object, a list of site objects, the
# particle that is moved, and the coordinates of the particle that is moved.
# It calculate the potential energy of this particle with all the other 
# particles in the system. This subroutine is called twice for one monte carlo
# move. It returns the potential energy of the particle for the state
# in atom.
#@nb.njit(nopython=True)
def atomic_pe(sim, atomx, atomy, atomz, particle):
    # Variables
    junk=np.float64(34)
    hL=np.float64(np.float64(sim.length)/2.0)
    N=np.int64(sim.N)
    length=np.float64(sim.length)
    
    # Zero out the potential energy accumulator
    u=np.float64(0.0)
    
    # Loop around the neighbors of the selected particle to get the
    # energy.
    for i in range(np.int64(sim.N)):
        if particle != i: # exclude particle i from itself
            dx=np.float64(atomx[i]-atomx[particle])
            dy=np.float64(atomy[i]-atomy[particle])
            dz=np.float64(atomz[i]-atomz[particle])
            
            # Apply minimum image convection
            if np.abs(dx)>hL:
                if dx < 0.0: dx=dx+length
                else: dx=dx-length
            if np.abs(dy)>hL:
                if dy < 0.0: dy=dy+length
                else: dy=dy-sim.length
            if np.abs(dz)>hL:
                if dz < 0.0: dz=dz+length
                else: dz=dz-length
            
            # Distance and energy calculation        
            dr=dx*dx+dy*dy+dz*dz
            if dr<sim.rc2:
                dr2=1.0/dr
                dr4=dr2*dr2
                dr6=dr4*dr2
                dr12=dr6*dr6
                u+=4.0*(dr12-dr6)
                
    return(u)