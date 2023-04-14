# visc is part of ljpy for Lennard Jones simulations.                       #
# Copyright (C) 2023 Thomas Allen Knotts IV - All Rights Reserved          	#
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
# visc.py                                                                 	#
#                                                                          	#
# Thomas A. Knotts IV                                                      	#
# Brigham Young University                                                 	#
# Department of Chemical Engineering                                       	#
# Provo, UT  84606                                                         	#
# Email: thomas.knotts@byu.edu                                             	#
# ========================================================================= #
# Version 1.0 - April 2023                                              	#
# ========================================================================= #

"""
This module is part of ljpy. It contains functionality to calculate the 
viscosity from the autocorrelation pressure tensor.
"""
# import relevant libraries
import numpy as np
from numba import njit

@njit
def ptensor(sim,atom,iprop):
    V=sim.length*sim.length*sim.length
    tensor=np.zeros(6,dtype=np.float64)
    for i in range(sim.N):
        tensor[0]+=atom[i].vx*atom[i].vx
        tensor[1]+=atom[i].vy*atom[i].vy
        tensor[2]+=atom[i].vz*atom[i].vz
        tensor[3]+=atom[i].vx*atom[i].vy
        tensor[4]+=atom[i].vx*atom[i].vz
        tensor[5]+=atom[i].vy*atom[i].vz
    
    for i in range(6):
        tensor[i]+=iprop.stress[i]
    
    return(tensor/V)
    
    
            
            

