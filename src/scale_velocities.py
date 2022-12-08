# scale_velocities is part of ljpy for Lennard Jones simulations.           #
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
# scale_velocities.py                                                   	#
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
This module has a funciton that scales the velocities to the desired 
temperature. This implements an isokinetic method which should only be used to
equilibrate an MD simulation to the desired average run temperature.
It should note be used in production steps as it doesn't generate the
correct NVT distribution of velocities.
"""

# Import relevant libraries
import numpy as np

# This function is passed the simulation object, a list of site objects and 
# the desired temperature.
# It scales the velocities to the desired temperature.
def scalevelocities(sim,atomvx, atomvy, atomvz, temp):
    scale=np.sqrt(sim.T/temp)
    for i in range(sim.N):
        atomvx[i]=atomvx[i]*scale;
        atomvy[i]=atomvy[i]*scale;
        atomvz[i]=atomvz[i]*scale;