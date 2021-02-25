# scale_delta is part of ljpy for Lennard Jones simulations.                #
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
# scale_delta.py                                                           	#
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
This module is part of ljpy. It adjusts the maximum displacement for the MC
simulation to obtain an acceptance ratio of 30%.
"""
# Import relevant libraries
import numpy as np

# This function takes a simulation object, and two props objects--one
# for the instantaneous properties and one for the average properties--
# and tries to scale the maximum displacement to acheive a 30% 
# acceptance ratio.

def scale_delta(sim,iprop,aprop):
    # Set the desired acceptance ratio
    dratio=0.3
    
    # Calculate the acceptance ratio since the last time 
    # scale_delta was called
    ratio=np.double(iprop.naccept)/np.double(iprop.ntry)
    
    # Increase the maximum displacement if the ratio is lower than
    # the desired value or decrease the ratio if it is higher
    if sim.dt < 2.0: # set a max on the displacement
        if ratio < dratio-0.02 or ratio > dratio+0.02:
            if ratio < dratio: sim.dt*=0.95
            if ratio > dratio: sim.dt*=1.05
    
    aprop.naccept+=iprop.naccept
    aprop.ntry+=iprop.ntry
    iprop.naccept=0
    iprop.ntry=0