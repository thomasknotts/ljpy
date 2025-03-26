# move is part of ljpy for Lennard Jones simulations.                       #
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
# move.py                                                        	        #
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
This module is part of ljpy. It has a function that generates a random
displacement on a random particle, calculates the new energy, and accepts
or rejects the new position according to the Metropolis criterion.  It is
the main propogation subroutine for an MC simulation.    
"""

# Import relevant libraries
import random
import numpy as np
from src.atomic_pe import atomic_pe
from src.forces import forces
from numba import njit

# This function accepts a simulation object, a list of site objects, the
# current state of the random number generator, and a property object.
# It returns True if the move is accepted. It returns False if the move is
# rejected.
@njit
def move(sim, atom, iprop):
    # Set the state of the random number generator
    
    # Select a random particle
    iprop.ntry+=1
    particle=random.randint(0, sim.N-1)
    
    # Save the current (old) position of the system in case
    # the move is not accepted 
    xold=atom[particle].x
    yold=atom[particle].y
    zold=atom[particle].z
                           
    # Propose a new move                        
    xnew=atom[particle].x + random.uniform(-1, 1)*sim.dt
    ynew=atom[particle].y + random.uniform(-1, 1)*sim.dt
    znew=atom[particle].z + random.uniform(-1, 1)*sim.dt

    # Apply periodic boundary conditions
    if xnew<0.0: 
        xnew+=sim.length
    elif xnew > sim.length:
        xnew-=sim.length
    
    if ynew<0.0: 
        ynew+=sim.length
    elif ynew > sim.length:
        ynew-=sim.length
    
    if znew<0.0: 
        znew+=sim.length
    elif znew > sim.length:
        znew-=sim.length
    
    # Calculate the energy of the proposed move
    atom[particle].x=xnew
    atom[particle].y=ynew
    atom[particle].z=znew
    penew=atomic_pe(sim, atom, particle)
    
    # Calculate the different in energy between the 
    # proposed and old state of the system
    de=penew-atom[particle].pe
    
    # Accept/Reject the move
    if random.uniform(0,1) < np.exp(-de/sim.T): # accept
        iprop.naccept+=1
        iprop.pe, iprop.virial=forces(sim, atom)
        iprop.pe2=iprop.pe*iprop.pe
        atom[particle].pe=penew # set "current" site pe to new pe 
        return(True)
    else: #reject
        # revert the positions to the "old" or previous state
        atom[particle].x=xold
        atom[particle].y=yold
        atom[particle].z=zold
        return(False)
    
    
    