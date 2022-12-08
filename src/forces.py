# forces is part of ljpy for Lennard Jones simulations.                     #
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
# forces.py                                                               	#
#                                                                          	#
# Thomas A. Knotts IV                                                      	#
# Brigham Young University                                                 	#
# Department of Chemical Engineering                                       	#
# Provo, UT  84606                                                         	#
# Email: thomas.knotts@byu.edu                                             	#
# ========================================================================= #
# Version 1.0 - February 2021                                              	#
# Version 2.0 - December 2022 Changed from atom class to arrays for numba. 	#
# ========================================================================= #

"""
This module is part of ljpy. This function calculates the energies and forces
between each Lennard Jones particle. It returns the potential energy
and virial for the state of the list of particles passed.
"""

# Import relevant libraries
import numpy as np
#from numba import jit

# This function is passed a simulation object and a list of site object.
# It returns the potential energy of the system and also assigns the 
# forces on each site. 
#@jit(nopython=True)
def forces(sim, atomx, atomy, atomz):
    # initialize the force arrays and system accumulators
    atomfx=np.zeros(sim.N)
    atomfy=np.zeros(sim.N)
    atomfz=np.zeros(sim.N)
    virial=0.0  # virial portion of pressure
    pe=0.0      # potential energy

    # Variables
    hL=np.float(sim.length*0.5)   # half the box length
    N=np.int(sim.N)
    
    # Zero out the force accumulators for each particle
    # for i in range(sim.N):
        # atomfx[i]=0.0
        # atomfy[i]=0.0
        # atomfz[i]=0.0
       
    # Calculate the forces by looping over all pairs of sites
    for i in range(N-1):
        for j in range(i+1, N):
            # Calculate the distance between sites i and j
            dx=atomx[i]-atomx[j]
            dy=atomy[i]-atomy[j]
            dz=atomz[i]-atomz[j]
            
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
            
            dr2=dx*dx+dy*dy+dz*dz
            
            # Calculate the energy and force for the pair
            if dr2 < sim.rc2: # apply cutoff
                d2=1.0/dr2
                d4=d2*d2
                d8=d4*d4
                d14=d8*d4*d2
                fr=48.0*(d14-0.5*d8)
                
                # components of forces
                atomfx[i]=atomfx[i]+fr*dx
                atomfy[i]=atomfy[i]+fr*dy
                atomfz[i]=atomfz[i]+fr*dz
                atomfx[j]=atomfx[j]-fr*dx
                atomfy[j]=atomfy[j]-fr*dy
                atomfz[j]=atomfz[j]-fr*dz   
                
                # virial and potential energy
                virial=virial+dr2*fr
                pe=pe+4.0*(d14-d8)*dr2

    return(np.float(pe), np.float(virial), atomfx, atomfy, atomfz)