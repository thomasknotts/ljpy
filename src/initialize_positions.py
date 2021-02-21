# initialize_positions is part of ljpy for Lennard Jones simulations.       #	 	                             #
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
# initialize_positions.py                                                  	#
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
This module is part of ljpy. It initializes the coordinates for the 
simulation. If the user specifies a coordinate file in the input file, the 
coordinates are read in from the file. Otherwise, they are generated on a 
lattice.
"""

# Import relevant libraries
import sys, os
import numpy as np
from src.ljpyclasses import site

# This function is passed a simulation object from the main program
# It returns a list of objects of type sites which is all the atoms (sites)
# in the system.
def initializepositions(sim):
    # initialize the atom list
    atom=[]
    
    # If the input file specificies "generate", then place the 
    # specified number of particles on a lattice.
    if sim.icoord == "generate" or sim.icoord == None:
        # zero out the counters
        particle = 0
        case = 0
        
        # Determine the number of particles on each needed on each side
        # of the box
        nlin=np.uint((sim.N/4.0)**(1.0/3.0))
        if nlin**3 < sim.N/4.0: nlin = nlin + 1 # add 1 if not a perfect cube
        
        # Calculate the length of one unit cell
        a=sim.length/nlin
        
        # Loop around all the particles and assign the positions to a lattice.
        for zdir in range(nlin):
            for ydir in range(nlin):
                for xdir in range(nlin):
                    for i in range(4):
                        if particle == sim.N: return(atom)
                        atom.append(site())
                        if case == 0:
                            atom[particle].x=0.0+xdir*a
                            atom[particle].y=0.0+ydir*a
                            atom[particle].z=0.0+zdir*a
                            case=1
                        elif case == 1:
                            atom[particle].x=0.0+xdir*a
                            atom[particle].y=0.5*a+ydir*a
                            atom[particle].z=0.5*a+zdir*a
                            case=2
                        elif case == 2:
                            atom[particle].x=0.5*a+xdir*a
                            atom[particle].y=0.0+ydir*a
                            atom[particle].z=0.5*a+zdir*a
                            case=3
                        else:
                            atom[particle].x=0.5*a+xdir*a
                            atom[particle].y=0.5*a+ydir*a
                            atom[particle].z=0.0+zdir*a
                            case=0
                        particle=particle + 1
    else:
        # Check if the specified input file exists.
        if not os.path.isfile(sim.icoord): 
            print("Input file \"" + sim.icoord +"\" does not exist.\n")
            sys.exit("Error: Specified coordinate file missing.")
        
