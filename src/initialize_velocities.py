# initialize_velocities is part of ljpy for Lennard Jones simulations.      #
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
# initialize_velocities.py                                                	#
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
This module is part of ljpy. It initializes the velocities for the 
simulation. If the user specifies a velocity file in the input file, the 
velocities are read in from the file. Otherwise, they are generated to give a
distribution of velocities with an average magnitude corresponding
to the desired temperature of the simulation.
"""

# Import relevant libraries
import sys, os, random
import numpy as np
from src.momentum_correct import zeromomentum
from src.kinetic import temperature
from src.scale_velocities import scalevelocities


# This function is passed a simulation object and a list of site objects
# from the main program. It sets the velocities of each particle in the
# site object.
def initializevelocities(sim,atom):
    # If the input file specificies "generate" or the vel keywork is
    # omitted, then randomly generate the velocities
    if sim.ivel == "generate" or sim.ivel == "":
        # Loop around the sites and assign random velocities from -1.0 to 1.0  
        for i in range(sim.N):
            atom[i].vx=random.uniform(-1, 1)
            atom[i].vy=random.uniform(-1, 1)
            atom[i].vz=random.uniform(-1, 1)
            
        # Zero out the linear momentum
        momentum_flag=zeromomentum(atom)
        if momentum_flag:
            sys.exit("The linear momentum could not be zeroed out when " +
                     "the velocities were initialized.")
        # Determine the temperature of the randomly-assigned velocities
        T=temperature(atom)
        
        # Scale the temperatures to the system temperature
        scalevelocities(sim,atom,T)
        
    # If a file with velocities is supplied, read the velocities.
    else:
        # Check if the specified input file exists.
        print("Reading velocities from ", sim.ivel)
        if not os.path.isfile(sim.ivel): 
            print("Input file \"" + sim.ivel +"\" does not exist.\n")
            sys.exit("Error: Specified velocity file missing.")
        
        # Open the file with the 
        fp=open(sim.ivel)
        
        # Define a counter for the number of particles
        particle = 0

        # Loop around the lines in the file and assign the velocities
        while True:
            line=fp.readline()
            # break if there is a blank line or the end of file
            if not line: break
            xyz = line.split()
            if len(xyz) != 3:
                sys.exit("There is a problem with the velocities for " +
                         "atom " + str(particle+1) + " in \"" + sim.ivel + 
                         "\"\n")
            else:
                try:
                    atom[particle].vx = np.float64(xyz[0])
                    atom[particle].vy = np.float64(xyz[1])
                    atom[particle].vz = np.float64(xyz[2])
                except ValueError:
                    sys.exit("There is a problem with the velocities for " +
                             "atom " + str(particle+1) + " in \"" + 
                             sim.ivel + "\"\n")
                particle=particle+1
        fp.close()    
        # Check to see if the number of velocities read in from 
        # the velocity file is the same as the number specified in 
        # the input file.
        if particle != sim.N:
            sys.exit("The number of velocities (" + str(particle) + ") " + 
                     "in \"" + sim.ivel + "\" is not equal to the number " +
                     "of atoms " + "(" + str(sim.N) + ") in \"" + 
                     sim.inputfile + "\"\n")