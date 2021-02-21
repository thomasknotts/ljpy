# ljpy performs NVT MC or NVE MD simulations of a Lennard Jones fluid.	 	#
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.    	#

# ========================================================================= #
# ljpy.py 	                                                             	#
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
This program performs NVT Monte Carlo or NVE Molecular Dynamics simulations 
of a Lennard Jones fluid. It outputs instanteous property values and block
averages and errors according to input specifications. It can also output a
.trr file for movies and calculate the radial distribution function. It is 
written for pedagogical purposes and thus favors readability over speed. 
It is run with the following command.

  python ljpy.py <inputfilename> <outputfilename>

<inputfilename> is the name of the input file 
And example could be MCN500T85R9.input which would do and NVT MC simulations
of 500 particles at a dimensionless temperature of 0.85 and a dimensionless
density of 0.9.

<outputfilename> is the desired name of the output file
An example could be MCN500T85R9.output

The structure of the input file is a keyword followed by another entry or
multiple entries. An example input file is included with the program. Please
also refer to the included documentation.
"""

# Import relevant libraries
#import numpy as np
import sys

#from ljpyclasses import simulation, site, props
#from src.ljpyclasses import site
from src.read_input import readinput
from src.initialize_positions import initializepositions
from src.initialize_files import initializefiles
        
# ========================================================================= #
# Check the command line arguments for the input and output file names.     #
# ========================================================================= #
if len(sys.argv) != 3:
    print("This program requires two command line arguments: the path to",
          "the input file\nand the path to the output file.\n")
    sys.exit("Error: Invalid or missing command line arguments.")

# ========================================================================= #
# Read the input file and store the simulation parameters to an object.     #
# ========================================================================= #
sim=readinput(sys.argv)       

# ========================================================================= #
# Initialize or read in positions.                                          #
# ========================================================================= #
atom=initializepositions(sim)

# ========================================================================= #
# Initialize the output files.                                              #
# ========================================================================= #
initializefiles(sim,atom)
    



print("sim.N = ", sim.N)
    

    


    

print(atom[17].x)                       