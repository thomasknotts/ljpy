# nvemd is part of ljpy for Lennard Jones simulations.                      #
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
# nvemd.py                                                                	#
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
This module is part of ljpy. It is the main driver for the NVE MD simulations.
"""

# Import relevant libraries
from src.forces import forces
from src.kinetic import ke_and_T
from src.verlet import verlet1, verlet2
from src.ljpyclasses import props
import numpy as np

def nvemd(sim, atom):
    # Set variables
    rescale_freq=100
    
    # Create objects for the instanteous and average properties
    iprop=props()
    aprop=props()
    
    # Determine the initial properties (Iteration 0) and write to file.
    iprop.pe, iprop.virial = forces(sim, atom)
    iprop.ke, iprop.T = ke_and_T(atom)
    P=sim.rho*iprop.T + 1.0/3.0/sim.length**3.0*iprop.virial + sim.ptail
    fp=open(sim.outputfile, "a")
    fp.write("{:<13}    {:13.6f}    {:13.6f}    {:13.6f}    {:13.6f}    " \
             "{:13.6f}    {:13.6f}    {:13.6f}\n" \
             .format(0, iprop.T, iprop.T, P, P, iprop.ke/sim.N, \
                     iprop.pe/sim.N + sim.utail, \
                     (iprop.ke + iprop.pe)/sim.N + sim.utail))
    fp.close()

    # Perform equilibration steps
    # During equilibration, the velocities are rescaled periodically
    # to the set point temperature. After equilibration, during production,
    # the velocities are no longer rescaled.
    for i in range(1,np.int(sim.eq+1)):
        verlet1(sim, atom) # first half of velocity verlet algorithm
        iprop.pe, iprop.virial = forces(sim, atom) # calculate the forces
        verlet2(sim, atom) # second half of velocity verlet algorithm
        iprop.ke, iprop.T = ke_and_T(atom) # kinetic and potential energy
        
        # Accumulate the properties
        aprop.pe=aprop.pe + iprop.pe
        aprop.ke=aprop.ke + iprop.ke
        aprop.T=aprop.T + iprop.T
        aprop.virial=aprop.virial + iprop.virial
        
        # Output instantaneous properties at the interval
        # specified in the input file.
        if i%sim.output == 0:
            P=sim.rho*iprop.T + 1.0/3.0/sim.length**3.0*iprop.virial + \
              sim.ptail
            Pave=sim.rho*aprop.T/i + 1.0/3.0/sim.length**3.0*aprop.virial/i + \
              sim.ptail 
            fp=open(sim.outputfile, "a")
            fp.write("{:<13}    {:13.6f}    {:13.6f}    {:13.6f}    {:13.6f}    " \
                     "{:13.6f}    {:13.6f}    {:13.6f}\n" \
                     .format(i, iprop.T, aprop.T/i, P, Pave, iprop.ke/sim.N, \
                             iprop.pe/sim.N + sim.utail, \
                             (iprop.ke + iprop.pe)/sim.N + sim.utail))
            fp.close()
            print("Equilibration Step " + str(i) + "\n")
            
        # Rescale the velocities to achieve the temperature specified
        # in the input file. This is only done during the equilibration
        # steps of MD simulations.
        if i%rescale_freq == 0: scalevelocities(sim, atom, aprop.T/i)
    
    print("pe =",iprop.pe,", virial =",iprop.virial,", ke =",iprop.ke,", T =",iprop.T)