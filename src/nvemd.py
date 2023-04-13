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
from src.scale_velocities import scalevelocities
import numpy as np
import src.dhist as dh
from src.rdf import rdf_accumulate
from src.finalize_file import finalizefile

def nvemd(sim, atom):
    # Set variables
    rescale_freq=100
    
    # Create objects for the instanteous and average properties
    iprop=props()
    aprop=props()
    
    # Determine the initial properties (Iteration 0) and write to file.
    iprop.pe, iprop.virial, iprop.stress = forces(sim, atom)
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
        iprop.pe, iprop.virial, iprop.stress = forces(sim, atom) # calculate the forces
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
        
    # Reset the accumulators for the production steps
    aprop.pe=0.0
    aprop.ke=0.0
    aprop.T=0.0
    aprop.virial=0.0
    for i in range(sim.N):
        atom[i].dx=0.0
        atom[i].dy=0.0
        atom[i].dz=0.0
    
    # Initialize the radial distribution function histogram
    Nrdfcalls=0
    if sim.rdf:
        rdfh=dh.hist(sim.rdfmin, sim.rdfmax,sim.rdfN)
    else:
        rdfh=dh.hist(0.8, 4.0, 100) # this is the default

    # Initialize the pressure tensor correlation function histogram
    if sim.visc:
        ttxxh=dh.hist(0, sim.visc, np.int64(sim.visc/sim.dt))
        ttyyh=dh.hist(0, sim.visc, np.int64(sim.visc/sim.dt))
        ttzzh=dh.hist(0, sim.visc, np.int64(sim.visc/sim.dt))
        ttxyh=dh.hist(0, sim.visc, np.int64(sim.visc/sim.dt))
        ttxzh=dh.hist(0, sim.visc, np.int64(sim.visc/sim.dt))
        ttyzh=dh.hist(0, sim.visc, np.int64(sim.visc/sim.dt))

    # Perform the production steps
    # During production, accumulate all the properties.
    for i in range(1,np.int(sim.pr+1)):
        verlet1(sim, atom) # first half of velocity verlet algorithm
        iprop.pe, iprop.virial, iprop.stress = forces(sim, atom) # calculate the forces
        verlet2(sim, atom) # second half of velocity verlet algorithm
        iprop.ke, iprop.T = ke_and_T(atom) # kinetic and potential energy

        # Accumulate the properties
        aprop.pe+=iprop.pe
        aprop.ke+=iprop.ke
        aprop.T+=iprop.T
        aprop.virial+=iprop.virial
        aprop.pe2+=iprop.pe*iprop.pe
        aprop.stress+=iprop.stress
        
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
            print("Production Step " + str(i) + "\n")
            
        # Accumulate the radial distribution function
        if sim.rdf:
            if i%sim.rdf == 0:
                Nrdfcalls+=1
                rdf_accumulate(sim, atom, rdfh)
        
    # Finalize the output file
    finalizefile(sim, atom, aprop, rdfh, Nrdfcalls)

        
    #print("pe = %.2f virial = %.3f ke = %.2f T = %.4f" % (iprop.pe,iprop.virial,iprop.ke,iprop.T))