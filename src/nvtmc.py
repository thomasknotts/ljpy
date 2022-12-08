# nvtmc is part of ljpy for Lennard Jones simulations.                      #
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
# nvtmc.py                                                      	        #
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
This subroutine is the main driver for NVT MC simulations.  
"""

# Import relevant libraries
from src.atomic_pe import atomic_pe
from src.move import move
from src.forces import forces
from src.ljpyclasses import props
from src.scale_delta import scale_delta
from src.rdf import rdf_accumulate
from src.finalize_file import finalizefile
import src.dhist as dh
import numpy as np

def nvtmc(sim, atomx, atomy, atomz): 
    # Variables
    atompe=np.zeros(sim.N,np.float64)
    freq_scale_delta=1 # frequency to scale the maximum displacement
    
    # Create objects for the instanteous and average properties
    iprop=props()
    aprop=props()
    
    # Determine the initial properties (Iteration 0) and write to file.
    iprop.pe, iprop.virial, atomfx, atomfy, atomfz = forces(sim, atomx, atomy, atomz) # calculate the forces
    P=sim.rho*sim.T + 1.0/3.0/sim.length**3.0*iprop.virial + sim.ptail
    fp=open(sim.outputfile, "a")
    fp.write("{:<13}    {:13.6f}    {:13.6f}    {:13.6f}\n" \
             .format(0, P, P, iprop.pe/sim.N + sim.utail))
    fp.close()
    
    # Initialize the potential energy of each site
    # These are the "old" or "current" energies needed
    # to calculate the change in energy between the current state
    # and a proposed state (move).
    for i in range(sim.N):
        atompe[i]=atomic_pe(sim, atomx, atomy, atomz, i)
        
    
    # Perform the equilibration steps
    # Each step proposes sim.N moves (one Monte Carlo "sweep").
    for i in range(1, np.int(sim.eq+1)):
        for j in range(sim.N): # This loop performs sim.N moves per step
            # Propose and accept or reject a move
            move(sim, atomx, atomy, atomz, atompe, iprop)
            
            # Accumulate the properties for the move
            aprop.pe+=iprop.pe
            aprop.pe2+=iprop.pe2
            aprop.virial+=iprop.virial
        
        # Output equilibration progress at the interval specified
        # in the input file
        if i%sim.output == 0:
            P=sim.rho*sim.T + 1.0/3.0/sim.length**3.0*iprop.virial + sim.ptail
            Pave=sim.rho*sim.T + \
                 1.0/3.0/sim.length**3.0*aprop.virial/i/sim.N + sim.ptail
            fp=open(sim.outputfile, "a")
            fp.write("{:<13}    {:13.6f}    {:13.6f}    {:13.6f}\n" \
                     .format(i, P, Pave, iprop.pe/sim.N + sim.utail))
            fp.close()
            print("Equilibration Step " + str(i) + "\n")
        
        # Scale delta to obtain desired acceptance of moves
        if i%freq_scale_delta == 0: scale_delta(sim,iprop,aprop)
        
    # Reset accumulators for production steps
    iprop.ntry=0
    iprop.naccept=0
    aprop.ntry=0
    aprop.naccept=0
    aprop.pe=0.0
    aprop.pe2=0.0
    aprop.virial=0.0
        
    # Initialize the radial distribution function histogram
    if sim.rdf:
        rdfh=dh.hist(sim.rdfmin, sim.rdfmax,sim.rdfN)
        Nrdfcalls=0
    
    # Perform the production steps
    # During production, accumulate all the properties.    
    for i in range(1, np.int(sim.pr+1)):
        for j in range(sim.N): # This loop performs sim.N moves per step
            # Propose and accept or reject a move
            move(sim, atomx, atomy, atomz, atompe, iprop)
            
            # Accumulate the properties for the move
            aprop.pe+=iprop.pe
            aprop.pe2+=iprop.pe2
            aprop.virial+=iprop.virial
        
        # Accumulate the radial distribution function
        if sim.rdf:
            if i%sim.rdf == 0:
                Nrdfcalls+=1
                rdf_accumulate(sim, rdfh, atomx, atomy, atomz)
       
        # Output production progress at the interval specified
        # in the input file
        if i%sim.output == 0:
            P=sim.rho*sim.T + 1.0/3.0/sim.length**3.0*iprop.virial + sim.ptail
            Pave=sim.rho*sim.T + \
                 1.0/3.0/sim.length**3.0*aprop.virial/i/sim.N + sim.ptail
            fp=open(sim.outputfile, "a")
            fp.write("{:<13}    {:13.6f}    {:13.6f}    {:13.6f}\n" \
                     .format(i, P, Pave, iprop.pe/sim.N + sim.utail))
            fp.close()
            print("Production Step " + str(i) + "\n")
        
        # Scale delta to obtain desired acceptance of moves
        if i%100 == 0: scale_delta(sim,iprop,aprop)
        
    # Finalize the output file after all equilibration and production    
    # steps are finished.  This calculates and write the averages to the 
    # output file.
    dummy=np.zeros(sim.N) # need dummay arrays for next function
    finalizefile(sim, aprop, rdfh, Nrdfcalls, atomx, atomy, atomz,  \
                                              dummy, dummy, dummy, \
                                              dummy, dummy, dummy)    
    
    