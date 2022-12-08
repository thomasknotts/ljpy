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
# Version 2.0 - December 2022 Changed from atom class to arrays for numba. 	#
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

def nvemd(sim, atomx,  atomy,  atomz, \
               atomvx, atomvy, atomvz):

    # initialize arrays for displacements
    atomdx=np.zeros(sim.N)    # x displacement for diffusion
    atomdy=np.zeros(sim.N)    # y displacement for diffusion
    atomdz=np.zeros(sim.N)    # z displacement for diffusion
    #atomdr2=np.zeros(sim.N)   # MSD accumulator for diffusion

    # Set variables
    rescale_freq=100
    
    # Create objects for the instanteous and average properties
    iprop=props()
    aprop=props()
    
    # Determine the initial properties (Iteration 0) and write to file.
    iprop.pe, iprop.virial, atomfx, atomfy, atomfz = forces(sim, atomx, atomy, atomz)
    iprop.ke, iprop.T = ke_and_T(atomvx, atomvy, atomvz)
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
        verlet1(sim, atomx,  atomy,  atomz,  \
                     atomvx, atomvy, atomvz, \
                     atomfx, atomfy, atomfz, \
                     atomdx, atomdy, atomdz) # first half of velocity verlet algorithm
        iprop.pe, iprop.virial, atomfx, atomfy, atomfz = forces(sim, atomx, atomy, atomz) # calculate the forces
        verlet2(sim, atomvx, atomvy, atomvz, atomfx, atomfy, atomfz) # second half of velocity verlet algorithm
        iprop.ke, iprop.T = ke_and_T(atomvx, atomvy, atomvz) # kinetic and potential energy
        
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
        if i%rescale_freq == 0: scalevelocities(sim, atomvx, atomvy, atomvz, aprop.T/i)
        
    # Reset the accumulators for the production steps
    aprop.pe=0.0
    aprop.ke=0.0
    aprop.T=0.0
    aprop.virial=0.0
    atomdx=np.zeros(sim.N)
    atomdy=np.zeros(sim.N)
    atomdz=np.zeros(sim.N)
    #for i in range(sim.N):
    #    atomdx[i]=0.0
    #    atomdy[i]=0.0
    #    atomdz[i]=0.0
    
    # Initialize the radial distribution function histogram
    if sim.rdf:
        rdfh=dh.hist(sim.rdfmin, sim.rdfmax,sim.rdfN)
        Nrdfcalls=0
        
    # Perform the production steps
    # During production, accumulate all the properties.
    for i in range(1,np.int(sim.pr+1)):
        verlet1(sim, atomx,  atomy,  atomz,  \
                     atomvx, atomvy, atomvz, \
                     atomfx, atomfy, atomfz, \
                     atomdx, atomdy, atomdz) # first half of velocity verlet algorithm
        iprop.pe, iprop.virial, atomfx, atomfy, atomfz = forces(sim, atomx, atomy, atomz) # calculate the forces
        verlet2(sim, atomvx, atomvy, atomvz, atomfx, atomfy, atomfz) # second half of velocity verlet algorithm
        iprop.ke, iprop.T = ke_and_T(atomvx, atomvy, atomvz) # kinetic and potential energy

        
        # Accumulate the properties
        aprop.pe+=iprop.pe
        aprop.ke+=iprop.ke
        aprop.T+=iprop.T
        aprop.virial+=iprop.virial
        aprop.pe2+=iprop.pe*iprop.pe
        
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
                rdf_accumulate(sim, rdfh, atomx, atomy, atomz)
        
    # Finalize the output file
    finalizefile(sim, aprop, rdfh, Nrdfcalls, atomx,  atomy,  atomz,  \
                                              atomvx, atomvy, atomvz, \
                                              atomdx, atomdy, atomdz)

    print("pe =",iprop.pe,", virial =",iprop.virial,", ke =",iprop.ke,", T =",iprop.T)