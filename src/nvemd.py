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
from pickletools import TAKEN_FROM_ARGUMENT1
from src.forces import forces
from src.kinetic import ke_and_T
from src.verlet import verlet1, verlet2
from src.ljpyclasses import props
from src.scale_velocities import scalevelocities
import numpy as np
import src.dhist as dh
from src.rdf import rdf_accumulate
from src.finalize_file import finalizefile
from numba.typed import List
from src.visc import ptensor, tautaucorr_accumulate, traceless, tautaucorr_finalize

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
    # tautaucorr[0:5] hold the xx, yy, zz, xy, xz, yz components
    # tautaucorr[6] holds the number of times each taucorr[0:5] was incremented
    # so that an average may be determined. The histogram bin values are 
    # the time.
    if sim.visc:
        tautaucorr=List()
        time0=0.0
        tau0=ptensor(sim,atom,iprop)
        tau0traceless=traceless(tau0)
        for i in range(7):
            tautaucorr.append(dh.hist(0, sim.visc, np.int64(sim.visc/sim.dt)))
    
    # Perform the production steps
    # During production, accumulate all the properties.
    for i in range(1,np.int(sim.pr+1)):
        verlet1(sim, atom) # first half of velocity verlet algorithm
        iprop.pe, iprop.virial, iprop.stress = forces(sim, atom) # calculate the forces
        verlet2(sim, atom) # second half of velocity verlet algorithm
        iprop.ke, iprop.T = ke_and_T(atom) # kinetic and potential energy
        iprop.press = ptensor(sim,atom,iprop)
        # Accumulate the properties
        aprop.pe+=iprop.pe
        aprop.ke+=iprop.ke
        aprop.T+=iprop.T
        aprop.virial+=iprop.virial
        aprop.pe2+=iprop.pe*iprop.pe
        aprop.stress+=iprop.stress
        aprop.press+=iprop.press
        
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
        
        # Accumulate the pressure tensor correlation function.
        # First calculate 
        if sim.visc:
            time = i*sim.dt - time0
            tau = ptensor(sim,atom,iprop)
            tau=traceless(tau) # make the tensor traceless
            tautau0 = tau0*tau
            tautaucorr_accumulate(time, tautau0, tautaucorr)
            if time==sim.visc: # this resets the zero in time
                time0=time
                tau0=tau

        
    # Finalize the output file
    finalizefile(sim, atom, aprop, rdfh, Nrdfcalls)

    # Write the tau*tau0 correlation file
    for i in range(6): # get the average tautaucorr across the simulation
        tautaucorr[i]=dh.div(tautaucorr[i], tautaucorr[6]) 
        finaltautau=tautaucorr_finalize(tautaucorr)
        fn=sim.outputfile.split('.')[0] + "tautau" + str(i) + "." + sim.outputfile.split('.')[1]
        fp=open(fn,"w")
        for j in range(tautaucorr[i].N): 
            fp.write("{:>10}\t{:13.6f}\t{:13.6f}\n".format(j+1, \
                     tautaucorr[i].range[j], tautaucorr[i].bin[j]))
        fp.close()
    fn=sim.outputfile.split('.')[0] + "tautauave" + "." + sim.outputfile.split('.')[1]
    fp=open(fn,"w")
    for j in range(tautaucorr[0].N): 
        fp.write("{:>10}\t{:13.6f}\t{:13.6f}\n".format(j+1, \
                  tautaucorr[0].range[j], finaltautau[j]))
    fp.close()

        
    #print("pe = %.2f virial = %.3f ke = %.2f T = %.4f" % (iprop.pe,iprop.virial,iprop.ke,iprop.T))