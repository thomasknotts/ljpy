# finalize_file is part of ljpy for Lennard Jones simulations.              #
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
# finalize_file.py                                                         	#
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
This module is part of ljpy. It writes the final properties of the system
in the output file.
"""

# Import relevant libraries
from src.rdf import rdf_finalize

def finalizefile(sim, aprop, rdfh, rdfcalls, atomx,  atomy,  atomz,  \
                                             atomvx, atomvy, atomvz, \
                                             atomdx, atomdy, atomdz):
    # Variables
    pr=sim.pr
    N=sim.N
    
    # Calculate simple averages
    pe=aprop.pe/pr
    pe2=aprop.pe2/pr
    virial=aprop.virial/pr
    if sim.method == "md":
        ke=aprop.ke/pr
        T=aprop.T/pr
    else:
        # Divide by N here because each production step
        # is composed of N trials (1 MC sweep)
        pe/=N
        pe2/=N
        virial/=N
        T=sim.T
    
    # Calculate heat capacity and pressure
    P = sim.rho*T + 1.0/3.0/sim.length**3.0*virial + sim.ptail
    if sim.method == "mc": cv=(pe2 - pe*pe)/(T*T)/N + 3.0/2.0  # nvt expression
    else: cv = 3.0/2.0/(1 - 2.0/3.0*(pe2 - pe*pe)/N/(T*T))     # nve expression
   
    # Calculate the diffusivity from the MSD
    # This is zero for mc simulations
    Dmsd = 0.0
    for i in range(N):
        Dmsd+=atomdx[i]*atomdx[i] + atomdy[i]*atomdy[i] + \
              atomdz[i]*atomdz[i]
    Dmsd=Dmsd/pr/N/6.0/sim.dt
    
    # Write the data to file
    fp=open(sim.outputfile, "a")
    
    fp.write("\n    ***FINAL POSITIONS, XYZ Format***\n")
    fp.write(str(N) + "\nYou can copy these coordinates to a file to " +
             "open in a viewer.\n")
    for i in range(N):
        fp.write("C\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atomx[i], \
                                                            atomy[i], \
                                                            atomz[i]))

    if sim.method == "md":
        fp.write("\n         ***FINAL VELOCITIES***\n");
        for i in range(N):
            fp.write("\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atomvx[i], \
                                                               atomvy[i], \
                                                               atomvz[i]))

    if sim.rdf:
        fp.write("\n***Radial Distribution Function***\n\n");
        rdf_finalize(sim, rdfh, rdfcalls)
        rdfh.write(fp)


    if sim.pr > 0:
        fp.write("\n***Simulation Averages***\n\n")
        fp.write("Temperature:            {:10.6f}\n".format(T))
        fp.write("Pressure:               {:10.6f}\n".format(P))
        fp.write("Heat Capacity:          {:10.6f}\n".format(cv))
        fp.write("Potential Energy:       {:10.6f}\n".format(pe/N + sim.utail))
        if sim.method == "md":
            fp.write("Kinetic Energy:         {:10.6f}\n".format(ke/N))
            fp.write("Total Energy:           {:10.6f}\n".format((ke+pe) \
                                                                 /N+sim.utail))
            fp.write("Diffusivity             {:10.6f}\n".format(Dmsd))
        if sim.method == "mc":
            if aprop.ntry != 0:
                fp.write("MC Moves Accepted:      {:10.6f}\n" \
                         .format(aprop.naccept/aprop.ntry))
                fp.write("Final Max Displacment:  {:10.6f}\n" \
                         .format(sim.dt))
    else:
        fp.write("\nNo productions steps were specified, so simulation " +
                 "averages were not calculated.\n\n")
        
    fp.close()