# initialize_files is part of ljpy for Lennard Jones simulations.           #
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
# initialize_files.py                                                    	#
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
This module is part of ljpy. It writes the initial parts of the output file
such as the simulations parameters, the initial positions, and the
initial velocities 
"""

def initializefiles(sim,atom):
    
    # Initialize the movie files.
    if sim.movie:
        # Initilaize the .trr file.
        fi=open(sim.moviefile,"w")
        fi.close()
        
        # Write the .xyz file needed for loading the .trr file.
        fi=open(sim.moviefile.split(".")[0]+".xyz","w")
        fi.write(str(sim.N)+"\nLoad this file in VMD before the .trr file\n")
        for i in range(sim.N):
            fi.write("C\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atom[i].x, \
                                                                atom[i].y, \
                                                                atom[i].z))
        fi.close()
        
    # Write the header output file including the simulation parameters,
    # the initial, and the initial velocities.
    fi=open(sim.outputfile,"w")
    fi.write(str(sim.method).upper() + " simulation of " + str(sim.N) + 
             " LJ Particles at T*={:.4f}".format(sim.T) + " and " +
             "rho*={:.4f}\n\n".format(sim.rho))
    fi.write("Input File:         " + sim.inputfile + "\n")
    fi.write("Output File:        " + sim.outputfile + "\n\n")

    if sim.seedkeyvalue == "generate":
        fi.write("The random number seed was generated from the system " +
                 "clock.\n")
    elif sim.seedkeyvalue == "specified":
        fi.write("The random number seed was specified in the input file.\n")
    else:
        fi.write("The seed for the random number generator was not " +
                 "specified. The default value was used.\n");
    fi.write("Random Number Seed: " + str(sim.seed) + "\n\n")
    fi.write("\n    ***Input Parameters***\n");
    fi.write("sim         " + sim.method + "\n")
    fi.write("N           " + str(sim.N) + "\n")
    fi.write("temp        " + str(sim.T) + "\n")
    fi.write("rho         " + str(sim.rho) + "\n")
    fi.write("esteps      " + str(sim.eq) + "\n")
    fi.write("psteps      " + str(sim.pr) + "\n")
    fi.write("rcut        " + str(sim.rc) + "\n")
    fi.write("dt          " + str(sim.dt) + "\n")
    if sim.output != 0:
        fi.write("output      " + str(sim.output) + "\n")
    if sim.seedkeyvalue != None:
        if sim.seedkeyvalue == "generate":
            fi.write("seed        generate\n")
        else:
            fi.write("seed        " + str(sim.seed) + "\n")
    if sim.icoord != None:
        fi.write("coord       " + sim.icoord + "\n")
    if sim.ivel != None:
        fi.write("vel         " + sim.ivel + "\n")
#    if sim.moviefile != None:
#        fi.write("movie       " + sim.moviefile + "  " + str(sim.movie) + "\n")
    if sim.rdf != 0:
        fi.write("rdf         " + str(sim.rdfmin) + "  " +
                 str(sim.rdfmax) + "  " + str(sim.rdfN) +"  " +
                 str(sim.rdf) + "\n")
    fi.write("\n")

    fi.write("Energy Tail Correction:    {:.8f}\n".format(sim.utail))
    fi.write("Pressure Tail Correction:  {:.8f}\n".format(sim.ptail))

    fi.write("\n    ***INITIAL POSITIONS, XYZ Format***\n")
    fi.write(str(sim.N) + "\nYou can copy these coordinates to a file to " +
             "open in a viewer.\n")
    for i in range(sim.N):
        fi.write("C\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atom[i].x, \
                                                            atom[i].y, \
                                                            atom[i].z))

    if sim.method == "md":
        fi.write("\n         ***INITIAL VELOCITIES***\n");
        for i in range(sim.N):
            fi.write("\t{:13.6f}\t{:13.6f}\t{:13.6f}\n".format(atom[i].vx, \
                                                               atom[i].vy, \
                                                               atom[i].vz))
        fi.write("\n\nIteration                T              T Ave.       " +
                 "       P              P Ave.            KE               " +
                 "PE               TE\n\n")
    else: fi.write("\n\nIteration                P              P Ave. " +
                       "             PE\n\n") 
    fi.close()
    
    
    
    
    
    
    