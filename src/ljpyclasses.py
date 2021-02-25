# ljpyclasses is part of ljpy for Lennard Jones simulations.                #
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
# ljpyclasses.py                                                         	#
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
This module is part of ljpy. It defines the classes needed for simulations. 
Three classes are defined.

    site:       a class that holds the information for each site (atom) in the 
                system such as x,y,z, position, x,y,z velocities, etc.
    simualtion: a class that holds the information for the simulation as read
                from the input file such as number of particles, temperature
                density, etc.
    props:      a class that holds the properties of the simulation such as 
                potential energy, kinetic energy, pressure, etc.
"""
# The class for each site in the system
class site:
    def __init__(self):
        self.x=0.0      # x position
        self.y=0.0      # y position
        self.z=0.0      # z position
        self.vx=0.0     # x velocity
        self.vy=0.0     # y velocity
        self.vz=0.0     # z velocity
        self.fx=0.0     # x force
        self.fy=0.0     # y force
        self.fz=0.0     # z force
        self.dx=0.0     # x displacement for diffusion (MD)
        self.dy=0.0     # y displacement for diffusion (MD)
        self.dz=0.0     # z displayement for diffusion (MD)   
        self.dr2=0.0    # MSD accumulator for diffusion (MD)
        self.pe=0.0     # potential energy of site (MC)
        
# The class to hold the simulation information
class simulation:
    def __init__(self):
        self.method = None      # simulation method (md or mc)
        self.T = 0.0            # temperature [T*]
        self.rho = 0.0          # density [rho*]
        self.N=0                # number of particles
        self.eq=0               # equilibration steps
        self.pr=0               # production steps
        self.itrr=0             # interval for saving trajectories
        self.rc=0.0             # cutoff radius [r*]
        self.rc2=0.0            # square of cutoff radius [r*]
        self.dt=0.0             # time step for md or max delta for mc
        self.icoord=None        # filename for initial coordinates
        self.ivel=None          # filename for initial velocities
        self.inputfile=None     # name of input file
        self.outputfile=None    # name of output file
        self.length=0.0         # length of simulation blox
        self.output=0           # interval for ouput of instan. props.
        self.movie=0            # interval for movie frames
        self.moviefile=None     # name of movie file
        self.utail=0.0          # tail correction to energy
        self.ptail=0.0          # tail correction to pressure
        self.seed=-1            # seed to the random number generator
        self.seedkeyvalue=None  # key value for seed
        self.rdfmin=0.0         # minimum r value for rdf
        self.rdfmax=0.0         # maximum r value for rdf
        self.rdfN=0             # number of bins for rdf
        self.rdf=0              # frequency to accumulate the rdf
        
# The class to hold the simulation properties
class props:
    def __init__(self):
        self.ke=0.0             # kinetic energy
        self.pe=0.0             # potential energy
        self.pe2=0.0            # squared potential energy (for Cv)
        self.T=0.0              # temperature
        self.virial=0.0         # virial for pressure
        self.naccept=0          # number of mc moves accepted
        self.ntry=0             # number of mc moves tried
        self.Nhist=0            # number of times accumulated