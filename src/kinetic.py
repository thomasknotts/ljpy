# kinetic is part of ljpy for Lennard Jones simulations.                    #	 	                             #
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
# kinetic.py                                                            	#
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
This module is part of ljpy. It has three functions. One returns the
kinetic energy of the system, one returns the temperature of the system, and
one returns both the kinetic energy and the temperature.
"""

# This function is passed a list of site objects.
# It returns the kinetic energy of the system.
def kinetic_energy(atom):
    # Determine the number of particles
    N=len(atom)
    
    # Zero out the accumulator for the energy
    ke=0.0
    # Loop around the farticls to calculate the kinetic energy
    for i in range(N):
        v2=atom[i].vx*atom[i].vx + atom[i].vy*atom[i].vy + \
           atom[i].vz*atom[i].vz
        ke=ke+0.5*v2
    
    return(ke)

# This function is passed the a list of site objects.
# It returns the temperature of the system.
def temperature(atom):
    # Determine the number of particles
    N=len(atom)
    # Calculate the temperature from the kinetic energy
    T=2.0/3.0/N*kinetic_energy(atom)
    return(T)

# This function is passed the a list of site objects.
# It returns a tuple with the kinetic energy and the temperature of
# the system.
def ke_and_T(atom):
    # Determine the number of particles
    N=len(atom)
    # Calculate the kinetic energy of the system.
    ke=kinetic_energy
    # Calculate the temperature from the kinetic energy
    T=2.0/3.0/N*ke
    return(ke,T)



