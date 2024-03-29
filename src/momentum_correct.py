# momentum_correct is part of ljpy for Lennard Jones simulations.           #
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
# momentum_correct.py                                                   	#
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
This module is part of ljpy. It contains two functions. The first checks to
see if the linear momentrum of the box is zero. Thd second adjusts the 
velocities to ensure that the linear momentum is zero.
"""

# This function is passed a list of site objects.
# It returns 0 if the linear momentum is zero (or close to zero)
# It returns 1 if the linear momentum is not zero
def checkmomentum(atom):
    # Determine the number of particles
    N=len(atom)
    
    # Zero out the momentum counters
    vcumx=0.0
    vcumy=0.0
    vcumz=0.0
    
    # Loop around the atoms in the system to determine the momentum
    # The dimensionless mass is equal to 1, so the momentum of each 
    # atom is equal to its velocity.
    for i in range(N):
        vcumx=vcumx+atom[i].vx
        vcumy=vcumy+atom[i].vy
        vcumz=vcumz+atom[i].vz
    if(vcumx+vcumy+vcumz)< 10.0**-10.0: return(0)
    else: return(1)


# This function is passed a list of site objectss.
# It attemps to zero out the linear momentum.
# It returns 0 if the linear momentum is zero (or close to zero)
# It returns 1 if the linear momentum is not zero    
def zeromomentum(atom):
    # Determine the number of particles
    N=len(atom)
    
    # Zero out the momentum counters
    vcumx=0.0
    vcumy=0.0
    vcumz=0.0
    
    # Loop around the atoms in the system to determine the momentum
    # The dimensionless mass is equal to 1, so the momentum of each 
    # atom is equal to its velocity.
    for i in range(N):
        vcumx=vcumx+atom[i].vx
        vcumy=vcumy+atom[i].vy
        vcumz=vcumz+atom[i].vz
    vcumx=vcumx/N
    vcumy=vcumy/N
    vcumz=vcumz/N
    
    for i in range(N):
        atom[i].vx=atom[i].vx-vcumx    
        atom[i].vy=atom[i].vy-vcumy
        atom[i].vz=atom[i].vz-vcumz
        
    return(checkmomentum(atom))
    