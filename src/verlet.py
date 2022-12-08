# verlet is part of ljpy for Lennard Jones simulations.                     #
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
# verlet.py                                                                	#
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
This module is part of ljpy. It contains two functions to integrate the 
equations of motion using the velocity verlet algorithm.
"""
# Import relevant libraries
from numba import njit
# This function is passed a simulation object and a list of site objects.
# It is the first needed to use the velocity verlet algorithm. It uses
# the data at time step t to update the positions to the next time step and
# the velocities to the next half time step.
@njit
def verlet1(sim, atom):
    
    for i in range(sim.N):
        # Update the positions to a full time step.
        dx=sim.dt*atom[i].vx+sim.dt*sim.dt*atom[i].fx/2.0
        dy=sim.dt*atom[i].vy+sim.dt*sim.dt*atom[i].fy/2.0
        dz=sim.dt*atom[i].vz+sim.dt*sim.dt*atom[i].fz/2.0
        atom[i].x=atom[i].x+dx
        atom[i].y=atom[i].y+dy
        atom[i].z=atom[i].z+dz
        
        # Update the displacement accumulators for diffusivity
        atom[i].dx=atom[i].dx+dx
        atom[i].dy=atom[i].dy+dy
        atom[i].dz=atom[i].dz+dz
        
        # Apply periodic boundary conditions
        if atom[i].x < 0.0:          atom[i].x=atom[i].x+sim.length
        elif atom[i].x > sim.length: atom[i].x=atom[i].x-sim.length
        if atom[i].y < 0.0:          atom[i].y=atom[i].y+sim.length
        elif atom[i].y > sim.length: atom[i].y=atom[i].y-sim.length
        if atom[i].z < 0.0:          atom[i].z=atom[i].z+sim.length
        elif atom[i].z > sim.length: atom[i].z=atom[i].z-sim.length
        
        # Update the velocities to half a time step
        atom[i].vx=atom[i].vx+sim.dt*atom[i].fx/2.0
        atom[i].vy=atom[i].vy+sim.dt*atom[i].fy/2.0
        atom[i].vz=atom[i].vz+sim.dt*atom[i].fz/2.0
        
# This function is passed a simulation object and a list of site objects.
# It is the second function needed to use the velocity verlet       
# algorithm.  It updates the velocites from the half time step to the 
# full time step.
@njit
def verlet2(sim, atom):
    for i in range(sim.N):
        atom[i].vx=atom[i].vx+sim.dt*atom[i].fx/2.0
        atom[i].vy=atom[i].vy+sim.dt*atom[i].fy/2.0
        atom[i].vz=atom[i].vz+sim.dt*atom[i].fz/2.0       