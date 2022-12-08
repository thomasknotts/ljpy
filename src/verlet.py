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
# Version 2.0 - December 2022 Changed from atom class to arrays for numba. 	#
# ========================================================================= #

"""
This module is part of ljpy. It contains two functions to integrate the 
equations of motion using the velocity verlet algorithm.
"""

# This function is passed a simulation object.
# It is the first needed to use the velocity verlet algorithm. It uses
# the data at time step t to update the positions to the next time step and
# the velocities to the next half time step.
def verlet1(sim, atomx,  atomy,  atomz,  \
                 atomvx, atomvy, atomvz, \
                 atomfx, atomfy, atomfz, \
                 atomdx, atomdy, atomdz):
    
    for i in range(sim.N):
        # Update the positions to a full time step.
        dx=sim.dt*atomvx[i]+sim.dt*sim.dt*atomfx[i]/2.0
        dy=sim.dt*atomvy[i]+sim.dt*sim.dt*atomfy[i]/2.0
        dz=sim.dt*atomvz[i]+sim.dt*sim.dt*atomfz[i]/2.0
        atomx[i]=atomx[i]+dx
        atomy[i]=atomy[i]+dy
        atomz[i]=atomz[i]+dz
        
        # Update the displacement accumulators for diffusivity
        atomdx[i]=atomdx[i]+dx
        atomdy[i]=atomdy[i]+dy
        atomdz[i]=atomdz[i]+dz
        
        # Apply periodic boundary conditions
        if atomx[i] < 0.0:          atomx[i]=atomx[i]+sim.length
        elif atomx[i] > sim.length: atomx[i]=atomx[i]-sim.length
        if atomy[i] < 0.0:          atomy[i]=atomy[i]+sim.length
        elif atomy[i] > sim.length: atomy[i]=atomy[i]-sim.length
        if atomz[i] < 0.0:          atomz[i]=atomz[i]+sim.length
        elif atomz[i] > sim.length: atomz[i]=atomz[i]-sim.length
        
        # Update the velocities to half a time step
        atomvx[i]=atomvx[i]+sim.dt*atomfx[i]/2.0
        atomvy[i]=atomvy[i]+sim.dt*atomfy[i]/2.0
        atomvz[i]=atomvz[i]+sim.dt*atomfz[i]/2.0
        
# This function is passed a simulation object.
# It is the second function needed to use the velocity verlet       
# algorithm.  It updates the velocites from the half time step to the 
# full time step.
def verlet2(sim, atomvx, atomvy, atomvz, atomfx, atomfy, atomfz):
    for i in range(sim.N):
        atomvx[i]=atomvx[i]+sim.dt*atomfx[i]/2.0
        atomvy[i]=atomvy[i]+sim.dt*atomfy[i]/2.0
        atomvz[i]=atomvz[i]+sim.dt*atomfz[i]/2.0       