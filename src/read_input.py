# read_input is part of ljpy for Lennard Jones simulations.                 #	 	                             #
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
# read_input.py                                                           	#
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
This module is part of ljpy. It has functions to read and process the input 
script with the parameters for the simulation such as number of particles,
temperature, density, etc. The input script is in the form of key words
followed by values. The script first reads the entire file, removes comments
and white space, and creates a python dictionary of keywords and values. Then, 
the values are assigned to an object which holds all the information for the 
simulation. Checks are performed to ensure all required key words are present
and the valid values for the keywords are given. 
"""
# Import relevant libraries
import numpy as np
import sys, os
from src.ljpyclasses import simulation


# This function is passed sys.argv from the main program
# It returns an object of type simulation with the values
# read from the input file.
def readinput(args):
    # check to see if the input files exists
    if not os.path.isfile(args[1]): 
        print("Input file \"" + args[1] +"\" does not exist.\n")
        sys.exit("Error: Input file missing.")
        
    # Parse the input file
    fi=open(args[1])        # open the file
    content=fi.readlines()  # read the file into a variable
    fi.close()              # close the file
    params={}               # make a dicitonary to hold the keywords and values
    
    
    for line in content:          # interate through each line of text in file
        linetext=line.strip()     # get rid of whitespace on each end of line
        if not linetext: continue # skip empty lines
        # Remove the end of line comment (everything after '#') and
        # and split the lines at all spaces
        linetext=linetext.split('#',1)[0].split()
        if not linetext: continue # skip a line that was only comments
        # Separate the line into the key and the value pair
        key=linetext[0]
        val=linetext[1:]
        # Place the key and value into the dictionary
        params[key]=val
    
    # Assign the input parameter values to the simulation object.
    # This requires changing some values from strings to numbers such as int,
    # double, float, etc. Also, each keyword is checked for proper type and 
    # the presenence of required keywords are checked. The program terminates
    # if any required keywords are missing or if invalid values are provided.
    
    sim=simulation()     # create the sim object that will hold the information
    
    # -------- sim keyword -------- #
    methodget=params.get('sim')
    if not methodget: sys.exit("The required keyword \"sim\" is missing in " +
                               "the input file.")
    if methodget[0] =='md' or methodget == 'mc': sim.method=params['sim'][0]
    else: sys.exit("The value of keyword \"sim\" in the input file must be " +
                   "either \"md\" or \"mc\".\n")
    
    # --------- N keyword --------- #    
    Nget=params.get('N')
    if not Nget: sys.exit("The required keyword \"N\" is missing in the " + 
                          "input file.")
    try:
        sim.N=np.ulonglong(params['N'][0])
    except ValueError:
        sys.exit("The value of required keyword \"N\" in the input file " +
                 "is not a valid integer.\n")
    
    # ------- temp keyword -------- #  
    Tget=params.get('temp')
    if not Tget: sys.exit("The required keyword \"temp\" is missing in the " + 
                          "input file.")
    try:
        sim.T=np.float(params['temp'][0])
    except ValueError:
        sys.exit("The value of required keyword \"temp\" in the input file " +
                 "is not a valid number.\n")
    
    # -------- rho keyword -------- # 
    rhoget=params.get('rho')
    if not rhoget: sys.exit("The required keyword \"rho\" is missing in the " + 
                            "input file.")
    try:
        sim.rho=np.float(params['rho'][0])
    except ValueError:
        sys.exit("The value of required keyword \"rho\" in the input file " +
                 "is not a valid number.\n")
    
    # ------- esteps keyword ------- # 
    eqget=params.get('esteps')
    if not eqget: sys.exit("The required keyword \"esteps\" is missing in " +
                           "the input file.")
    try:
        sim.eq=np.ulonglong(params['esteps'][0])
    except ValueError:
        sys.exit("The value of required keyword \"esteps\" in the input " +
                 "file is not a valid integer.\n")
    
    # ------- psteps keyword ------- # 
    prget=params.get('psteps')
    if not prget: sys.exit("The required keyword \"psteps\" is missing in " +
                           "the input file.")
    try:
        sim.pr=np.ulonglong(params['psteps'][0])
    except ValueError:
        sys.exit("The value of required keyword \"psteps\" in the input " + 
                 "file is not a valid integer.\n")
    
    # ------- rcut keyword -------- # 
    rcget=params.get('rcut')
    if not rcget: sys.exit("The required keyword \"rcut\" is missing in the " + 
                           "input file.")
    try:
        sim.rc=np.float(params['rcut'][0])
    except ValueError:
        sys.exit("The value of required keyword \"rcut\" in the input file " +
                 "is not a valid number.\n")
    
    sim.rc2=sim.rc*sim.rc # square of rcut
    
    # -------- dt keyword --------- #
    dtget=params.get('dt')
    if not dtget: sys.exit("The required keyword \"dt\" is missing in the " + 
                           "input file.")
    try:
        sim.dt=np.float(params['dt'][0])
    except ValueError:
        sys.exit("The value of required keyword \"dt\" in the input file " +
                 "is not a valid number.\n")
    
    # ------- coord keyword ------- #
    icoordget=params.get('coord')
    if icoordget: 
        sim.icoord=params['coord'][0]
    
    # -------- vel keyword -------- #
    ivelget=params.get('vel')
    if ivelget: 
        sim.ivel=params['vel'][0]
    
    # ------ output keyword ------- #
    outputget=params.get('output')
    if outputget: 
        try:
            sim.output=np.ulonglong(params['output'][0])
        except ValueError:
            sys.exit("The value of keyword \"output\" in the input " +
                     "file is not a valid integer greater than zero.\n")
        if sim.output == 0:
             sys.exit("The value of keyword \"output\" in the input " +
                     "file must be an integer greater than zero.\n")           
        
    # ------- seed keyword -------- #
    seedget=params.get('seed')
    if seedget: 
        if params['seed'][0] =='generate': sim.seedkeyvalue=params['seed'][0]
        else:
            try: 
                sim.seed=np.longlong(params['seed'][0])
                sim.seedkeyvalue="specified"
            except ValueError:
                sys.exit("The value of keyword \"seed\" in the " +
                         "input file is not a valid negative integer.\n")
        if sim.seed >=0:
            sys.exit("The value of keyword \"seed\" in the input file " +
                     "must be a negative integer.\n")
    
    # ------- movie keyword ------- #
    movieget=params.get('movie')
    if movieget: 
        try:  
            sim.movie=np.uint(params['movie'][1])
        except (ValueError, IndexError):
            sys.exit("The value for the interval of keyword \"movie\" in " +
                     "the input file is missing or is not a valid integer.\n")
        if sim.movie == 0:
            sys.exit("The value for the interval of keyword \"movie\" in " +
                     "the input file must be an integer greater than zero.\n")            
        try:
            sim.moviefile=params['movie'][0]
        except ValueError:
            sys.exit("The name for the movie file is either missing in the " +
                     "the input file or is not valid.\n")                

    # -------- rdf keyword -------- #
    rdfget=params.get('rdf')
    if rdfget:
        if len(params['rdf']) != 4: 
            sys.exit("The rdf keyword must be followed by four inputs:" +
                     "\n- the minimum " +
                     "value of r\n- the maximum value of r\n- the number of " +
                     "bins\n- the frequency at which to accumulate the " +
                     "histograms\n")
        try:  
            sim.rdfmin=np.float(params['rdf'][0])
            sim.rdfmax=np.float(params['rdf'][1])
            sim.rdfN=np.uint(params['rdf'][2])      
            sim.rdf=np.uint(params['rdf'][3])
        except ValueError:
            sys.exit("One or more of the parameters for the radial\n" +
                     "distribution function " +
                     "(rdf) in the input file are incorrect.\n" +
                     "The rdf keyword must be followed by four inputs:" +
                     "\n- the minimum " +
                     "value of r\n- the maximum value of r\n- the number of " +
                     "bins\n- the frequency at which to accumulate the " +
                     "histograms\n")
        if sim.rdfmin < 0.0:
            sys.exit("The minimum r value of the radial distribution\n"+
                     "function must be >= 0. (rdf Input File Error)")
        if sim.rdfmax < 0.0:
            sys.exit("The maximum r value of the radial distribution\n"+
                     "function must be >= 0. (rdf Input File Error)")
        if sim.rdfmax < sim.rdfmin:
            sys.exit("The maximum r value of the radial distribution\n"+
                     "function must be greater than the mininimum " +
                     "r value.\n(rdf Input File Error)")
        if sim.rdf < 1:
            sys.exit("The interval for keyword \"rdf\" must be "+
                     "an integer greater than zero.")            
            
    sim.inputfile=args[1]
    sim.outputfile=args[2]
    sim.length = np.double(sim.N/sim.rho)**(1.0/3.0)
    sim.utail = 8.0 / 3.0*np.pi*sim.rho*(1.0 / 3.0 * sim.rc**(-9.0) -
                                         sim.rc**(-3.0))
    sim.ptail = 16.0 / 3.0*np.pi*sim.rho*sim.rho*(2.0 / 3.0 * sim.rc**(-9.0) -
                                                  sim.rc**(-3.0))
    
    return sim





