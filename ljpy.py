# ljpy performs NVT MC or NVE MD simulations of a Lennard Jones fluid.	 	#
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
# along with ljmcmd.  If not, see <http://www.gnu.org/licenses/>.          	#

# ========================================================================= #
# ljpy.py 	                                                             	#
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
This program performs NVT Monte Carlo or NVE Molecular Dynamics simulations 
of a Lennard Jones fluid. It outputs instanteous property values and block
averages and errors according to input specifications. It can also output a
.trr file for movies and calculate the radial distribution function. It is 
written for pedagogical purposes and thus favors readability over speed. 
It is run with the following command.

  python ljpy.py <inputfilename> <outputfilename>

<inputfilename> is the name of the input file 
And example could be MCN500T85R9.input which would do and NVT MC simulations
of 500 particles at a dimensionless temperature of 0.85 and a dimensionless
density of 0.9.

<outputfilename> is the desired name of the output file
An example could be MCN500T85R9.output

The structure of the input file is a keyword followed by another entry or
multiple entries. An example input file is included with the program. Please
also refer to the included documentation.
"""

# Import relevant libraries
import numpy as np
import sys, os

# Create a class for each site in the system
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
        
# Create a class to hold the simulation information
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
        self.seedkeyvalue=-1    # key value for seed
        self.rdfmin=0.0         # minimum r value for rdf
        self.rdfmax=0.0         # maximum r value for rdf
        self.rdfN=0             # number of bins for rdf
        self.rdf=0              # frequency to accumulate the rdf
        
# Create a class to hold the simulation properties
class props:
    def __init__(self):
        self.ke=0.0             # kinetic energy
        self.pe=0.0             # potential energy
        self.pe2=0.0            # squared potential energy (for Cv)
        self.T=0.0              # temperature
        self.virial=0.0         # virial for pressure
        self.naccept=0          # number of mc moves accepted
        self.ntrys=0            # number of mc moves tried
        self.Nhist=0            # number of times accumulated
        
# ========================================================================= #
# Check the command line arguments for the input and output file names.     #
# ========================================================================= #


if len(sys.argv) != 3:
    print("This program requires two command line arguments: the path to",
          "the input file\nand the path to the output file.\n")
    sys.exit("Error: Input arguments missing.")

# ========================================================================= #
# Read the input file.                                                      #
# ========================================================================= #
       
# check to see if the input files exists
if not os.path.isfile(sys.argv[1]): 
    print("Input file \"" + sys.argv[1] +"\" does not exist.\n")
    sys.exit("Error: Input file missing.")
    
# Parse the input file
fi=open(sys.argv[1])    # open the file
content=fi.readlines()  # read the file into a variable
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

# Assign the input parameter values. This requires changing some values from
# strings to numbers such as int, double, float, etc. Also, for each 
# keyword is checked for proper type.
sim=simulation()     # create the sim object that will hold the information

methodget=params.get('sim') # sim keyword
if not methodget: sys.exit("The required keyword \"sim\" is missing in the " + 
                        "input file.")
if methodget[0] =='md' or methodget == 'mc': sim.method=params['sim'][0]
else: sys.exit("The value of keyword \"sim\" in the input file must be " +
               "either \"md\" or \"mc\".\n")
    
Nget=params.get('N') # N keyword
if not Nget: sys.exit("The required keyword \"N\" is missing in the " + 
                      "input file.")
try:
    sim.N=np.ulonglong(params['N'][0])
except ValueError:
    sys.exit("The value of required keyword \"N\" in the input file " +
             "is not a valid integer.\n")

Tget=params.get('temp') # temp keyword
if not Tget: sys.exit("The required keyword \"temp\" is missing in the " + 
                      "input file.")
try:
    sim.T=np.float(params['temp'][0])
except ValueError:
    sys.exit("The value of required keyword \"temp\" in the input file " +
             "is not a valid number.\n")

rhoget=params.get('rho') # rho keyword
if not rhoget: sys.exit("The required keyword \"rho\" is missing in the " + 
                        "input file.")
try:
    sim.rho=np.float(params['rho'][0])
except ValueError:
    sys.exit("The value of required keyword \"rho\" in the input file " +
             "is not a valid number.\n")

eqget=params.get('esteps') # esteps keyword
if not eqget: sys.exit("The required keyword \"esteps\" is missing in the " + 
                       "input file.")
try:
    sim.eq=np.ulonglong(params['esteps'][0])
except ValueError:
    sys.exit("The value of required keyword \"esteps\" in the input file " +
             "is not a valid integer.\n")

prget=params.get('psteps') # psteps keyword
if not prget: sys.exit("The required keyword \"psteps\" is missing in the " + 
                       "input file.")
try:
    sim.pr=np.ulonglong(params['psteps'][0])
except ValueError:
    sys.exit("The value of required keyword \"psteps\" in the input file " +
             "is not a valid integer.\n")
    
rcget=params.get('rcut') # rcut keyword
if not rcget: sys.exit("The required keyword \"rcut\" is missing in the " + 
                       "input file.")
try:
    sim.rc=np.float(params['rcut'][0])
except ValueError:
    sys.exit("The value of required keyword \"rcut\" in the input file " +
             "is not a valid number.\n")

sim.rc2=sim.rc*sim.rc # square of rcut

dtget=params.get('dt') # dt keyword
if not dtget: sys.exit("The required keyword \"dt\" is missing in the " + 
                       "input file.")
try:
    sim.dt=np.float(params['dt'][0])
except ValueError:
    sys.exit("The value of required keyword \"dt\" in the input file " +
             "is not a valid number.\n")

icoordget=params.get('coord') # coord keyword
if not icoordget: sys.exit("The required keyword \"coord\", or its value, " +
                           "is missing in the input file.")
if icoordget[0] =="" or not icoordget[0]:
    sys.exit("The value of keyword \"coord\" in the input file is " +
             "missing.\n")
else: sim.icoord=params['coord'][0]

ivelget=params.get('vel') # vel keyword
if not ivelget: sys.exit("The required keyword \"vel\", or its value, " +
                         "is missing in the input file.")
if ivelget[0] =="" or not ivelget[0]:
    sys.exit("The value of keyword \"vel\" in the input file is " +
             "missing.\n")
else: sim.ivel=params['vel'][0]

outputget=params.get('output') # output keyword
if not outputget: sys.exit("The required keyword \"output\" is missing " +
                           "in the input file.")
try:
    sim.output=np.ulonglong(params['output'][0])
except ValueError:
    sys.exit("The value of required keyword \"output\" in the input file " +
             "is not a valid integer.\n")
    
seedget=params.get('seed') # seed keyword
if not seedget: sys.exit("The required keyword \"seed\" is missing " +
                         "in the input file.")
if params['seed'][0] =='generate': sim.seedkeyvalue=params['seed'][0]
else:
    try: 
        sim.seed=np.longlong(params['seed'][0])
        sim.seedkeyvalue="specified"
    except ValueError:
        sys.exit("The value of required keyword \"seed\" in the input file " +
                 "is not a valid negative integer.\n")
if sim.seed >=0:
    sys.exit("The value of required keyword \"seed\" in the input file " +
             "must be a negative integer.\n")


sim.movie=np.uint(params['movie'][1])
sim.moviefile=params['movie'][0]
sim.inputfile=sys.argv[1]
sim.outputfile=sys.argv[2]
print("new = ", sim.output)
print(sim.moviefile)






print(sim.N)
    

    


    
atom=[]
atom.append(site())
print(atom[0].x)                       