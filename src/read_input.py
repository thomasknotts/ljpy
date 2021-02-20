import numpy as np
import sys, os
from src.ljpyclasses import simulation

def readinput(args):
    # check to see if the input files exists
    if not os.path.isfile(args[1]): 
        print("Input file \"" + args[1] +"\" does not exist.\n")
        sys.exit("Error: Input file missing.")
        
    # Parse the input file
    fi=open(args[1])    # open the file
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
    sim.inputfile=args[1]
    sim.outputfile=args[2]

    return sim





