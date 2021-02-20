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