import numpy as np
import sys
import shutil

from fargo_setup import ShearingBox, Output, FargoSetup
from polydust import Polydust, SizeDistribution
import run_script

# Name of the setup. Any existing setup with the same name will be
# overwritten!
setup_name = 'DVSI_gas_2D_Z05_eta05_a5_tse3_hg4x5'
setup_name = 'test'
folder = 'DVSI/'
folder = 'IRS48/'


sys.path.insert(0, '/Users/jmatthijsse/Documents/PhD_Project2/createsetup/psi_stratfield')

# FARGO directory: /path/to/fargo/public
# If None, it is assumed that FARGO was installed as a submodule.
##############################################################################
fargo_dir = '/Users/jmatthijsse/Documents/PhD_project2/fargo3d'
output_direc = '@outputs/'+setup_name

# If we want to copy the base FARGO directory, new folder to run on DelftBlue
##############################################################################
# dir_fargo_base = '/Users/jmatthijsse/Documents/PhD_project2/delftblue/fargo3d'
# fargo_dir =  '/Users/jmatthijsse/Documents/PhD_project2/delftblue/'+folder+setup_name+ '/fargo3d'
# shutil.copytree(dir_fargo_base, fargo_dir, dirs_exist_ok=False)
# # add run_script to path
# run_script.write_runscript(setup_name, time=[119, 59, 59], n_tasks=128)
# shutil.copy('run.sub', fargo_dir + '/run.sub')
# output_direc = './'

##############################################################################
# Define how output: how much and when
#
# Follow the FARGO terminology: dt is the amount of time between lightweight
# output, and Ninterm is the amount of lightweight outputs between full data
# dumps. Therefore, the code will perform a dump every dt*Ninterm. The total
# simulation time is dt*Ntot.
##############################################################################

# Data output directory, to be put in FARGO par file

# Output object: how many dumps, how much total time
output = Output(output_direc, dt=0.1, Ninterm=10, Ntot=5000)

##############################################################################
# Gas and dust properties
##############################################################################


# Gas density at the midplain
gas_density = 1.0

# Total metalicity
dust_density = 0.05
# dust_density = 0.0425

##############################################################################
# Define the Stokes number distribution
#
# We can either have a discrete number of Stokes numbers, or a continuous
# distribution between a minimum and a maximum Stokes number.
###########################################################################

# Number of dust species. In the case of a continuous size distribution, this
# is the number of Gauss-Legendre nodes
n_dust = 1

# For discrete multifluid: list all Stokes numbers
stokes_range = [1e-3]
# stokes_range = [3.5e-1]

# For continuum: min and max Stokes numbers
# stokes_range = [1e-3, 0.1]


# Discrete sizes: list of dust densities
size_dist = [1.0]
# Continuous size distribution
# size_dist = SizeDistribution(stokes_range)


##############################################################################
# Define the shearing box
##############################################################################
Hg = 1.0

Lx = Hg*0.5           # 'radial' box size
Lz = Hg*4.0              # vertical box size default:4Hg

Nz = int(2*Lz*1024)      # vertical number of grid points default: 2048 pixels/Hg
Nz = int(2*Lz*1024/8)      # vertical number of grid points default: 2048 pixels/Hg

Nx = int(Nz*Lx/Lz)       # vertical number of grid points

# If 1D box, set Nx = 1, Lx = 2*Nx*Lz/Nz
ONED = False
if ONED:
    Nx = 1
    Lx = 2*Nx*Lz/Nz

viscous_alpha = 0.0

# Create polydust object
pd = Polydust(n_dust=n_dust,
              stokes_range=stokes_range,
              dust_density=dust_density,
              gas_density=gas_density,
              alpha_viscosity_gas=viscous_alpha,
              height_box=Lz,
              size_distribution=size_dist,
              gauss_legendre=True,
              discrete_equilibrium=False,
              pressure_gradient=True)  #if pressure gradient is True, eta=1.0 else eta=0.0

pd.alpha_viscosity_dust = 1e-5# Set to None if no dust viscosity is used
pd.alpha_viscosity_dust = 1e-5# Set to None if no dust viscosity is used

shearing_box = ShearingBox(dims=[0, Lx, Lz], mesh_size=[1, Nx, Nz])

# Shear parameter, set to None for no shear
shear =  0.05 # default in Lehmann & Lin 2023 is 0.05:

################################################################################
# Define the perturbation
##############################################################################
f_wn = 1e-3
kx = np.pi/0.125 # radial wavenumber
kz = np.pi/2     # vertical wavenumber
kx_min = np.pi/(2*Lx)
kz_min = np.pi/(2*Lz)
kx_max = np.pi*Nx/(2*Lx) * (1/16) #cutt-off the smalllest wavelengths
kz_max = np.pi*Nz/(2*Lz) * (1/16) #cutt-off the smalllest wavelengths

# if not ONED:
#     if kx < kz/(2*shear):
#         raise Exception('kx must be larger than kz/(2*shear) for a VSI setup!')  # check if kx is large enough for VSI
#     if 2*np.pi/kx <= 16*2*Lx/Nx:
#         raise Exception('wave is not resolved in radial direction!')  # check if wave is resolved in radial direction
#     if 2*np.pi/kz <= 16*2*Lz/Nz:
#         raise Exception('wave is not resolved in vertical direction!')  # check if wave is resolved in vertical direction

# Single mode perturbation, set to None for no perturbation and float for white noise perturbation
single_mode = [f_wn, kx, kz]
multiple_mode = [f_wn, [kx_min, kx_max], [kz_min, kz_max], 256] # amplitude, [kx_min, kx_max], [kz_min, kz_max], number of modes
if type(single_mode) == float:
    print('White noise perturbation: amplitude = {}'.format(single_mode))
else:
    print('Single mode perturbation: amplitude = {}, kx = {}, kz = {}'.format(single_mode[0], single_mode[1], single_mode[2]))

# If you want to use a tracer fluid, 
# tracer = [Stokes tracer, Mass fraction (m_tracer/m_gas), H_tracer] default: [1e-6, 1e-4, 1e-2]
# If you do not want to use a tracer, set tracer = None
# tracer = [[1e-6, 1e-2, 1e-2], [1e-2, 1e-2, 1e-2]]
# tracer = [[3.73e-4, 7.5e-3, 1e-1], [3.5e-1, 0.045, 1e-2]]



##############################################################################
# Write the FARGO setup files
##############################################################################

try:
    setup = FargoSetup(setup_name, fargo_dir=fargo_dir)
except:
    raise

# If aceleration_gas is True, the gas will feel an accelration:
# either due to the velocity difference between the gas and dust in the equilibrium state,
# or to enforce the shear factor which leads an azimuthal velocity gradient qz.

# If dust=True, the dust will be included in the setup.
# If dust=False, the dust will not be included in the setup and dust related parameters
# will oly be used to calculate the acceleration field on the gas.
'''
Arguments create setup   
polydust: PolyDust object
shearing_box: ShearingBox object
output: Output object
cfl (optional): Courant number to use. Defaults to 0.44.
pertubrbation (optional): white noise or single mode pertibution. str for WN or [amp, Kx, Kz] for single mode. Defaults None
acceleration_gas (optional): boolian, add z dependent acceleration from dust feedback in equilibrium. Defaults False
dust (optional): boolian, adds dust to simulation if True. Defaults True
equlibrium (optional): boolian, uses table with z dependent equlibrium values if true. Defaults False
shear (optional): str, qz constant vertical azymuthal shear. Defaults None
tracer (optional): list ([Stokes tracer, Mass fraction (m_tracer/m_gas), H_tracer]) of a tracer particle, no tracer if None. Defaults None'''
setup.create(pd, shearing_box, output, cfl=0.44,
             equilibrium=True, dust=False,
             perturbation=f_wn,
             acceleration_gas=True,
            #  tracer=tracer,
             # shear=shear,
             # eta=0.05,
             # gauss=True,
             nubound=1e-4
             )