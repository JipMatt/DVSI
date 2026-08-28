import shutil
import numpy as np
import pandas as pd
import os

import initial_conditions
import fargo_opt
import fargo_par
import fargo_boundary
import fargo_boundary2
import fargo_boundary3
import fargo_boundary_outflow
import fargo_boundary_fixed


from psi_stratified.equilibrium import Equilibrium
from psi_stratified.stokesdensity import StokesDensity

class Output:
    '''Class governing FARGO output

    Args:
        output_dir: Output directory
        dt: Time between fine-grain outputs
        Ninterm: Number of fine-grain outputs per full output
        Ntot: Total number of full outputs
    '''
    def __init__(self, output_dir, dt, Ninterm, Ntot):
        self.output_dir = output_dir
        self.dt = dt
        self.Ninterm = Ninterm
        self.Ntot = Ntot

class ShearingBox:
    '''Class governing shearing box parameters

    Args:
        dims: Dimensions of the box, should have 3 elements.
        mesh_size: grid points for each dimension, should have 3 elements.
    '''
    def __init__(self, dims, mesh_size):
        self.dims = dims
        self.mesh_size = mesh_size

    @classmethod
    def from_output_direc(cls, direc):
        y = np.loadtxt(direc+'/domain_y.dat')
        ny = len(y) - 7
        z = np.loadtxt(direc+'/domain_z.dat')
        nz = len(z) - 7

        Ly = -2*y[3]
        Lz = -2*z[3]

        return cls([0, Ly, Lz],[1, ny, nz])

class FargoSetup:
    '''Class for creating a new PSI FARGO setup

    Args:
        setup_name: Name of setup to be created
        fargo_dir (optional): Path to public FARGO directory when not using vanila submodule.
    '''
    def __init__(self, setup_name, fargo_dir=None):
        self.setup_name = setup_name

        # Path of current file, should be /path/to/psi_fargo/psi
        self.psi_dir = os.path.dirname(os.path.abspath(__file__))

        self.fargo_dir = fargo_dir

    def create(self, polydust, shearing_box, output, cfl=0.44, eta=0.05,
               perturbation=None, acceleration_gas=False, dust=True, equilibrium=False, shear=None, tracer=None, gauss=None, nubound=None):
        '''Create FARGO setup

        Args:
            polydust: PolyDust object
            shearing_box: ShearingBox object
            output: Output object
            cfl (optional): Courant number to use. Defaults to 0.44.
            pertubrbation (optional): white noise or single mode pertibution. str for WN or [amp, Kx, Kz] for single mode. Defaults None
            acceleration_gas (optional): boolian, add z dependent acceleration from dust feedback in equilibrium. Defaults False
            dust (optional): boolian, adds dust to simulation if True. Defaults True
            equilibrium (optional): boolian, uses table with z dependent equlibrium values if true. Defaults False
            shear (optional): str, qz constant vertical azymuthal shear. Defaults None
            tracer (optional): list ([Stokes tracer, Mass fraction (m_tracer/m_gas), H_tracer]) of a tracer particle, no tracer if None. Defaults None

        '''
        # Going to create a .opt file, a .par file and initial conditions
        self.polydust = polydust
        created_files = [self.setup_name + '.opt',
                         self.setup_name + '.par',
                         'condinit.c']
        fargo_opt.write_opt_file(self.setup_name,
                                polydust,
                                acceleration_gas=acceleration_gas,
                                dust=dust,
                                tracer=tracer,
                                shear=shear,
                                nubound=nubound)
        fargo_par.write_par_file(self.setup_name,
                                polydust,
                                shearing_box,
                                output,
                                eta=eta,
                                cfl=cfl,
                                dust=dust,
                                tracer=tracer,
                                shear=shear,
                                nubound=nubound)

        # Write initial conditions
        nghz = 3
        if equilibrium:
            df, dens, gas_vx, gas_vy = self.initial(polydust, shearing_box, eta=eta, gauss=gauss)
            path_equilibrium = './setups/' + self.setup_name
            # created_files.extend(\
            # fargo_boundary.write_boundary_files(polydust, dens=dens[:6], v=[[gas_vx[3],gas_vx[-3]],
            #                                    [gas_vy[3], gas_vy[-3]]], dust=dust,
            #                                    tracer=tracer, shear=shear))
            eq_bound = [dens, gas_vx, gas_vy]
            created_files.extend(\
            fargo_boundary2.write_boundary_files(polydust, nghz, eq_bound=eq_bound, dust=dust,
                                    tracer=tracer, shear=shear))
        else: 
            path_equilibrium = None
            stokes, sigma, v = polydust.initial_conditions()
            created_files.extend(\
            fargo_boundary.write_boundary_files(polydust, v=[[[v[0],v[0]],
                                               [v[1], v[1]]]], dust=dust,
                                               tracer=tracer, shear=shear))
        initial_conditions.write_condinit_file(polydust, perturbation=perturbation,
                                              Nz=shearing_box.mesh_size[2],
                                              equilibrium=path_equilibrium,
                                              acceleration_gas=acceleration_gas,
                                              dust=dust, tracer=tracer, shear=shear)

        # Setup directory
        output_dir = self.psi_dir + '/../public/setups/' + self.setup_name
        # If not using submodule FARGO
        if self.fargo_dir is not None:
            output_dir = self.fargo_dir + '/setups/' + self.setup_name

        # Create setup directory if not exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        else:
            print("Warning: setup directory exists; contents will be overwritten")

        # Move all created files to setup directory
        for f in created_files:
            os.replace(f, output_dir + '/' + f)
        

        shutil.copy(self.psi_dir + '/fargo_src/substep1_x.c', output_dir + '/' + 'substep1_x.c')
        shutil.copy(self.psi_dir + '/fargo_src/substep1_y.c', output_dir + '/' + 'substep1_y.c')
        shutil.copy(self.psi_dir + '/fargo_src/substep1_z.c', output_dir + '/' + 'substep1_z.c')
        shutil.copy(self.psi_dir + '/fargo_src/substep1_y.c', output_dir + '/' + 'substep1_y.c')
        shutil.copy(self.psi_dir + '/fargo_src/ymin_bound_0.c', output_dir + '/' + 'ymin_bound_0.c')
        shutil.copy(self.psi_dir + '/fargo_src/ymax_bound_0.c', output_dir + '/' + 'ymax_bound_0.c')
        shutil.copy(self.psi_dir + '/fargo_src/monitor.c', output_dir + '/' + 'monitor.c')
        shutil.copy(self.psi_dir + '/fargo_src/change_arch.c', output_dir + '/' + 'change_arch.c')
        shutil.copy(self.psi_dir + '/fargo_src/dust_diffusion_coefficients.c', output_dir + '/' + 'dust_diffusion_coefficients.c')
        shutil.copy(self.psi_dir + '/fargo_src/LowTasks.c', output_dir + '/' + 'LowTasks.c')
        shutil.copy(self.psi_dir + '/fargo_src/stretch.c', output_dir + '/' + 'stretch.c')
        shutil.copy(self.psi_dir + '/fargo_src/visctensor_cart.c', output_dir + '/' + 'visctensor_cart.c')
        # Remove collision_kernel.h, global.h and define.h from src directory
        # and copy them from psi_src directory
        os.remove(self.fargo_dir + '/src/collision_kernel.h')
        os.remove(self.fargo_dir + '/src/global.h')
        os.remove(self.fargo_dir + '/src/define.h')
        shutil.copy(self.psi_dir + '/fargo_src/global.h', self.fargo_dir + '/src/global.h')
        shutil.copy(self.psi_dir + '/fargo_src/define.h', self.fargo_dir + '/src/define.h')
        shutil.copy(self.psi_dir + '/fargo_src/collision_kernel.h', self.fargo_dir + '/src/collision_kernel.h')
        if equilibrium:
            df.to_csv(self.fargo_dir + '/setups/' + self.setup_name + '/initial.csv', index=False)
        print('Setup created: ' + output_dir)


    def initial(self, polydust, shearing_box, nghz = 3, eta=0.05, gauss=None):
        '''Create dictionary with initial conditions that are verticaly stratified
        Args:
            polydust: PolyDust object
            shearing_box: ShearingBox object
            shear: Shear parameter, if None, no shear is applied.
        Returns:
            df: DataFrame with initial conditions dependent on z
            vel_x: Gas velocity in radial-direction
            vel_y: Gas velocity in azimuthal-direction
        '''
        n_coll = 100
        zmin = np.linspace(-polydust.height_box, polydust.height_box, shearing_box.mesh_size[2], endpoint=False)
        z = zmin+0.5*np.diff(zmin)[0]
        # Monodisperse

        z = np.pad(z, nghz, mode='reflect', reflect_type='odd')
        zmin = np.pad(zmin, (nghz,nghz+1), mode='reflect', reflect_type='odd')
        dict = {}

        stokes_density = StokesDensity(polydust.stokes_range)
        metallicity = polydust.dust_density/polydust.gas_density
        if polydust.alpha_viscosity_dust == 0:
            tau, sigma, v = polydust.initial_conditions()
            vgx = v[0]*eta
            vgy = v[1]*eta
            vgz = v[2]

            vdx = v[3::3]*eta
            vdy = v[4::3]*eta
            vdz = v[5::3]
            lrho = []
            rho = 1
            lrhop = []
            rho = 1
            for i in range(int(len(z)/2)+1, len(z)):
                rhonew = rho*(2 - (zmin[i+1]-zmin[i])*z[i])/(2 + (zmin[i]-zmin[i-1])*z[i])
                lrhop.append(rhonew)
                rho = rhonew
            rho = lrhop[0]
            lrhom = []
            for i in range(int(len(z)/2)+1, 0,-1):
                rhoold = rho*(2 + (zmin[i]-zmin[i-1])*z[i])/(2 - (zmin[i+1]-zmin[i])*z[i])
                lrhom.append(rhoold)
                rho = rhoold

            dict['rho_g'] = np.exp(-0.5*z**2)
            dict['vx_g'] = np.ones_like(z)*vgx#np.pad(gas_vx, nghz, mode='edge')
            dict['vy_g'] = np.ones_like(z)*vgy#np.pad(gas_vy, nghz, mode='edge')
            dict['vz_g'] = np.ones_like(z)*vgz#np.pad(gas_vz, nghz, mode='edge')
            if gauss:
                dict['vy_g'] += 0.5*np.exp(-0.5*z**2)
                dict['ax_g'] = -np.exp(-0.5*z**2)#np.pad(np.sum(sigma*(dust_ux-np.tile(gas_vx, (polydust.N,1)))/polydust.dust_nodes(), axis=0)/rhog, nghz, mode='edge')
                dict['ay_g'] = np.zeros_like(z)
                dict['az_g'] = np.zeros_like(z)
            for i in range(int(polydust.N)):
                dict['rho_d{}'.format(i+1)] = sigma[i]*np.append(lrhom[::-1], lrhop)#np.pad(sigma[i,:], nghz, mode='edge')
                dict['ux_d{}'.format(i+1)] = vdx[i]*np.ones_like(z)#np.pad(dust_ux[i,:], nghz, mode='edge')
                dict['uy_d{}'.format(i+1)] = vdy[i]*np.ones_like(z)#np.pad(dust_uy[i,:], nghz, mode='edge')
                dict['uz_d{}'.format(i+1)] = vdz[i]*np.ones_like(z)#np.pad(dust_uz[i,:], nghz, mode='edge')

            
        else:
            alpha = polydust.alpha_viscosity_dust if polydust.alpha_viscosity_dust is not None else polydust.alpha_viscosity_gas
            equilibrium = Equilibrium(metallicity, stokes_density, polydust.N, alpha)
            if polydust.alpha_viscosity_dust == None and polydust.alpha_viscosity_gas != 0:
                neglect_gas_viscosity = False
            else:
                neglect_gas_viscosity = True
            equilibrium.solve_horizontal_velocities(n_coll, neglect_gas_viscosity=neglect_gas_viscosity)
            eta = eta if polydust.pressure_gradient else 0.0
            rhog, sigma, dg_ratio, dust_rho, gas_vx, gas_vy, gas_vz, dust_ux, dust_uy, dust_uz = equilibrium.get_state(z, eta=eta)
            rhog_min, sigma_min, dg_ratio_min, dust_rho_min, gas_vx_min, gas_vy_min, gas_vz_min, dust_ux_min, dust_uy_min, dust_uz_min = equilibrium.get_state(zmin[:-1], eta=eta)
            sigma = np.where(sigma < 1e-11, 1e-11, sigma)
            sigma_min = np.where(sigma_min < 1e-11, 1e-11, sigma_min)
            dict['rho_g'] = rhog#rhog, nghz, mode='reflect', reflect_type='odd')
            dict['vx_g'] = gas_vx#np.pad(gas_vx, nghz, mode='edge')
            dict['vy_g'] = gas_vy#np.pad(gas_vy, nghz, mode='edge')
            dict['vz_g'] = gas_vz_min#np.pad(gas_vz, nghz, mode='edge')

            dict['ax_g'] = np.sum(sigma*(dust_ux-np.tile(gas_vx, (polydust.N,1)))/polydust.dust_nodes(), axis=0)/rhog#np.pad(np.sum(sigma*(dust_ux-np.tile(gas_vx, (polydust.N,1)))/polydust.dust_nodes(), axis=0)/rhog, nghz, mode='edge')
            dict['ay_g'] = np.sum(sigma*(dust_uy-np.tile(gas_vy, (polydust.N,1)))/polydust.dust_nodes(), axis=0)/rhog#np.pad(np.sum(sigma*(dust_uy-np.tile(gas_vy, (polydust.N,1)))/polydust.dust_nodes(), axis=0)/rhog, nghz, mode='edge')
            dict['az_g'] = np.sum(sigma_min*(dust_uz_min-np.tile(gas_vz_min, (polydust.N,1)))/polydust.dust_nodes(), axis=0)/rhog_min#np.pad(np.sum(sigma*(dust_uz-np.tile(gas_vz, (polydust.N,1)))/polydust.dust_nodes(), axis=0)/rhog, nghz, mode='edge')

            for i in range(int(polydust.N)):
                dict['rho_d{}'.format(i+1)] = sigma[i,:]#np.pad(sigma[i,:], nghz, mode='edge')
                dict['ux_d{}'.format(i+1)] = dust_ux[i,:]#np.pad(dust_ux[i,:], nghz, mode='edge')
                dict['uy_d{}'.format(i+1)] = dust_uy[i,:]#np.pad(dust_uy[i,:], nghz, mode='edge')
                dict['uz_d{}'.format(i+1)] = dust_uz_min[i,:]#np.pad(dust_uz[i,:], nghz, mode='edge')
        dict['z'] = z
        dict['zmin'] = zmin[:-1]
        df = pd.DataFrame(dict)
        limit  = int(nghz*2)
        return df, np.append(dict['rho_g'][:limit],dict['rho_g'][-limit:]), np.append(dict['vx_g'][:limit],dict['vx_g'][-limit:]), np.append(dict['vy_g'][:limit], dict['vy_g'][-limit:])

    # End of fargo_setup.py