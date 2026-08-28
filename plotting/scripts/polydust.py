import numpy as np
from scipy.integrate import quad
import os
from scipy.special import roots_legendre, erf



class SizeDistribution():
    '''Class holding a size (or stopping time) distribution

    Args:
        stokes_range: minimum and maximum Stokes number
    '''
    def __init__(self, stokes_range):
        self.taumin = stokes_range[0]
        if len(stokes_range) > 1:
            self.taumax = stokes_range[1]
        else:
            self.taumax = self.taumin

    def sigma(self, x):
        '''MRN size distribution, normalized to unity'''
        return 0.5/(np.sqrt(x)*(np.sqrt(self.taumax) - np.sqrt(self.taumin)))
        
    def sigma38(self, x):
        '''MRN size distribution, normalized to unity'''
        return 4/(5*(self.taumax**(4/5) - self.taumin**(4/5)))*x**(-1/5)
    
    def sigma32(self, x):
        '''MRN size distribution, normalized to unity'''
        return 1/(5*(self.taumax**(1/5) - self.taumin**(1/5)))*x**(-4/5)

class Polydust():
    def __init__(self, n_dust, stokes_range,
                 dust_density, gas_density,
                 alpha_viscosity_gas,
                 height_box,
                 size_distribution=None,
                 gauss_legendre=True,
                 discrete_equilibrium=False,
                 pressure_gradient=False):
        '''Class for polydisperse dust component in FARGO3D simulations.

        Args:
            n_dust: Number of collocation points ('dust fluids')
            stokes_range: Range of Stokes numbers to consider. In case of discrete sizes, it can be a list of Stokes numbers.
            dust_density: background dust density
            gas_density: background gas density
            alpha_viscosity: the diffusion factor
            height_box: vertical domain of the shearing box
            size_distribution (optional): Either a SizeDistribution object, or a list of dust densities, in which case the numer of elements should be equal to the number of elements of stokes_range. Defaults to None, which is a continuous MRN size distribution.
            gauss_legendre (optional): Flag whether to use Gauss-Legendre collocation points. If False, use equidistant nodes in log(tau) space. Defaults to True.
            discrete_equilibrium (optional): In case of equidistant nodes, force a discrete equilibrium a la Krapp et al. and Chao-Chin & Zhu. Defaults to False.
        '''

        self.N = n_dust
        self.stokes_range = np.asarray(stokes_range)
        self.dust_density = dust_density
        self.gas_density = gas_density
        self.alpha_viscosity_gas = alpha_viscosity_gas
        self.alpha_viscosity_dust = None  # Default is gas viscosity
        self.height_box = height_box
        self.gauss_legendre = gauss_legendre
        self.discrete_equilibrium = discrete_equilibrium
        self.size_distribution = size_distribution
        self.pressure_gradient = pressure_gradient

        # Set the size density function
        self.sigma = None
        if isinstance(self.size_distribution, SizeDistribution):
            # Continuous size distribution
            print('Continuous size distribution with {} nodes'.format(self.N))
            if self.gauss_legendre is True:
                print('Using Gauss-Legendre quadrature')
            else:
                print('Using equidistant nodes')
                if self.discrete_equilibrium is True:
                    print('Using discrete equilibrium')
            self.sigma = self.size_distribution.sigma
        else:
            if self.size_distribution is None:
                # MRN, normalized to unity
                print('Continuous MRN size distribution with {} nodes'.format(self.N))
                if self.gauss_legendre is True:
                    print('Using Gauss-Legendre quadrature')
                else:
                    print('Using equidistant nodes')
                    if self.discrete_equilibrium is True:
                        print('Using discrete equilibrium')
                self.size_distribution = SizeDistribution(stokes_range)
                self.sigma = self.size_distribution.sigma
            else:
                print('Discrete dust sizes with {} dust species'.format(self.N))
                self.N = len(stokes_range)
                # Discrete sizes passed, size should be the same as stokes_range
                self.sigma = np.asarray(self.size_distribution)
                # Normalize to total dust density unity
                self.sigma = self.sigma/np.sum(self.sigma)

    @classmethod
    def from_output_direc(cls, direc):
        n_dust = len([name for name in os.listdir(direc) if name.endswith('dens0_0.dat')]) - 1

        print('Number of dust species:', n_dust)

        # NOTE: Stokes range and gas and dust densities not read!
        return cls(n_dust, [0.01, 0.1], 3.0, 1.0)
    
    def beta(self, st):
        '''dust continuos solution'''
        # return (1-(1-4*st*st)**0.5)/(2*st)
        return st

    def D(self, st):
        '''dust difussion coeficcient'''
        # return (1+st+4*st*st)*(1+st*st)**(-2)*self.alpha_viscosity
        if self.alpha_viscosity_dust is None:
            # Use gas viscosity
            return  self.alpha_viscosity_gas
        else:
            # Use dust viscosity
            return self.alpha_viscosity_dust
        
    def Hd(self, st):
        '''Scale height dust'''
        # return np.sqrt(self.D(st)/(self.D(st)+st))
        return np.sqrt(self.D(st)/(st*(1+st*st)))
    
    def norm(self, st):
        ''''normilzation factor'''
        return (np.sqrt(2*np.pi)*erf(self.height_box/np.sqrt(2)))/(np.sqrt(2*np.pi)*self.Hd(st)*erf(self.height_box/(np.sqrt(2)*self.Hd(st))))

    def nodes_and_weights(self):
        '''Return Gauss-Legendre nodes and weights'''
        return roots_legendre(self.N)

    def dust_nodes(self):
        '''Return collocation points'''

        if callable(self.sigma):
            # Continuous size distribution
            if self.gauss_legendre is True:
                # Gauss-Legendre collocation points
                xi, weights = self.nodes_and_weights()

                xi = np.asarray(xi)
                weights = np.asarray(weights)

                if len(self.stokes_range) > 1:
                    q = self.stokes_range[1]/self.stokes_range[0]
                else:
                    q = 1.0

                # Stopping time nodes
                tau = self.stokes_range[0]*np.power(q, 0.5*(xi + 1))
            else:
                # Constant spacing in log tau space
                xi1 = np.log(self.stokes_range[0])
                xi2 = np.log(self.stokes_range[1])

                tau = np.exp(np.linspace(xi1, xi2, self.N))

                xi_edge = np.linspace(xi1, xi2, self.N + 1)
                tau = xi_edge + 0.5*(xi_edge[1] - xi_edge[0])
                tau = np.exp(tau)[0:-1]

        else:
            # Discrete dust sizes
            tau = self.stokes_range

        return tau

    def dust_densities(self):
        '''Return densities at collocation points'''
        if callable(self.sigma):
            # Continuous size distribution
            if self.gauss_legendre is True:
                # Gauss-Legende collocation points
                xi, weights = roots_legendre(self.N)

                xi = np.asarray(xi)
                weights = np.asarray(weights)

                if len(self.stokes_range) > 1:
                    q = self.stokes_range[1]/self.stokes_range[0]
                else:
                    q = 1.0

                # Stopping time nodes
                tau = self.stokes_range[0]*np.power(q, 0.5*(xi + 1))

                # 'density' at nodes
                dens_dust = 0.5*np.log(q)*weights*tau*self.dust_density*self.sigma(tau)
            else:
                # Constant spacing in log tau space
                xi1 = np.log(self.stokes_range[0])
                xi2 = np.log(self.stokes_range[1])

                xi_edge = np.linspace(xi1, xi2, self.N + 1)
                dxi = xi_edge[1] - xi_edge[0]
                tau = np.exp(xi_edge + 0.5*dxi)[0:-1]

                dens_dust = self.dust_density*self.sigma(tau)*dxi*tau

                # Set dust density so mass in size bin is exact MRN
                if (self.discrete_equilibrium is True):
                    t_l = np.exp(0.5*xi_edge)[0:-1]
                    t_r = np.exp(0.5*np.roll(xi_edge, -1))[0:-1]
                    dens_dust = self.dust_density*(t_r - t_l)/(t_r[-1] - t_l[0])

                    # Make sure self.sigma is no longer callable: discrete!
                    self.sigma = dens_dust/self.dust_density
                    self.stokes_range = tau

        else:
            # Discrete dust sizes
            tau = self.stokes_range


            # Normalize so total dust_density is correct
            dens_dust = self.sigma*self.dust_density

        print('Total dust density: ', np.sum(dens_dust))

        return tau, dens_dust


    def dust_densities_midplane(self):
        '''Return densities at collocation points'''
        if callable(self.sigma):
            # Continuous size distribution
            if self.gauss_legendre is True:
                # Gauss-Legende collocation points
                xi, weights = roots_legendre(self.N)

                xi = np.asarray(xi)
                weights = np.asarray(weights)

                if len(self.stokes_range) > 1:
                    q = self.stokes_range[1]/self.stokes_range[0]
                else:
                    q = 1.0

                # Stopping time nodes
                tau = self.stokes_range[0]*np.power(q, 0.5*(xi + 1))

                dens_dust = 0.5*np.log(q)*weights*tau*self.sigma(tau)*self.norm(tau)
            else:
                # Constant spacing in log tau space
                xi1 = np.log(self.stokes_range[0])
                xi2 = np.log(self.stokes_range[1])

                xi_edge = np.linspace(xi1, xi2, self.N + 1)
                dxi = xi_edge[1] - xi_edge[0]
                tau = np.exp(xi_edge + 0.5*dxi)[0:-1]
                dens_dust = self.dust_density*self.sigma(tau)*dxi*tau*self.norm(tau)

                # Set dust density so mass in size bin is exact MRN
                if (self.discrete_equilibrium is True):
                    t_l = np.exp(0.5*xi_edge)[0:-1]
                    t_r = np.exp(0.5*np.roll(xi_edge, -1))[0:-1]
                    dens_dust = self.dust_density*(t_r - t_l)/(t_r[-1] - t_l[0])*self.norm(tau)

                    # Make sure self.sigma is no longer callable: discrete!

                    self.sigma = dens_dust/self.dust_density
                    self.stokes_range = tau

        else:
            # Discrete dust sizes
            tau = np.asarray(self.stokes_range)

            # Normalize so total dust_density is correct
            dens_dust = self.sigma*self.dust_density*self.norm(tau)
        return tau, dens_dust
        
            
    
    def gas_density_z(self, z, tau):
        #We assume that the Hg = 1, aspectratio = 1 and soundspeed = 1.
        lz = np.asarray(z)
        rho_g = []
        for z in lz:
            if callable(self.sigma):
                f = lambda x: self.D(np.exp(x))*self.norm(np.exp(x))*self.sigma(np.exp(x))\
                    *(np.exp(-self.beta(np.exp(x))*z*z/(2*self.D(np.exp(x))))-1)
                intergral = quad(f, np.log(self.stokes_range[0]),
                           np.log(self.stokes_range[1]))[0]
                rho_g.append(np.exp(-0.5*z*z)*np.exp(intergral*self.dust_density/self.gas_density))
            else:
               tau = np.asarray(self.stokes_range)
               rho_g.append(np.exp(-0.5*z*z)*np.exp(np.sum((self.D(tau)/tau)\
                            *self.dust_density*self.norm(tau)*self.sigma*(np.exp(-self.beta(tau)*z*z/(2*self.D(tau)))-1))))
        return np.asarray(rho_g)    

    def dust_density_z(self, z, gas_dens_z):
        #We assume that the Hg = 1, aspectratio = 1 and soundspeed = 1.
        lz = np.asarray(z)

        dust_dens_z = np.empty((len(lz), self.N))
        for i, z in enumerate(lz):
            if callable(self.sigma):
                if self.gauss_legendre is True:
                    # Gauss-Legende collocation points
                    xi, weights = roots_legendre(self.N)

                    xi = np.asarray(xi)
                    weights = np.asarray(weights)

                    if len(self.stokes_range) > 1:
                        q = self.stokes_range[1]/self.stokes_range[0]
                    else:
                        q = 1.0

                    # Stopping time nodes
                    tau = self.stokes_range[0]*np.power(q, 0.5*(xi + 1))

                    dust_dens_z[i] = gas_dens_z[i]*self.dust_density*0.5*np.log(q)*weights*tau\
                        *self.sigma(tau)*self.norm(tau)*(np.exp(-self.beta(tau)*z*z/(2*self.D(tau))))
                else:
                    # Constant spacing in log tau space
                    xi1 = np.log(self.stokes_range[0])
                    xi2 = np.log(self.stokes_range[1])

                    xi_edge = np.linspace(xi1, xi2, self.N + 1)
                    dxi = xi_edge[1] - xi_edge[0]
                    tau = np.exp(xi_edge + 0.5*dxi)[0:-1]
                    dust_dens_z[i] = gas_dens_z[i]*self.dust_density*self.sigma(tau)*dxi*tau\
                        *self.norm(tau)*(np.exp(-self.beta(tau)*z*z/(2*self.D(tau))))

                    # Set dust density so mass in size bin is exact MRN
                    if (self.discrete_equilibrium is True):
                        t_l = np.exp(0.5*xi_edge)[0:-1]
                        t_r = np.exp(0.5*np.roll(xi_edge, -1))[0:-1]
                        dust_dens_z[i] = gas_dens_z[i]*self.dust_density*(t_r - t_l)/(t_r[-1] - t_l[0])\
                            *self.norm(tau)*(np.exp(-self.beta(tau)*z*z/(2*self.D(tau))))

            
            else:
                # Discrete dust sizes
                tau = np.asarray(self.stokes_range)

                # Normalize so total dust_density is correct
                dust_dens_z[i] = gas_dens_z[i]*self.sigma*self.dust_density\
                    *self.norm(tau)*(np.exp(-self.beta(tau)*z*z/(2*self.D(tau))))
        return np.where(dust_dens_z > 1e-11, dust_dens_z, 1e-11)
    
    def equilibrium_velocities(self, tau):
        '''Return equilibrium dust and gas velocities

        Args:
           tau: list of stopping times
        '''
        if callable(self.sigma):
            # Continuous size distribution
            f = lambda x: np.exp(x)*self.sigma(np.exp(x))/(1.0 + np.exp(2*x))
            J0 = quad(f, np.log(self.stokes_range[0]),
                         np.log(self.stokes_range[1]))[0]
            f = lambda x: np.exp(2*x)*self.sigma(np.exp(x))/(1.0 + np.exp(2*x))
            J1 = quad(f, np.log(self.stokes_range[0]),
                         np.log(self.stokes_range[1]))[0]

            J0 = J0*self.dust_density/self.gas_density
            J1 = J1*self.dust_density/self.gas_density

            denom = (1 + J0)*(1 + J0) + J1*J1
            v = []
            if self.pressure_gradient:
                print(2*J1/denom)
                v.append(2*J1/denom)      # Gas vx
                v.append(-(1 + J0)/denom) # Gas vy
                v.append(0.0)             # Gas vz

                # Dust velocities
                for t in tau:
                    v.append(2*(J1 - t*(1 + J0))/((1 + t*t)*denom))
                    v.append(-(1 + J0 + t*J1)/((1 + t*t)*denom))
                    v.append(0.0)
            else:
                v.append(0.0) # Gas vx
                v.append(0.0) # Gas vy
                v.append(0.0) # Gas vz

                # Dust velocities
                for t in tau:
                    v.append(0.0)
                    v.append(0.0)
                    v.append(0.0)                
        else:
            # Discrete dust sizes, override argument
            tau = np.asarray(self.stokes_range)
            
            mu = self.dust_density/self.gas_density
            AN = mu*np.sum(self.sigma*tau/(1 + tau*tau))
            BN = 1.0 + mu*np.sum(self.sigma/(1 + tau*tau))

            v = []
            if self.pressure_gradient:
                v.append(2*AN/(AN*AN + BN*BN))   # Gas vx
                v.append(-BN/(AN*AN + BN*BN))    # Gas vy
                v.append(0.0)                    # Gas vz

                for t in tau:
                    v.append((v[0] + 2*t*v[1])/(1 + t*t))
                    v.append((v[1] - 0.5*t*v[0])/(1 + t*t))
                    v.append(0.0)
            else:
                v.append(0.0)   # Gas vx
                v.append(0.0)   # Gas vy
                v.append(0.0)   # Gas vz

                for t in tau:
                    v.append(0.0)
                    v.append(0.0)
                    v.append(0.0)                
            
        return np.asarray(v)


    def equilibrium_velocities_z(self, lz, tau, gas_dens_z, dust_dens_z):
        '''Return equilibrium dust and gas velocities

        Args:
           tau: list of stopping times
        '''
        v_z = np.empty((len(lz),int(3+self.N*3)))
        for i, z in enumerate(lz):
            if callable(self.sigma):
                # Continuous size distribution
                # this is the average sigma over the whole Z range at not simga at a specific height
                # There is no analytical solution for a stratifield vx and vy,
                #  we make a small correction with z dependent dust-gas-ratio
                f = lambda x: np.exp(x)*gas_dens_z[i]*self.sigma(np.exp(x))*self.norm(np.exp(x))\
                    *np.exp(-self.beta(np.exp(x))*z*z/(2*self.D(np.exp(x))))/(1.0 + np.exp(2*x))
                J0 = quad(f, np.log(self.stokes_range[0]),
                            np.log(self.stokes_range[1]))[0]
                
                f = lambda x: np.exp(2*x)*gas_dens_z[i]*self.sigma(np.exp(x))*self.norm(np.exp(x))\
                    *np.exp(-self.beta(np.exp(x))*z*z/(2*self.D(np.exp(x))))/(1.0 + np.exp(2*x))
                J1 = quad(f, np.log(self.stokes_range[0]),
                            np.log(self.stokes_range[1]))[0]

                J0 = J0*dust_dens_z.sum(1)[i]/gas_dens_z[i]
                J1 = J1*dust_dens_z.sum(1)[i]/gas_dens_z[i]
                mu = dust_dens_z.sum(1)[i]/gas_dens_z[i]
                v = []
                denom = (1 + J0)*(1 + J0) + J1*J1
                if self.pressure_gradient:
                    v.append(2*mu*(1+mu)**(-3)*0.1*z)      # Gas vx
                    v.append(-1/(1+mu)) # Gas vy
                    v.append(0.0)             # Gas vz

                    # Dust velocities
                    for t in tau:
                        v.append(2*mu*(1+mu)**(-3)*t*z)
                        v.append(-1/(1+mu))
                        v.append(-self.beta(t)*z)
                    
                else:
                    v.append(0.0) # Gas vx
                    v.append(0.0) # Gas vy
                    v.append(-self.beta(t)*z) # Gas vz

                    # Dust velocities
                    for t in tau:
                        v.append(0.0)
                        v.append(0.0)
                        v.append(0.0)
                v_z[i] = np.asarray(v)         
            else:
                # Discrete dust sizes, override argument
                tau = np.asarray(self.stokes_range)

                # sigma = dust_dens_z[i]/dust_dens_z.sum(1)[i]
                mu = dust_dens_z.sum(1)[i]/gas_dens_z[i]

                # AN = mu*np.sum(sigma*tau/(1 + tau*tau))
                # BN = 1.0 + mu*np.sum(sigma/(1 + tau*tau))

                v = []
                if self.pressure_gradient:
                    v.append(2*mu*(1+mu)**(-3)*0.1*z)   # Gas vx
                    v.append(-1/(1+mu))    # Gas vy
                    v.append(0.0)                    # Gas vz

                    for t in tau:
                        v.append(2*mu*(1+mu)**(-3)*t*z)
                        v.append(-1/(1+mu))
                        v.append(-self.beta(t)*z)
                else:
                    v.append(0.0)   # Gas vx
                    v.append(0.0)   # Gas vy
                    v.append(0.0)   # Gas vz

                    for t in tau:
                        v.append(0.0)
                        v.append(0.0)
                        v.append(-self.beta(t)*z)
                v_z[i] = np.asarray(v)         
             

        return np.asarray(v_z)
    
    def initial_conditions(self):
        '''Return FARGO initial conditions (equilibrium)'''

        # Collocation points and densities
        tau, dens = self.dust_densities()

        # Equilibrium velocities
        v = self.equilibrium_velocities(tau)

        return tau, dens, v

    def initial_conditions_z(self, z):
        '''Return FARGO initial conditions (equilibrium)'''

        # Collocation points and densities
        tau, dens = self.dust_densities_midplane()
        
        gas_dens_z = self.gas_density_z(z, tau)

        dust_dens_z = self.dust_density_z(z, gas_dens_z)

        # Equilibrium velocities
        v = self.equilibrium_velocities_z(z, tau, gas_dens_z, dust_dens_z)

        return tau, z, gas_dens_z, dust_dens_z, v