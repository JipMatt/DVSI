import numpy as np
import math
import sys, os

from scipy.special import roots_legendre
from scipy.interpolate import BarycentricInterpolator
from scipy.optimize import curve_fit 

dir = os.path.dirname(os.path.dirname(__file__))+'/plotting/code/'
sys.path.append(dir)

from psi_stratified.equilibrium import Equilibrium
from psi_stratified.stokesdensity import StokesDensity

def l_mass(pf,coord, dt, max_frame, fluid=0):
    dz = coord.z[1]-coord.z[0]
    pf.read(0)
    mass0 = sum(pf.Fluids[fluid].dens[:,0,0])
    t = []
    mass = []
    for i in range(max_frame):
        pf.read(i)
        t.append(i*dt)
        mass.append(sum(pf.Fluids[fluid].dens[:,0,0])*dz)
    return np.array(t), np.array(mass), mass0

def Hd(alpha, St):
    return np.sqrt(alpha/(St*(1+St**2)))

def gauss(z, rho0, h):
    return rho0*np.exp(-.5*np.power(z/h,2))

def num_gauss(z):
    zmin = z-0.5*np.diff(z)[0]
    zmin = np.append(zmin, max(z)+0.5*np.diff(z)[0])
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
    return np.append(lrhom[::-1], lrhop)

def n_gauss(h, zmax):
    return np.sqrt(2*np.pi)*h*math.erf(zmax/(np.sqrt(2)*h))

def sigma(x, taumin=1e-3, taumax=1e-1):
    '''MRN size distribution, normalized to unity'''
    return 0.5/(np.sqrt(x)*(np.sqrt(taumax) - np.sqrt(taumin)))

def anytical(ndust=1, tau=[1e-5],  hg=4, metallicity=.01, eta=0.05, alpha_viscosity=1e-6, n_coll=100):
    z = np.linspace(-hg, hg, int(2*hg*1024))
    stokes_density = StokesDensity(tau)

    equilibrium = Equilibrium(metallicity, stokes_density, ndust, alpha_viscosity)
    equilibrium.solve_horizontal_velocities(n_coll, neglect_gas_viscosity=False)
    rhog, sigma, dg_ratio, dust_rho, gas_vx, gas_vy, gas_vz, dust_ux, dust_uy, dust_uz = equilibrium.get_state(z, eta=eta)

    dens_dust = equilibrium.weights[:, np.newaxis]*equilibrium.stokes_numbers[:, np.newaxis]*sigma
    dens_dust = np.where(dens_dust < 1e-11, 1e-11, dens_dust)
    
    print('Ndust: {}, tau:{}, Zmax:{}Hg, Z:{}, alpha:{:.1e}'.format(ndust,tau,hg,metallicity, alpha_viscosity))
    print('Mdust: {}, Mgas: {} in domain'.format(np.sum(dust_rho*np.diff(z)[0]), np.sum(rhog*np.diff(z)[0])))
    print('Dust-to-gas ratio in domain:', np.sum(dust_rho)/np.sum(rhog))
    print('Peak dust density:', np.max(dust_rho))
    print('dz: {}\n'.format(np.diff(z)[0]))
    
    ax_gas = np.sum(sigma * (dust_ux - gas_vx) / equilibrium.stokes_numbers[:, np.newaxis], axis=0) / rhog
    ay_gas = np.sum(sigma * (dust_uy - gas_vy) / equilibrium.stokes_numbers[:, np.newaxis], axis=0) / rhog
    az_gas = np.sum(sigma * (dust_uz - gas_vz) / equilibrium.stokes_numbers[:, np.newaxis], axis=0) / rhog


    #azymutal shear
    fwhm_min_gas = np.min(z[rhog > np.nanmax(rhog)/2])
    fwhm_max_gas = np.max(z[rhog > np.nanmax(rhog)/2])
    fwhm_min_vx = np.min(z[gas_vx > np.nanmax(gas_vx)/2])
    fwhm_max_vx = np.max(z[gas_vx > np.nanmax(gas_vx)/2])
    fwhm_min_vy = np.min(z[gas_vy+1 > np.nanmax(gas_vy+1)/2])
    fwhm_max_vy = np.max(z[gas_vy+1 > np.nanmax(gas_vy+1)/2])
    fwhm_min_dust = np.min(z[dust_rho > np.nanmax(dust_rho)/2])
    fwhm_max_dust = np.max(z[dust_rho > np.nanmax(dust_rho)/2])
    dust_ux0 = dust_ux[0,:]-dust_ux[0,0]
    try:
        fwhm_min_ux = np.min(z[dust_ux0 > np.nanmax(dust_ux0)/2])
        fwhm_max_ux = np.max(z[dust_ux0 > np.nanmax(dust_ux0)/2])
    except:
        fwhm_min_ux = np.min(z)
        fwhm_max_ux = np.max(z)



    print('sigma gas: {}'.format((fwhm_max_gas-fwhm_min_gas)/2.355))
    print('sigma vx: {}'.format((fwhm_max_vx-fwhm_min_vx)/2.355))
    print('sigma vy: {} \n'.format((fwhm_max_vy-fwhm_min_vy)/2.355))
    print('sigma dust: {}'.format((fwhm_max_dust-fwhm_min_dust)/2.355))
    print('sigma ux: {} \n'.format((fwhm_max_ux-fwhm_min_ux)/2.355))

    zsplit=(z[1:]+z[:-1])/2

    z_dmax_vx = zsplit[np.argmax(np.diff(gas_vx)/np.diff(z))]
    z_dmin_vx = zsplit[np.argmin(np.diff(gas_vx)/np.diff(z))]
    z_dmax_ux = zsplit[np.argmax(np.diff(dust_ux[0,:])/np.diff(z))]
    z_dmin_ux = zsplit[np.argmin(np.diff(dust_ux[0,:])/np.diff(z))]
    z_dmax_vy = zsplit[np.argmax(np.diff(gas_vy)/np.diff(z))]
    z_dmin_vy = zsplit[np.argmin(np.diff(gas_vy)/np.diff(z))]
    z_dmin_hmax = np.max(zsplit[(np.diff(gas_vy)/np.diff(z)) < np.min(np.diff(gas_vy)/np.diff(z))/2])
    z_dmin_hmin = np.min(zsplit[(np.diff(gas_vy)/np.diff(z)) < np.min(np.diff(gas_vy)/np.diff(z))/2])
    z_dmax_hmax = np.max(zsplit[(np.diff(gas_vy)/np.diff(z)) > np.max(np.diff(gas_vy)/np.diff(z))/2])
    z_dmax_hmin = np.min(zsplit[(np.diff(gas_vy)/np.diff(z)) > np.max(np.diff(gas_vy)/np.diff(z))/2])
    z_dmax_ax = zsplit[np.argmax(np.diff(ax_gas)/np.diff(z))]
    z_dmin_ax = zsplit[np.argmin(np.diff(ax_gas)/np.diff(z))]
    z_dmax_ay = zsplit[np.argmax(np.diff(ay_gas)/np.diff(z))]
    z_dmin_ay = zsplit[np.argmin(np.diff(ay_gas)/np.diff(z))]

    print('peak dvx/dz at z={} is {}'.format((z_dmin_vx-z_dmax_vx)/2, np.max(np.diff(gas_vx)/np.diff(z))))
    print('peak dvy/dz at z={} is {}, and has a FWHM of {} or sigma of {}'.format((z_dmin_vy-z_dmax_vy)/2, np.max(np.diff(gas_vy)/np.diff(z)),((z_dmax_hmax-z_dmax_hmin)+(z_dmin_hmax-z_dmin_hmin))/2,((z_dmax_hmax-z_dmax_hmin)+(z_dmin_hmax-z_dmin_hmin))/(4*np.sqrt(2*np.log(2)))))
    print('peak dax/dz at z={} is {}'.format((z_dmax_ax-z_dmin_ax)/2, np.max(np.diff(ax_gas)/np.diff(z))))
    print('peak day/dz at z={} is {}'.format((z_dmin_ay-z_dmax_ay)/2, np.max(np.diff(ay_gas)/np.diff(z))))
    print('peak dux/dz at z={} is {}\n'.format((z_dmin_ux-z_dmax_ux)/2, np.max(np.diff(dust_ux[0,:])/np.diff(z))))

    return z, rhog, gas_vx, gas_vy, gas_vz, ax_gas, ay_gas, az_gas, dust_rho, dust_ux, dust_uy, dust_uz, sigma, dens_dust, equilibrium.stokes_numbers, equilibrium.weights

def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='x',anchorpad=0,**kw):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw)) 
                    for text,color in zip(list_of_strings,list_of_colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=(0.2, -0.09),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',rotation=90,**kw)) 
                     for text,color in zip(list_of_strings[::-1],list_of_colors) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(-0.10, 0.2), 
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)
    
def sizedistributie(sigma, stokes, ndust):
    xi, weights = roots_legendre(ndust)
    xi = np.asarray(xi)
    weights = np.asarray(weights)
    x=xi
    xi = np.linspace(-1, 1, 100)
    res = BarycentricInterpolator(x, sigma)(xi)
    logts = np.exp(BarycentricInterpolator(x, np.log(stokes))(xi))
    return logts, res


def reynolds(number, pf, coord, shear=0, v=np.array([0,0,0]), minmean=True):
    '''
    This function calculates the reneyolds stress tensor R_xz and R_xy R_zy R_zx.
    Input
    -----
    pf: polydust.object
    coord: array, from polydust.object coordinates of the simulation box
    number: int, framenumber
    '''

    vn = []
    for i in range(3):
        if not isinstance(v[i], (np.ndarray)):
            vn.append(v[i]*np.ones_like(coord.z))
        else:
            vn.append(v[i])
    v = np.array(vn)

    dx = np.diff(coord.x)[0]
    dz = np.diff(coord.z)[0]
    pf.read(number)
    meshX, meshZ = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')

    dens = pf.Fluids[0].dens[:,0,:]
    velx = pf.Fluids[0].velx[:,0,:]-v[0][:,None]
    vely = pf.Fluids[0].vely[:,0,:]-v[1][:,None]+ 1.5*meshX - shear*meshZ
    velz = pf.Fluids[0].velz[:,0,:]-v[2][:,None]

    if minmean:
        delta_velx = velx - np.nanmean(velx, axis=1, keepdims=True)
        delta_vely = vely - np.nanmean(vely, axis=1, keepdims=True)
        delta_velz = velz - np.nanmean(velz, axis=1, keepdims=True)

    R_xz_2d = dens * velx * delta_velz*dx*dz
    R_xy_2d = dens * velx * delta_vely*dx*dz
    R_zy_2d = dens * velz * delta_vely*dx*dz
    

    zmass = np.nansum(dens*dx*dz, axis=1)
    R_xz_1dz = np.nansum(R_xz_2d, axis=1)/zmass
    R_xy_1dz = np.nansum(R_xy_2d, axis=1)/zmass
    R_zy_1dz = np.nansum(R_zy_2d, axis=1)/zmass
    E_x_1dz = np.nansum(.5* dens * velx *velx * dx * dz, axis=1)
    E_y_1dz = np.nansum(.5* dens * vely *vely * dx * dz, axis=1)
    E_z_1dz = np.nansum(.5* dens * velz *velz * dx * dz, axis=1)
    E_1dz = np.nansum(.5* dens * (velx *velx + vely * vely + velz * velz) * dx * dz, axis=1)

    R_xz = np.nanmean(R_xz_1dz)
    R_xy = np.nanmean(R_xy_1dz)
    R_zy = np.nanmean(R_zy_1dz)
    E_x = np.nansum(E_x_1dz)
    E_y = np.nansum(E_y_1dz)
    E_z = np.nansum(E_z_1dz)
    E = np.nansum(E_1dz)

    return {'reyn':{'2D':{'R_xz': R_xz_2d, 'R_xy': R_xy_2d, 'R_zy': R_zy_2d},
                    '1D':{'R_xz': R_xz_1dz, 'R_xy': R_xy_1dz, 'R_zy': R_zy_1dz},
                  'mean':{'R_xz': R_xz, 'R_xy': R_xy, 'R_zy': R_zy}},
            'Ekin':{'1D':{'E_x': E_x_1dz,'E_y': E_y_1dz, 'E_z': E_z_1dz, 'E': E_1dz},
                   'tot':{'E_x': E_x,'E_y': E_y, 'E_z': E_z, 'E': E}}}


def density_structure(fluid, number, pf, coord,  crit_dens=1e-8):
    '''
    This function calculates the radial strucure dependent on Z from the FARGO3D simulation
    of the density for the gas or dust.
    Input
    -----
    number: int, framenumber
    dt: float, timestep between frames
    fname: str, directory name to FAGRO3D files (exmpl: /folder/name/to/simulations/gas or /fname/dust1)
    x,z: array, coordinates simulations box
    '''

    z = coord.z
    
    if fluid == 'tracer':
        fluid = -1
    pf.read(number)
    if fluid =='sum':
        dens = pf.dust_density()[:,0,:]
    elif fluid=='gas' or fluid==0:
        dens = pf.Fluids[0].dens[:,0,:]
    elif fluid < 0 or fluid == 'tracer':
        dens = pf.Fluids[abs(fluid)+pf.n_dust].dens[:,0,:]
    else:
        dens = pf.Fluids[fluid].dens[:,0,:]
    mu = dens/pf.Fluids[0].dens[:,0,:]

    dens = np.where(dens < crit_dens, np.nan, dens)
    mu = np.where(mu < crit_dens, np.nan, mu)


    dens_r = np.nanmean(dens, axis=1)
    std_r = np.nanstd(dens, axis=1)
    max_r =  np.nanmax(dens,axis=1)

    mu_r = np.nanmean(mu, axis=1)
    std_mu_r = np.nanstd(mu, axis=1)
    max_mu_r =  np.nanmax(mu,axis=1)

    dens_norm = dens / dens_r[:,None]
    std_norm_r = np.nanstd(dens_norm,axis=1)
    max_norm_r = np.nanmax(dens_norm,axis=1)

    fwhm_min = np.min(z[dens_r > np.nanmax(dens_r)/2])
    fwhm_max = np.max(z[dens_r > np.nanmax(dens_r)/2])
    try:
        parameters, covariance = curve_fit(gauss, z, dens_r, nan_policy='omit')
    except: 
        parameters = np.array([np.nan, np.nan])
    if parameters[1] < 0: parameters[1]=np.nan

    fwhm_min_mu = np.min(z[mu_r > np.nanmax(mu_r)/2])
    fwhm_max_mu = np.max(z[mu_r > np.nanmax(mu_r)/2])
    try:
        parameters_mu, covariance_mu = curve_fit(gauss, z, mu_r, nan_policy='omit')
    except:
        parameters_mu = np.array([np.nan, np.nan])
    if parameters_mu[1] <0: parameters_mu[1]=np.nan

    index_mid = abs(z-0).argmin()
    index_hg1 = abs(z-1).argmin()
    index_hg2 = abs(z+1).argmin()
    index_hd1 = abs(z-parameters[1]).argmin()
    index_hd2 = abs(z+parameters[1]).argmin()
    index_max = np.nanargmax(abs(max_r))


    index_hmu1 = abs(z-parameters_mu[1]).argmin()
    index_hmu2 = abs(z+parameters_mu[1]).argmin()
    index_mu_max = np.nanargmax(abs(max_mu_r))

    dens_structure = {
    'norm': { 'std': {'tot': np.nanstd(dens_norm),
                      'mid': np.nanstd(dens_norm[index_mid,:]),
                       'hg': (np.nanstd(dens_norm[index_hg1,:]) + np.nanstd(dens_norm[index_hg2,:]))/2,
                       'hd': (np.nanstd(dens_norm[index_hd1,:]) + np.nanstd(dens_norm[index_hd2,:]))/2,
                    'h_mu': (np.nanstd(dens_norm[index_hmu1,:]) + np.nanstd(dens_norm[index_hmu2,:]))/2,
                   'max_mu':  np.nanstd(dens_norm[index_mu_max,:])},
              'max': {'tot': np.nanmax(dens_norm),
                      'mid': np.nanmax(dens_norm[index_mid,:]),
                       'hg': (np.nanmax(dens_norm[index_hg1,:]) + np.nanmax(dens_norm[index_hg2,:]))/2,
                       'hd': (np.nanmax(dens_norm[index_hd1,:]) + np.nanmax(dens_norm[index_hd2,:]))/2,
                     'h_mu': (np.nanmax(dens_norm[index_hmu1,:]) + np.nanmax(dens_norm[index_hmu2,:]))/2,
                   'max_mu':  np.nanmax(dens_norm[index_mu_max,:])}},
      'mu': {'mean': {'tot': np.nanmean(mu),
                      'mid': np.nanmean(mu[index_mid,:]),
                       'hg': (np.nanmean(mu[index_hg1,:]) + np.nanmean(mu[index_hg2,:]))/2,
                       'hd': (np.nanmean(mu[index_hd1,:]) + np.nanmean(mu[index_hd2,:]))/2,
                     'h_mu': (np.nanmean(mu[index_hmu1,:]) + np.nanmean(mu[index_hmu2,:]))/2,
                   'max_mu':  np.nanmean(mu[index_mu_max,:])},
              'std': {'tot': np.nanstd(mu),
                      'mid': np.nanstd(mu[index_mid,:]),
                       'hg': (np.nanstd(mu[index_hg1,:]) + np.nanstd(mu[index_hg2,:]))/2,
                       'hd': (np.nanstd(mu[index_hd1,:]) + np.nanstd(mu[index_hd2,:]))/2,
                    'h_mu': (np.nanstd(mu[index_hmu1,:]) + np.nanstd(mu[index_hmu2,:]))/2,
                   'max_mu':  np.nanstd(mu[index_mu_max,:])},
              'max': {'tot': np.nanmax(mu),
                      'mid': np.nanmax(mu[index_mid,:]),
                       'hg': (np.nanmax(mu[index_hg1,:]) + np.nanmax(mu[index_hg2,:]))/2,
                       'hd': (np.nanmax(mu[index_hd1,:]) + np.nanmax(mu[index_hd2,:]))/2,
                     'h_mu': (np.nanmax(mu[index_hmu1,:]) + np.nanmax(mu[index_hmu2,:]))/2,
                   'max_mu':  np.nanmax(mu[index_mu_max,:])},
               'loc': {'sigma': parameters_mu[1],
                     'fwhm': (fwhm_max_mu-fwhm_min_mu)/2,
                      'max': z[index_mu_max]}},
    'dens': {'mean': {'tot': np.nanmean(dens),
                      'mid': np.nanmean(dens[index_mid,:]),
                       'hg': (np.nanmean(dens[index_hg1,:]) + np.nanmean(dens[index_hg2,:]))/2,
                       'hd': (np.nanmean(dens[index_hd1,:]) + np.nanmean(dens[index_hd2,:]))/2,
                     'h_mu': (np.nanmean(dens[index_hmu1,:]) + np.nanmean(dens[index_hmu2,:]))/2,
                   'max_mu':  np.nanmean(dens[index_mu_max,:])},
              'std': {'tot': np.nanstd(dens),
                      'mid': np.nanstd(dens[index_mid,:]),
                       'hg': (np.nanstd(dens[index_hg1,:]) + np.nanstd(dens[index_hg2,:]))/2,
                       'hd': (np.nanstd(dens[index_hd1,:]) + np.nanstd(dens[index_hd2,:]))/2,
                     'h_mu': (np.nanstd(dens[index_hmu1,:]) + np.nanstd(dens[index_hmu2,:]))/2,
                   'max_mu':  np.nanstd(dens[index_mu_max,:])},
              'max': {'tot': np.nanmax(dens),
                      'mid': np.nanmax(dens[index_mid,:]),
                       'hg': (np.nanmax(dens[index_hg1,:]) + np.nanmax(dens[index_hg2,:]))/2,
                       'hd': (np.nanmax(dens[index_hd1,:]) + np.nanmax(dens[index_hd2,:]))/2,
                     'h_mu': (np.nanmax(dens[index_hmu1,:]) + np.nanmax(dens[index_hmu2,:]))/2,
                   'max_mu':  np.nanmax(dens[index_mu_max,:])},
            'loc': {'sigma': parameters[1],
                     'fwhm': (fwhm_max-fwhm_min)/2, 
                      'max': z[index_max]}}}


    return dens_structure, (dens_r, std_r, max_r), (mu_r, std_mu_r, max_mu_r), (std_norm_r, max_norm_r)

def radial_power_spectrum(frame, pf, coord, direc=None, v=np.array([0,0,0]), shear=0):
    """
    Compute the 1D (isotropic) radial power spectrum P(k) from a 2D field.
    
    Parameters
    ----------
    input : (Nx, Ny) array
        Real-space field sampled on a uniform grid.
    coord : object
        Contains grid information (nx, ny, dx, dy).

    Returns
    -------
    k_centers : (Nbins,) array
        Bin centers in angular wavenumber units [rad / unit].
    Pk : (Nbins,) array
        Azimuthal average of |F(k)|^2 per bin (numpy 'ortho' normalization).
        If you need a specific physical normalization, see the note below.
    n_modes : (Nbins,) array
        Number of Fourier modes contributing to each bin.
    """

    vn = []
    for i in range(3):
        if not isinstance(v[i], (np.ndarray)):
            vn.append(v[i]*np.ones_like(coord.z))
        else:
            vn.append(v[i])
    v = np.array(vn)

    # 2D FFT and power
    X, Z = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')
    pf.read(frame)
    velx = pf.Fluids[0].velx[:,0,:]-v[0][:,None]
    vely = pf.Fluids[0].vely[:,0,:]-v[1][:,None]+1.5*X-shear*Z
    velz = pf.Fluids[0].velz[:,0,:]-v[2][:,None]

    Vx = np.fft.fft2(velx); Vy = np.fft.fft2(vely); Vz = np.fft.fft2(velz)

    kx = 2*np.pi*np.fft.fftfreq(coord.nx, d=coord.dx)
    ky = 2*np.pi*np.fft.fftfreq(coord.nz, d=coord.dz)
    KX, KY = np.meshgrid(kx, ky, indexing='xy')
    KR = np.hypot(KX, KY)
    A = (coord.dx*coord.dz)/(coord.nx*coord.nz)
    if direc is None:
        E2D = 0.5*A*(np.abs(Vx)**2 + np.abs(Vy)**2 + np.abs(Vz)**2)
    elif direc == 'x':
        E2D = 0.5 * A * np.abs(Vx)**2
    elif direc == 'y':
        E2D = 0.5 * A * np.abs(Vy)**2
    elif direc == 'z':
        E2D = 0.5 * A * np.abs(Vz)**2
    # radial shells (Lx<Ly → use finer spacing 2π/Ly)
    dk = 2*np.pi / max(coord.lx, coord.lz)  # equals 2π/Ly if Lx<Ly
    kmax = np.pi/coord.dx
    edges = np.arange(0.0, kmax+dk, dk)
    mask = KR <= kmax
    KRv = KR[mask].ravel()
    Wv = E2D[mask].ravel()

    which = np.digitize(KRv, edges) - 1
    nb = len(edges)-1
    which = np.clip(which, 0, nb-1)
    Esum = np.bincount(which, weights=Wv, minlength=nb)
    Ncnt = np.bincount(which, minlength=nb)
    valid = Ncnt > 0
    kcent = 0.5*(edges[:-1] + edges[1:])[valid]
    Ek = Esum[valid] / dk

    return kcent, Ek, Ncnt

def max_amplitude(frame, pf, coord, order=10,  v=np.array([0,0,0]), shear=0):
    # 2D FFT (properly normalized)

    vn = []
    for i in range(3):
        if not isinstance(v[i], (np.ndarray)):
            vn.append(v[i]*np.ones_like(coord.z))
        else:
            vn.append(v[i])
    v = np.array(vn)

    X, Z = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')
    pf.read(frame)
    velx = pf.Fluids[0].velx[:,0,:]-v[0][:,None]
    vely = pf.Fluids[0].vely[:,0,:]-v[1][:,None]+1.5*X-shear*Z
    velz = pf.Fluids[0].velz[:,0,:]-v[2][:,None]

    vel = np.sqrt(velx**2+vely**2+velz**2)
    
    F_mag = np.abs(np.fft.fftshift(np.fft.fft2(vel))) / (coord.nx * coord.nz)
    Fx_mag = np.abs(np.fft.fftshift(np.fft.fft2(velx))) / (coord.nx * coord.nz)
    Fy_mag = np.abs(np.fft.fftshift(np.fft.fft2(vely))) / (coord.nx * coord.nz)
    Fz_mag = np.abs(np.fft.fftshift(np.fft.fft2(velz))) / (coord.nx * coord.nz)


    # freq centers (angular frequency)
    kx_c = np.fft.fftshift(np.fft.fftfreq(coord.nx, d=coord.dx)) * 2*np.pi
    kz_c = np.fft.fftshift(np.fft.fftfreq(coord.nz, d=coord.dz)) * 2*np.pi
    # Outline the actual (possibly skewed) quad for the max cell

    max_amp = {'vel': {}, 'velx': {}, 'vely': {}, 'velz': {}}
    for i in range(order):

        i_max, j_max = np.unravel_index(np.argmax(F_mag), F_mag.shape)
        ix_max, jx_max = np.unravel_index(np.argmax(Fx_mag), Fx_mag.shape)
        iy_max, jy_max = np.unravel_index(np.argmax(Fy_mag), Fy_mag.shape)
        iz_max, jz_max = np.unravel_index(np.argmax(Fz_mag), Fz_mag.shape)

        max_amp['vel']['max{}'.format(i)] = {'amp': np.max(F_mag), 'kx': np.asanyarray(kx_c)[j_max], 'kz': np.asanyarray(kz_c)[i_max]}
        max_amp['velx']['max{}'.format(i)] = {'amp': np.max(Fx_mag), 'kx': np.asanyarray(kx_c)[jx_max], 'kz': np.asanyarray(kz_c)[ix_max]}
        max_amp['vely']['max{}'.format(i)] = {'amp': np.max(Fy_mag), 'kx': np.asanyarray(kx_c)[jy_max], 'kz': np.asanyarray(kz_c)[iy_max]}
        max_amp['velz']['max{}'.format(i)] = {'amp': np.max(Fz_mag), 'kx': np.asanyarray(kx_c)[jz_max], 'kz': np.asanyarray(kz_c)[iz_max]}

        F_mag[i_max, j_max] = 0
        Fx_mag[ix_max, jx_max] = 0
        Fy_mag[iy_max, jy_max] = 0
        Fz_mag[iz_max, jz_max] = 0  


    return max_amp