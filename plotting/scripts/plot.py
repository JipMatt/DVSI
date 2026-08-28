
#import libaries
import sys, os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.patches import Polygon


import numpy as np
from scipy.optimize import curve_fit 

dir = os.path.dirname(os.path.dirname(__file__))+'/plotting/code/'
sys.path.append(dir)
import func as f

seismic = mpl.cm.get_cmap("seismic", 256)
red = colors.LinearSegmentedColormap.from_list(
    "red", seismic(np.linspace(0.5, 1, 256)))
blue = colors.LinearSegmentedColormap.from_list(
    "blue", seismic(np.linspace(0, .5, 256)))

def heatmap(fluid, fig, ax, number, dt, pf, coord, v=np.array([0.,0.,0.]), norm=False, mu=False, shear=0, rhomin=1e-2, zmin=1, vel_centered=False, cmap='viridis'):
    '''
    This function plots the heatmap of FARGO3D simulation
    of the density and x,y,z-velocity for the gas or dust.
    Input
    -----
    fig: matplotlib.object
    ax: matplotlib.object
    number: int, framenumber
    dt: float, timestep between frames
    fname: str, directory name to FAGRO3D files (exmpl: /folder/name/to/simulations/gas or /fname/)
    x,z: array, coordinates simulations box
    '''
    vn = []
    for i in range(3):
        if not isinstance(v[i], (np.ndarray)):
            vn.append(v[i]*np.ones_like(coord.z))
        else:
            vn.append(v[i])
    v = np.array(vn)

    z = coord.z
    domain = (z > -zmin) & (z < zmin)
    meshX, meshZ = np.meshgrid(coord.x, z, sparse=False, indexing='xy')
    keplershear =  1.5*meshX

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
    if mu:
        dens_mu = dens/pf.Fluids[0].dens[:,0,:]
    elif norm:
        dens = np.where(dens < 1e-8, np.nan, dens)
        dens_r =np.nanmean(dens, axis=1)
        dens = dens / dens_r[:,None]

    rho = ax[0].pcolormesh(coord.x, z[domain], dens[domain,:], cmap=cmap, norm=mpl.colors.LogNorm(vmin=float(max(rhomin, np.nanmin(dens[domain,:]))), vmax=float(np.nanmax(dens[domain,:]))),shading='auto', rasterized=True)
    ax[0].set_aspect('equal')
    xmin, xmax = ax[0].get_xlim()
    xmin_tick = round_to_nearest(xmin, 0.05)
    xmax_tick = round_to_nearest(xmax, 0.05)
    ax[0].set_xticks([xmin_tick, xmax_tick])
    if np.min(dens[domain,:]) < rhomin:
        cbar = fig.colorbar(rho, ax=ax[0],extend='min')
    else:    
        cbar = fig.colorbar(rho, ax=ax[0])   
    if fluid=='gas' or fluid==0:
        if norm:
            cbar.set_label(r'$\rho_{\mathrm{g}}/\left< \rho_{\mathrm{g}} \right>_\mathrm{R}$')
        else:
            cbar.set_label(r'$\rho_{\mathrm{g}}$')
    else:
        if norm:
            cbar.set_label(r'$\rho_{\mathrm{d}}/\left< \rho_{\mathrm{d}} \right>_\mathrm{R}$')
        else:
            cbar.set_label(r'$\rho_{\mathrm{d}}$')

    if mu:
        rho_mu = ax[1].pcolormesh(coord.x, z[domain], dens_mu[domain,:], cmap=cmap, norm=mpl.colors.LogNorm(vmin=float(max(rhomin, np.nanmin(dens_mu[domain,:]))), vmax=float(np.nanmax(dens_mu[domain,:]))),shading='auto', rasterized=True)
        ax[1].set_aspect('equal')
        xmin, xmax = ax[1].get_xlim()
        xmin_tick = round_to_nearest(xmin, 0.05)
        xmax_tick = round_to_nearest(xmax, 0.05)
        ax[1].set_xticks([xmin_tick,  xmax_tick])
        if np.min(dens_mu[domain,:]) < rhomin:
            cbar = fig.colorbar(rho_mu, ax=ax[1],extend='min')
        else:    
            cbar = fig.colorbar(rho_mu, ax=ax[1])  
        cbar.set_label(r'$\mu$')

    if fluid =='sum':
        velx = pf.dust_velx()[:,0,:]-v[0][:,None]
        vely = pf.dust_vely()[:,0,:]+keplershear-shear*meshZ-v[1][:,None]
        velz = pf.dust_velz()[:,0,:]-v[2][:,None]
        vel = np.sqrt(velx**2 + vely**2 + velz**2)
    elif fluid=='gas' or fluid==0:
        velx= pf.Fluids[0].velx[:,0,:]-v[0][:,None]
        vely = pf.Fluids[0].vely[:,0,:]+keplershear-shear*meshZ-v[1][:,None]
        velz = pf.Fluids[0].velz[:,0,:]-v[2][:,None]
        vel = np.sqrt(velx**2 + vely**2 + velz**2)
    elif fluid<0 or fluid == 'tracer':
        velx= pf.Fluids[abs(fluid)+pf.n_dust].velx[:,0,:]-v[0][:,None]
        vely = pf.Fluids[abs(fluid)+pf.n_dust].vely[:,0,:]+keplershear-shear*meshZ-v[1][:,None]
        velz = pf.Fluids[abs(fluid)+pf.n_dust].velz[:,0,:]-v[2][:,None]
        vel = np.sqrt(velx**2 + vely**2 + velz**2)
    else:
        velx= pf.Fluids[fluid].velx[:,0,:]-v[0][:,None]
        vely = pf.Fluids[fluid].vely[:,0,:]+keplershear-shear*meshZ-v[1][:,None]
        velz = pf.Fluids[fluid].velz[:,0,:]-v[2][:,None]
        vel = np.sqrt(velx**2 + vely**2 + velz**2)


    if norm or mu:
        v = ax[2].pcolormesh(coord.x, z[domain], vel[domain,:], cmap=plt.cm.gnuplot, norm=colors.Normalize(),shading='auto',  rasterized=True)
        ax[2].set_aspect('equal')
        xmin, xmax = ax[2].get_xlim()
        xmin_tick = round_to_nearest(xmin, 0.05)
        xmax_tick = round_to_nearest(xmax, 0.05)
        ax[2].set_xticks([xmin_tick, xmax_tick])
        cbar = fig.colorbar(v, ax=ax[2])
        if fluid=='gas' or fluid==0:
            cbar.set_label(r'$v_{\mathrm{g}}-v_{\mathrm{g}}^0$')
        else:
            cbar.set_label(r'$v_{\mathrm{d}}-v_{\mathrm{d}}^0$')

    else:
        if vel_centered:
            if np.max(velx[domain,:]) <= 0:
                v_x = ax[1].pcolormesh(coord.x, z[domain], velx[domain,:], cmap= blue, norm=mpl.colors.Normalize(vmin=np.min(velx[domain,:]), vmax=0),rasterized=True)
            elif np.min(velx[domain,:]) >= 0:
                v_x = ax[1].pcolormesh(coord.x, z[domain], velx[domain,:], cmap= red, norm=mpl.colors.Normalize(vmin=0, vmax=np.max(velx[domain,:])),rasterized=True)
            else:
                v_x = ax[1].pcolormesh(coord.x, z[domain], velx[domain,:], cmap='seismic', norm=mpl.colors.TwoSlopeNorm(vcenter=0, vmin=np.min(velx[domain,:]), vmax=np.max(velx[domain,:])), rasterized=True)
        else:
            if np.max(velx[domain,:]) <= 0:
                v_x = ax[1].pcolormesh(coord.x, z[domain], velx[domain,:], cmap= blue, norm=mpl.colors.Normalize(vmin=np.min(velx[domain,:]), vmax=np.max(velx[domain,:])),rasterized=True)
            elif np.min(velx[domain,:]) >= 0:
                v_x = ax[1].pcolormesh(coord.x, z[domain], velx[domain,:], cmap= red, norm=mpl.colors.Normalize(vmin=np.min(velx[domain,:]), vmax=np.max(velx[domain,:])),rasterized=True)
            else:
                v_x = ax[1].pcolormesh(coord.x, z[domain], velx[domain,:], cmap='seismic', norm=mpl.colors.Normalize(vmin=np.min(velx[domain,:]), vmax=np.max(velx[domain,:])), rasterized=True)
        ax[1].set_aspect('equal')
        xmin, xmax = ax[1].get_xlim()
        xmin_tick = round_to_nearest(xmin, 0.05)
        xmax_tick = round_to_nearest(xmax, 0.05)
        ax[1].set_xticks([xmin_tick, xmax_tick])
        cbar = fig.colorbar(v_x, ax=ax[1])
        if fluid=='gas' or fluid==0:
            cbar.set_label(r'$v_{\mathrm{g},x}-v_{\mathrm{g},x}^0$')
        else:
            cbar.set_label(r'$v_{\mathrm{d},x}-v_{\mathrm{d},x}^0$')

        if fluid =='sum':
            vely = pf.dust_vely()[:,0,:]-v[1][:,None]
        elif fluid=='gas' or fluid==0:
            vely= pf.Fluids[0].vely[:,0,:]-v[1][:,None]
        elif fluid<0 or fluid == 'tracer':
            vely= pf.Fluids[abs(fluid)+pf.n_dust].vely[:,0,:]-v[1][:,None]
        else:
            vely= pf.Fluids[fluid].vely[:,0,:]-v[1][:,None]
        vely += keplershear
        vely -= shear*meshZ
        if vel_centered:
            if np.max(vely[domain,:]) <= 0:
                v_y = ax[2].pcolormesh(coord.x, z[domain], vely[domain,:], cmap= blue, norm=mpl.colors.Normalize(vmin=np.min(vely[domain,:]), vmax=0),rasterized=True)
            elif np.min(vely[domain,:]) >= 0:
                v_y = ax[2].pcolormesh(coord.x, z[domain], vely[domain,:], cmap= red, norm=mpl.colors.Normalize(vmin=0, vmax=np.max(vely[domain,:])),rasterized=True)
            else:
                v_y = ax[2].pcolormesh(coord.x, z[domain], vely[domain,:], cmap='seismic', norm=mpl.colors.TwoSlopeNorm(vcenter=0, vmin=np.min(vely[domain,:]), vmax=np.max(vely[domain,:])), rasterized=True)
        else:
            if np.max(vely[domain,:]) <= 0:
                v_y = ax[2].pcolormesh(coord.x, z[domain], vely[domain,:], cmap= blue, norm=mpl.colors.Normalize(vmin=np.min(vely[domain,:]), vmax=np.max(vely[domain,:])),rasterized=True)
            elif np.min(vely[domain,:]) >= 0:
                v_y = ax[2].pcolormesh(coord.x, z[domain], vely[domain,:], cmap= red, norm=mpl.colors.Normalize(vmin=np.min(vely[domain,:]), vmax=np.max(vely[domain,:])),rasterized=True)
            else:
                v_y = ax[2].pcolormesh(coord.x, z[domain], vely[domain,:], cmap='seismic', norm=mpl.colors.Normalize(vmin=np.min(vely[domain,:]), vmax=np.max(vely[domain,:])), rasterized=True)
        ax[2].set_aspect('equal')
        xmin, xmax = ax[2].get_xlim()
        xmin_tick = round_to_nearest(xmin, 0.05)
        xmax_tick = round_to_nearest(xmax, 0.05)
        ax[2].set_xticks([xmin_tick, xmax_tick])
        cbar = fig.colorbar(v_y, ax=ax[2])
        if fluid=='gas' or fluid==0:
            cbar.set_label(r'$v_{\mathrm{g},y} - v_{\mathrm{g},y}^0$')
        else:
            cbar.set_label(r'$v_{\mathrm{g},d} - v_{\mathrm{g},d}^0$')

        # velz   
        if fluid =='sum':
            velz = pf.dust_velz()[:,0,:]-v[2][:,None]
        elif fluid=='gas' or fluid==0:
            velz= pf.Fluids[0].velz[:,0,:]-v[2][:,None]
        elif fluid<0 or fluid == 'tracer':
            velz= pf.Fluids[abs(fluid)+pf.n_dust].velz[:,0,:]-v[2][:,None]
        else:
            velz= pf.Fluids[fluid].velz[:,0,:]-v[2][:,None]
        if vel_centered:
            if np.max(velz[domain,:]) <= 0:
                v_z = ax[3].pcolormesh(coord.x, z[domain], velz[domain,:], cmap= blue, norm=mpl.colors.Normalize(vmin=np.min(velz[domain,:]), vmax=0),rasterized=True)
            elif np.min(velz[domain,:]) >= 0:
                v_z = ax[3].pcolormesh(coord.x, z[domain], velz[domain,:], cmap= red, norm=mpl.colors.Normalize(vmin=0, vmax=np.max(velz[domain,:])),rasterized=True)
            else:
                v_z = ax[3].pcolormesh(coord.x, z[domain], velz[domain,:], cmap='seismic', norm=mpl.colors.TwoSlopeNorm(vcenter=0, vmin=np.min(velz[domain,:]), vmax=np.max(velz[domain,:])), rasterized=True)
        else:
            if np.max(velz[domain,:]) <= 0:
                v_z = ax[3].pcolormesh(coord.x, z[domain], velz[domain,:], cmap= blue, norm=mpl.colors.Normalize(vmin=np.min(velz[domain,:]), vmax=np.max(velz[domain,:])),rasterized=True)
            elif np.min(velz[domain,:]) >= 0:
                v_z = ax[3].pcolormesh(coord.x, z[domain], velz[domain,:], cmap= red, norm=mpl.colors.Normalize(vmin=np.min(velz[domain,:]), vmax=np.max(velz[domain,:])),rasterized=True)
            else:
                v_z = ax[3].pcolormesh(coord.x, z[domain], velz[domain,:], cmap='seismic', norm=mpl.colors.Normalize(vmin=np.min(velz[domain,:]), vmax=np.max(velz[domain,:])), rasterized=True)
        ax[3].set_aspect('equal')
        xmin, xmax = ax[3].get_xlim()
        xmin_tick = round_to_nearest(xmin, 0.05)
        xmax_tick = round_to_nearest(xmax, 0.05)
        ax[3].set_xticks([xmin_tick, xmax_tick])
        cbar = fig.colorbar(v_z, ax=ax[3])
        if fluid=='gas'or fluid==0:
            cbar.set_label(r'$v_{\mathrm{g},z} - v_{\mathrm{g},z}^0$')
        else:
            cbar.set_label(r'$v_{\mathrm{d},z} - v_{\mathrm{d},z}^0$')

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel(r'x $[H_\mathrm{g}]$')
    plt.ylabel(r'z $[H_\mathrm{g}]$')
    if fluid =='gas' or fluid==0:
        plt.suptitle(r'Gas: $\Omega t:${:.2f}'.format(dt*number),fontsize=16)   
    elif fluid=='sum':
        plt.suptitle(r'Sum dust: $\Omega t:${:.2f}'.format(dt*number),fontsize=14)
    elif fluid<0 or fluid == 'tracer':
        plt.suptitle(r'Tracer{}: $\Omega t:${:.2f}'.format(abs(fluid),dt*number),fontsize=14)
    elif fluid > pf.n_dust:
        plt.suptitle(r'Tracer{}: $\Omega t:${:.2f}'.format(fluid-pf.n_dust,dt*number),fontsize=14)
    else:
        plt.suptitle(r'Dust{}: $\Omega t:${:.2f}'.format(fluid,dt*number),fontsize=14)
    return fig, ax

def radial_avg(fluid, fig, ax, number, dt, pf, coord, mu=False, v=np.array([0., 0., 0.]), shear=0, rhocrit=1e-6, scaleheight=True, color='tab:blue'):
    '''
    This function plots the radial average values of FARGO3D simulation
    of the density and x,y,z-velocity for the gas or dust.
    Input
    -----
    fig: matplotlib.object
    ax: matplotlib.object
    number: int, framenumber
    dt: float, timestep between frames
    fname: str, directory name to FAGRO3D files (exmpl: /folder/name/to/simulations/gas or /fname/dust1)
    x,z: array, coordinates simulations box
    '''

    vn = []
    for i in range(3):
        if not isinstance(v[i], (np.ndarray)):
            vn.append(v[i]*np.ones_like(coord.z))
        else:
            vn.append(v[i])
    v = np.array(vn)

    z = coord.z
    meshX, meshZ = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')
    keplershear =  1.5*meshX
    qz = shear*meshZ
    if fluid == 'tracer':
        fluid = -1
    pf.read(number)
    if fluid =='sum':
        dens = pf.dust_density()[:,0,:]
    elif fluid=='gas' or fluid==0:
        dens = pf.Fluids[0].dens[:,0,:]
    elif fluid < 0:
        dens = pf.Fluids[abs(fluid)+pf.n_dust].dens[:,0,:]
    else:
        dens = pf.Fluids[fluid].dens[:,0,:]
    dens = dens.mean(1)
    domain = dens>rhocrit
    if mu:
        dens_mu = dens/np.mean(pf.Fluids[0].dens[:,0,:],axis=1)
        domain = dens_mu>rhocrit

    try:
        fwhm_min = np.min(z[dens > np.max(dens)/2])
        fwhm_max = np.max(z[dens > np.max(dens)/2])
        parameters, covariance = curve_fit(f.gauss, z, dens)
    except:
        fwhm_min, fwhm_max = 0
        parameters = [0,0]
    if parameters[1] < 0: parameters[1] = np.nan
    
    h = ([parameters[1], (fwhm_max-fwhm_min)/2])

    ax[0].plot(z[domain], dens[domain],color=color)
    if scaleheight:
        ax[0].vlines([fwhm_min, fwhm_max], 0, max(dens[domain]*1.01), color='k', alpha=.75, linestyle='dashed')
        if (parameters[1] > 0) & (parameters[1] < np.max(z[domain])):
            ax[0].vlines([-parameters[1], +parameters[1]], 0, max(dens[domain]*1.01), color='k', alpha=.75,  linestyle='dotted')
    if fluid=='gas' or fluid==0:
        ax[0].set_ylabel(r'$\rho_\mathrm{g}$')
    else:
        ax[0].set_ylabel(r'$\rho_\mathrm{d}$')

    if mu:
        try:
            fwhm_min = np.min(z[dens_mu > np.max(dens_mu)/2])
            fwhm_max = np.max(z[dens_mu > np.max(dens_mu)/2])
            parameters, covariance = curve_fit(f.gauss, z, dens_mu)
        except:
            fwhm_min, fwhm_max = 0
            parameters = [0,0]
        if parameters[1] < 0: parameters[1] = np.nan

        ax[1].plot(z[domain], dens_mu[domain],color=color)
        if scaleheight:
            ax[1].vlines([fwhm_min, fwhm_max], 0, max(dens_mu[domain]*1.01), color='k', linestyle='dashed')
            if (parameters[1] > 0) & (parameters[1] < np.max(z[domain])):
                ax[1].vlines([-parameters[1], +parameters[1]], 0, max(dens_mu[domain]*1.01), color='k', linestyle='dotted')
            else:
                parameters[1] = np.nan
        ax[1].set_ylabel(r'$\mu$')
        h.append(parameters[1])
        h.append((fwhm_max-fwhm_min)/2)


    if mu is False:
        if fluid =='sum':
            velx = pf.dust_velx()[:,0,:]-v[0][:,None]
        elif fluid=='gas' or fluid==0:
            velx= pf.Fluids[0].velx[:,0,:]-v[0][:,None]
        else:
            velx= pf.Fluids[fluid].velx[:,0,:]-v[0][:,None]
        velx = velx.mean(1)

        ax[1].plot(z[domain], velx[domain],color=color)
        ax[1].set_ylabel(r'$v_{\mathrm{g},x} - v_{\mathrm{g},x}^0$')

        if fluid =='sum':
            vely = pf.dust_vely()[:,0,:]-v[1][:,None]
        elif fluid=='gas' or fluid==0:
            vely= pf.Fluids[0].vely[:,0,:]-v[1][:,None]
        else:
            vely= pf.Fluids[fluid].vely[:,0,:]-v[1][:,None]
        vely += keplershear
        vely -= qz
        vely = vely.mean(1)

        ax[2].plot(z[domain], vely[domain],color=color)
        ax[2].set_ylabel(r'$v_{\mathrm{g},y} - v_{\mathrm{g},y}^0$')


        if fluid =='sum':
            velz = pf.dust_velz()[:,0,:]-v[2][:,None]
        elif fluid=='gas' or fluid==0:
            velz= pf.Fluids[0].velz[:,0,:]-v[2][:,None]
        else:
            velz= pf.Fluids[fluid].velz[:,0,:]-v[2][:,None]
        velz = velz.mean(1)

        ax[3].plot(z[domain], velz[domain], color=color)
        ax[3].set_ylabel(r'$v_{\mathrm{g},z} - v_{\mathrm{g},z}^0$')

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    if fluid =='gas' or fluid==0:
        plt.suptitle(r'Gas: $\Omega t:${:.2f}'.format(dt*number),fontsize=14)   
    elif fluid=='sum':
        plt.suptitle(r'Sum dust: $\Omega t:${:.2f}'.format(dt*number),fontsize=14)
    elif fluid<0:
        plt.suptitle(r'Tracer{}: $\Omega t:${:.2f}'.format(abs(fluid),dt*number),fontsize=14)
    elif fluid > pf.n_dust:
        plt.suptitle(r'Tracer{}: $\Omega t:${:.2f}'.format(fluid-pf.n_dust,dt*number),fontsize=14)
    else:
        plt.suptitle(r'Dust{}: $\Omega t:${:.2f}'.format(fluid,dt*number),fontsize=14) 
    plt.xlabel(r'z $[H_\mathrm{g}]$')
    return fig, ax, h

def amplitude(fig, ax, frame, pf, coord, dt, direc=None, v=np.array([0,0,0]), shear=0, scale='log', postitive=True, pmin=1e-6):
    vn = []
    for i in range(3):
        if not isinstance(v[i], (np.ndarray)):
            vn.append(v[i]*np.ones_like(coord.z))
        else:
            vn.append(v[i])
    v = np.array(vn)
    

    # 2D FFT (properly normalized)
    X, Z = np.meshgrid(coord.x, coord.z, sparse=False, indexing='xy')
    pf.read(frame)
    velx = pf.Fluids[0].velx[:,0,:]-v[0][:,None]
    vely = pf.Fluids[0].vely[:,0,:]-v[1][:,None]+1.5*X-shear*Z
    velz = pf.Fluids[0].velz[:,0,:]-v[2][:,None]

    if direc is None:
        vel = np.sqrt(velx**2+vely**2+velz**2)
    elif direc == 'x':
        vel = velx
    elif direc == 'y':
        vel = vely
    elif direc == 'z':
        vel = velz
    
    F = np.fft.fftshift(np.fft.fft2(vel))
    F_mag = np.abs(F) / (coord.nx * coord.nz)

    # freq centers (angular frequency)
    kx_c = np.fft.fftshift(np.fft.fftfreq(coord.nx, d=coord.dx)) * 2*np.pi
    kz_c = np.fft.fftshift(np.fft.fftfreq(coord.nz, d=coord.dz)) * 2*np.pi

    kx_e = cell_edges(kx_c)
    kz_e = cell_edges(kz_c)

    if direc is None:
        vcm = ax[0].pcolormesh(X,Z,vel, cmap='gnuplot', norm=mpl.colors.Normalize(vmin=0, vmax=np.max(vel)),rasterized=True)
        label = r'$v/v^0$'
    else:
        if np.max(vel) <= 0:
            vcm = ax[0].pcolormesh(X, Z, vel, cmap= blue, norm=mpl.colors.Normalize(vmin=np.min(vel), vmax=np.max(vel)),rasterized=True)
        elif np.min(vel) >= 0:
            vcm = ax[0].pcolormesh(X, Z, vel, cmap= red, norm=mpl.colors.Normalize(vmin=np.min(vel), vmax=np.max(vel)),rasterized=True)
        else:
            vcm = ax[0].pcolormesh(X, Z, vel, cmap='seismic', norm=mpl.colors.TwoSlopeNorm(vcenter=0, vmin=np.min(vel), vmax=np.max(vel)), rasterized=True)

        label = r'$v_{'+direc+r'}/v^0$'
    ax[0].set_xlabel(r"$x \, [H_\mathrm{g}]$")
    ax[0].set_ylabel(r"$z \, [H_\mathrm{g}]$")
    xmin, xmax = ax[0].get_xlim()
    xmin_tick = round_to_nearest(xmin, 0.05)
    xmax_tick = round_to_nearest(xmax, 0.05)
    ax[0].set_xticks([xmin_tick, xmax_tick])
    ax[0].set_aspect('equal')
    cbar = fig.colorbar(vcm, ax=ax[0],label=label)

    # plot aplitude
    pcm = ax[1].pcolormesh(kx_e, kz_e, F_mag, shading='flat', norm=mpl.colors.LogNorm(vmin=float(max(pmin, np.nanmin(F_mag))), vmax=float(np.nanmax(F_mag))), rasterized=True)  # note transpose to match edges (x along axis 0)

    # Vertical/horizontal lines at ±linthresh
    ax[1].axvline(np.min(kx_e[kx_e>0]), color='w', linestyle='--', linewidth=0.75)
    ax[1].axvline(np.max(kx_e[kx_e<0]), color='w', linestyle='--', linewidth=0.75)

    ax[1].axhline(np.min(kz_e[kz_e>0]), color='w', linestyle='--', linewidth=0.75)
    ax[1].axhline(np.max(kz_e[kz_e<0]), color='w', linestyle='--', linewidth=0.75)
    
    # Outline the actual (possibly skewed) quad for the max cell
    i_max, j_max = np.unravel_index(np.argmax(F_mag), F_mag.shape)

    poly = Polygon(cell_corners(kx_e, kz_e, i_max, j_max),
               closed=True, fill=False, edgecolor='red', linewidth=1)
    ax[1].add_patch(poly)
    if direc is None:
        f_mag2 = F_mag.copy()
        f_mag2[i_max, j_max] = 0
        # Outline the actual (possibly skewed) quad for the max cell
        i_max2, j_max2 = np.unravel_index(np.argmax(f_mag2), f_mag2.shape)

        poly2 = Polygon(cell_corners(kx_e, kz_e, i_max2, j_max2),
                closed=True, fill=False, edgecolor='orange', linewidth=1)
        ax[1].add_patch(poly2)
    
    if scale == 'log':
        linth_x = np.min(kx_e[kx_e > 0])
        linth_z = np.min(kz_e[kz_e > 0])
        ax[1].set_xscale('symlog', linthresh=linth_x, linscale=.1)
        ax[1].set_yscale('symlog', linthresh=linth_z, linscale=.1)
        if np.max(kx_e)<1e3:
            ax[1].set_xticks([-1e2, -1e1,0,1e1, 1e2])
            ax[1].set_yticks([-1e2, -1e1,-1e0, 0,1e0,1e1, 1e2])
        else:
            ax[1].set_xticks([-1e3, -1e2, -1e1, 0, 1e1, 1e2, 1e3])
            ax[1].set_yticks([-1e3, -1e2, -1e1, -1e0, 0, 1e0, 1e1, 1e2, 1e3])


        # Reduce crowding of tick labels
    ax[1].set_xlabel(r"$k_x$ [rad / $H_\mathrm{g}$]")
    ax[1].set_ylabel(r"$k_z$ [rad / $H_\mathrm{g}$]")
    ax[1].set_box_aspect(1)

    if postitive:
        ax[1].set_xlim(0,np.max(kx_e))
        ax[1].set_ylim(0,np.max(kz_e))

    label2 = r'$|F_k$('+label+r')$| / (N_x N_z)$'
    if np.nanmin(F_mag) < pmin:
        cbar = fig.colorbar(pcm, ax=ax[1], label=label2, extend='min' )
    else:
        cbar = fig.colorbar(pcm, ax=ax[1], label=label2)

    # Compute radial spectrum
    k_centers, Ek, Nct = f.radial_power_spectrum(frame, pf, coord, direc, v, shear)

    # Plot
    ax[2].loglog(k_centers, Ek, marker='.')
    ax[2].set_xlabel(r'$k$ [rad / $H_\mathrm{g}$]')
    if direc is None:
        label = r'$E(k)$'
    else:
        label = r'$E_{'+direc+r'}(k)$'
    ax[2].set_ylabel(label)
    ax[2].set_box_aspect(1/3)


    if direc is None:
        fig.suptitle(r'Gas Velocity: $\Omega t:${:.2f}'.format(dt*frame),fontsize=14)
    elif direc == 'x':
        fig.suptitle(r'Gas Radial Velocity: $\Omega t:${:.2f}'.format(dt*frame),fontsize=14)
    elif direc == 'y':
        fig.suptitle(r'Gas Azimuthal Velocity: $\Omega t:${:.2f}'.format(dt*frame),fontsize=14)
    elif direc == 'z':
        fig.suptitle(r'Gas Vertical Velocity: $\Omega t:${:.2f}'.format(dt*frame),fontsize=14)
    return fig, ax



def round_to_nearest(val, base=0.05):
    return round(val / base) * base

# --- build cell edges from centers so pcolormesh draws proper pixels ---
def cell_edges(grid):
    mid = 0.5 * (grid[1:] + grid[:-1])
    left = grid[0] - (mid[0] - grid[0])
    right = grid[-1] + (grid[-1] - mid[-1])
    return np.concatenate(([left], mid, [right]))

def cell_corners(X, Y, i, j):
    """
    Return the 4 vertices of cell (i, j) in data coords, counterclockwise.
    Works with 1D edges or 2D grids.
    """
    X = np.asanyarray(X); Y = np.asanyarray(Y)
    if X.ndim == 1 and Y.ndim == 1:  # 1D edges
        return [(X[j],   Y[i]),
                (X[j+1], Y[i]),
                (X[j+1], Y[i+1]),
                (X[j],   Y[i+1])]
    else:  # 2D grids
        return [(X[i, j],     Y[i, j]),
                (X[i, j+1],   Y[i, j+1]),
                (X[i+1, j+1], Y[i+1, j+1]),
                (X[i+1, j],   Y[i+1, j])]