import numpy as np
import pandas as pd

def mean_z(interval, pf, coord):
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

    z = coord.z
    x = coord.x
    
    meshX, meshZ = np.meshgrid(x, z, sparse=False, indexing='xy')
    shear =  1.5*meshX
    nghz = 3
    dict = {}
    for number in range(interval[0], interval[1]):
        pf.read(number)
        dens = pf.Fluids[0].dens[:,0,:]
        dens = dens.mean(1)
        velx= pf.Fluids[0].velx[:,0,:]
        velx = velx.mean(1)
        vely= pf.Fluids[0].vely[:,0,:]
        vely += shear
        vely = vely.mean(1)
        velz= pf.Fluids[0].velz[:,0,:]
        velz = velz.mean(1)
        if number == interval[0]:
            dict['rho_g'] = np.pad(dens, nghz, mode='edge')
            dict['vx_g'] = np.pad(velz, nghz, mode='edge')
            dict['vy_g'] = np.pad(vely, nghz, mode='edge')
            dict['vz_g'] = np.pad(velz, nghz, mode='edge')
            dict['ax_g'] = np.zeros_like(np.pad(velx, nghz, mode='edge'))
            dict['ay_g'] = np.zeros_like(np.pad(vely, nghz, mode='edge'))
            dict['az_g'] = np.zeros_like(np.pad(velz, nghz, mode='edge'))
        else:
            dict['rho_g'] += np.pad(dens, nghz, mode='edge')
            dict['vx_g'] += np.pad(velx, nghz, mode='edge')
            dict['vy_g'] += np.pad(vely, nghz, mode='edge')
            dict['vz_g'] += np.pad(velz, nghz, mode='edge')
            dict['ax_g'] += np.zeros_like(np.pad(velx, nghz, mode='edge'))
            dict['ay_g'] += np.zeros_like(np.pad(vely, nghz, mode='edge'))
            dict['az_g'] += np.zeros_like(np.pad(velz, nghz, mode='edge'))

        for i in range(int(pf.n_dust)):
            dens = pf.Fluids[i+1].dens[:,0,:]
            dens = dens.mean(1)
            velx= pf.Fluids[i+1].velx[:,0,:]
            velx = velx.mean(1)
            vely= pf.Fluids[i+1].vely[:,0,:]
            vely +=shear
            vely = vely.mean(1)
            velz= pf.Fluids[i+1].velz[:,0,:]
            velz = velz.mean(1)
            if number == interval[0]:
                dict['rho_d{}'.format(i+1)] = np.pad(dens, nghz, mode='edge')
                dict['ux_d{}'.format(i+1)] = np.pad(velx, nghz, mode='edge')
                dict['uy_d{}'.format(i+1)] = np.pad(vely, nghz, mode='edge')
                dict['uz_d{}'.format(i+1)] = np.pad(velz, nghz, mode='edge')
            else:
                dict['rho_d{}'.format(i+1)] += np.pad(dens, nghz, mode='edge')
                dict['ux_d{}'.format(i+1)] += np.pad(velx, nghz, mode='edge')
                dict['uy_d{}'.format(i+1)] += np.pad(vely, nghz, mode='edge')
                dict['uz_d{}'.format(i+1)] += np.pad(velz, nghz, mode='edge')
    df = pd.DataFrame(dict)

    return df/(interval[1]-interval[0])