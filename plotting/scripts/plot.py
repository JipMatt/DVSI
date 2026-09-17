
#import libaries
import sys, os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import Divider, Size
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter




import numpy as np
from scipy.optimize import curve_fit 

dir = os.path.dirname(os.path.dirname(__file__))+'/plotting/code/'
sys.path.append(dir)
import func as f

seismic = mpl.pyplot.get_cmap("seismic", 256)
red = colors.LinearSegmentedColormap.from_list(
    "red", seismic(np.linspace(0.5, 1, 256)))
blue = colors.LinearSegmentedColormap.from_list(
    "blue", seismic(np.linspace(0, .5, 256)))

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

def add_side_profile(
    fig,
    ax_parent,
    profile,
    y,
    profile_gap=0.025,
    profile_width=0.14,
    xlabel=r"$q_z$",
    show_xlabels=False,
):
    bbox = ax_parent.get_position()
    fig_width, _ = fig.get_size_inches()

    gap_norm = profile_gap / fig_width
    width_norm = profile_width / fig_width

    ax_side = fig.add_axes(
        [
            bbox.x1 + gap_norm,
            bbox.y0,
            width_norm,
            bbox.height,
        ],
        sharey=ax_parent,
    )

    ax_side.plot(
        profile,
        y,
        color="black",
        linewidth=0.8,
    )

    ax_side.axvline(
        0.0,
        color="0.5",
        linestyle="--",
        linewidth=0.75,
        zorder=0,
    )

    # Fixed q_z range and tick positions
    ax_side.set_xlim(-0.1, 0.1)
    ax_side.set_xticks([-0.1, 0.0, 0.1])

    # No y-axis labels or ticks
    ax_side.tick_params(
        axis="y",
        left=True,
        right=True,
        labelleft=False,
        labelright=False,
    )

    # Tick marks on both top and bottom for both rows
    ax_side.tick_params(
        axis="x",
        top=True,
        bottom=True,
        labeltop=False,
        labelbottom=show_xlabels,
        labelsize=7,
        # length=2.0,
        # width=0.8,
        pad=1.5,
    )

    if show_xlabels:
        ax_side.set_xticklabels(
            [r"$-0.1$", r"$0$", r"$0.1$"]
        )
        ax_side.set_xlabel(
            xlabel,
            fontsize=9,
            labelpad=1.5,
        )
    else:
        # No axis label and no numerical labels on top row
        ax_side.set_xlabel("")
        ax_side.set_xticklabels([])

    return ax_side


def multi_heatmaps(
    dens, coordx, coordz,
    *,
    width_plot=7.08661,
    max_panels=6,

    # layout / margins (inches)
    gap_x=0.08,
    gap_y=0.10,
    left=0.33, right=0.06,
    bottom=0.32, top=0.06,

    # colorbar geometry (inches)
    cbar_h=0.09,
    cbar_pad_top=0.025,
    cbar_pad_bottom=0.12,

    # reserve extra clearance specifically for the shared xlabel in 2-row layout
    xlabel_clearance=0.10,

    # appearance
    cmap="viridis",
    rasterized=True,

    # color scaling
    cbar_mode="per_panel", #'per_panel', 'first_individual_rest_shared', 'shared'
    vmin_floor=1e-4,
    vmax_floor=1e-2,
    symmetric=False,
    use_lognorm=True,
    vmin=None,
    vmax=None,

    # labels
    xlabel=r"$x$",
    ylabel=r"$z$",
    omegat=None,
    loc_omegat=[0.45, 1],
    cbar_labels=None,
    lwidth=1,

    # per-panel corner text
    panel_text=None,
    panel_text_kwargs=None,

    # ticks
    x_tick_step=0.5,
    y_tick_step=0.5,
    x_tick_fmt="%.1f",
    y_tick_fmt="%.1f",
    show_inner_yticks=False,
    show_inner_xticks=False,

    # layout rule for the original automatic layout
    two_row_aspect_threshold=2.0,

    # layout mode
    #   "auto"            -> original behavior
    #   "3large_2small"   -> 5-panel paper layout:
    #                         [ large ][ large ]
    #                         [ large ][small][small]
    layout="auto",
):
    """
    Plot N=1..max_panels heatmaps with fixed figure width and minimal whitespace.

    Parameters
    ----------
    dens : sequence of 2D arrays
        Heatmap data. For each panel i, dens[i] should have shape
        (len(coordz_i), len(coordx_i)) when coordinates are cell centers,
        or one fewer element than coordinate lengths when coordinates are edges.

    coordx, coordz : 1D array or sequence of 1D arrays
        May be a single coordinate array shared by all panels, or a list/tuple
        of length N containing per-panel coordinates.

    layout : {"auto", "3large_2small"}
        "auto" keeps the original one-row/two-row layout logic.

        "3large_2small" requires exactly five panels and creates

            [ panel 0 large ] [ panel 1 large ]
            [ panel 2 large ] [ panel 3 small ][ panel 4 small ]

        Panel widths are derived from the *physical x-domain extents*, not from
        the number of grid cells. A single physical scale (inches per coordinate
        unit) is used for x and z, so dx == dz visually when
        ax.set_aspect("equal") is used. This is intended for the resolution /
        half-domain layout in the paper.

        For this special layout, `gap_x` is treated as the total horizontal
        gutter budget per row: the top row has one gap of `gap_x`, while the
        bottom row has two gaps of `gap_x/2`. This lets both rows line up at the
        same left/right edges while preserving exact physical panel widths.

    vmin, vmax, vmin_floor, vmax_floor
        Each may be None, a scalar, a supported string keyword, or a list/tuple.

    Supported string keywords
    -------------------------
    "auto"  : let the function determine the bound from the data
    "min"   : use the data minimum
    "max"   : use the data maximum
    "floor" : accepted for floor parameters
    """

    dens = list(dens)
    N = len(dens)
    if not (1 <= N <= max_panels):
        raise ValueError(f"Expected 1..{max_panels} heatmaps; got {N}.")

    valid_layouts = {"auto", "3large_2small"}
    if layout not in valid_layouts:
        raise ValueError(f"`layout` must be one of {valid_layouts}; got {layout!r}.")

    if layout == "3large_2small" and N != 5:
        raise ValueError("layout='3large_2small' requires exactly 5 heatmaps.")

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def _as_per_panel_coords(c, name):
        if isinstance(c, (list, tuple)) and len(c) == N:
            return list(c)
        if isinstance(c, np.ndarray) and c.dtype == object and c.ndim == 1 and len(c) == N:
            return list(c)
        return [c] * N

    def _coords_equal(a, b):
        a = np.asarray(a)
        b = np.asarray(b)
        return a.shape == b.shape and np.array_equal(a, b)

    def _broadcast_norm_param(x, name, cbar_mode, N):
        """Broadcast norm-related parameters to a list of length N."""
        if not isinstance(x, (list, tuple, np.ndarray)) or isinstance(x, str):
            return [x] * N

        if isinstance(x, np.ndarray):
            if x.ndim != 1:
                raise ValueError(f"`{name}` array must be 1D.")
            vals = list(x)
        else:
            vals = list(x)

        L = len(vals)

        if cbar_mode == "per_panel":
            if L != N:
                raise ValueError(
                    f"With cbar_mode='per_panel', `{name}` must have length N={N}; got {L}."
                )
            return vals

        if cbar_mode == "shared":
            if L == 0:
                raise ValueError(f"`{name}` list cannot be empty.")
            if L == 1:
                return vals * N
            if L == N:
                return vals
            raise ValueError(
                f"With cbar_mode='shared', `{name}` must be a scalar/str/None, "
                f"or a list of length 1{' or N' if N > 1 else ''}; got {L}."
            )

        if cbar_mode == "first_individual_rest_shared":
            if L == 0:
                raise ValueError(f"`{name}` list cannot be empty.")
            if L == 1:
                return vals * N
            if L == 2:
                return [vals[0]] + [vals[1]] * (N - 1)
            if L == N:
                return vals
            raise ValueError(
                f"With cbar_mode='first_individual_rest_shared', `{name}` must be "
                f"a scalar/str/None, or a list of length 1, 2, or N={N}; got {L}."
            )

        raise ValueError(f"Unknown cbar_mode={cbar_mode!r}")

    def _broadcast_bool_param(x, name, cbar_mode, N):
        """Broadcast boolean norm flags to a list of length N."""
        if isinstance(x, (bool, np.bool_)):
            return [bool(x)] * N

        if isinstance(x, np.ndarray):
            if x.ndim != 1:
                raise ValueError(f"`{name}` array must be 1D.")
            vals = list(x)
        elif isinstance(x, (list, tuple)):
            vals = list(x)
        else:
            raise TypeError(f"`{name}` must be a bool or a list/tuple/1D ndarray of bools.")

        if not all(isinstance(v, (bool, np.bool_)) for v in vals):
            raise TypeError(f"All entries of `{name}` must be booleans.")

        vals = [bool(v) for v in vals]
        L = len(vals)

        if cbar_mode == "per_panel":
            if L != N:
                raise ValueError(
                    f"With cbar_mode='per_panel', `{name}` must have length N={N}; got {L}."
                )
            return vals

        if cbar_mode == "shared":
            if L == 0:
                raise ValueError(f"`{name}` list cannot be empty.")
            if L == 1:
                return vals * N
            if L == N:
                return vals
            raise ValueError(
                f"With cbar_mode='shared', `{name}` must be a bool, "
                f"or a list of length 1{' or N' if N > 1 else ''}; got {L}."
            )

        if cbar_mode == "first_individual_rest_shared":
            if L == 0:
                raise ValueError(f"`{name}` list cannot be empty.")
            if L == 1:
                return vals * N
            if L == 2:
                return [vals[0]] + [vals[1]] * (N - 1)
            if L == N:
                return vals
            raise ValueError(
                f"With cbar_mode='first_individual_rest_shared', `{name}` must be a bool, "
                f"or a list of length 1, 2, or N={N}; got {L}."
            )

        raise ValueError(f"Unknown cbar_mode={cbar_mode!r}")

    def _finite_stats(d):
        d = np.asarray(d, dtype=float)
        finite = d[np.isfinite(d)]
        if finite.size == 0:
            return {
                "finite": finite,
                "finite_pos": finite,
                "has_finite": False,
                "has_pos": False,
                "min": None,
                "max": None,
                "pos_min": None,
                "pos_max": None,
            }

        finite_pos = finite[finite > 0]
        return {
            "finite": finite,
            "finite_pos": finite_pos,
            "has_finite": True,
            "has_pos": finite_pos.size > 0,
            "min": float(np.nanmin(finite)),
            "max": float(np.nanmax(finite)),
            "pos_min": float(np.nanmin(finite_pos)) if finite_pos.size else None,
            "pos_max": float(np.nanmax(finite_pos)) if finite_pos.size else None,
        }

    def _resolve_bound(value, *, stats, which, use_lognorm_i, default=None):
        if value is None:
            return default

        if isinstance(value, str):
            key = value.lower().strip()
            if key == "auto":
                return default
            if key == "min":
                if use_lognorm_i:
                    return stats["pos_min"] if stats["has_pos"] else default
                return stats["min"] if stats["has_finite"] else default
            if key == "max":
                if use_lognorm_i:
                    return stats["pos_max"] if stats["has_pos"] else default
                return stats["max"] if stats["has_finite"] else default
            raise ValueError(
                f"Unsupported string for `{which}`: {value!r}. "
                f"Use None, number, 'auto', 'min', or 'max'."
            )

        return float(value)

    def _resolve_floor(value, *, stats, which, use_lognorm_i, default=None):
        if value is None:
            return None

        if isinstance(value, str):
            key = value.lower().strip()
            if key in {"auto", "floor"}:
                return default
            if key == "min":
                if use_lognorm_i:
                    return stats["pos_min"] if stats["has_pos"] else default
                return stats["min"] if stats["has_finite"] else default
            if key == "max":
                if use_lognorm_i:
                    return stats["pos_max"] if stats["has_pos"] else default
                return stats["max"] if stats["has_finite"] else default
            raise ValueError(
                f"Unsupported string for `{which}`: {value!r}. "
                f"Use None, number, 'auto', 'floor', 'min', or 'max'."
            )

        return float(value)

    def make_norm_for_data(
        d,
        *,
        vmin_i,
        vmax_i,
        vmin_floor_i,
        vmax_floor_i,
        symmetric_i,
        use_lognorm_i,
    ):
        d = np.asarray(d, dtype=float)
        stats = _finite_stats(d)

        if symmetric_i and use_lognorm_i:
            raise ValueError(
                "`symmetric=True` and `use_lognorm=True` are incompatible for the same panel."
            )

        if symmetric_i:
            if not stats["has_finite"]:
                vv_auto = 1.0
            else:
                vv_auto = float(np.nanmax(np.abs(stats["finite"])))

            vv = _resolve_bound(
                vmax_i,
                stats=stats,
                which="vmax",
                use_lognorm_i=False,
                default=vv_auto,
            )
            v0 = _resolve_bound(
                vmin_i,
                stats=stats,
                which="vmin",
                use_lognorm_i=False,
                default=-vv,
            )

            if v0 is None:
                v0 = -vv
            if vv is None:
                vv = abs(v0)

            return mpl.colors.Normalize(vmin=float(v0), vmax=float(vv)), "neither"

        if use_lognorm_i:
            if not stats["has_pos"]:
                raise ValueError("LogNorm requested but data has no positive finite values.")

            dmin = stats["pos_min"]
            dmax = stats["pos_max"]

            vmax_floor_res = _resolve_floor(
                vmax_floor_i,
                stats=stats,
                which="vmax_floor",
                use_lognorm_i=True,
                default=vmax_floor_i,
            )
            vmin_floor_res = _resolve_floor(
                vmin_floor_i,
                stats=stats,
                which="vmin_floor",
                use_lognorm_i=True,
                default=vmin_floor_i,
            )

            vv = _resolve_bound(
                vmax_i,
                stats=stats,
                which="vmax",
                use_lognorm_i=True,
                default=dmax,
            )
            if vmax_floor_res is not None:
                vv = max(vv, float(vmax_floor_res))

            if vmin_i is None or (
                isinstance(vmin_i, str) and vmin_i.lower().strip() == "auto"
            ):
                if vmin_floor_res is None:
                    v0 = dmin
                else:
                    v0 = max(float(vmin_floor_res), dmin)
            else:
                v0 = _resolve_bound(
                    vmin_i,
                    stats=stats,
                    which="vmin",
                    use_lognorm_i=True,
                    default=dmin,
                )

            below = dmin < v0
            above = dmax > vv
            if below and above:
                ext = "both"
            elif below:
                ext = "min"
            elif above:
                ext = "max"
            else:
                ext = "neither"

            v0 = max(float(v0), np.finfo(float).tiny)
            vv = max(float(vv), v0 * (1 + 1e-12))

            return mpl.colors.LogNorm(vmin=v0, vmax=vv), ext

        if not stats["has_finite"]:
            v0 = _resolve_bound(
                vmin_i,
                stats=stats,
                which="vmin",
                use_lognorm_i=False,
                default=0.0,
            )
            vv = _resolve_bound(
                vmax_i,
                stats=stats,
                which="vmax",
                use_lognorm_i=False,
                default=1.0,
            )
            return mpl.colors.Normalize(vmin=float(v0), vmax=float(vv)), "neither"

        v0 = _resolve_bound(
            vmin_i,
            stats=stats,
            which="vmin",
            use_lognorm_i=False,
            default=stats["min"],
        )
        vv = _resolve_bound(
            vmax_i,
            stats=stats,
            which="vmax",
            use_lognorm_i=False,
            default=stats["max"],
        )

        if v0 == vv:
            vv = v0 + 1e-12

        return mpl.colors.Normalize(vmin=float(v0), vmax=float(vv)), "neither"

    def _setup_ticks(ax):
        ax.xaxis.set_major_locator(MultipleLocator(x_tick_step))
        ax.yaxis.set_major_locator(MultipleLocator(y_tick_step))
        ax.xaxis.set_major_formatter(FormatStrFormatter(x_tick_fmt))
        ax.yaxis.set_major_formatter(FormatStrFormatter(y_tick_fmt))

    def add_panel_text(ax, text):
        if text is None or text == "":
            return
        kw = dict(panel_text_kwargs)
        kw["transform"] = ax.transAxes
        ax.text(0.02, 0.98, text, **kw)

    def _plot_span(coord, n_cells, name):
        """
        Physical span represented by coordinates.

        If len(coord) == n_cells, coordinates are treated as cell centers and
        half a cell is added at both ends. If len(coord) == n_cells + 1,
        coordinates are treated as cell edges.
        """
        c = np.asarray(coord, dtype=float)
        if c.ndim != 1:
            raise ValueError(f"{name} must be 1D.")
        if c.size < 2:
            raise ValueError(f"{name} must contain at least two coordinates.")

        if c.size == n_cells + 1:
            return float(abs(c[-1] - c[0]))

        if c.size == n_cells:
            left_half = 0.5 * abs(c[1] - c[0])
            right_half = 0.5 * abs(c[-1] - c[-2])
            return float(abs(c[-1] - c[0]) + left_half + right_half)

        # pcolormesh will generally reject incompatible dimensions too, but this
        # message is more useful for the special physical-layout calculation.
        raise ValueError(
            f"{name} has length {c.size}, but the corresponding data dimension "
            f"has {n_cells} cells. Expected {n_cells} centers or {n_cells + 1} edges."
        )

    # ------------------------------------------------------------------
    # coords
    # ------------------------------------------------------------------
    coordx_list = _as_per_panel_coords(coordx, "coordx")
    coordz_list = _as_per_panel_coords(coordz, "coordz")

    can_share = True
    for i in range(1, N):
        if (
            not _coords_equal(coordx_list[0], coordx_list[i])
            or not _coords_equal(coordz_list[0], coordz_list[i])
        ):
            can_share = False
            break

    # ------------------------------------------------------------------
    # colorbar mode
    # ------------------------------------------------------------------
    valid_cbar_modes = {"per_panel", "shared", "first_individual_rest_shared"}
    if cbar_mode not in valid_cbar_modes:
        raise ValueError(f"`cbar_mode` must be one of {valid_cbar_modes}; got {cbar_mode!r}")

    # ------------------------------------------------------------------
    # colormaps
    # ------------------------------------------------------------------
    if isinstance(cmap, (list, tuple)):
        cmaps = list(cmap)

        if len(cmaps) == N:
            pass
        elif cbar_mode == "first_individual_rest_shared" and len(cmaps) == 2:
            cmaps = [cmaps[0]] + [cmaps[1]] * (N - 1)
        elif cbar_mode == "shared" and len(cmaps) == 1:
            cmaps = cmaps * N
        else:
            raise ValueError(
                f"If `cmap` is a list/tuple, it must have length N={N}"
                f"{' or length 2 when cbar_mode=first_individual_rest_shared' if cbar_mode == 'first_individual_rest_shared' else ''}; "
                f"got {len(cmaps)}."
            )
    else:
        cmaps = [cmap] * N

    # ------------------------------------------------------------------
    # panel text
    # ------------------------------------------------------------------
    if panel_text is None:
        panel_texts = [None] * N
    elif isinstance(panel_text, str):
        panel_texts = [panel_text] * N
    else:
        panel_texts = list(panel_text)
        if len(panel_texts) != N:
            raise ValueError(
                f"If `panel_text` is a list/tuple, it must have length N={N}; "
                f"got {len(panel_texts)}."
            )

    if panel_text_kwargs is None:
        panel_text_kwargs = dict(
            color="white",
            ha="left",
            va="top",
            transform=None,
            bbox=dict(
                boxstyle="round,pad=0.15",
                facecolor="black",
                alpha=0.35,
                linewidth=0.0,
            ),
        )

    # ------------------------------------------------------------------
    # cbar labels
    # ------------------------------------------------------------------
    if cbar_mode == "per_panel":
        if cbar_labels is None:
            cbar_labels = [""] * N
        elif isinstance(cbar_labels, str):
            cbar_labels = [cbar_labels] * N
        else:
            cbar_labels = list(cbar_labels)
        if len(cbar_labels) != N:
            raise ValueError(
                "With cbar_mode='per_panel', cbar_labels must be None, a str, "
                "or a list of N strings."
            )

    elif cbar_mode == "shared":
        if isinstance(cbar_labels, (list, tuple)):
            cbar_labels = cbar_labels[0] if len(cbar_labels) else ""
        elif cbar_labels is None:
            cbar_labels = ""

    elif cbar_mode == "first_individual_rest_shared":
        if cbar_labels is None:
            cbar_labels = ["", ""]
        elif isinstance(cbar_labels, str):
            cbar_labels = [cbar_labels, cbar_labels]
        else:
            cbar_labels = list(cbar_labels)
            if len(cbar_labels) != 2:
                raise ValueError(
                    "With cbar_mode='first_individual_rest_shared', "
                    "cbar_labels must be None, a str, or a 2-item list: "
                    "[first_panel_label, shared_rest_label]."
                )

    # ------------------------------------------------------------------
    # broadcast norm parameters per panel
    # ------------------------------------------------------------------
    vmin_list = _broadcast_norm_param(vmin, "vmin", cbar_mode, N)
    vmax_list = _broadcast_norm_param(vmax, "vmax", cbar_mode, N)
    vmin_floor_list = _broadcast_norm_param(vmin_floor, "vmin_floor", cbar_mode, N)
    vmax_floor_list = _broadcast_norm_param(vmax_floor, "vmax_floor", cbar_mode, N)
    symmetric_list = _broadcast_bool_param(symmetric, "symmetric", cbar_mode, N)
    use_lognorm_list = _broadcast_bool_param(use_lognorm, "use_lognorm", cbar_mode, N)

    # ------------------------------------------------------------------
    # norms
    # ------------------------------------------------------------------
    if cbar_mode == "per_panel":
        norms_ext = [
            make_norm_for_data(
                dens[i],
                vmin_i=vmin_list[i],
                vmax_i=vmax_list[i],
                vmin_floor_i=vmin_floor_list[i],
                vmax_floor_i=vmax_floor_list[i],
                symmetric_i=symmetric_list[i],
                use_lognorm_i=use_lognorm_list[i],
            )
            for i in range(N)
        ]

    elif cbar_mode == "shared":
        all_data = np.concatenate([np.ravel(np.asarray(d)) for d in dens])
        all_data = all_data[np.isfinite(all_data)]
        if all_data.size == 0:
            raise ValueError("All dens arrays are empty/non-finite.")

        global_norm, global_ext = make_norm_for_data(
            all_data,
            vmin_i=vmin_list[0],
            vmax_i=vmax_list[0],
            vmin_floor_i=vmin_floor_list[0],
            vmax_floor_i=vmax_floor_list[0],
            symmetric_i=symmetric_list[0],
            use_lognorm_i=use_lognorm_list[0],
        )
        norms_ext = [(global_norm, global_ext)] * N

    elif cbar_mode == "first_individual_rest_shared":
        if N == 1:
            norms_ext = [
                make_norm_for_data(
                    dens[0],
                    vmin_i=vmin_list[0],
                    vmax_i=vmax_list[0],
                    vmin_floor_i=vmin_floor_list[0],
                    vmax_floor_i=vmax_floor_list[0],
                    symmetric_i=symmetric_list[0],
                    use_lognorm_i=use_lognorm_list[0],
                )
            ]
        else:
            first_norm_ext = make_norm_for_data(
                dens[0],
                vmin_i=vmin_list[0],
                vmax_i=vmax_list[0],
                vmin_floor_i=vmin_floor_list[0],
                vmax_floor_i=vmax_floor_list[0],
                symmetric_i=symmetric_list[0],
                use_lognorm_i=use_lognorm_list[0],
            )

            rest_data = np.concatenate([np.ravel(np.asarray(d)) for d in dens[1:]])
            rest_data = rest_data[np.isfinite(rest_data)]
            if rest_data.size == 0:
                raise ValueError("Heatmaps 1..N-1 are empty/non-finite.")

            rest_norm_ext = make_norm_for_data(
                rest_data,
                vmin_i=vmin_list[1],
                vmax_i=vmax_list[1],
                vmin_floor_i=vmin_floor_list[1],
                vmax_floor_i=vmax_floor_list[1],
                symmetric_i=symmetric_list[1],
                use_lognorm_i=use_lognorm_list[1],
            )
            norms_ext = [first_norm_ext] + [rest_norm_ext] * (N - 1)

    # ==================================================================
    # SPECIAL LAYOUT: 3 large + 2 small
    # ==================================================================
    if layout == "3large_2small":
        # Physical plotting spans, accounting for center-vs-edge coordinates.
        xspan = []
        zspan = []
        for i in range(N):
            d = np.asarray(dens[i])
            if d.ndim != 2:
                raise ValueError(f"dens[{i}] must be 2D; got shape {d.shape}.")

            xspan.append(_plot_span(coordx_list[i], d.shape[1], f"coordx[{i}]"))
            zspan.append(_plot_span(coordz_list[i], d.shape[0], f"coordz[{i}]"))

        xspan = np.asarray(xspan, dtype=float)
        zspan = np.asarray(zspan, dtype=float)

        if np.any(xspan <= 0) or np.any(zspan <= 0):
            raise ValueError("All coordinate spans must be positive.")

        # Same z extent -> same physical panel height.
        if not np.allclose(zspan, zspan[0], rtol=1e-3, atol=1e-12):
            raise ValueError(
                "layout='3large_2small' assumes all five panels have the same "
                "physical z extent so they can have the same height. "
                f"Got z spans {zspan}."
            )

        # For the intended layout, the total x-domain represented by each row
        # should be the same: x0+x1 == x2+x3+x4.
        top_domain = xspan[0] + xspan[1]
        bottom_domain = xspan[2] + xspan[3] + xspan[4]

        if not np.isclose(top_domain, bottom_domain, rtol=1e-3, atol=1e-12):
            raise ValueError(
                "For layout='3large_2small', the total physical x extent of the "
                "top and bottom rows must match so the outer edges line up. "
                f"Got top={top_domain:g}, bottom={bottom_domain:g}."
            )

        # One global physical scale: inches per coordinate unit.
        # `gap_x` is the total gutter budget in each row.
        W_data = width_plot - left - right - gap_x
        if W_data <= 0:
            raise ValueError("width_plot too small for requested margins/gaps.")

        scale = W_data / top_domain
        widths = scale * xspan
        h_ax = scale * zspan[0]

        # Horizontal positions.
        # Top row: one full gap_x.
        # Bottom row: two half-gaps, so total gutter width remains gap_x.
        x0 = left
        x1 = x0 + widths[0] + gap_x

        gap_bottom = 0.5 * gap_x
        x2 = left
        x3 = x2 + widths[2] + gap_bottom
        x4 = x3 + widths[3] + gap_bottom

        right_edge_top = x1 + widths[1]
        right_edge_bottom = x4 + widths[4]
        expected_right_edge = width_plot - right

        if not (
            np.isclose(right_edge_top, expected_right_edge, rtol=0, atol=1e-10)
            and np.isclose(right_edge_bottom, expected_right_edge, rtol=0, atol=1e-10)
        ):
            raise RuntimeError("Internal geometry error in layout='3large_2small'.")

        # Vertical positions / total figure height.
        if cbar_mode == "per_panel":
            y_cbar_bot = bottom + xlabel_clearance
            y_bot = y_cbar_bot + cbar_h + cbar_pad_bottom
            y_top = y_bot + h_ax + gap_y
            y_cbar_top = y_top + h_ax + cbar_pad_top
            height_fig = y_cbar_top + cbar_h + top
        else:
            y_bot = bottom + xlabel_clearance
            y_top = y_bot + h_ax + gap_y
            y_cbar_top = y_top + h_ax + cbar_pad_top
            height_fig = y_cbar_top + cbar_h + top
            y_cbar_bot = None

        fig = plt.figure(figsize=(width_plot, height_fig))

        def _inch_rect(x, y, w, h):
            return [
                x / width_plot,
                y / height_fig,
                w / width_plot,
                h / height_fig,
            ]

        panel_positions = [
            (x0, y_top, widths[0], h_ax),  # 0 large
            (x1, y_top, widths[1], h_ax),  # 1 large
            (x2, y_bot, widths[2], h_ax),  # 2 large
            (x3, y_bot, widths[3], h_ax),  # 3 small
            (x4, y_bot, widths[4], h_ax),  # 4 small
        ]

        axes = []
        ims = []

        for i, (xp, yp, wp, hp) in enumerate(panel_positions):
            ax = fig.add_axes(_inch_rect(xp, yp, wp, hp))
            _setup_ticks(ax)

            norm, _ = norms_ext[i]
            im = ax.pcolormesh(
                coordx_list[i],
                coordz_list[i],
                dens[i],
                cmap=cmaps[i],
                norm=norm,
                rasterized=rasterized,
                shading="auto",
            )

            # Because panel box sizes were derived from the physical coordinate
            # extents using one global scale, this should not need to shrink the
            # axes. It simply enforces one rendered unit in x == one in z.
            ax.set_aspect("equal", adjustable="box")

            add_panel_text(ax, panel_texts[i])
            axes.append(ax)
            ims.append(im)

        # Hide inner tick labels.
        if not show_inner_yticks:
            # Only panel 0 and panel 2 are in the left column.
            for i in (1, 3, 4):
                axes[i].tick_params(labelleft=False)

        if not show_inner_xticks:
            # Only top-row x labels are hidden.
            axes[0].tick_params(labelbottom=False)
            axes[1].tick_params(labelbottom=False)

        # --------------------------------------------------------------
        # colorbars for special layout
        # --------------------------------------------------------------
        cbs = []

        if cbar_mode == "shared":
            cax = fig.add_axes(
                _inch_rect(
                    left,
                    y_cbar_top,
                    width_plot - left - right,
                    cbar_h,
                )
            )
            _, extend = norms_ext[0]
            cb = fig.colorbar(
                ims[-1],
                cax=cax,
                orientation="horizontal",
                extend=extend,
            )
            cax.xaxis.set_ticks_position("top")
            cax.xaxis.set_label_position("top")
            cb.set_label(cbar_labels)
            cbs = [cb]

        elif cbar_mode == "first_individual_rest_shared":
            # Panel 0 gets its own colorbar.
            cax0 = fig.add_axes(
                _inch_rect(x0, y_cbar_top, widths[0], cbar_h)
            )
            _, extend0 = norms_ext[0]
            cb0 = fig.colorbar(
                ims[0],
                cax=cax0,
                orientation="horizontal",
                extend=extend0,
            )
            cax0.xaxis.set_ticks_position("top")
            cax0.xaxis.set_label_position("top")
            cb0.set_label(cbar_labels[0])
            cbs.append(cb0)

            # Panels 1..4 share the second colorbar. Place it over the right
            # half of the figure, matching the original function's convention.
            caxr = fig.add_axes(
                _inch_rect(
                    x1,
                    y_cbar_top,
                    expected_right_edge - x1,
                    cbar_h,
                )
            )
            _, extendr = norms_ext[1]
            cb_rest = fig.colorbar(
                ims[1],
                cax=caxr,
                orientation="horizontal",
                extend=extendr,
            )
            caxr.xaxis.set_ticks_position("top")
            caxr.xaxis.set_label_position("top")
            cb_rest.set_label(cbar_labels[1])
            cbs.append(cb_rest)

        elif cbar_mode == "per_panel":
            # Top-row colorbars above panels 0 and 1.
            for idx in (0, 1):
                xp, _, wp, _ = panel_positions[idx]
                cax = fig.add_axes(_inch_rect(xp, y_cbar_top, wp, cbar_h))
                _, extend = norms_ext[idx]
                cb = fig.colorbar(
                    ims[idx],
                    cax=cax,
                    orientation="horizontal",
                    extend=extend,
                )
                cax.xaxis.set_ticks_position("top")
                cax.xaxis.set_label_position("top")
                cb.set_label(cbar_labels[idx])
                cbs.append(cb)

            # Bottom-row colorbars below panels 2, 3 and 4.
            for idx in (2, 3, 4):
                xp, _, wp, _ = panel_positions[idx]
                cax = fig.add_axes(_inch_rect(xp, y_cbar_bot, wp, cbar_h))
                _, extend = norms_ext[idx]
                cb = fig.colorbar(
                    ims[idx],
                    cax=cax,
                    orientation="horizontal",
                    extend=extend,
                )
                cax.xaxis.set_ticks_position("bottom")
                cax.xaxis.set_label_position("bottom")
                cb.set_label(cbar_labels[idx])
                cbs.append(cb)

        # --------------------------------------------------------------
        # shared figure labels
        # --------------------------------------------------------------
        for ax in axes:
            ax.set_xlabel("")
            ax.set_ylabel("")

        y_xlab = (bottom + 0.50 * xlabel_clearance) / height_fig
        fig.text(
            0.5 + (left * -0.05) / width_plot,
            y_xlab,
            xlabel,
            ha="center",
            va="center",
        )

        heatmaps_center_y = y_bot + h_ax + 0.5 * gap_y
        fig.text(
            (left * 0.20) / width_plot,
            heatmaps_center_y / height_fig,
            ylabel,
            ha="center",
            va="center",
            rotation="vertical",
        )

        if omegat is not None:
            fig.text(
                loc_omegat[0] + (left * 0.40) / width_plot,
                loc_omegat[1] - (top * 0.5) / height_fig,
                omegat,
                ha="center",
                va="center",
            )

        return fig, axes, cbs

    # ==================================================================
    # ORIGINAL / AUTOMATIC LAYOUT
    # ==================================================================

    # ------------------------------------------------------------------
    # panel sizing
    # ------------------------------------------------------------------
    aspects = []
    for i in range(N):
        cx = np.asarray(coordx_list[i])
        cz = np.asarray(coordz_list[i])
        if cx.ndim != 1 or cz.ndim != 1:
            raise ValueError("coordx/coordz must be 1D arrays (or lists of 1D arrays).")
        aspects.append(len(cz) / max(1, len(cx)))
    data_aspect = float(max(aspects))

    use_two_rows = (N % 2 == 0) and (N > 1) and (data_aspect < two_row_aspect_threshold)
    ncols_panels = (N // 2) if use_two_rows else N

    # ------------------------------------------------------------------
    # geometry
    # ------------------------------------------------------------------
    W_axes_total = width_plot - left - right - (ncols_panels - 1) * gap_x
    if W_axes_total <= 0:
        raise ValueError("width_plot too small for requested margins/gaps.")

    w_ax = W_axes_total / ncols_panels
    h_ax = w_ax * data_aspect

    if not use_two_rows:
        vbands = [bottom, h_ax, cbar_pad_top, cbar_h, top]
        vnames = ["bottom", "heat_top", "pad_top", "cbar_top", "top"]
    else:
        if cbar_mode == "per_panel":
            vbands = [
                bottom,
                xlabel_clearance,
                cbar_h,
                cbar_pad_bottom,
                h_ax,
                gap_y,
                h_ax,
                cbar_pad_top,
                cbar_h,
                top,
            ]
            vnames = [
                "bottom",
                "xlabel_clear",
                "cbar_bot",
                "pad_bot",
                "heat_bot",
                "gap_y",
                "heat_top",
                "pad_top",
                "cbar_top",
                "top",
            ]
        else:
            vbands = [
                bottom,
                xlabel_clearance,
                h_ax,
                gap_y,
                h_ax,
                cbar_pad_top,
                cbar_h,
                top,
            ]
            vnames = [
                "bottom",
                "xlabel_clear",
                "heat_bot",
                "gap_y",
                "heat_top",
                "pad_top",
                "cbar_top",
                "top",
            ]

    height_fig = float(sum(vbands))
    fig = plt.figure(figsize=(width_plot, height_fig))

    horiz = [Size.Fixed(left)]
    for j in range(ncols_panels):
        horiz.append(Size.Fixed(w_ax))
        if j < ncols_panels - 1:
            horiz.append(Size.Fixed(gap_x))
    horiz.append(Size.Fixed(right))

    divider = Divider(
        fig,
        (0, 0, 1, 1),
        horiz,
        [Size.Fixed(v) for v in vbands],
        aspect=False,
    )
    name_to_ny = {name: i for i, name in enumerate(vnames)}

    def nx_for_col(j):
        return 1 + 2 * j

    # ------------------------------------------------------------------
    # axes + heatmaps
    # ------------------------------------------------------------------
    axes, ims = [], []
    ax_ref = None

    if not use_two_rows:
        top_row_indices = list(range(N))
        bot_row_indices = []
    else:
        top_row_indices = list(range(ncols_panels))
        bot_row_indices = list(range(ncols_panels, N))

    ny_top_heat = name_to_ny["heat_top"]
    ny_bot_heat = name_to_ny.get("heat_bot", None)

    # Top row
    for j, idx in enumerate(top_row_indices):
        nx = nx_for_col(j)

        if (ax_ref is None) or (not can_share):
            ax = fig.add_axes(
                divider.get_position(),
                axes_locator=divider.new_locator(nx=nx, ny=ny_top_heat),
            )
            _setup_ticks(ax)
            if ax_ref is None and can_share:
                ax_ref = ax
        else:
            ax = fig.add_axes(
                divider.get_position(),
                axes_locator=divider.new_locator(nx=nx, ny=ny_top_heat),
                sharex=ax_ref,
                sharey=ax_ref,
            )

        norm, _ = norms_ext[idx]
        cx = coordx_list[idx]
        cz = coordz_list[idx]
        im = ax.pcolormesh(
            cx,
            cz,
            dens[idx],
            cmap=cmaps[idx],
            norm=norm,
            rasterized=rasterized,
            shading="auto",
        )
        ax.set_aspect("equal")

        add_panel_text(ax, panel_texts[idx])

        axes.append(ax)
        ims.append(im)

    # Bottom row
    for j, idx in enumerate(bot_row_indices):
        nx = nx_for_col(j)

        if can_share:
            ax = fig.add_axes(
                divider.get_position(),
                axes_locator=divider.new_locator(nx=nx, ny=ny_bot_heat),
                sharex=ax_ref,
                sharey=ax_ref,
            )
        else:
            ax = fig.add_axes(
                divider.get_position(),
                axes_locator=divider.new_locator(nx=nx, ny=ny_bot_heat),
            )
            _setup_ticks(ax)

        norm, _ = norms_ext[idx]
        cx = coordx_list[idx]
        cz = coordz_list[idx]
        im = ax.pcolormesh(
            cx,
            cz,
            dens[idx],
            cmap=cmaps[idx],
            norm=norm,
            rasterized=rasterized,
            shading="auto",
        )
        ax.set_aspect("equal")
        add_panel_text(ax, panel_texts[idx])

        axes.append(ax)
        ims.append(im)

    # ------------------------------------------------------------------
    # tick label hiding
    # ------------------------------------------------------------------
    if not show_inner_yticks:
        for i_ax, ax in enumerate(axes):
            col = i_ax if not use_two_rows else (i_ax % ncols_panels)
            if col != 0:
                ax.tick_params(labelleft=False)

    if use_two_rows and (not show_inner_xticks):
        for ax in axes[:ncols_panels]:
            ax.tick_params(labelbottom=False)

    # ------------------------------------------------------------------
    # colorbars
    # ------------------------------------------------------------------
    cbs = []

    if cbar_mode == "shared":
        ny_cbar = name_to_ny["cbar_top"]
        nx0 = nx_for_col(0)
        nx1 = nx_for_col(ncols_panels - 1) + 1
        cax = fig.add_axes(
            divider.get_position(),
            axes_locator=divider.new_locator(nx=nx0, nx1=nx1, ny=ny_cbar),
        )
        _, extend = norms_ext[0]
        cb = fig.colorbar(ims[-1], cax=cax, orientation="horizontal", extend=extend)
        cax.xaxis.set_ticks_position("top")
        cax.xaxis.set_label_position("top")
        cb.set_label(cbar_labels)
        cbs = [cb]

    elif cbar_mode == "per_panel":
        ny_cbar_top = name_to_ny["cbar_top"]
        for j, idx in enumerate(top_row_indices):
            nx = nx_for_col(j)
            cax = fig.add_axes(
                divider.get_position(),
                axes_locator=divider.new_locator(nx=nx, ny=ny_cbar_top),
            )
            _, extend = norms_ext[idx]
            cb = fig.colorbar(
                ims[idx],
                cax=cax,
                orientation="horizontal",
                extend=extend,
            )
            cax.xaxis.set_ticks_position("top")
            cax.xaxis.set_label_position("top")
            cb.set_label(cbar_labels[idx])
            cbs.append(cb)

        if use_two_rows:
            ny_cbar_bot = name_to_ny["cbar_bot"]
            for j, idx in enumerate(bot_row_indices):
                nx = nx_for_col(j)
                cax = fig.add_axes(
                    divider.get_position(),
                    axes_locator=divider.new_locator(nx=nx, ny=ny_cbar_bot),
                )
                _, extend = norms_ext[idx]
                cb = fig.colorbar(
                    ims[idx],
                    cax=cax,
                    orientation="horizontal",
                    extend=extend,
                )
                cax.xaxis.set_ticks_position("bottom")
                cax.xaxis.set_label_position("bottom")
                cb.set_label(cbar_labels[idx])
                cbs.append(cb)

    elif cbar_mode == "first_individual_rest_shared":
        ny_cbar = name_to_ny["cbar_top"]

        cax0 = fig.add_axes(
            divider.get_position(),
            axes_locator=divider.new_locator(nx=nx_for_col(0), ny=ny_cbar),
        )
        _, extend0 = norms_ext[0]
        cb0 = fig.colorbar(ims[0], cax=cax0, orientation="horizontal", extend=extend0)
        cax0.xaxis.set_ticks_position("top")
        cax0.xaxis.set_label_position("top")
        cb0.set_label(cbar_labels[0])
        cbs.append(cb0)

        if N > 1:
            if (not use_two_rows) and N > 2:
                nx0 = nx_for_col(1)
                nx1 = nx_for_col(ncols_panels - 1) + 1
                caxr = fig.add_axes(
                    divider.get_position(),
                    axes_locator=divider.new_locator(nx=nx0, nx1=nx1, ny=ny_cbar),
                )
            elif (not use_two_rows) and N == 2:
                caxr = fig.add_axes(
                    divider.get_position(),
                    axes_locator=divider.new_locator(nx=nx_for_col(1), ny=ny_cbar),
                )
            elif use_two_rows and ncols_panels > 1:
                nx0 = nx_for_col(1)
                nx1 = nx_for_col(ncols_panels - 1) + 1
                caxr = fig.add_axes(
                    divider.get_position(),
                    axes_locator=divider.new_locator(nx=nx0, nx1=nx1, ny=ny_cbar),
                )
            else:
                caxr = fig.add_axes(
                    divider.get_position(),
                    axes_locator=divider.new_locator(nx=nx_for_col(0), ny=ny_cbar),
                )

            _, extendr = norms_ext[1]
            cb_rest = fig.colorbar(
                ims[1],
                cax=caxr,
                orientation="horizontal",
                extend=extendr,
            )
            caxr.xaxis.set_ticks_position("top")
            caxr.xaxis.set_label_position("top")
            cb_rest.set_label(cbar_labels[1])
            cbs.append(cb_rest)

    # ------------------------------------------------------------------
    # shared fig labels
    # ------------------------------------------------------------------
    for ax in axes:
        ax.set_xlabel("")
        ax.set_ylabel("")

    if use_two_rows:
        y_xlab = (bottom + 0.50 * xlabel_clearance) / height_fig
    else:
        y_xlab = (bottom * 0.45) / height_fig

    fig.text(
        0.5 + (left * -0.05) / width_plot,
        y_xlab,
        xlabel,
        ha="center",
        va="center",
    )
    fig.text(
        (left * 0.20) / width_plot,
        0.5 + y_xlab,
        ylabel,
        ha="center",
        va="center",
        rotation="vertical",
    )

    if omegat is not None:
        fig.text(
            loc_omegat[0] + (left * 0.40) / width_plot,
            loc_omegat[1] - (top * 0.5) / height_fig,
            omegat,
            ha="center",
            va="center",
        )

    return fig, axes, cbs
