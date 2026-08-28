from matplotlib import lines
from common import *

def add_option(lines, option):
    '''Add option to opt file content'''
    add_line(lines, 'FARGO_OPT +=  -D'+ option + '\n')

def write_opt_file(setup_name, pd, acceleration_gas=False, dust=True, tracer=None, shear=None, nubound=None):
    '''Write FARGO .opt file

    Args:
        setup_name: Name of the setup being created
        pd: Polydust object
    '''
    n_dust = pd.N

    if tracer is not None:
        n_tracer = len(tracer)
        print('n_tracer', n_tracer)
        if dust:
            lines = ['# PSI FARGO2D setup using ' + str(n_dust) + ' dust fluids and a tracer fluid\n\n']
        else:
            lines = ['# PSI FARGO2D setup using a tracer fluid but no dust\n\n']
    else:
        if dust:
            lines = ['# PSI FARGO2D setup using ' + str(n_dust) + ' dust fluids\n\n']
        else:
            lines = ['# PSI FARGO2D setup using no dust fluids or tracer fluids\n\n']

    line = 'FLUIDS := 0'
    if dust:
        for i in range(1, n_dust + 1):
            line = line + ' ' + str(i)
        if tracer is not None:
            for i in range(1, n_tracer+1):
                line = line + ' ' + str(n_dust + i)
    else:
        if tracer is not None:
            for i in range(1, n_tracer + 1):
                line = line + ' ' + str(i)
    add_line(lines, line +'\n')
    if dust:
        if tracer is not None:
            add_line(lines, 'NFLUIDS = ' + str(n_dust + n_tracer + 1) + '\n')
        else:
            add_line(lines, 'NFLUIDS = ' + str(n_dust + 1) + '\n')
    else:
        if tracer is not None:
            add_line(lines, 'NFLUIDS = ' + str(n_tracer + 1) + '\n')
        else:
            add_line(lines, 'NFLUIDS = 1\n')
    add_line(lines, 'FARGO_OPT += -DNFLUIDS=${NFLUIDS}\n\n')

    add_option(lines, 'X')
    add_option(lines, 'Y')
    add_option(lines, 'Z')

    add_option(lines, 'CARTESIAN')
    add_option(lines, 'SHEARINGBOX')
    add_option(lines, 'SHEARINGBC')
    add_option(lines, 'ISOTHERMAL')
    add_option(lines, 'VISCOSITY')
    add_option(lines, 'NOSUBSTEP2')
    if nubound is not None:
        add_option(lines, 'DIFFUSIONBOUND')
    if shear is not None:
        add_option(lines, 'VERTICALSHEAR')


    if acceleration_gas:
        add_option(lines, 'ACCELERATIONGAS')
    if dust or tracer is not None:
        add_option(lines, 'DRAGFORCE')
    add_option(lines, 'CONSTANTSTOKESNUMBER')
    if pd.pressure_gradient is False:
        add_option(lines, 'PRESSUREGRADIENT')
    add_option(lines, 'STRATIFIELD')
    add_option(lines, 'FLOOR')
    add_option(lines, 'HARDBOUNDARIES')
    #add_option(lines, 'COLLISIONPREDICTOR')
    if dust:
        add_option(lines, 'DUSTDIFFUSION')
        add_option(lines, 'DUSTVISCOSITY')

    add_line(lines, 'MONITOR_SCALAR = MASS | EKTOT | EKY | EKZ | RHOMAX | RRHO | RTOT | RYY | RZZ | RENYX | RENZX | RYX | RYZ | RZX | MOM_X | MOM_Y | MOM_Z | \n')
    add_line(lines, 'MONITOR_Z_RAW = MASS | EKTOT | EKY | EKZ | RRHO | RTOT | RYY | RZZ | RENYX | RENZX | RYX | RYZ | RZX | MOM_X | MOM_Y | MOM_Z | \n')

    add_line(lines, '\n#Cuda blocks\n')
    add_line(lines, 'ifeq (${GPU}, 1)\n')
    add_line(lines, 'FARGO_OPT += -DBLOCK_X=16\n')
    add_line(lines, 'FARGO_OPT += -DBLOCK_Y=16\n')
    add_line(lines, 'FARGO_OPT += -DBLOCK_Z=1\n')
    add_line(lines, 'endif\n')

    fname = setup_name + '.opt'
    write_file(fname, lines)

    return fname
