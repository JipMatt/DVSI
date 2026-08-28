from common import *

def write_runscript(setup_name,time=[3,59,59], n_tasks=48):
    lines = ['#!/bin/sh\n',
             '#\n',
             '#SBATCH --job-name={}\n'.format(setup_name),
             '#SBATCH --partition=compute\n',
             '#SBATCH --time={:02d}:{:02d}:{:02d} #max 5 days HH:MM:SS\n'.format(time[0], time[1], time[2]),
             '#SBATCH --ntasks={:d} #max 48 nodes\n'.format(n_tasks),
             '#SBATCH --cpus-per-task=1\n',
             '#SBATCH --mem-per-cpu=1G\n',
             '#SBATCH --account=research-ae-spe\n\n',
             'module load 2024r1\n',
             'module load openmpi\n',
             'srun ./fargo3d setups/{}/{}.par\n'.format(setup_name, setup_name),
             '#for restarting a job\n',
             '#srun ./fargo3d -S n setups/{}/{}.par #n last saved timeframe]\n'.format(setup_name, setup_name)]
    write_file('run.sub', lines)

