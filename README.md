# Code underlying the manuscript "Dust-driven vertical shear instability in an isothermal shearing box"

This repository contains the setup scripts required for the project titled [Dust-driven Vertical Shear Instability](https://github.com/JipMatt/DVSI), which is documented in the accompanying research paper. 

This reposetory is openly avialble on `4TU.ResearchData` under the name: [Code underlying the manuscript "Dust-driven vertical shear instabilty in an isothermal shearing box"](https://doi.org/10.4121/0ab6e524-0c51-4c40-986e-fd5bca604953)

The following sections describe the structure of the repository, how to install dependencies, and how to set up the environment to use the scripts and data.

## Table of Contents

- [Project Structure](#project-structure)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Data Description](#data-description)
- [Setup Instructions](#setup-instructions)
- [Usage](#usage)
- [Paper](#paper)
- [License](#license)
- [Contact](#contact)

## Project Structure

```plaintext
|-- DVSI_I.ipynb        # Script generates the plots of manuscript I from the raw data and monitor files.
|-- DVSI_II.ipynb       # Script generates the plots of manuscript II from the raw data and monitor files.
|
|-- scripts/    
    |-- vistools.py         # Class to read in raw data in .dat format
    |-- polydust.py         # Class holding a size (or stopping time) distribution
    |-- func.py             # Script that with fuction used to manipulate the raw data and monitoring files.
    |-- plot.py             # Script with universal plotting funtions used to generating the plots
|
|-- README.md           # This file
|-- requirements.txt    # List of Python dependencies
|-- environment.yml     # Conda environment file (if using Conda)
```

## Prerequisites

Before running the scripts, ensure that you have the following software installed:

- **Python 3.x** (or the version specified in the `requirements.txt` or `environment.yml`)
- **FARGO3d** [FARGO3D](`https://github.com/FARGO3D/fargo3d`)
- **`pip` or `conda`**
- **Git** (for cloning a repository [`psitools`](https://github.com/psitools/psitools-public), ['spectral_solver'](https://github.com/SijmeJan/spectral-solvers), and ['psi_stratified'](https://github.com/SijmeJan/psi-stratified))


## Installation

### Step 1: Clone the Repository

If you haven't cloned the repository, you can do so using the following command:

```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
```

### Step 2: Install Dependencies

Using `pip`:

```bash
pip install -r requirements.txt
```

Alternatively, using `conda`:

```bash
conda env create -f environment.yml
conda activate your-env-name
```

Other repositories have be cloned though git psitools, spectral-slover, and psi-stratified

```bash
git clone https://github.com/psitools/psitools-public.git
pip install -e  ~/path/to/psitools-public

git clone https://github.com/SijmeJan/spectral-solvers.git
pip install -e ~/path/to/spectral-solvers

git clone https://github.com/SijmeJan/psi-stratified.git
pip install -e ~/path/to/psi-stratified
```

## Data Description

The data that support the findings of this study are openly available in `4TU.ResearchData` under the name: [Data underlying the manuscript "Dust-driven vertical shear instability in an isothermal shearing box"](http://doi.org/10.4121/571e09ec-10db-452d-b541-b7d86f573e50).

### Directory Structure

```plaintext
|-- /data                   # Directory containing raw and processed data
    |-- /[run_name]           
        |-- /raw_data       # Raw datasets (output of the code FARGO3D)
        |-- /monitor        # Monitoring files (output of the code FAGRO3D)
        |-- /setup          # Setup files for FARGO3D
|-- /plots                  # Plots that appear in the coresponding articles
|-- README.md               # This file
```

The data is stored in `/data` directory. Subdivided in the run_names, these folders are futher subdivided into thee subdirectories:

1. **Raw Data (`/[run_name]/raw_data`)**: These folders contain the raw output of the specific snapshots from FARGO3D used in the corresponding paper. Each dataset has been named according to the convention of FARGO3D`[dataset_name].dat`. Although, due to space constraints only specific snapshots relevant to the paper are stored. 

2. **Setup Files (`/[run_name]/setup`)**: These folders contain the FARGO3D setup files. These files can be reproduced by running the `psi/createsetup.py` to make a setup and run this setup with FARGO3D will reproduce the raw data.

3. **Monitor Files (`/data/[run_name]/monitor`])**: The folders contain the monitor files that exist of scaler or 1D timeseries of the FARGO3D runs are in a subfolders `[run_name]/monitor/[fluidtype]/` where `fluidtype` can be `gas`, `dust` or `tracer`.


### File Descriptions

- `data/[run_name]/[dataset_name].dat`: All the outputs of FARGO3D for the different fluids at a specific snapshot.
- `data/[run_name]/monitor/[fluidtype]/[scaler]_1d_Z_raw.dat`: 1D timeseries per horizonal pixex of the FARGO3D runs used for the final analysis.
- `data/[run_name]/monitor/[fluidtype]/[scaler].dat`: Timeseries of `scaler` of the FARGO3D runs used for the final analysis.

## Usage

Running the analysis script you reproduce the plot in the paper.

Example:

```bash
# Example command to run an analysis
python scripts/DVSI_I.ipynb --input data/[run_name]/[dataset_name].dat  --output name_figure.png
```

### Important Notes

- Ensure the paths to datasets are correctly defined in the scripts.
- Modify the parameters inside the scripts (e.g., hyperparameters, input file path) if necessary.

## Paper

    Dust-driven vertical shear instabilty in an isothermal shearing box: I. Turbulence and gas dynamics without dust feedback
    Matthijsse, Jip; Aly, Hossam; Paardekooper, Sijme-Jan; Nelson, Richard p.
    Astronomy & Astrophysics Journal,
    DOI: http://doi.org/10.1051/0004-6361/202662433 
    SciX: 
    arXiv: https://arxiv.org/abs/2610.12648

    Dust-driven vertical shear instabilty in an isothermal shearing box: II. Vertical mixing of dust and the stabilising effect of dust feedback
    Matthijsse, Jip; Aly, Hossam; Paardekooper, Sijme-Jan; van der Marel, Nienke; Nelson, Richard p.
    Submitted to Astronomy & Astrophysics Journal (Not yet accepted),
    DOI: 
    SciX: 
    arXiv:

## License

This project is licensed under the Apache 2.0 License. See the [LICENSE](https://www.apache.org/licenses/LICENSE-2.0) file for details.

## Contact

* Jip Matthijsse
* Hossam Aly
* Sijme-Jan Paardekooper
* Richard P. Nelson
* Nienke van der Marel

For any inquiries or questions, please contact:

- **Author Name**: Jip Matthijsse
- **Email**:       j.p.matthijsse@tudelft.nl
- **Institution**: Technical University Delft
 