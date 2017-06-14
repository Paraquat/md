# md
A classical MD code to test concepts and implementation for use in CONQUEST, in particular, extended-Lagrangian schemes for
thermostatting, barostatting and constraining.

## Progress

MD can run classical molecular dynamics simulations using Lennard-Jones and Morse pair potentials at present. Available 
ensembles are NVE with velocity Verlet integration, NVT with velocity rescaling or Nose-Hoover chains with the MTK
implmenetation. NPT is currently under development using the the Parrinello-Rahman barostat and Nose-Hoover Chains in the
MTTK implementation.

## Installation

Requirements:
* BLAS/LAPACK
* Fortran compiler (I use gfortran)
* Python 3 for analysis scripts
* Scipy
* Matplotlib
* VMD for visualisation

Simply run "make" in the libzamaan directory and molecular_dynamics directories, or run "build.sh" in the molecular dynamics
directory. It may be necessary to change some paths in the Makefile.inc files, but it runs on my Macbook Pro.

## Usage

Once compiled, run the main program using the "md" executable. The program will look for the following three input files in 
the working directory (examples are in the "examples" directory):
1. cell.in - the atomic configuration and cell
2. pp.in - pair potential parameters
3. md.in - molecular dynamics paramters

It will pipe a load of text to stdout, and additionally generate the following output files:
1. stats.out - thermodynamic data, including energies, temperature, pressure
2. dump.out - detailed data dump, including atomic configuration, cell parameters, stress tensor, etc for each step
3. trajectory.xsf - the trajectory in .xsf format, can be visualised using VMD.

The md_analysis script will look for the above outputs, and will generate:
1. stats.pdf - a plot of the thermodynamic data
2. stress.pdf - plot of the stress tensor elements and pressure
3. vacf.pdf - velocity autocorrelation functions
4. vdistr.pdf - velocity distribution
5. msd.pdf - mean squared displacement

A VMD script (view.vmd) is included with the examples; trajectory.xsf can be visualised by typing,

```
vmd -e view.vmd
```

Assuming the vmd executable is in your $PATH.
