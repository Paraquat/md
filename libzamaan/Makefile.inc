# Library paths
LPATH = -L/Users/zamaan/Software/lib
IPATH = -I/Users/zamaan/Software/include

# Uncomment for debugging
DEBUG = -O0 -g -fbounds-check

# Optimisations
#OPT = -O3

# Fortran compiler
FC = gfortran $(OPT) $(DEBUG)
