# Library paths
LPATH = -L/home/zamaan/Software/lib -L../libzamaan
IPATH = -I/home/zamaan/Software/include -I../libzamaan

# Uncomment for debugging
DEBUG = -O0 -g -fbounds-check

# Optimisations
#OPT = -O3

# Fortran compiler
FC = gfortran $(OPT) $(DEBUG)
