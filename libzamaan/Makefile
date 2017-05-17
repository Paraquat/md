include Makefile.inc

.DEFAULT_GOAL = libzamaan.a

LDFLAGS = $(LPATH) $(IPATH) -lfftw3 -llapack -lblas
FCFLAGS = -fbackslash # gfortran ONLY!
SRCS = datatypes.f90 constants.f90 rng.f90 vector.f90 progress.f90 grft.f90 data.f90 cell.f90 trajectory.f90 pairdist.f90
OBJS = datatypes.o constants.o rng.o vector.o progress.o grft.o data.o cell.o trajectory.o pairdist.o

.SUFFIXES: $(SUFFIXES) .f90
.PHONY: clean

libzamaan.a: $(OBJS)
	ar rc libzamaan.a $(OBJS)

clean:
	rm -f $(OBJS) *.mod libzamaan.a

.f90.o:
	$(FC) -c $< $(LDFLAGS)

%.o: %.mod

pairdist.o: datatypes.o constants.o progress.o cell.o data.o grft.o pairdist.f90
trajectory.o: datatypes.o cell.o trajectory.f90
cell.o: datatypes.o cell.f90 vector.o
data.o: datatypes.o data.f90
grft.o: datatypes.o grft.f90
progress.o: datatypes.o progress.f90
vector.o: datatypes.o constants.o vector.f90
rng.o: datatypes.o rng.f90
constants.o: datatypes.o constants.f90
