include Makefile.inc

.DEFAULT_GOAL = libzamaan.a

LDFLAGS = $(LPATH) $(IPATH) -llapack -lblas
FCFLAGS = -fbackslash # gfortran ONLY!
SRCS = datatypes.f90 constants.f90 maths.f90 rng.f90 vector.f90 progress.f90 grft.f90 data.f90 cell.f90 pairdist.f90 pairpotential.f90
OBJS = datatypes.o constants.o rng.o maths.o vector.o progress.o grft.o data.o cell.o pairdist.o pairpotential.o

.SUFFIXES: $(SUFFIXES) .f90
.PHONY: clean

libzamaan.a: $(OBJS)
	ar rc libzamaan.a $(OBJS)

clean:
	rm -f $(OBJS) *.mod libzamaan.a

.f90.o:
	$(FC) -c $< $(LDFLAGS)

%.o: %.mod

pairpotential.o: datatypes.o constants.o pairpotential.f90
pairdist.o: datatypes.o constants.o progress.o cell.o data.o grft.o pairdist.f90
cell.o: datatypes.o cell.f90 vector.o
data.o: datatypes.o data.f90
grft.o: datatypes.o grft.f90
progress.o: datatypes.o progress.f90
vector.o: datatypes.o constants.o vector.f90
rng.o: datatypes.o rng.f90
maths.o: datatypes.o constants.o maths.f90
constants.o: datatypes.o constants.f90
