include Makefile.inc
.DEFAULT_GOAL = md
PROG = md

#FCFLAGS = -fbackslash
LDFLAGS = $(LPATH) $(IPATH) -lzamaan  -lfftw3 -llapack -lblas
SRCS = input_module.f90 thermostat_module.f90 barostat_module.f90 md_module.f90 cmdline.f90 main.f90
OBJS = input_module.o thermostat_module.o barostat_module.o md_module.o cmdline.o main.o
DEBUG = -O0 -g -fbounds-check

.SUFFIXES: $(SUFFIXES) .f90
.PHONY: clean

FC = gfortran $(DEBUG) $(FCFLAGS)

all: $(PROG)

clean:
	rm -f $(PROG) $(OBJS) *.mod

$(PROG): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

.f90.o:
	$(FC) -c $< $(LDFLAGS)

%.o: %.mod

input_module.o: input_module.f90
thermostat_module.o: thermostat_module.f90 ../libzamaan/libzamaan.a
barostat_module.o: barostat_module.f90 ../libzamaan/libzamaan.a
md_module.o: md_module.f90 thermostat_module.f90 barostat_module.f90 ../libzamaan/libzamaan.a
cmdline.o: cmdline.f90
main.o: cmdline.o md_module.f90 main.f90 ../libzamaan/libzamaan.a
