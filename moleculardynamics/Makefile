include Makefile.inc 
.DEFAULT_GOAL = md
PROG = md

#FCFLAGS = -fbackslash
LDFLAGS = $(LPATH) $(IPATH) -lzamaan  -lfftw3 -llapack -lblas
SRCS = cmdline.f90 main.f90 md_module.f90
OBJS = cmdline.o main.o md_module.o
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

main.o: cmdline.o md_module.f90 main.f90 ../libzamaan/libzamaan.a
cmdline.o: cmdline.f90
md_module.f90: md_module.o ../libzamaan/libzamaan.a
