include Makefile.inc 
.DEFAULT_GOAL = rdf
PROG = rdf

#FCFLAGS = -fbackslash
LDFLAGS = $(LPATH) $(IPATH) -lzamaan -llapack -lblas
SRCS = cmdline.f90 main.f90
OBJS = cmdline.o main.o
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

main.o: cmdline.o main.f90 ../libzamaan/libzamaan.a
cmdline.o: cmdline.f90
