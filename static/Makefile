include Makefile.inc
.DEFAULT_GOAL = static
PROG = static

#FCFLAGS = -fbackslash
LDFLAGS = $(LPATH) $(IPATH) -lzamaan  -lfftw3 -llapack -lblas
SRCS = input_module.f90 cmdline.f90 static_module.f90 main.f90
OBJS = input_module.o cmdline.o  static_module.o main.o
DEBUG = -O0 -g -fbounds-check

.SUFFIXES: $(SUFFIXES) .f90
.PHONY: clean

FC = gfortran $(DEBUG) $(FCFLAGS)

all: $(PROG)

clean:
	rm -rf $(PROG) $(OBJS) *.mod $(PROG).dSYM

$(PROG): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

.f90.o:
	$(FC) -c $< $(LDFLAGS)

%.o: %.mod

$(OBJS): $(SRCS)
input_module.o: input_module.f90
cmdline.o: cmdline.f90
static_module.o: static_module.f90 input_module.o ../libzamaan/libzamaan.a
main.o: cmdline.o input_module.o static_module.o main.f90 ../libzamaan/libzamaan.a
