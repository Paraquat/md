include Makefile.inc
.DEFAULT_GOAL = md
PROG = md

#FCFLAGS = -fbackslash
LDFLAGS = $(LPATH) $(IPATH) -lzamaan -llapack -lblas
SRCS = input_module.f90 md_control.f90 md_model.f90 md_module.f90 cmdline.f90 main.f90
OBJS = input_module.o md_control.o md_model.o md_module.o cmdline.o main.o

.SUFFIXES: $(SUFFIXES) .f90
.PHONY: clean

all: $(PROG)

clean:
	rm -f $(PROG) $(OBJS) *.mod

$(PROG): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

.f90.o:
	$(FC) -c $< $(LDFLAGS)

%.o: %.mod

input_module.o: input_module.f90
md_control.o: md_control.f90 input_module.o ../libzamaan/libzamaan.a
md_model.o: md_model.f90 md_control.o input_module.o ../libzamaan/libzamaan.a
md_module.o: md_module.f90 md_control.o input_module.o md_model.o md_control.f90
cmdline.o: cmdline.f90
main.o: cmdline.o md_module.f90 main.f90 ../libzamaan/libzamaan.a
