
OBJS1 = numtype.o cheby.o ho1d10.o rk4step.o

PROG1 = ho1d

OBJS2 = numtype.o cheby.o ho1d10_energy.o rk4step_energy.o

PROG2 = ho2d

F95 = gfortran

F95FLAGS = -O3 -funroll-loops -fexternal-blas

LIBS = -framework vecLib

LDFLAGS = $(LIBS)

all: $(PROG1) $(PROG2)

$(PROG1): $(OBJS1)
	$(F95) $(LDFLAGS) -o $@ $(OBJS1) 

$(PROG2): $(OBJS2)
	$(F95) $(LDFLAGS) -o $@ $(OBJS2) 

clean:
	rm -f $(PROG1) *.{o,mod}
	rm -f $(PROG2) *.{o,mod}

.SUFFIXES: $(SUFFIXES) .f95

.f95.o:
	$(F95) $(F95FLAGS) -c $<

