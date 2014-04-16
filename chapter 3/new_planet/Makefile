
OBJS1 =	numtype.o planet.o rk4step.o 

PROG1 =	plan

F95 = gfortran

F95FLAGS = -O3 -funroll-loops -fexternal-blas

LDFLAGS = -framework vecLib

all: $(PROG1) 

$(PROG1): $(OBJS1)
	$(F95) $(LDFLAGS) -o $@ $(OBJS1) 

clean:
	rm -f $(PROG1)  *.{o,mod}

.SUFFIXES: $(SUFFIXES) .f95 

.f95.o:
	$(F95) $(F95FLAGS) -c $<

