
OBJS =	numtype.o int_mod.o d01b.o int.o

PROG =	int

F95 = gfortran

F95FLAGS = -O3 -funroll-loops -fexternal-blas

LDFLAGS = -framework vecLib

all: $(PROG) 

$(PROG): $(OBJS)
	$(F95) $(LDFLAGS) -o $@ $(OBJS) 

clean:
	rm -f $(PROG)  *.{o,mod}

.SUFFIXES: $(SUFFIXES) .f95 

.f95.o:
	$(F95) $(F95FLAGS) -c $<

