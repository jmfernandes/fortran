
OBJS =	numtype.o mtest_papp_version.o

PROG =	mat

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

