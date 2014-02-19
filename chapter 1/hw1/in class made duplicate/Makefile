
OBJS = numtype.o scattering_dup.o

PROG = scattering_dup

F95 = gfortran

F95FLAGS = -O3 -funroll-loops -fexternal-blas

LIBS = -framework vecLib

LDFLAGS = $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F95) $(LDFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(PROG) *.{o,mod}


.SUFFIXES: $(SUFFIXES) .f95

.f95.o:
	$(F95) $(F95FlAGS) -c $<