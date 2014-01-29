
OBJS = numtype.o ball.o

PROG = ball

F95 = gfortran

F95FLAGS = -O3 -funroll-loops -march=native -fexternal-blas

LIBS = -framework vecLib

LDFLAGS = $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F95) $(LDFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(PROG) *,{0,mod}


.SUFFICES: $(SUFFICES) .f95

.f95.o:
	$(F95) $(F95FlAGS) -c $<