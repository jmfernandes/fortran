
PROG1 =	pend


OBJS1 =	numtype.o adpend_cop.o rkf45step_cop.o

F95 = gfortran

F95FLAGS = -O3  

LDFLAGS =  

all: $(PROG1) 

$(PROG1): $(OBJS1)
	$(F95) $(LDFLAGS) -o $@ $(OBJS1) 

clean:
	rm -f $(PROG1)  *.{o,mod}

.SUFFIXES: $(SUFFIXES) .f95 

.f95.o:
	$(F95) $(F95FLAGS) -c $<

