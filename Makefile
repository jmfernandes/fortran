
OBJS = numtype.o ball.o

PROG = ball

F95 = gfortran

F95FLAGS = -O3 -funroll-loops -march=native -fexternal-blas

LIBS = -framework vecLib