
CC = gcc
CFLAGS = -O2 -c
OBJS = file1.o file2.o
SRCS = file1.c file2.c

all : executable
.c.o :
        $(CC) $(CFLAGS) $<
executable : $(OBJS)
        $(CC) -o $@ $(OBJS)CC = gcc
depend :
        makedepend -I. $(SRC)
% cat file1.c
#include "file1.h"
main() {}
% cat file2.c
#include "file2.h"
#include "file1.h"
% makedepend -I. file1.c file2.c
% cat Makefile
CC = gcc
CFLAGS = -O2 -c
OBJS = file1.o file2.o
SRCS = file1.c file2.c

all : executable
.c.o :
        $(CC) $(CFLAGS) $<
executable : $(OBJS)
        $(CC) -o $@ $(OBJS)CC = gcc
depend :
        makedepend -I. $(SRC)

# DO NOT DELETE

file1.o: ./file1.h
file2.o: ./file2.h ./file1.h



/home/vassago/Tools/petsc-3.5.3/linux-c-debug-complex/bin/mpicc -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -g3 -O0 -O0 -Wl,-rpath,/home/vassago/Tools/slepc-3.5.3/linux-c-debug-complex/lib -L/home/vassago/Tools/slepc-3.5.3/linux-c-debug-complex/lib -lslepc       -L/home/vassago/Tools/petsc-3.5.3/linux-c-debug-complex/lib  -lpetsc -llapack -lblas -lX11 -lpthread -lssl -lcrypto -lm -L/usr/lib/gcc/x86_64-unknown-linux-gnu/4.9.2 -lmpichf90 -lgfortran -lm -Wl,-rpath,/home/vassago/Tools/petsc-3.5.3/linux-c-debug-complex/lib -lgfortran -lm -lquadmath -lm -lmpichcxx -lstdc++ -L/home/vassago/Tools/petsc-3.5.3/linux-c-debug-complex/lib -L/usr/lib/gcc/x86_64-unknown-linux-gnu/4.9.2 -ldl -Wl,-rpath,/home/vassago/Tools/petsc-3.5.3/linux-c-debug-complex/lib -lmpich -lopa -lmpl -lrt -lpthread -lgcc_s -ldl   -L. -lglsa -o hyperh main.c -I./Libs  -I./Solvers/Utils -I./Solvers/Arnoldi -I./Tuning

