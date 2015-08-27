/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
// #include "mpi.h"
// #include "Libs/real2complex.h"
#include "Libs/mpi_lsa_com.h"
// #include <stdio.h>
// #include "petsc.h"
// #include "petscvec.h"
// #include "petscmat.h"
#include "slepceps.h"
#include <unistd.h>
#include <sys/stat.h>
#include <time.h>
#include "./Solvers/GMRES/gmres.h"
#include "./Solvers/Arnoldi/arnoldi.h"
#include "./Libs/read_matrix.h"
#include "./Solvers/LS/lsqr.h"
#include "./Solvers/Father/father.h"

