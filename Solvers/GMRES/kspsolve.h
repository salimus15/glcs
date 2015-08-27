/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#ifndef MY_KSP_SOLVE_H
#define MY_KSP_SOLVE_H
#include "petsc.h"
#include "../../Libs/mpi_lsa_com.h"
#include "gmres_solve.h"


PetscErrorCode MyKSPSolve(KSP ksp,Vec b,Vec x,com_lsa * com);


#endif
