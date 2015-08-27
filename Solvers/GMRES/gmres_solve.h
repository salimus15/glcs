/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#ifndef GMRES_SOLVE_H
#define GMRES_SOLVE_H

#include "petsc.h"
#include "gmres_cycle.h"
#include "gmres_precond.h"
#include "../../Libs/mpi_lsa_com.h"
#include "fgmresimpl.h"

//static PetscErrorCode MyKSPFGMRESResidual(KSP ksp);
PetscErrorCode MyKSPSolve_FGMRES(KSP ksp,com_lsa * com);

#endif
