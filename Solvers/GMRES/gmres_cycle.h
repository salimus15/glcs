/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#ifndef GMRES_CYCLE_H
#define GMRES_CYCLE_H

/* just one include necessary */
#include "fgmresimpl.h"
#include "petsc.h"
#include "petsc/private/kspimpl.h"


PetscErrorCode MyKSPFGMRESCycle(PetscInt *itcount,KSP ksp);

PetscErrorCode MyKSPFGMRESGetNewVectors(KSP,PetscInt);
PetscErrorCode MyKSPFGMRESUpdateHessenberg(KSP,PetscInt,PetscBool,PetscReal *);
PetscErrorCode MyKSPFGMRESBuildSoln(PetscScalar*,Vec,Vec,KSP,PetscInt);

#endif
