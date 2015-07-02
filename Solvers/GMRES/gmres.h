/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#ifndef GMRES_H
#define GMRES_H

//#include "fgmresimpl.h"
#include "petsc.h"
#include "petscvec.h"
#include "petscksp.h"
#include "mpi.h"
#include "../../Libs/mpi_lsa_com.h"
#include "gmres_solve.h"
#include "kspsolve.h"


PetscErrorCode launchGMRES(com_lsa * com, Vec * b, Mat * A);

#endif

