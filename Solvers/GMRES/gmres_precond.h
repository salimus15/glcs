/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#ifndef GMRES_PRECOND_H
#define GMRES_PRECOND_H

#include "fgmresimpl.h"
// #include "../Arnoldi/arnoldi.h"
#include "petsc.h"
#include "../Utils/lib.h"
#include "../../Libs/mpi_lsa_com.h"
#include <unistd.h>

#ifndef EIGEN_ALL
#define EIGEN_ALL 100
#endif

#ifndef LS_POWER
#define LS_POWER 15
#endif

#ifndef LS_LATENCY
#define LS_LATENCY 2
#endif

PetscErrorCode GmresLSAPrecond(com_lsa * com, KSP ksp);

#endif