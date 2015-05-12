/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#ifndef ARNOLDI_H
#define ARNOLDI_H 

#include "slepceps.h"
#include <unistd.h>
#include "../../Libs/mpi_lsa.h"
#include "../../Libs/mpi_lsa_com.h"
#include "../../Libs/data_rw.h"

//#include "epsimpl.h" // not supposed to do this.... but that's just to overcome the slepc problem of setting the initial vector


#ifndef EIGEN_ALL
#define EIGEN_ALL 20
#endif

PetscErrorCode Arnoldi(com_lsa * com, Mat * A, Vec * v);

#endif