/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "petsc.h"
#include "petscpc.h"
#include "petscmat.h"
#include "petscksp.h"
#include <stdlib.h>
#include "../Utils/lib.h"

#define at(a, i, j) (a)[j][i]

PetscErrorCode LSPrecond(PetscReal a_ell, PetscReal d_ell, PetscReal c_ell,
			PetscScalar * eta, PetscScalar * alpha, PetscScalar * beta, 
			PetscScalar * delta, PetscScalar * c, PetscScalar * d, 
			PetscInt * mu, PetscInt * nb_eigen, PetscInt * min_eigen, PetscInt * nb_eigen_all);
