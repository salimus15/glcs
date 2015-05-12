/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#ifndef CONVHULL_H_

#include "petsc.h"
#include <math.h>

int convhull(PetscScalar * ab, PetscScalar * c, PetscScalar * d, PetscInt n, PetscInt * ne, PetscInt offset, PetscInt mu);


#endif