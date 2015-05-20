/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "read_matrix.h"


PetscErrorCode read_matrix_vector(Mat * A, Vec * v, int * communicator){
	char filea[PETSC_MAX_PATH_LEN];
	char fileb[PETSC_MAX_PATH_LEN];
	char err[PETSC_MAX_PATH_LEN];
	PetscErrorCode ierr;
	PetscBool flaga,flagb;
	PetscViewer fd;
	PetscInt size,sizea;
	PetscScalar scal;


	ierr=PetscOptionsGetString(PETSC_NULL,"-mfile",filea,PETSC_MAX_PATH_LEN-1,&flaga);CHKERRQ(ierr);
	if (!flaga) {
		sprintf(err,"Error : mfile is not properly set -> %s\n",filea);
		SETERRQ(*communicator,(PetscErrorCode) 83,err);
	}

	/* read matrix file */
	PetscPrintf(PETSC_COMM_WORLD,"Loading Matrix : %s\n",filea);
//	PetscPrintf(PETSC_COMM_WORLD,"I am right here \n");
	ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,filea,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	//PetscPrintf(PETSC_COMM_WORLD," PetscViewerBinaryOpen gone good \n");
	ierr = MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
/*	PetscPrintf(PETSC_COMM_WORLD," STEP 1 IS OK\n");*/
//	ierr = MatSetFromOptions(*A);CHKERRQ(ierr);	   ierr =  MatSetType(*A,MATSEQAIJ); CHKERRQ(ierr);
//	PetscPrintf(PETSC_COMM_WORLD," STEP 2 IS OK \n");
	ierr=MatLoad(*A,fd);CHKERRQ(ierr);
	//PetscPrintf(PETSC_COMM_WORLD,"MatLoad gone OK\n");
	ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ierr=MatGetSize(*A,&size,&sizea);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Loaded Matrix of size : %d %d\n",size,sizea);

	ierr=PetscOptionsGetString(PETSC_NULL,"-vfile",fileb,PETSC_MAX_PATH_LEN-1,&flagb);CHKERRQ(ierr);
	if (!flagb) {
		sprintf(err,"Error : vfile is not properly set\n");
		//SETERRQ((PetscErrorCode) 83,err);
	}

	if (!flagb) {
		/* the user did not provide a vector, so generate it*/
		PetscPrintf(PETSC_COMM_WORLD,"Generating Vector \n");
		VecCreate(PETSC_COMM_WORLD,v);		// we create a PetscObject Vec 
		VecSetSizes(*v,PETSC_DECIDE,size); 	// we set it size to size according to size of the matrix A
		VecSetFromOptions(*v);				// setting the object Vec 
		scal=1.0/size;						
		VecSet(*v,scal); 					// we set all the elements of the vector to the single value "scal = 1/size"
		PetscPrintf(PETSC_COMM_WORLD,"Generated Vector of size : %d\n",size);
	} else {
		PetscPrintf(PETSC_COMM_WORLD,"Loading Vector : %s\n",fileb);
		ierr=PetscViewerBinaryOpen(PETSC_COMM_WORLD,fileb,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		VecCreate(PETSC_COMM_WORLD,v);
		ierr=VecLoad(*v,fd);CHKERRQ(ierr);
		ierr=PetscViewerDestroy(&fd);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD,"Loaded Vector of size : %d\n",size);
	}

	return 0;
}
