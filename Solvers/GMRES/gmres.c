/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "gmres.h"

PetscErrorCode launchGMRES(com_lsa * com, Vec * b, Mat * A){
	PetscErrorCode ierr;
	KSP ksp;
	Vec x;
	KSPConvergedReason reason;
	   PetscPrintf(PETSC_COMM_WORLD,"$}### GMRES mr rank : %d my group : %d my color : %d  I send to %d  and receive fom %d\n",com->rank_world, com->com_group,com->color_group,com->in_com,com->out_com);
	PetscPrintf(com->com_group,"#} GMRES Creating and setting vector x\n");

	ierr = VecDuplicate(*b, &x);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
	ierr = VecSet(x, 0.); CHKERRQ(ierr);


	PetscPrintf(com->com_group,"#} GMRES Creating KSP\n");

	ierr = KSPCreate(com->com_group, &ksp);CHKERRQ(ierr);

	ierr = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);

	ierr = KSPSetOperators(ksp, *A, *A);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

	PetscPrintf(com->com_group,"#} GMRES Solving...\n");
/* *******************************************************************
******************************************************************** */

	ierr = MyKSPSolve(ksp, *b, x,com); CHKERRQ(ierr);
	
/* **********************************************************************
********************************************************************** */
	PetscPrintf(com->com_group,"#} MyKSPSolve succesfully ended \n");
	KSPGetConvergedReason(ksp,&reason);
	PetscPrintf(com->com_group,"#} KSPGetConvergedReason  succesfully ended \n");
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	PetscPrintf(com->com_group,"#} KSPDestroy  succesfully ended \n");
	ierr = VecDestroy(&x);CHKERRQ(ierr);

	PetscPrintf(com->com_group,"#} GMRES Linear system Solved with reason : %D \n",reason);


	int exit_type=666;
	mpi_lsa_com_type_send(com,&exit_type);
	PetscPrintf(com->com_group,"#} GMRES Sending Exit message\n");


  return ierr;
}

