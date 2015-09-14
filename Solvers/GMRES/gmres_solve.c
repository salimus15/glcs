/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "gmres_solve.h"

#undef __FUNCT__
#define __FUNCT__ "MyKSPFGMRESResidual"
static PetscErrorCode MyKSPFGMRESResidual(KSP ksp)
{
  KSP_FGMRES     *fgmres = (KSP_FGMRES*)(ksp->data);
  Mat            Amat,Pmat;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PCGetOperators(ksp->pc,&Amat,&Pmat);CHKERRQ(ierr);

  /* put A*x into VEC_TEMP */
  ierr = KSP_MatMult(ksp,Amat,ksp->vec_sol,VEC_TEMP);CHKERRQ(ierr);
  /* now put residual (-A*x + f) into vec_vv(0) */
  ierr = VecWAXPY(VEC_VV(0),-1.0,VEC_TEMP,ksp->vec_rhs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MyKSPSolve_FGMRES"
PetscErrorCode MyKSPSolve_FGMRES(KSP ksp,com_lsa * com)
{
  PetscErrorCode ierr;
  PetscInt       cycle_its = 0; /* iterations done in a call to KSPFGMRESCycle */
  KSP_FGMRES     *fgmres   = (KSP_FGMRES*)ksp->data;
  PetscBool      diagonalscale, flag;
  Vec 		  vec_tmp, vec_t;
  int 		  type = -1, taille =0;
  Mat 		  Amat,Pmat;
  PetscInt nbr_prec =0;
  
  PetscFunctionBegin;
	
	ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_k_param",&nbr_prec,&flag);CHKERRQ(ierr);
	
	ierr = PCGetDiagonalScale(ksp->pc,&diagonalscale);CHKERRQ(ierr);
  if (diagonalscale) SETERRQ1(PetscObjectComm((PetscObject)ksp),PETSC_ERR_SUP,"Krylov method %s does not support diagonal scaling",((PetscObject)ksp)->type_name);

  ierr     = PetscObjectSAWsTakeAccess((PetscObject)ksp);CHKERRQ(ierr);
  ksp->its = 0;
  ierr     = PetscObjectSAWsGrantAccess((PetscObject)ksp);CHKERRQ(ierr);

  /* Compute the initial (NOT preconditioned) residual */
  if (!ksp->guess_zero) {
    ierr = MyKSPFGMRESResidual(ksp);CHKERRQ(ierr);
  } else { /* guess is 0 so residual is F (which is in ksp->vec_rhs) */
    ierr = VecCopy(ksp->vec_rhs,VEC_VV(0));CHKERRQ(ierr);
  }
  /* now the residual is in VEC_VV(0) - which is what
     KSPFGMRESCycle expects... */

  ierr = MyKSPFGMRESCycle(&cycle_its,ksp);CHKERRQ(ierr);
  while (!ksp->reason) {
  	  if(nbr_prec > 0){
		 if(!GmresLSAPrecond(com,ksp)){
		   PetscPrintf(PETSC_COMM_WORLD,"Préconditionnement LSA à %d itérations\n",ksp->its);
		  	nbr_prec--; 
		 }
		}
    ierr = MyKSPFGMRESResidual(ksp);CHKERRQ(ierr);
    if (ksp->its >= ksp->max_it) break;
    ierr = MyKSPFGMRESCycle(&cycle_its,ksp);CHKERRQ(ierr);
    
/*    if(!mpi_lsa_com_type_recv(com,&type)){*/
/*	  if(type==911){*/
		//   PetscPrintf(PETSC_COMM_WORLD,"$} GMRES received SOS message !!!!!!!!!!!'\n");
     /*  this is some deprecated especialy for GMRES due to the fact it is expensive and generates a copy the solution vector 
     but i have ot the choice for now and it is performed only when Anoldi has 0 converged values So lets do it
	*/
 		ierr = VecDuplicate(ksp->vec_sol,&vec_tmp);CHKERRQ(ierr);
// 		ierr = VecDuplicate(ksp->vec_sol,&vec_t);CHKERRQ(ierr);
//		ierr = VecCopy(ksp->vec_sol,vec_tmp);CHKERRQ(ierr);
// 		ierr = PCGetOperators(ksp->pc,&Amat,&Pmat);CHKERRQ(ierr);
//		KSPBuildSolution(ksp,vec_tmp,NULL);
		ierr=KSPBuildResidual(ksp,NULL,NULL,&vec_tmp);CHKERRQ(ierr);
//		ierr=MatMult(Amat,vec_tmp,vec_t);CHKERRQ(ierr);
//		ierr=VecAYPX(vec_t,-1.,ksp->vec_rhs);CHKERRQ(ierr);
		
		
		/* Step of Sending  to Arnoldi Under some conditions*/
		ierr=VecGetSize(vec_tmp,&taille);CHKERRQ(ierr);
		PetscPrintf(PETSC_COMM_WORLD," GMRES  THERE IS %d DATA TO SEND \n",taille);
		mpi_lsa_com_vec_send(com,&vec_tmp);
		   
/*	  }*/
/*	}*/
  }
  /* mark lack of convergence */
  if (ksp->its >= ksp->max_it && !ksp->reason) ksp->reason = KSP_DIVERGED_ITS;
  PetscFunctionReturn(0);
}

// #undef __FUNCT__
// #define __FUNCT__ "MyKSPSolve_FGMRES"
// PetscErrorCode MyKSPSolve_FGMRES(KSP ksp,com_lsa * com)
// {
//   PetscErrorCode ierr;
//   PetscInt       cycle_its = 0; /* iterations done in a call to FGMREScycle */
//   KSP_FGMRES     *fgmres = (KSP_FGMRES *)ksp->data;
//   PetscBool     diagonalscale;


//   PetscFunctionBegin;

//   ierr    = PCDiagonalScale(ksp->pc,&diagonalscale);CHKERRQ(ierr);
//   if (diagonalscale) SETERRQ1(PETSC_ERR_SUP,"Krylov method %s does not support diagonal scaling",((PetscObject)ksp)->type_name);
//   if (ksp->normtype != KSP_NORM_UNPRECONDITIONED) SETERRQ(PETSC_ERR_ARG_WRONGSTATE,"Can only use FGMRES with unpreconditioned residual (it is coded with right preconditioning)");

//   ierr = PetscObjectTakeAccess(ksp);CHKERRQ(ierr);
//   ksp->its = 0;
//   ierr = PetscObjectGrantAccess(ksp);CHKERRQ(ierr);


//   /* allocate preconditioning data */

//   /* Compute the initial (NOT preconditioned) residual */
//   if (!ksp->guess_zero) {
//     ierr = MyFGMRESResidual(ksp);CHKERRQ(ierr);
//   } else { /* guess is 0 so residual is F (which is in ksp->vec_rhs) */
//     ierr = VecCopy(ksp->vec_rhs,VEC_VV(0));CHKERRQ(ierr);
//   }
//   /* now the residual is in VEC_VV(0) - which is what
//      FGMREScycle expects... */

//   ierr    = MyFGMREScycle(&cycle_its,ksp);CHKERRQ(ierr);

//   while (!ksp->reason) {

// 	if(!GmresLSAPrecond(com,ksp))
// 	  PetscPrintf(PETSC_COMM_WORLD,"Préconditionnement LSA à %d itérations\n",ksp->its);


//   	ierr     = MyFGMRESResidual(ksp);CHKERRQ(ierr);
//   	if (ksp->its >= ksp->max_it) break;


//   	ierr     = GMREScycle(&cycle_its,ksp);CHKERRQ(ierr);


// // 	#ifdef DEBUG
// // 	PetscPrintf(PETSC_COMM_WORLD,"#} GMRES Sending residual to Arnoldi\n");
// // 	#endif
// // 	mpi_lsa_com_vec_send(com,&(ksp->vec_sol));
// 	}
//   /* mark lack of convergence */
//   if (ksp->its >= ksp->max_it && !ksp->reason) ksp->reason = KSP_DIVERGED_ITS;

//   if(ksp->heur_orth==HRAND)
//     PetscRandomDestroy(myRand);


//   PetscFunctionReturn(0);
// }

