/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "gmres_precond.h"



static int latency_count=0;

PetscErrorCode GmresLSAPrecond(com_lsa * com, KSP ksp)
{
  KSP_FGMRES     *fgmres = (KSP_FGMRES *)(ksp->data);
  Mat            Amat,Pmat;
  PetscErrorCode ierr;
  Vec r0_tmp,w0_tmp,w1_tmp,x_tmp,w_1_tmp,r1_tmp,sol_tmp,vec_tmp;
  PetscScalar data_tmp[EIGEN_ALL*2*3];
  PetscInt size_data,ls_power,latency,hang,timing;
  PetscInt nols;
  PetscBool flag;
  int i,j;

  PetscScalar alpha;
  PetscScalar * eta,*delta,*beta;
  PetscReal norm;


  PetscFunctionBegin;

  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_power",&ls_power,&flag);CHKERRQ(ierr);
  if(!flag) ls_power=LS_POWER;
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_latency",&latency,&flag);CHKERRQ(ierr);
  if(!flag) latency=LS_LATENCY;
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_m_hang",&hang,&flag);CHKERRQ(ierr);
  if(!flag) hang=ksp->max_it;
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_nols",&nols,&flag);CHKERRQ(ierr);
  if(!flag) nols=1;
  ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_timing",&timing,&flag);CHKERRQ(ierr);
  if(!flag) timing=60;

  latency_count++;

  if((latency_count%latency==0  && ksp->its>0) || ksp->its==2){
    #ifdef DEBUG
    printf("#} %d GMRESLSPrecond at %d its must wait %d seconds (soit %d minutes et %d secondes) before using LSQR\n",com->rank_world, ksp->its,timing,timing/60,timing%60);
    #endif
    sleep(timing);
  }



  /* received something ? */
  if(nols==0||mpi_lsa_com_array_recv(com,&size_data,data_tmp)){
    return 1;
  }
//  #ifdef DEBUG
  else
    printf("#}%d GMRESLSPrecond Received data from LSQR of size %d\n", com->rank_world,(PetscInt)data_tmp[0]);
 // #endif


  /* is data consistent ? */
  if((((PetscInt)data_tmp[0])*3+2)!=size_data){
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond data is not consistent : size mpi %d size data*3+2 %d %e data\n",
		size_data,(((PetscInt)data_tmp[0])*3+2),((PetscReal)data_tmp[0]));

    for(i=0;i<size_data;i++)
	PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond data[%d]=%e\n",i,(PetscReal)data_tmp[i]);

     return 1;
  }

  PetscMalloc(sizeof(PetscScalar)*(PetscInt)data_tmp[0],&eta);
  PetscMalloc(sizeof(PetscScalar)*(PetscInt)data_tmp[0],&beta);
  PetscMalloc(sizeof(PetscScalar)*(PetscInt)data_tmp[0],&delta);



  ierr=PetscMemcpy(&alpha,&data_tmp[1],1*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr=PetscMemcpy(eta,&data_tmp[2],((PetscInt)data_tmp[0])*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr=PetscMemcpy(beta,&data_tmp[2+((PetscInt)data_tmp[0])],((PetscInt)data_tmp[0])*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr=PetscMemcpy(delta,&data_tmp[2+2*((PetscInt)data_tmp[0])],((PetscInt)data_tmp[0])*sizeof(PetscScalar));CHKERRQ(ierr);


  /*get operators and vector data*/
  ierr = PCGetOperators(ksp->pc,&Amat,&Pmat);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&x_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&r0_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&w_1_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&w1_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&w0_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&r1_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&sol_tmp);CHKERRQ(ierr);
  ierr=VecDuplicate(ksp->vec_sol,&vec_tmp);CHKERRQ(ierr);


  ierr=VecSet(x_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(r0_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(w0_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(w1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(w_1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(r1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecSet(vec_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  ierr=VecCopy(ksp->vec_sol,sol_tmp);CHKERRQ(ierr);
  /* proceed to computation*/

    #ifdef DEBUG
    VecNorm(sol_tmp,NORM_2,&norm);
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond Sol norm %e\n",norm);
    #endif
    #ifdef DEBUGDATA
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB );
    VecView(sol_tmp,PETSC_VIEWER_STDOUT_WORLD);
    #endif

      #ifdef DEBUG
        for(i=0;i<(PetscInt)data_tmp[0];i++){
    	  PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond ne[%d] eta[%d]=%e %e beta[%d]=%e %e delta[%d]=%e %e\n",i,
		      i,PetscRealPart(eta[i]),PetscImaginaryPart(eta[i]),
		      i,PetscRealPart(beta[i]),PetscImaginaryPart(beta[i]),
		      i,PetscRealPart(delta[i]),PetscImaginaryPart(delta[i]));

	}
  #endif

  for(j=0;j<ls_power;j++){
    /* r0 = b-Ax*/
    /* put A*x into VEC_TEMP */
    ierr = MatMult(Amat,sol_tmp,vec_tmp);CHKERRQ(ierr);
    /* now put residual (-A*x + f) into vec_vv(0) */
    ierr = VecWAXPY(r0_tmp,-1.0,vec_tmp,ksp->vec_rhs);CHKERRQ(ierr);

    #ifdef DEBUG
    VecNorm(r0_tmp,NORM_2,&norm);
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond [%d] r0 norm %e\n",j,norm);
    #endif

    /* r0 = w0*/
    ierr=VecCopy(r0_tmp,w0_tmp);CHKERRQ(ierr);
     ierr=VecCopy(w0_tmp,x_tmp);CHKERRQ(ierr);
//     ierr=VecSet(x_tmp,0.0);CHKERRQ(ierr);
    /* x = eta[0] * w0 */
//     ierr=VecSet(x_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
//     ierr=VecAXPY(x_tmp,eta[0],w0_tmp);CHKERRQ(ierr);
    ierr=VecScale(x_tmp,eta[0]);CHKERRQ(ierr);
    ierr=VecSet(w_1_tmp,(PetscScalar)0.0);CHKERRQ(ierr);
  #ifdef DEBUGDATA
    VecNorm(x_tmp,NORM_2,&norm);
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond Temp [%d] x_tmp 1 norm %e (%e)\n",j,i,norm,PetscRealPart(eta[0]));
#endif
    /* depending ton the ls polynom size (PetscInt)data_tmp[0] */
    for(i=0;i<(PetscInt)data_tmp[0]-1;i++){
      /* w1=-alpha*w0 - delta[i]*w_1 ((  y = alpha x + delta y. )) (Vec y,PetscScalar alpha,PetscScalar beta,Vec x)*/
      ierr=VecCopy(w_1_tmp,w1_tmp);CHKERRQ(ierr);
      ierr=VecAXPBY(w1_tmp,-alpha,-(PetscScalar)delta[i],w0_tmp);CHKERRQ(ierr);
      #ifdef DEBUGDATA
      VecNorm(w1_tmp,NORM_2,&norm);
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond Temp [%d][%d] w1_tmp 1 norm %e (%e %e)\n",j,i,norm,PetscRealPart(-alpha),PetscRealPart(-delta[i]));
#endif
      /* w1 = w1 - A*w0 */
      ierr = MatMult(Amat,w0_tmp,vec_tmp);CHKERRQ(ierr);
      /* y = alpha x + y.  VecAXPY(Vec y,PetscScalar alpha,Vec x)*/
      ierr = VecAXPY(w1_tmp,1.0,vec_tmp);CHKERRQ(ierr);
      #ifdef DEBUGDATA
      VecNorm(w1_tmp,NORM_2,&norm);
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond Temp [%d][%d] w1_tmp 2 norm %e\n",j,i,norm);
#endif
      /* w1 = w1/beta[i] */
      ierr=VecScale(w1_tmp,1/beta[i]);CHKERRQ(ierr);
      #ifdef DEBUGDATA
      VecNorm(w1_tmp,NORM_2,&norm);
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond Temp [%d][%d] w1_tmp 3 norm %e (%e)\n",j,i,norm,PetscRealPart(1/beta[i]));
#endif
      /* w_1 = w0 */
      ierr=VecCopy(w0_tmp,w_1_tmp);CHKERRQ(ierr);

      /* w0 = w1 */
      ierr=VecCopy(w1_tmp,w0_tmp);CHKERRQ(ierr);

      /* x = x + (w0 * eta[i] ) */
      ierr=VecAXPY(x_tmp,eta[i+1],w0_tmp);CHKERRQ(ierr);
      #ifdef DEBUGDATA
      VecNorm(x_tmp,NORM_2,&norm);
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond Temp [%d][%d] x_tmp norm %e (%e)\n",j,i,norm,PetscRealPart(delta[i]));
#endif
    }

    /* x1= x1+x*/
    ierr=VecAXPY(sol_tmp,1.0,x_tmp);CHKERRQ(ierr);

    /* put A*x into VEC_TEMP */
    ierr = MatMult(Amat,sol_tmp,vec_tmp);CHKERRQ(ierr);
    /* now put residual (-A*x + f) into vec_vv(0) */
    ierr = VecWAXPY(r1_tmp,-1.0,vec_tmp,ksp->vec_rhs);CHKERRQ(ierr);

    /* compute norm and see if it's below epsilon */
    VecNorm(r1_tmp,NORM_2,&norm);

    #ifdef DEBUG
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond norm %e\n",norm);
    #endif
    #ifdef DEBUGDATA
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB );
    VecView(sol_tmp,PETSC_VIEWER_STDOUT_WORLD);
    #endif

    if(norm>epsilon()){
      PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond norm %e is over epsilon %e\n",norm,epsilon());
// FIXME:       PetscFunctionReturn(2);
    }

  }



  if(latency_count%latency!=0 || ksp->its==3){
    #ifdef DEBUG
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond must wait before using LSQR %d/%d mod %d=%d\n",ksp->its,fgmres->max_k,latency,(ksp->its/fgmres->max_k)%latency);
    #endif
     return 1;
  }

  #ifdef DEBUG
    PetscPrintf(PETSC_COMM_WORLD,"#} GMRESLSPrecond using vector for preconditionning (power = %d)\n",ls_power);
  #endif

  if(nols!=0) ierr=VecCopy(sol_tmp,ksp->vec_sol);CHKERRQ(ierr);

  ierr=PetscFree(eta);CHKERRQ(ierr);
  ierr=PetscFree(beta);CHKERRQ(ierr);
  ierr=PetscFree(delta);CHKERRQ(ierr);

  return 0;
}













