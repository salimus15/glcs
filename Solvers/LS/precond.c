/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "precond.h"


PetscErrorCode LSPrecond(PetscReal a_ell, PetscReal d_ell, PetscReal c_ell,
		PetscScalar * eta, PetscScalar * alpha, PetscScalar * beta,
		PetscScalar * delta, PetscScalar * c, PetscScalar * d,
		PetscInt * mu, PetscInt * nb_eigen, PetscInt * min_eigen, PetscInt * nb_eigen_all){

	PetscScalar ** gamma, * fact_tmp, ** mm_tmp;
	PetscScalar * res;
	Mat MM, F,fact;
	KSP ksplsqr, kspchol;
	PC  pcchol, pclsqr;
	Vec rhs,soln;
	int i,j,k,nu;
	PetscErrorCode ierr;


	/* allocate work array gamma and mm_tmp, init it to 0 */
	ierr=PetscMalloc(((*nb_eigen_all)+3)*sizeof(PetscScalar*),&gamma);CHKERRQ(ierr);
	for(i=0;i<(*nb_eigen_all)+3;i++){
		ierr=PetscMalloc(((*nb_eigen_all)+3)*sizeof(PetscScalar),&(gamma[i]));CHKERRQ(ierr);
			for(j=0;j<(*nb_eigen_all)+3;j++)
				gamma[i][j]=(PetscScalar)0.0;
	}

	ierr=PetscMalloc(((*nb_eigen_all)+1)*sizeof(PetscScalar*),&mm_tmp);CHKERRQ(ierr);
	for(i=0;i<(*nb_eigen_all)+1;i++){
		ierr=PetscMalloc(((*nb_eigen_all)+1)*sizeof(PetscScalar),&mm_tmp[i]);CHKERRQ(ierr);
		for(j=0;j<(*nb_eigen_all)+1;j++)
		  mm_tmp[i][j]=(PetscScalar)0.0;
	}



	for(i=0;i<*nb_eigen_all;i++){
	  beta[i]=(PetscScalar)0.0;
	  delta[i]=(PetscScalar)0.0;
	}

	/* begin computations */
	beta[0]=(PetscScalar)a_ell/2.0;
 	 *nb_eigen=*nb_eigen_all;
	for(i=0; i<*nb_eigen_all;i++){
		delta[i+1]=((PetscScalar)(d_ell))/(4*beta[i]);
		beta[i+1]=(PetscScalar)(a_ell) - delta[i+1];
		if(cabs(beta[i+1])<epsilon()){
			*nb_eigen=i;
			if(*nb_eigen<*min_eigen)
				return 0;
		}
	}
	*alpha=c_ell;

/*	#ifdef DEBUG*/
/*	for(i=0;i<*nb_eigen_all;i++){*/
/*	  PetscPrintf(PETSC_COMM_WORLD,"@} LSQR beta[%d]=%e delta[%d]=%e \n",i,(PetscReal)beta[i],i,(PetscReal)delta[i]);*/
/*	}*/
/*	PetscPrintf(PETSC_COMM_WORLD,"@} LSQR => COMPUTING M MATRIX\n");*/
/*	#endif*/



	/* computation of gamma */
	for(nu=0;nu<*mu;nu++){

		for(i=0;i<*nb_eigen_all+1;i++)
			  for(j=0;j<*nb_eigen_all+1;j++)
			      gamma[i][j]=(PetscScalar)0.0+PETSC_i*0.0;

		gamma[1][1]=(PetscScalar)1.0+PETSC_i*0.0;
		for(j=1; j<=*nb_eigen_all;j++){
			for(i=1;i<=j+1;i++){


				gamma[i][j+1]=(PetscReal)(
				  (PetscReal)(PetscRealPart(d[nu])/2.0*(PetscRealPart(gamma[i+1][j])+PetscRealPart(gamma[i-1][j]))
				- (PetscReal)PetscImaginaryPart(d[nu])/2*(PetscReal)(PetscImaginaryPart(gamma[i+1][j])+PetscImaginaryPart(gamma[i-1][j]))
				+ (PetscReal)(PetscRealPart(c[nu])-*alpha)*(PetscReal)PetscRealPart(gamma[i][j])
				-(PetscReal) PetscImaginaryPart(c[nu])*(PetscReal)PetscImaginaryPart(gamma[i][j])
				- (PetscReal)delta[j-1]*PetscRealPart(gamma[i][j-1]))/(PetscReal) beta[j-1])+PETSC_i*PetscImaginaryPart(gamma[i][j+1]);


				gamma[i][j+1]=(PetscScalar)PetscRealPart(gamma[i][j+1])+PETSC_i*(
				(PetscReal)(PetscRealPart(d[nu])/2.0*(PetscReal)(PetscImaginaryPart(gamma[i+1][j])+PetscImaginaryPart(gamma[i-1][j]))
				+ (PetscReal)PetscImaginaryPart(d[nu])/2.0*((PetscReal)PetscRealPart(gamma[i+1][j])+PetscRealPart(gamma[i-1][j]))
				+ (PetscReal)(PetscRealPart(c[nu])-*alpha)*(PetscReal)PetscImaginaryPart(gamma[i][j])
				+ (PetscReal)PetscImaginaryPart(c[nu])*(PetscReal)PetscRealPart(gamma[i][j])
				- (PetscReal)delta[j-1]*(PetscReal)PetscImaginaryPart(gamma[i][j-1]))/ (PetscReal)beta[j-1]);

			}

			gamma[0][j+1] = gamma[2][j+1];


		}

		/* computation of MM */
		for(j=0;j<=*nb_eigen_all;j++){
			for(i = 0; i <=j; i++){

				mm_tmp[i][j]=(PetscReal)mm_tmp[i][j]+4.*((PetscReal)(PetscRealPart(gamma[1][j+1])*PetscRealPart(gamma[1][i+1])
							    +PetscImaginaryPart(gamma[1][j+1])*PetscImaginaryPart(gamma[1][i+1])));

				if(i>1){
					for(k=1;k<i;k++){
						mm_tmp[i][j]=(PetscReal)mm_tmp[i][j]+2.*((PetscReal)(PetscRealPart(gamma[k+1][j+1])*PetscRealPart(gamma[k+1][i+1])
						+PetscImaginaryPart(gamma[k+1][j+1])*PetscImaginaryPart(gamma[k+1][i+1])));
					}
				}
			}
		}
	}




	MatCreateSeqDense(PETSC_COMM_WORLD,(*nb_eigen_all)+1,(*nb_eigen_all)+1,PETSC_NULL,&MM);
	MatSetFromOptions(MM);

	/* remplissage de la partie triangulaire inférieure */
	for(j=0;j<=*nb_eigen_all;j++){
		for(i=0;i<=j;i++){
			MatSetValue(MM,i,j,(PetscReal)mm_tmp[i][j],INSERT_VALUES);
			MatSetValue(MM,j,i,(PetscReal)mm_tmp[i][j],INSERT_VALUES);
		}
	}

	MatAssemblyBegin(MM,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(MM,MAT_FINAL_ASSEMBLY);

/*	#ifdef DEBUGDATA*/
/*	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB );*/
/* 	MatView(MM,PETSC_VIEWER_STDOUT_WORLD);*/
/*	#endif*/

/*	#ifdef DEBUG*/
/*	PetscPrintf(PETSC_COMM_WORLD,"@} LSQR => PROCESSING TO CHOLESKY FACTORIZATION\n");*/
/*	#endif*/

	/*proceed to factorization*/
	ierr=KSPCreate(PETSC_COMM_WORLD,&kspchol);CHKERRQ(ierr);
	ierr=KSPSetOperators(kspchol,MM,MM);CHKERRQ(ierr);
	ierr=KSPSetType(kspchol,KSPPREONLY);CHKERRQ(ierr);
	ierr=KSPGetPC(kspchol,&pcchol);CHKERRQ(ierr);
	ierr=PCSetType(pcchol,PCCHOLESKY);CHKERRQ(ierr);
	ierr=KSPSetUp(kspchol);CHKERRQ(ierr);

/*	#ifdef DEBUG*/
/*	PetscPrintf(PETSC_COMM_WORLD,"@} LSQR => CHOLESKY DONE, EXTRACTING FACTORIZED MATRIX\n");*/
/*	#endif*/
	/*get factor matrix*/
	MatCreateSeqDense(PETSC_COMM_WORLD,(*nb_eigen_all)+1,(*nb_eigen_all)+1,PETSC_NULL,&fact);
	// MatSetFromOptions(fact);
	MatAssemblyBegin(fact,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(fact,MAT_FINAL_ASSEMBLY);

	ierr=PCFactorGetMatrix(pcchol,&fact); CHKERRQ(ierr);
	#ifdef DEBUGDATA
	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB );
 	MatView(fact,PETSC_VIEWER_STDOUT_WORLD);
	#endif

/*	#ifdef DEBUG*/
/*	PetscPrintf(PETSC_COMM_WORLD,"@} LSQR => CREATING F MATRIX\n");*/
/*	#endif*/

	/*Create the matrix operator that will be used in the QR factorization*/
	MatCreateSeqDense(PETSC_COMM_WORLD,(*nb_eigen)+1,(*nb_eigen),PETSC_NULL,&F);

/*	#ifdef DEBUG*/
/*	for(j=0;j<*nb_eigen;j++){*/
/*	  PetscPrintf(PETSC_COMM_WORLD,"@} LSQR beta[%d]=%e delta[%d]=%e alpha=%e\n",*/
/*	    j,(PetscReal)beta[j],j,(PetscReal)delta[j],(PetscReal)*alpha*/
/*	  );*/
/*	}*/
/*	#endif*/

	/*get matrix array, fact is of dense format so wouldn't be a problem for addressing*/
	ierr=MatSeqAIJGetArray(fact,&fact_tmp);CHKERRQ(ierr);
	MatSetValue(F,0,0,(PetscReal)(*alpha)*(PetscReal)fact_tmp[0]+(PetscReal)beta[0]*(PetscReal)fact_tmp[1],INSERT_VALUES);
	MatSetValue(F,1,0,(PetscReal)beta[0]*(PetscReal)fact_tmp[1+((*nb_eigen_all)+1)],INSERT_VALUES);

	for(j=1;j<*nb_eigen;j++){
	  for(i=0;i<j;i++){
	    ierr=MatSetValue(F,i,j,(PetscReal)delta[j-1]*(PetscReal)fact_tmp[(j-1)+i*((*nb_eigen_all)+1)]
			      +(PetscReal)(*alpha)*(PetscReal)fact_tmp[j+i*((*nb_eigen_all)+1)]
			      +(PetscReal)(beta[j])*(PetscReal)fact_tmp[(j+1)+i*((*nb_eigen_all)+1)],INSERT_VALUES);CHKERRQ(ierr);
	  }

	  ierr=MatSetValue(F,j,j,(PetscReal)(*alpha)*(PetscReal)fact_tmp[j+j*((*nb_eigen_all)+1)]
			    +(PetscReal)beta[j]*(PetscReal)fact_tmp[(j+1)+j*((*nb_eigen_all)+1)],INSERT_VALUES);CHKERRQ(ierr);
	  if(j+1<=*nb_eigen){
	    ierr=MatSetValue(F,j+1,j,(PetscReal)beta[j]*(PetscReal)fact_tmp[(j+1)+(j+1)*((*nb_eigen_all)+1)],INSERT_VALUES);CHKERRQ(ierr);
	  }
	}
	/* set the vectors*/
	MatGetVecs(F,&soln,&rhs);

	/*set the solution to zero*/
	ierr=VecSet(soln,(PetscScalar)0.0);CHKERRQ(ierr);
 	ierr=VecSet(rhs,(PetscScalar)0.0);CHKERRQ(ierr);

	/* rhs[0] must be setted to beta*/
	ierr=VecSetValue(rhs,0,(PetscReal)fact_tmp[0],INSERT_VALUES);CHKERRQ(ierr);

	ierr = VecAssemblyBegin(soln);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(soln);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(rhs);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(rhs);CHKERRQ(ierr);

	/*no longer need to access factored matrix, restore it*/
	ierr=MatSeqAIJRestoreArray(fact,&fact_tmp);CHKERRQ(ierr);

	/*assemble F for processing*/
	MatAssemblyBegin(F,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(F,MAT_FINAL_ASSEMBLY);

	#ifdef DEBUG
	int mlsizex,mlsizey,mgsizex,mgsizey;
	int vssizel,vssizeg,vrsizel,vrsizeg;
	MatGetSize(F,&mgsizex,&mgsizey);
	MatGetLocalSize(F,&mlsizex,&mlsizey);
	VecGetSize(soln,&vssizeg);
	VecGetSize(rhs,&vrsizeg);
	VecGetLocalSize(soln,&vssizel);
	VecGetLocalSize(rhs,&vrsizel);

	PetscPrintf(PETSC_COMM_WORLD,"@} LSQR Sizes Mat(%d ; %d) = (%d ; %d) Vec sol = (%d) (%d) Vec rhs = (%d) (%d)\n",
	    mgsizex,mgsizey,mlsizex,mlsizey,vssizeg,vssizel,vrsizeg,vrsizel
	);
	#endif
	#ifdef DEBUGDATA
 	PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_MATLAB );
 	MatView(F,PETSC_VIEWER_STDOUT_WORLD);
	#endif

	#ifdef DEBUG
	PetscPrintf(PETSC_COMM_WORLD,"@} LSQR => COMPUTE LSQR\n");
	#endif
	/*create the lsqr solver context and set it up*/
	ierr=KSPCreate(PETSC_COMM_WORLD,&ksplsqr);CHKERRQ(ierr);
	ierr=KSPSetOperators(ksplsqr,F,F);CHKERRQ(ierr);
	ierr=KSPGetPC(ksplsqr,&pclsqr);CHKERRQ(ierr);
	ierr=PCSetType(pclsqr,PCNONE);CHKERRQ(ierr);
	ierr=KSPSetType(ksplsqr,KSPLSQR);CHKERRQ(ierr);
	ierr=KSPSetInitialGuessNonzero(ksplsqr,PETSC_TRUE);CHKERRQ(ierr);
	ierr=KSPSetUp(ksplsqr);CHKERRQ(ierr);


	/* now we are ready to kick some ass and chew some bubble gum
	   unfortunately i'm all out of gum */
	ierr = KSPSolve(ksplsqr, rhs, soln); CHKERRQ(ierr);
	int its;
	KSPGetIterationNumber(ksplsqr,&its);
	#ifdef DEBUG
	PetscPrintf(PETSC_COMM_WORLD,"@} LSQR => LSQR COMPUTED\n");
	PetscPrintf(PETSC_COMM_WORLD,"@} LSQR Number of LSQR iterations = %d\n",its);
	#endif
	/* extract solution elements and place them into eta*/
	ierr = VecGetArray(soln,&res);

 	for(i=0;i<*nb_eigen_all;i++){
	  eta[i]=(PetscScalar)res[i];

	  #ifdef DEBUG
 	  PetscPrintf(PETSC_COMM_WORLD,"@} %d LSQR eta[%d]=%e beta[%d]=%e delta[%d]=%e \n",i,(PetscReal)eta[i],i,(PetscReal)beta[i],i,(PetscReal)delta[i]);
	  #endif
	}

	ierr = VecRestoreArray(soln,&res);


	KSPDestroy(&kspchol);
	KSPDestroy(&ksplsqr);
	VecDestroy(&rhs);
	VecDestroy(&soln);
	MatDestroy(&F);
	MatDestroy(&MM);

	/* free work array gamma */
	for(i=0;i<(*nb_eigen_all)+2;i++)
		PetscFree(gamma[i]);
	PetscFree(gamma);

	for(i=0;i<(*nb_eigen_all)+1;i++)
		PetscFree(mm_tmp[i]);
	PetscFree(mm_tmp);

	return 0;
}










