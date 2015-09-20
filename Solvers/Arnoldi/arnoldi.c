/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "arnoldi.h"

/* Compute cyclicly eigenvalue */
PetscErrorCode Arnoldi(com_lsa * com, Mat * A, Vec  *v){
	EPS eps; /* eigensolver context */
	char  load_path[PETSC_MAX_PATH_LEN],export_path[PETSC_MAX_PATH_LEN];
	PetscInt end,first,validated;
	PetscErrorCode ierr;
	/* eigenvalues number is set to 100, can be changed if needed
	   we choosed to fix it because mallocs weren't working properly */
	PetscScalar eigenvalues[1000], ei, er;
	PetscReal re,im,vnorm;
	PetscInt eigen_nb,j,i,size,one=1, taille;
	Vec initialv,nullv,*vs;
	PetscBool flag,data_load,data_export,continuous_export,load_any;
	int exit_type=0, counter = 0, l;
	int sos_type = 911;
	PetscScalar *vecteur_initial;
	PetscViewer viewer;	


	PetscBool need_new_init = PETSC_FALSE, exit = PETSC_FALSE;

	sprintf(load_path,"./arnoldi.bin");
	sprintf(export_path,"./arnoldi.bin");
	
	PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
//	 PetscViewerSetType(viewer,PETSCVIEWERBINARY);
	// if (skippheader) { PetscViewerBinarySetSkipHeader(viewer,PETSC_TRUE); }
//	 PetscViewerFileSetMode(viewer,FILE_MODE_WRITE);
//	 PetscViewerBinarySetUseMPIIO(viewer,PETSC_TRUE); 
// 	 PetscViewerFileSetName(viewer,"arnoldidbg.txt");

	/* create the eigensolver */
	ierr=EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
	/* set the matrix operator */
	ierr=EPSSetOperators(eps,*A,PETSC_NULL);
	/* set options */
	ierr=EPSSetType(eps,EPSARNOLDI);
	ierr=EPSSetFromOptions(eps);CHKERRQ(ierr);
	
	/* duplicate vector properties */
	ierr=VecDuplicate(*v,&initialv);CHKERRQ(ierr);
	ierr=VecDuplicate(*v,&nullv);CHKERRQ(ierr);
	ierr=VecSet(nullv,(PetscScalar)0.0);CHKERRQ(ierr);
/*	ierr=VecSet(initialv,(PetscScalar)1.0);CHKERRQ(ierr);*/
	
	ierr=VecSetRandom(initialv,PETSC_NULL);//initialize initial vector to random
	ierr=VecGetSize(initialv,&size);CHKERRQ(ierr);

	ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_eigen",&eigen_nb,&flag);CHKERRQ(ierr);
	if(!flag) eigen_nb=EIGEN_ALL;
/*	ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_arnoldi_load",load_path,PETSC_MAX_PATH_LEN,&data_load);CHKERRQ(ierr);*/
/*	ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_arnoldi_export",export_path,PETSC_MAX_PATH_LEN,&data_export);CHKERRQ(ierr);*/

/*	ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_arnoldi_load_any",&load_any);CHKERRQ(ierr);*/
/*	ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_arnoldi_cexport",&continuous_export);CHKERRQ(ierr);*/

/*	if(load_any) PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi loading default data file\n");*/
/*	PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi path in= %s out= %s\n",load_path,export_path);*/

	PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi allocating buffer of %d for invariant subspace\n",eigen_nb*2);
	vs=malloc(size*sizeof(Vec));
	for(i=0;i<size;i++){
		ierr=VecDuplicate(*v,&vs[i]);CHKERRQ(ierr);
	}
	vecteur_initial = malloc(size * sizeof(PetscScalar));


	end=0;
	first=1;
	validated=1;
	while(!end){
		/*check if the program need to exit */
		if(exit == PETSC_TRUE)
			break;
		/* check if we received an exit message from Father*/
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
		  if(exit_type==666){
		    end=1;
		    PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi Sending Exit message\n");
		    mpi_lsa_com_type_send(com,&exit_type);
		    break;
		  }
		}
		/* check if we received some data from GMRES		*/
/*		if(!mpi_lsa_com_vec_recv(com, &initialv)){*/
/*				printf(" =========   I RECEIVED SOME DATA FROM GMRES ============\n");*/
/*		}*/
		    

		for(j=0;j<eigen_nb;j++){
			eigenvalues[j]=(PetscScalar)0.0;
		}

/*		//FIXME: refactoriser les if suivants + flags file read, c'est très très moche*/
/*		if(data_load&&load_any){*/
/*		  load_any=PETSC_FALSE;*/
/*		  data_load=PETSC_TRUE;*/
/*		}*/
		  ierr = VecAssemblyBegin(initialv);CHKERRQ(ierr);
  		  ierr = VecAssemblyEnd(initialv);CHKERRQ(ierr);
/*		ierr = VecView(initialv,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);*/
	
/*		ierr=VecGetSize(initialv,&taille);CHKERRQ(ierr);*/
/*		PetscPrintf(PETSC_COMM_WORLD,"==== > OUR INITIALV IS OF SIZE %d\n", taille);*/
/*		ierr=VecGetArray(initialv, &vecteur_initial);CHKERRQ(ierr);*/
/*		for (i = 0; i < taille; i++)*/
/*			PetscPrintf(PETSC_COMM_WORLD,"==== > initialv[%d] = %e\n", i, vecteur_initial[i]);*/
/*	*/
	
/*		if(!(data_load^=load_any)){*/
		  ierr=EPSSetInitialSpace(eps,1,&initialv);CHKERRQ(ierr);
/*		} else {*/
/*				PetscPrintf(PETSC_COMM_WORLD,"==== > I AM LOADING DATA FROM FILE\n");*/
/*				PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi Reading file %s\n",load_path);*/
/*				ierr=readBinaryVecArray(load_path,(int*)one,&initialv);CHKERRQ(ierr);*/
/*				data_load=PETSC_FALSE;*/
/*				load_any=PETSC_FALSE;*/
/*				ierr=EPSSetInitialSpace(eps,1,&initialv);CHKERRQ(ierr);*/
/*				PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi Has Read file %s\n",load_path);*/
/*		}*/
		#ifdef DEBUG
/*	  	printf("*} Arnoldi  initial vector initiated\n");*/
/*		ierr=VecGetSize(initialv,&taille);CHKERRQ(ierr);*/
/*		PetscPrintf(PETSC_COMM_WORLD,"==== > OUR INITIALV IS OF SIZE %d\n", taille);*/
/*  		vecteur_initial = realloc(vecteur_initial,taille);  			*/
/*  		ierr=VecGetArray(initialv, &vecteur_initial);CHKERRQ(ierr);*/
/*		for (i = 0; i < taille; i++)*/
/*			PetscPrintf(PETSC_COMM_WORLD,"==== > initialv[%d] = %e\n", i, vecteur_initial[i]);*/
/*		ierr= VecRestoreArray(initialv, &vecteur_initial);CHKERRQ(ierr);*/
/*		*/
	
		 #endif 	

/*		#ifdef DEBUG*/
/*		PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi Launching EPSSolve\n");*/
/*		#endif*/
//		PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi errest = %e \n",eps->errest );
	
/*		ierr = VecView(initialv,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);*/
		/* compute eigenvalues */

//		 ierr=EPSSetUp(eps);CHKERRQ(ierr);
		ierr=EPSSolve(eps);CHKERRQ(ierr);
/*				ierr = VecView(initialv,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);*/
/*		ierr=VecGetSize(initialv,&taille);CHKERRQ(ierr);*/
/*		PetscPrintf(PETSC_COMM_WORLD,"==== > OUR INITIALV IS OF SIZE %d\n", taille);*/
/*  		vecteur_initial = realloc(vecteur_initial,taille);  			*/
/*  		ierr=VecGetArray(initialv, &vecteur_initial);CHKERRQ(ierr);*/
/*		for (i = 0; i < taille; i++)*/
/*			PetscPrintf(PETSC_COMM_WORLD,"==== > initialv[%d] = %e\n", i, vecteur_initial[i]);*/
/*		ierr= VecRestoreArray(initialv, &vecteur_initial);CHKERRQ(ierr);*/
		
		
//		 ierr= EPSVectorsView(eps,viewer);
				/*construct new initial vector*/
			ierr=EPSGetInvariantSubspace(eps, vs);CHKERRQ(ierr);
		++counter;
//		PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi Achieved %d Iterations of EPSSolve\n", counter);
/*		#ifdef DEBUG*/
/*		PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi Achieved EPSSolve\n");*/
/*		#endif*/

		/* get the number of guessed eigenvalues */
		ierr=EPSGetConverged(eps,&eigen_nb);CHKERRQ(ierr);

 		#ifdef DEBUG
		PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi %d converged eigenvalues\n",eigen_nb);
 		#endif

		/* send them */
		for(j=0;j<eigen_nb;j++){
			//EPSGetValue(eps,j,&er,&ei);
			//ierr = EPSGetEigenpair(eps,j,&er,&ei,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
			ierr = EPSGetEigenvalue(eps,j,&er,&ei);CHKERRQ(ierr);
			#ifdef PETSC_USE_COMPLEX
			  re=PetscRealPart(er);
			  im=PetscImaginaryPart(er);
			#else
			  re=er;
			  im=ei;
			#endif

			eigenvalues[j]=(PetscScalar)re+PETSC_i*(PetscScalar)im;


//	 		#ifdef DEBUG
				if(im!=0.0)
				  PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi %d/%d val : %e %e\n",j,eigen_nb,re,im);
				else
				  PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi  %d/%d val : %e\n",j,eigen_nb,er);
//			#endif

		}
		
/*		ierr=VecGetSize(initialv,&taille);CHKERRQ(ierr);*/
/*		PetscPrintf(PETSC_COMM_WORLD,"==== > OUR INITIALV IS OF SIZE %d\n", taille);*/
/*  		vecteur_initial = realloc(vecteur_initial,taille);  			*/
/*  		ierr=VecGetArray(initialv, &vecteur_initial);CHKERRQ(ierr);*/
/*		for (i = 0; i < taille; i++)*/
/*			PetscPrintf(PETSC_COMM_WORLD,"==== > initialv[%d] = %e\n", i, vecteur_initial[i]);*/
/*		ierr= VecRestoreArray(initialv, &vecteur_initial);CHKERRQ(ierr);*/
		
		if( eigen_nb != 0){
/*		#ifdef DEBUG*/
/*		PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi  Sending to LS\n");*/
/*		#endif*/
		/* send the data array */
		mpi_lsa_com_array_send(com, &eigen_nb, eigenvalues);
		
			/*construct new initial vector*/
/*			ierr=EPSGetInvariantSubspace(eps, vs);CHKERRQ(ierr);*/
			ierr=VecCopy(vs[0],initialv);CHKERRQ(ierr);
		
		
			for(j=1;j<eigen_nb;j++){
				ierr=VecAYPX(initialv,(PetscScalar)1.0,vs[j]);
			}
		
			ierr=VecNorm(initialv,NORM_2,&vnorm);CHKERRQ(ierr);
			ierr=VecAYPX(initialv,(PetscScalar)(1.0/vnorm),nullv);CHKERRQ(ierr);
				
	  		if(continuous_export){
/*		  		ierr=writeBinaryVecArray(data_export?export_path:"./arnoldi.bin", 1, &initialv);*/
			}
			if(!mpi_lsa_com_type_recv(com,&exit_type)){
   			if(exit_type==666){
		 			end=1;
		 			PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi Sending Exit message\n");

		 			mpi_lsa_com_type_send(com,&exit_type);
		 			need_new_init = PETSC_FALSE;
		 			exit = PETSC_TRUE;
					 break;
  				}	
			}
  		}else{
  			need_new_init = PETSC_TRUE;
			PetscPrintf(PETSC_COMM_WORLD, "!!! Arnoldi has not converged so we change the initial vector !!!\n"); 
 			while(need_new_init){
				// this was my first try to solve the poblem but it doesn't work better 
				/*ierr=VecSetRandom(initialv,PETSC_NULL);CHKERRQ(ier);*/
				// Now the best to do i think is to develop a kind of help from GMRES 
				// Arnoldi when no convergence observed will send a msg to GMRES like an SOS 
				//and GMRES will send a vector wich will be used to generate a new initial vector that make arnoldi converge **I HOPE **
  				//need_new_init = PETSC_TRUE
				// mpi_lsa_com_type_send(com,&sos_type); // here we send the message 
  				//	 ierr=VecDuplicate(initialv,&vec_tmp_receive);
  				//PetscPrintf(PETSC_COMM_WORLD, "!!! Arnoldi has not converged so we change the initial vector !!!\n");
  				 /* check if there's an incoming message */
				 if(!mpi_lsa_com_vec_recv(com, &initialv)){
					printf(" =========   I RECEIVED SOME DATA FROM GMRES ============\n");
					need_new_init = PETSC_FALSE;
		     	 }else{
		     			if(!mpi_lsa_com_type_recv(com,&exit_type)){
			      			if(exit_type==666){
				    			end=1;
				    			PetscPrintf(PETSC_COMM_WORLD,"*} Arnoldi Sending Exit message\n");

				    			mpi_lsa_com_type_send(com,&exit_type);
				    			need_new_init = PETSC_FALSE;
				    			exit = PETSC_TRUE;
				   			 break;
				  				}
						}
      		 }
			//goto checking;
		     	//return 1;
		     }	
		     if(exit == PETSC_TRUE)
		     		break;		 
  		}
  		
  		// i will check it later 	
		
	}


/*	if(data_export){*/
/*	  ierr=writeBinaryVecArray(export_path, 1, &initialv);*/
/*	}*/


	for(i=0;i<eigen_nb;i++){
		ierr=VecDestroy(&(vs[i]));CHKERRQ(ierr);
	}
	
/*	ierr=PetscFree(vs);CHKERRQ(ierr);*/

	/* and destroy the eps */
	ierr=EPSDestroy(&eps);CHKERRQ(ierr);
	ierr=VecDestroy(&initialv);CHKERRQ(ierr);
	ierr=VecDestroy(&nullv);CHKERRQ(ierr);

	return 0;
}
