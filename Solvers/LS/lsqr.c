/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "lsqr.h"

PetscErrorCode LSQR(com_lsa * com, int * vector_size){
	/* variables */
	PetscInt end,cumul,eigen_received,eigen_total,eigen_max;
	char  load_path[PETSC_MAX_PATH_LEN],export_path[PETSC_MAX_PATH_LEN];
	int i,info,exit_type=0;
	PetscBool flag,data_load,data_export,continuous_export,data_load_any;
	PetscScalar * data,*eigen_cumul,*eigen_tri,*d,*c;
	PetscReal a_ell,c_ell,d_ell,d_reel;
	int data_size;
	PetscInt chsign;
	PetscInt mu1,mu2,mu,result_array_size;
	PetscScalar alpha,*eta,*beta,*delta,scalar_tmp;
	PetscInt ls_eigen_min, ls_eigen; // use default values
	PetscErrorCode ierr;
 	PetscScalar * result_array,*data_buffer;/*[EIGEN_ALL*3+2]*/
	FILE * ftest = NULL;
	int descriptor;
	sprintf(load_path,"./lsqr.bin");
	sprintf(export_path,"./lsqr.bin");
	
	if((ftest = fopen((char *)load_path, "r")) == NULL)
	{
		ierr=PetscBinaryOpen(load_path,FILE_MODE_WRITE,&descriptor);CHKERRQ(ierr);
		ierr=PetscBinaryClose(descriptor);CHKERRQ(ierr);
	}else fclose(ftest);

	/* check if there is arguments for ls */
	ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_eigen_min",&ls_eigen_min,&flag);CHKERRQ(ierr);
	if(!flag) ls_eigen_min=EIGEN_MIN;
	ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_eigen",&ls_eigen,&flag);CHKERRQ(ierr);
	if(!flag) ls_eigen=EIGEN_ALL;
	/* check the number of eigenvalues that one will receive from arnoldi */
  	ierr=PetscOptionsGetInt(PETSC_NULL,"-ksp_ls_k_param",&eigen_max,&flag);CHKERRQ(ierr);
	if(!flag)eigen_max=ls_eigen;

	ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_ls_load",load_path,PETSC_MAX_PATH_LEN,&data_load);CHKERRQ(ierr);
	ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_ls_load_any",&data_load_any);CHKERRQ(ierr);
	ierr=PetscOptionsGetString(PETSC_NULL,"-ksp_ls_export",export_path,PETSC_MAX_PATH_LEN,&data_export);CHKERRQ(ierr);

	ierr=PetscOptionsHasName(PETSC_NULL,"-ksp_ls_cexport",&continuous_export);CHKERRQ(ierr);

/*// 	#ifdef DEBUG*/
/*	    PetscPrintf(PETSC_COMM_WORLD,"@} LSQR args : eigen min %d eigen all %d param k %d\n",ls_eigen_min,ls_eigen,eigen_max);*/
/*	    PetscPrintf(PETSC_COMM_WORLD,"@} LSQR args : loading data from file : %s\n",load_path);*/
/*	    PetscPrintf(PETSC_COMM_WORLD,"@} LSQR args : exporting data to file : %s\n",export_path);*/
/*// 	#endif*/

	/* allocations */
	ierr=PetscMalloc((*vector_size)*sizeof(PetscScalar),&eigen_tri);
	ierr=PetscMalloc((*vector_size)*sizeof(PetscScalar),&eigen_cumul);
	ierr=PetscMalloc(((*vector_size)+1)*sizeof(PetscScalar),&d);
	ierr=PetscMalloc(((*vector_size)+1)*sizeof(PetscScalar),&c);
	ierr=PetscMalloc((*vector_size)*sizeof(PetscScalar),&data);
	ierr=PetscMalloc(eigen_max*sizeof(PetscScalar),&data_buffer);

	/* data that will be sended to GMRES for it's preconditionning step */
	ierr=PetscMalloc((eigen_max+1)*sizeof(PetscScalar),&eta);
	ierr=PetscMalloc((eigen_max+1)*sizeof(PetscScalar),&beta);
	ierr=PetscMalloc((eigen_max+1)*sizeof(PetscScalar),&delta);
	ierr=PetscMalloc((eigen_max*3+2)*sizeof(PetscScalar),&result_array);
	result_array_size=2+3*eigen_max;

	for(i=0;i<(*vector_size);i++){
	  eigen_tri[i]=(PetscScalar)0.0;
	  eigen_cumul[i]=(PetscScalar)0.0;
	}

	for(i=0;i<(*vector_size)+1;i++){
	  d[i]=(PetscScalar)0.0;
	  c[i]=(PetscScalar)0.0;
	}

	cumul=0;
	eigen_received=0;
	eigen_total=0;
	end=0;
	eigen_total=0;
	ls_eigen=0;
	while(!end){
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
		  if(exit_type==666){
		    end=1;
		    PetscPrintf(PETSC_COMM_WORLD,"@} LSQR Receiving Exit message\n");
		    break;
		  }
		}

		/*in any case clear data array*/
		for(i=0;i<eigen_max;i++)
		  data[i]=(PetscScalar)0.0+PETSC_i*(PetscScalar)0.0;
		/* if received something */
		if(!mpi_lsa_com_array_recv(com, &data_size,data) || data_load || data_load_any){
			/* we received data or load it depending on the flags (for first step only*/

			if(data_load&&data_load_any){
			  data_load_any=PETSC_FALSE;
			  data_load=PETSC_TRUE;
			}
			if(!(data_load^=data_load_any)){
/*				#ifdef DEBUG*/
/*				for(i=0;i<data_size;i++){*/
/*					if(PetscImaginaryPart(data[i])!=0.0)*/
/*					  printf("@} LSQR RECEIVED %d/%d : %e %e\n",i,data_size,PetscRealPart(data[i]),PetscImaginaryPart(data[i]));*/
/*					else*/
/*					  printf("@} LSQR RECEIVED Real %d/%d : %e\n",i,data_size,PetscRealPart(data[i]));*/
/*				}*/
/*				#endif*/

				/* first we gonna remove some non-needed values */
				epurer(data,&data_size);

/*				#ifdef DEBUG*/
/*					printf("@} LSQR Epurer : retained %d values\n",data_size);*/
/*					for(i=0;i<data_size;i++){*/
/*					   printf("@} LSQR -> retained %e %e\n",PetscRealPart(data[i]),PetscImaginaryPart(data[i]));*/
/*					}*/
/*				#endif*/

				/* add them to the accumulated eigenvalues */
				/* if full renew full eigenvalues */
				if(eigen_total+data_size>*vector_size) eigen_total=0;

				/* select eigenvalues */
				for(i=0;i<data_size;i++){
					eigen_cumul[eigen_total+i]=data[i];
				}
				eigen_total+=data_size;


				if(cumul<eigen_total) cumul=eigen_total;

				for(i=0;i<cumul;i++){
					eigen_tri[i]=eigen_cumul[i];
				}
			} else {
/*			  	PetscPrintf(PETSC_COMM_WORLD,"@} LSQR Reading file %s\n",load_path);*/
				ierr=readBinaryScalarArray(load_path,&cumul, eigen_tri);CHKERRQ(ierr);
				data_load=PETSC_FALSE;
				data_load_any=PETSC_FALSE;
				data_size=cumul;
/*				PetscPrintf(PETSC_COMM_WORLD,"@} LSQR Has Read %s\n",load_path);*/
			}



			eigen_received+=data_size;
/*			#ifdef DEBUG*/
/*				printf("@} LSQR -> COMPARISON : r : %d  m : %d ()\n",eigen_received,ls_eigen_min);*/
/*			#endif*/
			/* if we didn't received enough eigenvalues */
			if(eigen_received<ls_eigen_min) continue;
			else {
/*				#ifdef DEBUG*/
/*				PetscPrintf(PETSC_COMM_WORLD,"@} LSQR Enough good eigenvalues to compute the ellipse\n");*/
/*				#endif*/
				eigen_received=0;

				tri(eigen_tri,cumul,&chsign);

				mu1=0;
				mu2=0;

/*				#ifdef DEBUG*/
/*				printf("@} LSQR tri = %d %d\n",cumul,chsign);*/
/*				for(i=0;i<cumul;i++){*/
/*					#ifdef PETSC_USE_COMPLEX*/
/*					printf("@} LSQR Values %d/%d-%d = %e %e\n",i,cumul,chsign,PetscRealPart(eigen_tri[i]),PetscImaginaryPart(eigen_tri[i]));*/
/*					#else*/
/*					printf("@} LSQR Values %d/%d-%d = %e \n",i,cumul,chsign,PetscRealPart(eigen_tri[i]));*/
/*					#endif*/
/*				}*/
/*				#endif*/

				/* convex hull computation */
				if(chsign>0){
					convhull(eigen_tri, c, d, chsign, &mu1, 0, 0);
/*					#ifdef DEBUG*/
/*					printf("@} LSQR convhul negatif chsigne %d cumul %d mu1 %d\n",chsign,cumul,mu1);*/
/*					#endif*/
				}
				if(chsign<cumul){


					convhull(eigen_tri, c, d, cumul-chsign, &mu2, chsign, mu1);
/*					#ifdef DEBUG*/
/*					printf("@} LSQR convhul positif chsigne %d cumul %d mu1 %d mu2 %d\n",chsign,cumul,mu1,mu2);*/
/*					#endif*/
				}
				mu=mu1+mu2;
/*				#ifdef DEBUG*/
/*					printf("@} LSQR cumul-chsign=%d mu2=%d chsign=%d mu=%d\n",cumul-chsign, mu2, chsign, mu1);*/
/*					for(i=0;i<cumul;i++){*/
/*					    printf("@} LSQR Values c[%d]=%e %e\n",i,PetscRealPart(c[i]),PetscImaginaryPart(c[i]));*/
/*					}*/
/*					for(i=0;i<cumul;i++){*/
/*					    printf("@} LSQR Values d[%d]=%e %e\n",i,PetscRealPart(d[i]),PetscImaginaryPart(d[i]));*/
/*					}*/
/*				#endif*/
/*				#ifdef DEBUG*/
/*					PetscPrintf(PETSC_COMM_WORLD,"@} LSQR Params Ellipse %d %d %d\n",mu, mu1,mu2);*/
/*				#endif*/

				/* Ellipse computation */
				ierr=ellipse(c,  d, mu+1, mu1, &c_ell, &a_ell, &d_ell, &d_reel, &info);CHKERRQ(ierr);

/*				#ifdef DEBUG*/
/*					PetscPrintf(PETSC_COMM_WORLD,"@} LSQR Ellipse computed %e %e %e %e\n",c_ell, a_ell, d_ell, d_reel);*/
/*				#endif*/

				if(fabs(d_ell)<epsilon()) d_ell = 1.;

				if(fabs(a_ell)<epsilon())
					ls_eigen=0;
				else {
					LSPrecond(a_ell, d_ell,c_ell,eta, &alpha, beta,
					  delta, c, d,&mu, &ls_eigen, &ls_eigen_min, &eigen_max);
				}

			}
		}



		if(ls_eigen>1){
/*			#ifdef DEBUG*/
/*			PetscPrintf(PETSC_COMM_WORLD,"@} LSQR Sending parameters to Father %d\n",ls_eigen);*/
/*			#endif*/

			/* place the computed results inside the array */
			scalar_tmp=(PetscScalar)ls_eigen;
			for(i=0;i<ls_eigen;i++){
/*			    PetscPrintf(PETSC_COMM_WORLD,"@} %d LSQR eta[%d]=%e %e beta[%d]=%e %e delta[%d]=%e %e\n",ls_eigen,*/
/*			      i,PetscRealPart(eta[i]),PetscImaginaryPart(eta[i]),*/
/*			      i,PetscRealPart(beta[i]),PetscImaginaryPart(beta[i]),*/
/*			      i,PetscRealPart(delta[i]),PetscImaginaryPart(delta[i]));*/
			}

			ierr=PetscMemcpy(&result_array[0],&scalar_tmp,1*sizeof(PetscScalar));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[1],&alpha,1*sizeof(PetscScalar));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[2],eta,ls_eigen*sizeof(PetscScalar));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[2+ls_eigen],beta,ls_eigen*sizeof(PetscScalar));CHKERRQ(ierr);
			ierr=PetscMemcpy(&result_array[2+2*ls_eigen],delta,ls_eigen*sizeof(PetscScalar));CHKERRQ(ierr);

			if(continuous_export){
			  ierr=writeBinaryScalarArray(data_export?export_path:"./lsqr.bin", cumul, eigen_tri);
			}

			/* and send it */
			result_array_size=2+3*ls_eigen;
			mpi_lsa_com_array_send(com, &result_array_size,result_array);
		}
		if(ls_eigen>1){
		  ls_eigen=0;
		}

	}

	if(data_export){
		ierr=writeBinaryScalarArray(export_path, cumul, eigen_tri);
	}

	/* Free the arrays */
	PetscFree(eigen_tri);
	PetscFree(eigen_cumul);
	PetscFree(d);
	PetscFree(c);
	PetscFree(data);
	PetscFree(eta);
	PetscFree(beta);
	PetscFree(delta);
 	ierr=PetscFree(result_array);CHKERRQ(ierr);

	return 0;
}

