/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "father.h"

/* Compute cyclicly eigenvalues, we create and destroy the EPS context inside the loop. 
   This is not the best thing to do because it may decrease performances due to setup overhead */
PetscErrorCode Father(com_lsa * com, Vec * v){
	Vec vec_tmp_receive, vec_tmp;
	PetscInt end;
	PetscErrorCode ierr;	
	int exit_type=0, i ;	
	end=0;
	MPI_Status status;
	int flag;
	Vec vecteur_initial;
        
	int taille = 1;
	
		
	ierr=VecDuplicate(*v,&vec_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(*v,&vec_tmp_receive);CHKERRQ(ierr);
	//vecteur_initial = malloc(taille * sizeof(PetscScalar));
		//setting_out_vec_sizes( com, v);
	
	while(!end){
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
		  if(exit_type==666){
		    end=1;
		    PetscPrintf(PETSC_COMM_WORLD,"$} Father Sending Exit message\n");

		 	// On purge le buffer d'abord avant d'envoyer le message de fin	    
		    for(i=0;i<com->out_vec_sended;i++){
		
				MPI_Test(&(com->vec_requests[i]),&flag,&status);
				/* if not cancel it */
				if(!flag){					
					/* we update the vector to send to the latest version */
				 	MPI_Cancel(&(com->vec_requests[i]));				
				}
		    }
		    	mpi_lsa_com_type_send(com,&exit_type);
		  		break;
		  }
	   }
	
		/* check if there's an incoming message */
		if(!mpi_lsa_com_vec_recv(com, &vec_tmp_receive)){
			
/*			if(!mpi_lsa_com_vec_recv_validate(com, &vec_tmp_receive)){*/
		
		
/*		ierr=VecGetSize(vec_tmp_receive,&taille);CHKERRQ(ierr);*/
/*		PetscPrintf(PETSC_COMM_WORLD,"==== > Father OUR INITIALV IS OF SIZE %d\n", taille);*/
/*  		vecteur_initial = realloc(vecteur_initial,taille);  			*/
/*  		ierr=VecGetArray(vec_tmp_receive, &vecteur_initial);CHKERRQ(ierr);*/
/*		for (i = 0; i < taille; i++)*/
/*			PetscPrintf(PETSC_COMM_WORLD,"==== > father[%d] = %e\n", i, vecteur_initial[i]);*/
	
/*		mpi_lsa_com_array_send(com, &taille, vecteur_initial);*/
/*		ierr= VecRestoreArray(vec_tmp_receive, &vecteur_initial);CHKERRQ(ierr);*/

					
// 				vec_tmp=vec_tmp_receive;
			if(com->in_received == com->in_number){
				ierr=VecCopy(vec_tmp_receive,vec_tmp);CHKERRQ(ierr);
				mpi_lsa_com_vec_send(com,&vec_tmp);
				com->in_received = 0;
			}
			
			
		}

/*		if(!mpi_lsa_receive_vec(com, &vecteur_initial)){*/
/*			mpi_lsa_send_vec(com, &vecteur_initial);*/
/*		}		*/
	}	
	return 0;
}
