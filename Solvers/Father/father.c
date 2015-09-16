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
	
	ierr=VecDuplicate(*v,&vec_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(*v,&vec_tmp_receive);CHKERRQ(ierr);
	
	while(!end){
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
		  if(exit_type==666){
		    end=1;
		    PetscPrintf(PETSC_COMM_WORLD,"$} Father Sending Exit message\n");

		 		    
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
			printf(" =========FATHER   I RECEIVED SOME DATA SO I AM GOING TO CHECK THEM ============\n");
			if(!mpi_lsa_com_vec_recv_validate(com, &vec_tmp_receive)){
				printf(" =========   I RECEIVED SOME DATA SO I AM GOING TO SEND THEM TO ARNOLDI ============\n");
				vec_tmp=vec_tmp_receive;
				mpi_lsa_com_vec_send(com,&vec_tmp);
			}
		}		
	}	
	return 0;
}
