/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "father.h"

/* Compute cyclicly eigenvalues, we create and destroy the EPS context inside the loop. 
   This is not the best thing to do because it may decrease performances due to setup overhead */
PetscErrorCode Father(com_lsa * com, Vec * v){
	Vec vec_tmp_receive, vec_tmp;
	PetscInt end;
	PetscErrorCode ierr;	
	int exit_type=0;	
	end=0;
 PetscPrintf(PETSC_COMM_WORLD,"$}### FATHER mr rank : %d my group : %d my color : %d  I send to %d  and receive fom %d\n",com->rank_world, com->com_group,com->color_group,com->in_com,com->out_com);
	ierr=VecDuplicate(*v,&vec_tmp);CHKERRQ(ierr);
	ierr=VecDuplicate(*v,&vec_tmp_receive);CHKERRQ(ierr);
	
	while(!end){
		if(!mpi_lsa_com_type_recv(com,&exit_type)){
		  if(exit_type==666){
		    end=1;
		    PetscPrintf(PETSC_COMM_WORLD,"$} Father Sending Exit message\n");

		    mpi_lsa_com_type_send(com,&exit_type);
		    break;
		  }
		}
		/* check if there's an incoming message */
		if(!mpi_lsa_com_vec_recv(com, &vec_tmp_receive)){
			 PetscPrintf(PETSC_COMM_WORLD,"$} ######### SOMETHING IS IN WE GONNA CHECK IT ##########\n");
			if(!mpi_lsa_com_vec_recv_validate(com, &vec_tmp_receive)){
				vec_tmp=vec_tmp_receive;
				mpi_lsa_com_vec_send(com,&vec_tmp);
			}
		}		
	}
	
	ierr = VecDestroy(&vec_tmp_receive); CHKERRQ(ierr);
	ierr = VecDestroy(&vec_tmp); CHKERRQ(ierr);
	
	
	
	return 0;
}
