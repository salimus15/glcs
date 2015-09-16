/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#include "mpi_lsa_com.h"

/***********************************************************************************
*	\brief Exchange information about the vector location within each intracommunicator to allow asynchronous communications 
*
*	
*
************************************************************************************/
int mpi_lsa_com_vecsize_init(com_lsa * com, Vec * v){
	PetscErrorCode ierr;
	int * size1,*size2;
	int number1,number2,i;
	MPI_Comm com1,com2;
	PetscInt global_size;
	PetscInt local_size;

	/* get vector size */
	ierr=VecGetLocalSize(*v,&local_size);CHKERRQ(ierr);
	ierr=VecGetSize(*v,&global_size);CHKERRQ(ierr);

	/* get communicator size */
	MPI_Comm_remote_size(com->out_com,&(com->out_number));
	MPI_Comm_remote_size(com->in_com,&(com->in_number));

	MPI_Barrier(MPI_COMM_WORLD);

	/* alloc arrays depending on the number of remote processors */
	com->out_size=malloc(sizeof(int)*(com->out_number));
	com->in_size=malloc(sizeof(int)*(com->in_number));

	/* because we need to select the proper communicator between the two that link
	 * each intracom we will copy it into a variable */
	if(com->color_group==GMRES_GR || com->color_group==ARNOLDI_GR){ // if many
		com1=com->in_com;
		com2=com->out_com;
		size1=com->in_size;
		size2=com->out_size;
		number1=com->in_number;
		number2=com->out_number;
	} else { // if one
		com1=com->out_com;
		com2=com->in_com;
		size1=com->out_size;
		size2=com->in_size;
		number1=com->out_number;
		number2=com->in_number;
	}
	com->vec_requests=malloc(sizeof(MPI_Request)*(com->out_number));
	com->type_requests=malloc(sizeof(MPI_Request)*(com->out_number));
	com->array_requests=malloc(sizeof(MPI_Request)*(com->out_number));

	com->vec_in_disp=malloc(sizeof(int)*(com->in_number));
	com->vec_out_disp=malloc(sizeof(int)*(com->out_number));

	com->out_sended=0;
	com->in_received=0;

	PetscMalloc(local_size*sizeof(PetscScalar),&(com->out_sended_buffer));
	PetscMalloc(local_size*sizeof(PetscScalar),&(com->in_received_buffer));

	PetscMalloc(local_size*sizeof(PetscScalar),&(com->array_out_sended_buffer));
	PetscMalloc(local_size*sizeof(PetscScalar),&(com->array_in_received_buffer));

	com->out_vec_sended=0;
	com->array_out_sended=0;


/*********************************************************************************** 
***************** 		TO SEE 		********************************************
***********************************************************************************/

	/* now process to the exchange */
	MPI_Allgather(&local_size,1,MPI_INT,size1,number1,MPI_INT,com1);
	MPI_Allgather(&local_size,1,MPI_INT,size2,number2,MPI_INT,com2);

	// displacements for each in or out node
	com->vec_in_disp[0]=0;
	for(i=1;i<com->in_number;i++){
		com->vec_in_disp[i]=com->vec_in_disp[i-1]+com->in_size[i-1];
	}

	com->vec_out_disp[0]=0;
	for(i=1;i<com->out_number;i++){
		com->vec_out_disp[i]=com->vec_out_disp[i-1]+com->out_size[i-1];

	}


	return 0;
}

/*********************************************************************************** 
*	\brief send data to one or many nodes through an intercommunicator 
*
*
*
************************************************************************************/
int mpi_lsa_com_vec_send(com_lsa * com, Vec * v){
	int i,incr,tmp_int,vsize,send_size;
	PetscScalar * array;
	PetscScalar tmp_global,tmp_local;
	PetscErrorCode ierr;
	MPI_Status status;
	int flag;


	/* check if previous requests where completed */
	for(i=0;i<com->out_vec_sended;i++){
		
		MPI_Test(&(com->vec_requests[i]),&flag,&status);
		/* if not cancel it */
		if(!flag){
			/* if no vector was sent by any of the group nodes*/
			tmp_local=(PetscScalar)com->out_vec_sended;
			tmp_global=(PetscScalar)0;

			/*now proceed to the communication*/
			// we comput the sum of tmp local of the group within each node of the group
			MPI_Allreduce(&tmp_local,&tmp_global,1,MPIU_SCALAR,MPI_SUM,com->com_group);
			// determines the size of com_group
			MPI_Comm_size(com->com_group,&vsize);

			/* we update the vector to send to the latest version */
			if(((int)(PetscRealPart(tmp_global)))==vsize){
				for(i=0;i<com->out_vec_sended;i++)
					MPI_Cancel(&(com->vec_requests[i]));
				com->out_vec_sended=0;
			}else{
				return 1;
			}
		}else{
			com->out_vec_sended--;
		}
	}

	MPI_Barrier(PETSC_COMM_WORLD);

	/* extract the vector array */
	// first we get the size of the array (stored in the local memory)
	ierr=VecGetLocalSize(*v,&vsize);CHKERRQ(ierr);
	// and now we get a pointer to that array  
	ierr=VecGetArray(*v,&array);CHKERRQ(ierr);

	// we copy the array extracted to the sended_buffer
	for(i=0;i<vsize;i++)
		(com->out_sended_buffer)[i]=(PetscScalar)array[i];
		
		
		
/*********************************************************************************** 
***************** 		TO SEE 		********************************************
***********************************************************************************/

	/* for each node in the out domain */
	for(i=0;i<com->out_number;i++){
		/* compute array displacement */
		if(i<1)
			incr=0;
		else
			incr=com->vec_out_disp[i];
			MPI_Comm_size(com->com_group,&tmp_int);

		send_size=com->out_size[i];
		if(vsize<com->out_size[i]) send_size=vsize;

		/* we send a portion of data */
		MPI_Isend(&(com->out_sended_buffer)[incr],send_size,MPIU_SCALAR,i,i,com->out_com,&(com->vec_requests[i]));
		com->out_vec_sended++;
	}


	/* restore array */
	ierr=VecRestoreArray(*v,&array);CHKERRQ(ierr);
	return 0;
}

/*********************************************************************************** 
*	\brief 
*
*
*
************************************************************************************/

// for the validaig funcion check on the bottom
int mpi_lsa_com_vec_recv(com_lsa * com, Vec * v){
	int flag,count,i;
	PetscInt size,tmp_int;
	PetscErrorCode ierr;
	MPI_Status status,statas;

	MPI_Comm_size(com->com_group,&tmp_int);

	/* first we check if there's data to receive */
	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,com->in_com,&flag,&status);
	/* did we received something ? */
	if(!flag){
		return 1; // didn't received anything
	}

	/* just verify coerancy of data */
	MPI_Get_count(&status,MPIU_SCALAR,&count);
	ierr=VecGetSize(*v,&size);CHKERRQ(ierr);

	//   if(count!=(int)size) {
	// 	printf("ERROR: count %d - size %d\n",count,size);
	// 	return -1; // ALERT !!! We got a real problem here
	// }

	/* get the array, must be done because we can't restore the array otherwise
	 * this operation is faster than setting the values of the array (i think).
	 * the array pointer is only pointing to the local vector memory of the node */
	MPI_Recv(&com->in_received_buffer[com->vec_in_disp[status.MPI_SOURCE]],count,MPIU_SCALAR,status.MPI_SOURCE,status.MPI_TAG,com->in_com,&statas);
	com->in_received++;
	
	if(!(mpi_lsa_com_vec_recv_validate(com, v, size))){
		for( i = 0; i < size; ++i){
			ierr=VecSetValue(*v, (PetscInt)i,com->in_received_buffer[i], INSERT_VALUES);CHKERRQ(ierr);
		}
		return 0;	
	}

	return 1;
}




int mpi_lsa_com_array_send(com_lsa * com, int * size, PetscScalar * data){
	MPI_Status status;
	int flag,i,vsize;
	PetscScalar tmp_global,tmp_local;

	/* check if previous requests where completed */
	for(i=0;i<com->array_out_sended;i++){
		MPI_Test(&(com->array_requests[i]),&flag,&status);
		/* if not cancel it */
		if(!flag){
			/* if no vector was sent by any of the group nodes*/
			tmp_local=(PetscScalar)com->array_out_sended;
			tmp_global=(PetscScalar)0;
			MPI_Allreduce(&tmp_local,&tmp_global,1,MPIU_SCALAR,MPI_SUM,com->com_group);

			MPI_Comm_size(com->com_group,&vsize);

			/* we update the vector to send to the latest version */
			if(((int)(PetscRealPart(tmp_global)))==vsize){
				for(i=0;i<com->array_out_sended;i++)
					MPI_Cancel(&(com->array_requests[i]));
				com->array_out_sended=0;
			}else{
				return 1;
			}
		}else{
			com->array_out_sended--;
		}
	}

	/*fill the send buffer with the data */
	for(i=0;i<*size;i++){
		com->array_out_sended_buffer[i]=(PetscScalar)data[i];
	}

	/* for each node in the out domain */
	for(i=0;i<com->out_number;i++){
		/* we send a portion of data */
		MPI_Isend(com->array_out_sended_buffer,*size,MPIU_SCALAR,i,i,com->out_com,&(com->array_requests[i]));
		com->out_vec_sended++;
	}
	printf(" \n\n there is %d message sended \n\n",*size);
	return 0;
}


int mpi_lsa_com_array_recv(com_lsa * com, int * size, PetscScalar * data){
	MPI_Status status;
	int flag;
	int my_com_size;

/*	printf(" %d : j'entre bien dans mpi_lsa_com_array_recv\n",com->rank_world);*/
	/* first we check if there's data to receive */
	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,com->in_com,&flag,&status);
	MPI_Get_count(&status,MPI_INT,size);
	/* did we received something ? */
	if(!flag || *size==1)
		return 1; // didn't received anything

	/* how large will be the array */
	MPI_Get_count(&status,MPIU_SCALAR,size);
//	if(size == 1){ size =0;	return 1;
	printf("Receive size=%d (%d scalars)\n",*size,*size/8);
	/* receive the data array */
//	PetscMalloc(sizeof(PetscScalar)*(PetscInt)size,&data);
	MPI_Recv(data,*size,MPIU_SCALAR,status.MPI_SOURCE,status.MPI_TAG,com->in_com,&status);
	return 0;
}



/* validate if a vector was entirely received */
// int mpi_lsa_com_array_recv_validate(com_lsa * com){
// 	int vector_chunks;
// 	int my_com_size;
// 	PetscScalar tmp_local,tmp_global,*array;
// 	PetscInt local_size,i;
// 	PetscErrorCode ierr;
//
// 	/* first we get the number of vector portions that the group should receive */
// 	MPI_Comm_size(com->com_group,&my_com_size);
// 	ierr=VecGetLocalSize(*v,&local_size);CHKERRQ(ierr);
//
// 	if(my_com_size==1) { /* node is alone in it's group */
// 		vector_chunks=com->in_number;
// 	} else { /* many nodes in the group */
// 		vector_chunks=my_com_size;
// 	}
//
// 	tmp_local=(PetscScalar)com->in_received;
// 	tmp_global=(PetscScalar)0;
// 	/* now computes how many data chunks we got */
// 	ierr=PetscGlobalSum(&tmp_local,&tmp_global,com->com_group);CHKERRQ(ierr);
//
// 	if(((int)(PetscRealPart(tmp_global)))==vector_chunks){ /* we got the chunks on all processors or from all processors */
// 		ierr=VecGetArray(*v,&array);CHKERRQ(ierr);
// 		for(i=0;i<local_size;i++)
// 			array[i]=com->in_received_buffer[i];
// 		ierr=VecRestoreArray(*v,&array);CHKERRQ(ierr);
// 		com->in_received=0;
// 		return 0;
// 	}
//
// 	return 1;
// }

int mpi_lsa_com_free(com_lsa * com){
	/* free arrays */
	free(com->out_size);
	free(com->in_size);
	free(com->vec_requests);
	free(com->type_requests);
	free(com->array_requests);
	free(com->vec_in_disp);
	free(com->vec_out_disp);
	/* and mpi communicators */
	MPI_Comm_free(&(com->out_com));
	MPI_Comm_free(&(com->in_com));
	MPI_Comm_free(&(com->com_group));

	return 0;
}





int mpi_lsa_com_type_recv(com_lsa * com, int * type){
	MPI_Status status;
	int flag,size;

	/* first we check if there's data to receive */
	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,com->in_com,&flag,&status);
	/* did we received something ? */
	MPI_Get_count(&status,MPI_INT,&size);
	if(!flag || size != 1)
		return 1; // didn't received anything
	printf("size type = %d\n",size);
	MPI_Recv(type,1,MPI_INT,status.MPI_SOURCE,status.MPI_TAG,com->in_com,&status);

	return 0;
}

int mpi_lsa_com_type_send(com_lsa * com, int * type){
	MPI_Status status;
	int flag,i;

    //FIXME : Ceci part d'une bonne idée mais fait plus de mal que de bien Donc à revoir
     
	/* check if previous requests where completed */
/*	for(i=0;i<com->out_sended;i++){*/
/*		MPI_Test(&(com->type_requests[i]),&flag,&status);*/
/*		/* if not cancel it */
/*		if(!flag)*/
/*			MPI_Cancel(&(com->type_requests[i]));*/
/*	}*/

	/* for each node in the out domain */
	for(i=0;i<com->out_number;i++){
		MPI_Isend(type,1,MPI_INT,i,i,com->out_com,&(com->type_requests[i]));
		com->out_sended++;
	}

	return 0;
}


int mpi_lsa_com_vec_recv_validate(com_lsa * com, Vec * v, int size){
	int vector_chunks;
	int my_com_size;
	PetscScalar  *array;
	int global_size,local_size,i;
	PetscErrorCode ierr;

	/* first we get the number of vector portions that the group should receive */
	MPI_Comm_size(com->com_group,&my_com_size);
	
	if(my_com_size == 1){ //this means that the process is alone in his group  
	   ierr=VecGetSize(*v,&local_size);CHKERRQ(ierr);
	   if(size == local_size){
	   	ierr=VecGetArray(*v,&array);CHKERRQ(ierr);
		for(i=0;i<local_size;i++)
			array[i]=com->in_received_buffer[i];
		ierr=VecRestoreArray(*v,&array);CHKERRQ(ierr);
		com->in_received=0;
		return 0;
	   } else return 1;
	}
	
	if(my_com_size > 1){
	   global_size=0;
	   
	   ierr=VecGetLocalSize(*v,&local_size);CHKERRQ(ierr);
	   MPI_Allreduce( &local_size, &global_size, 1, MPI_INT, MPI_SUM, com->com_group);
	   if(global_size == size){
	   	
	   	   ierr=VecGetSize(*v,&local_size);CHKERRQ(ierr);
		   if(size == local_size){
		        ierr=VecGetArray(*v,&array);CHKERRQ(ierr);
			   for(i=0;i<local_size;i++)
					array[i]=com->in_received_buffer[i];
			   ierr=VecRestoreArray(*v,&array);CHKERRQ(ierr);
			   com->in_received=0;
		   return 0;
		   }
		}else return 1; 	
	}
	
	
	return 1; //this happends when my_com_size has strenge value "less than 0" 
}

/*/* validate if a vector was entirely received */
/*int mpi_lsa_com_vec_recv_validate(com_lsa * com, Vec * v){*/
/*	int vector_chunks;*/
/*	int my_com_size;*/
/*	PetscScalar tmp_local,tmp_global,*array;*/
/*	PetscInt local_size,i;*/
/*	PetscErrorCode ierr;*/

/*	/* first we get the number of vector portions that the group should receive */
/*	MPI_Comm_size(com->com_group,&my_com_size);*/
/*	ierr=VecGetLocalSize(*v,&local_size);CHKERRQ(ierr);*/

/*	if(my_com_size==1) { /* node is alone in it's group */
/*		vector_chunks=com->in_number;*/
/*	} else { /* many nodes in the group */
/*		vector_chunks=my_com_size;*/
/*	}*/

/*	tmp_local=(PetscScalar)com->in_received;*/
/*	tmp_global=(PetscScalar)0;*/
/*	/* now computes how many data chunks we got */
/*	MPI_Allreduce(&tmp_local,&tmp_global,1,MPIU_SCALAR,MPI_SUM,com->com_group);*/

/*	if(((int)(PetscRealPart(tmp_global)))==vector_chunks){ /* we got the chunks on all processors or from all processors */
/*		ierr=VecGetArray(*v,&array);CHKERRQ(ierr);*/
/*		for(i=0;i<local_size;i++)*/
/*			array[i]=com->in_received_buffer[i];*/
/*		ierr=VecRestoreArray(*v,&array);CHKERRQ(ierr);*/
/*		com->in_received=0;*/
/*		return 0;*/
/*	}*/

/*	return 1;*/
/*}*/
