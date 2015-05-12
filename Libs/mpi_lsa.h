/*CREATED BY PIERRE-YVES AQUILANTI 2011*/
#ifndef MPI_LSA_H_
#define MPI_LSA_H_


#include "mpi.h"
#include "args_handler.h"
#include "petsc.h"

/* some constants */
#define GMRES_GR 0
#define FATHER_GR 1
#define ARNOLDI_GR 2
#define LS_GR 3


typedef struct _com_lsa{
	int com_world;

	union {
		int com[4];
		int gmres,father,arnoldi,ls;
	}group;

	union {
		int com[4];
		int gmres,father,arnoldi,ls;
	} inter;

	union {
		int com[4];
		int gmres,father,arnoldi,ls;
	} size;

	union {
		int com[4];
		int gmres,father,arnoldi,ls;
	} master;

	int size_world;
	int com_group;
	int rank_group;
	int color_group;

	int rank_world;

	int * in_size;
	int * out_size;

	int * vec_requests;
	int * type_requests;
	int * array_requests;

	int vec_requests_nb;
	int type_requests_nb;
	int array_requests_nb;

	int * vec_in_disp;
	int * vec_out_disp;

	int in_number;
	int out_number;

	int out_sended;
	int in_received;

	PetscScalar * out_sended_buffer;
	PetscScalar * in_received_buffer;


	PetscScalar * array_out_sended_buffer;
	PetscScalar * array_in_received_buffer;

	int out_vec_sended;
	int array_out_sended;

	int in_com;
	int out_com;

} com_lsa;

int mpi_lsa_init(int argc, char ** argv, com_lsa * com);

int mpi_lsa_create_groups(com_lsa * com);

int mpi_lsa_create_intercoms(com_lsa * com);

int mpi_lsa_print(char * s,com_lsa * com);

#endif

