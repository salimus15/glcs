#CREATED BY PIERRE-YVES AQUILANTI 2011*/
ALL: all

CFLAGS   =  -DEBUG  ${PETSC_CC_INCLUDES}
SOURCEC  =  mpi_lsa.c args_handler.c mpi_lsa_com.c read_matrix.c data_rw.c real2complex.c
SOURCEH  =  mpi_lsa.h args_handler.h mpi_lsa_com.h read_matrix.h data_rw.h real2complex.h
CLEANFILE = mpi_lsa.o args_handler.o mpi_lsa_com.o read_matrix.o data_rw.o real2complex.o

%.o: %.c
#	scorep ${CC} ${CFLAGS} -c $< -o ./$@
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

#cleanall: ${CLEANFILE}
#	${RM} ${GLSA_OBJ_DIR}/$<

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
