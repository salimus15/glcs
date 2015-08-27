##CREATED BY PIERRE-YVES AQUILANTI 2011*/
ALL: all

CFLAGS   = -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -DDEBUG 
SOURCEC  =  gmres_cycle.c gmres_solve.c kspsolve.c gmres.c gmres_precond.c
SOURCEH  =  gmres_cycle.h gmres_solve.h kspsolve.h gmres.h gmres_precond.h
CLEANFILE    = gmres_cycle.o gmres_solve.o kspsolve.o gmres.o gmres_precond.o

%.o: %.c
#	scorep ${CC} ${CFLAGS} -c $< -o ./$@
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

#cleanall: ${CLEANFILE}
#	${RM} ${GLSA_OBJ_DIR}/$<

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
