##CREATED BY PIERRE-YVES AQUILANTI 2011*/
ALL: all

CFLAGS   = -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include -DDEBUG 
SOURCEC  =  arnoldi.c
SOURCEH  =  arnoldi.h
CLEANFILE    = arnoldi.o

%.o: %.c
#	scorep ${CC} ${CFLAGS} -c $< -o ./$@
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

#cleanall: ${CLEANFILE}
#	${RM} ${GLSA_OBJ_DIR}/$<

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
