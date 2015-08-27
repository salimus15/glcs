#CREATED BY PIERRE-YVES AQUILANTI 2011*/
ALL: all

CFLAGS   = -DEBUG -I${PETSC_DIR}/include -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/${PETSC_ARCH}/lib/
SOURCEC  =  father.c
SOURCEH  =  father.h
CLEANFILE  = father.o

%.o: %.c
#	scorep ${CC} ${CFLAGS} -c $< -o ./$@
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

#cleanall: ${CLEANFILE}
#	${RM} ${GLSA_OBJ_DIR}/$<

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
