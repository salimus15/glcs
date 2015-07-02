#CREATED BY PIERRE-YVES AQUILANTI 2011*/
ALL: all

CFLAGS   =  -I${SLEPC_DIR}/include/private/ -I${SLEPC_DIR}/include -DDEBUG ${PETSC_CC_INCLUDES} -I${SLEPC_DIR}/${PETSC_ARCH}/include ${SLEPC_LIB}
SOURCEC  =  arnoldi.c
SOURCEH  =  arnoldi.h
CLEANFILE    = arnoldi.o

%.o: %.c
#	scorep ${CC} ${CFLAGS} -c $< -o ./$@
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

#cleanall: ${CLEANFILE}
#	${RM} ${GLSA_OBJ_DIR}/$<

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

