#CREATED BY PIERRE-YVES AQUILANTI 2011*/
ALL: all

CFLAGS   = -DEBUG -I${SLEPC_DIR}/include -I ${SLEPC_DIR}/${PETSC_ARCH}/include ${SLEPC_LIB} ${PETSC_LIB} -I${PETSC_DIR}/include 
SOURCEC  =  father.c
SOURCEH  =  father.h
CLEANFILE  = father.o

%.o: %.c
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

#cleanall: ${CLEANFILE}
#	${RM} ${GLSA_OBJ_DIR}/$<

include ${SLEPC_DIR}/conf/slepc_common
