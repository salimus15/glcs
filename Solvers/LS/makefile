#CREATED BY PIERRE-YVES AQUILANTI 2011*/
ALL: all

CFLAGS   =  -DDEBUG -I${PETSC_DIR}/include ${PETSC_LIB} -I${PETSC_CC_INCLUDES}
SOURCEC  =  lsqr.c precond.c
SOURCEH  =  lsqr.h precond.h
CLEANFILE    = lsqr.o precond.o

%.o: %.c
#	scorep ${CC} ${CFLAGS} -c $< -o ./$@
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

#cleanall: ${CLEANFILE}
#	${RM} ${GLSA_OBJ_DIR}/$<


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
