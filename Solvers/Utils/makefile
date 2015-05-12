#CREATED BY PIERRE-YVES AQUILANTI 2011*/
ALL: all

CFLAGS   = ${PETSC_CC_INCLUDES}
SOURCEC  =  ellipse.c convhul.c triep.c lib.c
SOURCEH  =  ellipse.h convhul.h triep.h lib.h
CLEANFILE    = ellipse.o convhul.o triep.o lib.o

%.o: %.c
	@${CC} ${CFLAGS} -c $< -o ./$@

all: ${CLEANFILE}

#cleanall: ${CLEANFILE}
#	${RM} ${GLSA_OBJ_DIR}/$<


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
