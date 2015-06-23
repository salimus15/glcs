#CREATED BY PIERRE-YVES AQUILANTI 2011*/
ALL:  blib exec

#compilation flags
DIRS    = Libs Tuning Solvers
EXEC    = hyperh
CFLAGS	= -O0
#-I${PETSC_DIR}/src/ksp/ksp/impls/gmres/fgmres/

CLEANFILES  = ${EXEC}
HFILES = -I./Libs  -I./Solvers/Utils -I./Solvers/Arnoldi -I./Tuning 
OFILES= ${wildcard ./Libs/*.o} ${wildcard ./Solvers/*/*.o} ${wildcard ./Tuning/*.o}

#matrices directory
MDIR=./data
#MDIR=/home/pyaquilanti/Data/Matrices

##################################################################
##################      Tuning  FLAGS      #######################
##################################################################

RESTART_MAX=  100
RESTART_MIN=  3
RESTART_INCREMENT = 3
RESTART_STRATEGY = none
RESTART_FULL_TYPE = none
ORTHOG_TYPE = 1
ORTHOG_SIZE = ${RESTART_MAX}
ORTHOG_HEURISTIC = 1
RESTART =  -ksp_gmres_restart_min ${RESTART_MIN} -restart_type ${RESTART_STRATEGY} -restart_type_full ${RESTART_FULL_TYPE}
ORTHOG = -iorthog_type ${ORTHOG_TYPE} -iorthog_size ${ORTHOG_SIZE} -iorthog_heurist ${ORTHOG_HEURISTIC}

##################################################################
##################      Solvers FLAGS      #######################
##################################################################

#debug options
# DEBUG_VALGRIND = valgrind --tool=memcheck -q
DEBUG_KSP_VIEW = -ksp_view
#gmres options
GMRES_PRECISION = 1e-20
GMRES_RESTART = ${RESTART_MAX}
GMRES_NB_NODES = 1
GMRES_MONITOR = -ksp_monitor_true_residual
GMRES_FLAGS = -ksp_rtol 1e-100 -ksp_divtol 1e1000 -ksp_max_it 10000 -pc_type none -ksp_atol ${GMRES_PRECISION} -ksp_gmres_restart ${GMRES_RESTART}\
		${GMRES_MONITOR} ${GMRES_VIEW} -lsa_gmres ${GMRES_NB_NODES} ${RESTART} ${ORTHOG}
#arnoldi options
ARNOLDI_PRECISION = 1e-5
ARNOLDI_NBEIGEN = 100
ARNOLDI_NB_NODES = 1
ARNOLDI_MONITOR = -eps_monitor
# ARNOLDI_LOAD_ANY = -ksp_arnoldi_load_any
ARNOLDI_FLAGS = -eps_type arnoldi -eps_true_residual -eps_largest_imaginary -eps_nev ${ARNOLDI_NBEIGEN} -eps_tol ${ARNOLDI_PRECISION} \
		${ARNOLDI_MONITOR} -lsa_arnoldi ${ARNOLDI_NB_NODES} -eps_max_it 5 -ksp_arnoldi_cexport ${ARNOLDI_LOAD_ANY}
#ls options
LS_POWER = 15
LS_POLY_APPL = 11
LS_LATENCY = 20
LS_PC_USE = 1
LS_HANG_IT = 2000
LS_HANG_TIME =  1
# LS_LOAD_ANY = -ksp_ls_load_any
LS_FLAGS = -ksp_ls_power ${LS_POWER} -ksp_ls_m_hang ${LS_HANG_IT} -ksp_ls_timing ${LS_HANG_TIME}  -ksp_ls_k_param ${LS_POLY_APPL} -ksp_ls_nopc ${LS_PC_USE} -ksp_ls_latency ${LS_LATENCY} -ksp_ls_cexport ${LS_LOAD_ANY}
#final flag composition
GLSA_FLAGS = ${DEBUGG} ${GMRES_FLAGS} ${ARNOLDI_FLAGS} ${LS_FLAGS} ${DEBUG_KSP_VIEW}
MPI_NODES = ${shell echo ${GMRES_NB_NODES}+${ARNOLDI_NB_NODES}+2 | bc}

##################################################################
##################   Compilation rules     #######################
##################################################################


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${SLEPC_DIR}/conf/slepc_common

# blib :
# 	-@echo "BEGINNING TO COMPILE LIBRARIES "
# 	-@echo "========================================="
# 	${OMAKE}  PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} ACTION=lib
# 	SHELL=/bin/bash

blib:
	@for a in $(DIRS); do \
	if [ -d $$a ]; then \
#	echo "processing folder $$a"; \
	$(MAKE) -C $$a; \
	fi; \
#	echo "Creating library from $$GLSA_OBJ_DIR"; \
#	${AR} rcs libglsa.a ${GLSA_OBJ_DIR}/*.o; \
	done;
	@echo "$$a Done!"

	# -@echo "Completed building libraries"
	# -@echo "========================================="

distclean :
	-@echo "Cleaning application and libraries"
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH}  PETSC_DIR=${PETSC_DIR} clean
	-${RM} ${OFILES} main.o
	-@echo "Finised cleaning application and libraries"
	-@echo "========================================="

rmat :
	-rm -drf matrices_tmp/*

exec: main.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	-@echo ${CLINKER}	
	-@echo " -----------------------------------------------------------"
	@${CLINKER} -g -v -o ${EXEC} main.o ${OFILES} ${HFILES} -I${PETSC_DIR}/include -L${SLEPC_LIB} -L${PETSC_DIR}/${PETSC_ARCH}/lib  -L.

effacer :
	-rm *.bin
	-rm *.o 
	-rm ./*/*.o
	-rm ./*/*/*.o


	-@echo "Completed building application"
	-@echo "========================================="

##################################################################
##################     Execution Rules     #######################
##################################################################
#valgrind --sigill-diagnostics=yes --show-below-main=yes --leak-check=full --show-leak-kinds=all
runl:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ${MDIR}/young4c.mtx_841x841_4089nnz \
	2>&1 | tee log.txt
	
runs:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ${MDIR}/mhd1280a.mtx_1280x1280_47906nnz \
	2>&1 | tee log.txt

runx:
	 ${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND}  ./hyperh ${GLSA_FLAGS} \
	-mfile ${MDIR}/waveguide3D.mtx_21036x21036_303468nnz  \
	-vfile ${MDIR}/waveguide3D_b.mtx_21036 \
	2>&1 | tee log.txt


runa:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh ${GLSA_FLAGS} \
	-mfile ${MDIR}/utm1700b/data/utm1700b.mtx_1700x1700_21509nnz.gz  \
	-vfile ${MDIR}/utm1700b/data/utm1700b_b.mtx_1700.gz \
	2>&1 | tee log.txt


runb:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ${MDIR}/pde225/data/pde225.mtx_225x225_1065nnz.gz \
	2>&1 | tee log.txt

runc:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ${MDIR}/utm300/data/utm300.mtx_300x300_3155nnz.gz \
	-vfile ${MDIR}/utm300/data/utm300_b.mtx_300.gz \
	2>&1 | tee log.txt

rund:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ${MDIR}/data/fs_183_1.mtx_183x183_1069nnz.gz  \
	2>&1 | tee log.txt

rune:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ${MDIR}/pde2961/data/pde2961.mtx_2961x2961_14585nnz.gz \
	2>&1 | tee log.txt

runf:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh ${GLSA_FLAGS} \
	-mfile ${MDIR}/epb0/data/epb0.mtx_1794x1794_7764nnz.gz \
	2>&1 | tee log.txt

rung:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh ${GLSA_FLAGS} \
	-mfile ${MDIR}/fs_541_1/data/fs_541_1.mtx_541x541_4285nnz.gz \
	2>&1 | tee log.txt


runh:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ${MDIR}/shermanACa/data/shermanACa.mtx_3432x3432_25220nnz.gz \
	2>&1 | tee log.txt


runi:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh ${GLSA_FLAGS} \
	-mfile ${MDIR}/utm3600/data/utm3060.mtx_3060x3060_42211nnz.gz \
	-vfile ${MDIR}/utm3600/data/utm3060_b.mtx_3060.gz \
	2>&1 | tee log.txt


runj:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile /home/perif/Data/datasets/utils/scripts/utm1700a.mtx_1700x1700_21313nnz.gz \
	-vfile /home/perif/Data/datasets/utils/scripts/utm1700a_b.mtx_1700.gz \
	2>&1 | tee log.txt



utm3:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
         -mfile ~/AutoSolverClone/Matrices/utm300.mtx_300x300_3155nnz.gz -vfile ~/AutoSolverClone/Matrices/utm300_b.mtx_300.gz | tee log.txt

utm7:
	 -@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
         -mfile ~/Data/utils/scripts/utm1700a.mtx_1700x1700_21313nnz.gz -vfile ~/Data/utils/scripts/utm1700a_b.mtx_1700.gz | tee log.txt


utm5:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
        -mfile ~/Data/Matrices/utm5940/data/utm5940.mtx_5940x5940_83842nnz.gz -vfile ~/Data/Matrices/utm5940/data/utm5940_b.mtx_5940.gz  | tee log.txt


dwa:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ~/Data/Matrices/dwa512/data/dwa512.mtx_512x512_2480nnz.gz | tee log.txt

dwb:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ~/Data/Matrices/dw256A/data/dw256A.mtx_512x512_2480nnz.gz | tee log.txt


pde:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ~/Data/Matrices/pde2961/data/pde2961.mtx_2961x2961_14585nnz.gz | tee log.txt

e2:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ~/Data/Matrices/epb2_a_confirmer/data/epb2.mtx_25228x25228_175027nnz.gz | tee log.txt

e1:
	-@${MPIEXEC} -np ${MPI_NODES} ${DEBUG_VALGRIND} ./hyperh  ${GLSA_FLAGS} \
	-mfile ~/Data/Matrices/epb1_a_confirmer/data/epb1.mtx_14734x14734_95053nnz.gz | tee log.txt


