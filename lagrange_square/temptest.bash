#!/bin/bash

###############################################################################
# Current directory
CURRENTDIR=`pwd`

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)
###############################################################################
## EDIT THIS
############

#PROGRAM="sq_po_semipara"
PROGRAM="sq_lgr"

RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 0 --ang 0 --rey 100 --noel 8 --itstimedir tempdir"

###############################################################################
#cd $OOMPH_ROOT_DIR && \
#./autogen_wipebuild_noselftest.sh --rebuild --jobs=4 && \
#cd $CURRENTDIR && \
#make $PROGRAM && \
#mpirun -np 2 ./$PROGRAM $RUNARGUMENTS


#cd $OOMPH_ROOT_DIR/src && make && make install && \
#cd $CURRENTDIR && \
#make $PROGRAM && mpirun -np 1 ./$PROGRAM $RUNARGUMENTS

cd $OOMPH_ROOT_DIR/src && make && make install && \
cd $CURRENTDIR && \
make $PROGRAM && mpirun -np 1 ./$PROGRAM $RUNARGUMENTS

###############################################################################

#cd $OOMPH_ROOT_DIR && \
#./autogen_wipebuild_noselftest.sh --rebuild --jobs=4 && \
#cd $CURRENTDIR && \
#make $PROGRAM && \
#mpirun -np 1 ./$PROGRAM --ns_solver 1 --visc 0 --ang 0 --rey 0 --noel 64 --itstimedir ray_temp

# I should get:
# RAYITS: 0    19 34    26.5(2)



