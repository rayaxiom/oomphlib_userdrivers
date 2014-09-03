#!/bin/bash

# Create the THISFILEBASE, this is the folder where the testing will be done.
THISFILE=$0 # This contains "./", which we do not want.
THISFILE=${THISFILE:2} # Gets rid of "./"
THISFILEBASE=${THISFILE%%.*} # Get rid of the extension (in this case, ".sh")
OUTFILE="${THISFILEBASE}.out"


SOLN_DIR="RESLT_SOLN_tmp"
TIME_DIR="RESLT_TIME_tmp"

rm -rf $SOLN_DIR $TIME_DIR
mkdir $SOLN_DIR
mkdir $TIME_DIR


DIST_PROB="--dist_prob"
PROB_ID="--prob_id 20"
DOC_SOLN="--doc_soln $SOLN_DIR"
#DOC_SOLN=""
VISC="--visc 1"
REY="--rey 100"
MAX_ITER="--max_solver_iter 100" # TO CHECK
ITSTIMEDIR="--itstimedir $TIME_DIR" # TO CHECK
#ITSTIMEDIR="" # TO CHECK
SOLVER_TYPE="--solver_type 2" # TO CHECK
DT="--dt -1.0"
TIME_START="--time_start 0.0"
TIME_END="--time_end 0.2"
NSPP="$DIST_PROB $PROB_ID $DOC_SOLN $VISC $REY $MAX_ITER $ITSTIMEDIR $SOLVER_TYPE $DT $TIME_START $TIME_END"

###########################################
PRINT_HYPRE="--print_hypre"

W_SOLVER="--w_solver 0"
NS_SOLVER="--ns_solver 1"
P_SOLVER="--p_solver 1"
#P_SOLVER="--p_solver 0"

FAMG_ITER="--f_amg_iter 1"
FAMG_SMITER="--f_amg_smiter 2"
FAMG_SSMOOTHER="--f_amg_sim_smoo 1"
#FAMG_CSMOOTHER="--f_amg_com_smoo 9" # Euclid
FAMG_DAMP="--f_amg_damp -1"
FAMG_STRN="--f_amg_str 0.668" # REMEMBER TO CHANGE THIS FOR SIMPLE/STRESS VISCOUS FORMS.
FAMG_COARSE="--f_amg_coarse 1" #RS - 1, Falgout - 6

F_SOLVER="--f_solver 96 $FAMG_ITER $FAMG_SMITER $FAMG_SSMOOTHER $FAMG_CSMOOTHER $FAMG_DAMP $FAMG_STRN $FAMG_COARSE"
#F_SOLVER="--f_solver 0"


LPH="$PRINT_HYPRE $W_SOLVER $NS_SOLVER $P_SOLVER $F_SOLVER"

#LPH="--w_solver 0 --ns_solver 0"
##########################################

NOEL="--noel 4"
CL="$NOEL"



###########################################
PARAM="$NSPP $LPH $CL"

echo $PARAM
PROGRAM="annular_wedge_threed_bk"

#make $PROGRAM && \
#mpirun -np 1 ./$PROGRAM $PARAM #&& \

make $PROGRAM && rm -rf ./$SOLN_DIR/*.dat && rm -rf ./$TIME_DIR/* && \
mpirun -np 1 ./$PROGRAM $PARAM #> $OUTFILE 2>&1 #&& \

#make $PROGRAM && rm -rf ./$SOLN_DIR/*.dat && rm -rf ./$TIME_DIR/* && \
#mpirun -np 1 ./$PROGRAM $PARAM #&& \

#diff RESLT/soln8.dat RESLT_A30_R200_N4/soln8.dat




