#!/bin/bash

PROGRAM="unstructured_three_d_fluid"

#PARAM1="--doc_soln RESLT --mesh_type 1 --mesh_area 0.4"

#PARAM2="--visc 0 --rey 100"
#./unstructured_three_d_fluid $PARAM1 $PARAM2

CURRENT_DIR=`pwd`

DIST_PROB="--dist_prob"
PROB_ID="--prob_id 0"
DOC_SOLN="--doc_soln tmp_soln"
VISC="--visc 0"
REY="--rey 100"
MAX_SOLVER_ITER="--max_solver_iter 1000"
ITSTIMEDIR="--itstimedir tmp_res_its"
SOLVER_TYPE="--solver_type 1"
DT="--dt 0.1"
#DT=""
TIME_START="--time_start 0"
TIME_END="--time_end 1"
MESH_TYPE="--mesh_type 1"

NSPP="$DIST_PROB $PROB_ID $DOC_SOLN $VISC $REY $MAX_SOLVER_ITER $ITSTIMEDIR $SOLVER_TYPE $DT $TIME_START $TIME_END $MESH_TYPE"

#######
PRINT_HYPRE="--print_hypre"
W_SOLVER="--w_solver 0"
NS_SOLVER="--ns_solver 1"
P_SOLVER="--p_solver 1"
#P_SOLVER="--p_solver 0"
F_SOLVER="--f_solver 96 --f_amg_iter 2 --f_amg_smiter 2 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.668 --f_amg_coarse 1"
#F_SOLVER="--f_solver 0"

LPH="$PRINT_HYPRE $W_SOLVER $NS_SOLVER $P_SOLVER $F_SOLVER"

OUTFILEBASE="sim_amg_0d"


#########
#AREA="8"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1
#
#########
#AREA="7"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1
#
#########
#AREA="6"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1
#
#########
#AREA="5"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1
#
#########
#AREA="4"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1
#
#########
#AREA="3"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1

#########
#AREA="2"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1
#
#########
#AREA="1"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1
#
#########
#AREA="09"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1
#
#########
#AREA="08"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1
#
#########
#AREA="07"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1
#
#########
#AREA="06"
#MESH_AREA="--mesh_area 0.$AREA"
#OUTFILE="$OUTFILEBASE$AREA"
#BL="$MESH_AREA"
#PARAM="$NSPP $LPH $BL"
#./$PROGRAM $PARAM > "$OUTFILE" 2>&1

########
AREA="05"
MESH_AREA="--mesh_area 0.$AREA"
OUTFILE="$OUTFILEBASE$AREA"
BL="$MESH_AREA"
PARAM="$NSPP $LPH $BL"
./$PROGRAM $PARAM > "$OUTFILE" 2>&1

########
AREA="04"
MESH_AREA="--mesh_area 0.$AREA"
OUTFILE="$OUTFILEBASE$AREA"
BL="$MESH_AREA"
PARAM="$NSPP $LPH $BL"
./$PROGRAM $PARAM > "$OUTFILE" 2>&1

########
AREA="03"
MESH_AREA="--mesh_area 0.$AREA"
OUTFILE="$OUTFILEBASE$AREA"
BL="$MESH_AREA"
PARAM="$NSPP $LPH $BL"
./$PROGRAM $PARAM > "$OUTFILE" 2>&1

########
AREA="02"
MESH_AREA="--mesh_area 0.$AREA"
OUTFILE="$OUTFILEBASE$AREA"
BL="$MESH_AREA"
PARAM="$NSPP $LPH $BL"
./$PROGRAM $PARAM > "$OUTFILE" 2>&1







