#!/bin/bash
 
set -u


RES_DIR="test2_ang67_1v2cRSGs_oomph"
mkdir -p $RES_DIR
PROGRAM="lgr_cube"

#################
TIMETYPE="--time_type 1"
SOLVERTYPE="--solver_type 1"
DISTPROB="--dist_prob"
MAXSOLVERITER="--max_solver_iter 1000"
DT=""
TIMESTART="--time_start 0.0"
TIMEEND="--time_end 1.0"
DOCSOLN=""
ITSTIMEDIR="--itstimedir $RES_DIR"

GENHELPER="$TIMETYPE $SOLVERTYPE $DISTPROB $MAXSOLVERITER "
GENHELPER+="$DT $TIMESTART $TIMEEND "
GENHELPER+="$DOCSOLN $ITSTIMEDIR"
#################

VISC="" ## CHANGE THIS --visc 0, 1 VARYING THIS #############################
REY="--rey 200"
NSHELPER="$VISC $REY"

#################
WSOLVER="--w_solver 0"
NSSOLVER="--ns_solver 1"
PSOLVER="--p_solver 0"
FSOLVER="--f_solver 1"
F_ITER="--f_amg_iter 1"
F_SMITER="--f_amg_smiter 2"
F_SMOOTHER="--f_amg_sim_smoo 1" ## CHANGE THIS
F_DAMP="--f_amg_damp -1"
F_COARSE="--f_amg_coarse 1"
F_PRINT="--print_f_hypre"

PRECHELPER="$WSOLVER $NSSOLVER $FSOLVER $PSOLVER "
PRECHELPER+="$F_ITER $F_SMITER $F_SMOOTHER $F_DAMP $F_COARSE "
PRECHELPER+="$F_PRINT"

#################

PROBID="--prob_id 0"
ANG="--ang 67"
NOEL=""  #NOEL="--noel 4" ## CHANGE THIS VARYING THIS #######################

PROBHELPER="$PROBID $ANG $NOEL "

#################

PARAM="$GENHELPER $NSHELPER $PRECHELPER $PROBHELPER"

VISC="--visc 0"

F_STRN="--f_amg_str 0.25"
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 4
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 6
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 8
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 10
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 12
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 14
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 16
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 18

VISC="--visc 1"
F_STRN="--f_amg_str 0.75"
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 4
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 6
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 8
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 10
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 12
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 14
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 16
mpirun -np 1 ./$PROGRAM $PARAM $F_STRN $VISC --noel 18





