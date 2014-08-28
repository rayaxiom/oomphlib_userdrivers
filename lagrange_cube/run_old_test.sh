#!/bin/bash

#--dist_prob --prob_id 21 --max_solver_iter 100 --solver_type 2 --dt 0.1 --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 96 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.668 --f_amg_coarse 1 --visc 0 --ang 30 --rey 200 --noel 4 --time_start 0.0 --time_end 1.0
PROGRAM="cube"


TIMETYPE="--time_type 2"
SOLVERTYPE="--solver_type 2"
DISTPROB="--dist_prob"
MAXSOLVERITER="--max_solver_iter 100"
DT="--dt 0.1"
TIMESTART="--time_start 0.0"
TIMEEND="--time_end 1.0"
DOCSOLN=""
ITSTIMEDIR=""

GENHELPER="$TIMETYPE $SOLVERTYPE $DISTPROB $MAXSOLVERITER "
GENHELPER+="$DT $TIMESTART $TIMEEND "
GENHELPER+="$DOCSOLN $ITSTIMEDIR"

NSHELPER=""

PRECHELPER=""

PROBHELPER="--prob_id 0 --ang 30 --noel 4 "



PARAM="$GENHELPER $NSHELPER $PRECHELPER $PROBHELPER"


make $PROGRAM && \
mpirun -np 1 ./$PROGRAM $PARAM | tee old_data && \
grep RAYITS old_data > old_data_its && \
diff old_data_its old_data_its_static
