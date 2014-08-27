#!/bin/bash

RES_DIR="2v2Strn0-75cRSEuclid_oomph"
mkdir -p $RES_DIR

PROGRAM="lgr_cube"
PARAM="--time_type 1 --solver_type 1 --dist_prob --max_solver_iter 100 --time_start 0.0 --time_end 1.0 --itstimedir $RES_DIR --rey 200 --prob_id 0 --ang 30 --w_solver 0 --ns_solver 1 --f_solver 1 --p_solver 0 --f_amg_iter 2 --f_amg_smiter 2 --f_amg_com_smoo 9 --f_amg_str 0.75 --f_amg_coarse 1 --print_f_hypre"

mpirun -np 1 ./$PROGRAM $PARAM --visc 0 --noel 4
mpirun -np 1 ./$PROGRAM $PARAM --visc 0 --noel 6
mpirun -np 1 ./$PROGRAM $PARAM --visc 0 --noel 8
mpirun -np 1 ./$PROGRAM $PARAM --visc 0 --noel 10
mpirun -np 1 ./$PROGRAM $PARAM --visc 0 --noel 12
mpirun -np 1 ./$PROGRAM $PARAM --visc 0 --noel 14
mpirun -np 1 ./$PROGRAM $PARAM --visc 0 --noel 16
mpirun -np 1 ./$PROGRAM $PARAM --visc 0 --noel 18

mpirun -np 1 ./$PROGRAM $PARAM --visc 1 --noel 4
mpirun -np 1 ./$PROGRAM $PARAM --visc 1 --noel 6
mpirun -np 1 ./$PROGRAM $PARAM --visc 1 --noel 8
mpirun -np 1 ./$PROGRAM $PARAM --visc 1 --noel 10
mpirun -np 1 ./$PROGRAM $PARAM --visc 1 --noel 12
mpirun -np 1 ./$PROGRAM $PARAM --visc 1 --noel 14
mpirun -np 1 ./$PROGRAM $PARAM --visc 1 --noel 16
mpirun -np 1 ./$PROGRAM $PARAM --visc 1 --noel 18
