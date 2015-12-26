#!/bin/bash
CURRENTDIR=$(pwd)
NSSRCDIR="/home/ray/oomphlib/mpi_debug_paranoid/src/navier_stokes"
PROGRAM="sq_lgr_mod_AMG_diag"
PARAM="--dist_prob --prob_id 11  --max_solver_iter 1000 --itstimedir res_its_time --solver_type 2 --rey 100 --visc 0 --print_hypre --w_solver 0 --ns_solver 1 --p_solver 0 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_coarse 1 --f_amg_sim_smoo 0 --f_amg_str 0.25 --f_amg_damp 0.1 --ang 30 --f_solver 9093 --noel 32"

#
#
#
# 9090 block diagonal
# 9091 upper triangular
# 9092 lower triangular
# 9093 full AMG



#cd $NSSRCDIR && make && make install && cd $CURRENTDIR && \
#  make sq_lgr && mpirun -np 1 ./sq_lgr $PARAM

make $PROGRAM && mpirun -np 1 ./$PROGRAM $PARAM


