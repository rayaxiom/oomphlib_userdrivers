#!/bin/bash

CURRDIR=`pwd`
BASEDIR="/home/ray/oomphlib/trunk"
PROG="sq_temp"


function runtest {
mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 110 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 0 --rey 100 --visc 0 --ang 30 --noel 16
mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 110 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 0 --rey 100 --visc 0 --ang 67 --noel 16

mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 110 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 0 --rey 100 --visc 1 --ang 30 --noel 16
mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 110 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 0 --rey 100 --visc 1 --ang 67 --noel 16



mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 110 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0 --rey 100 --visc 0 --ang 30 --noel 16
mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 110 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0 --rey 100 --visc 0 --ang 67 --noel 16

mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 110 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0 --rey 100 --visc 1 --ang 30 --noel 16
mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 110 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0 --rey 100 --visc 1 --ang 67 --noel 16



## AMG

#### Rey 100
mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 120 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --f_solver 96 --p_solver 96 --p_amg_coarse 1 --p_amg_str 0.25 --p_amg_sim_smoo 0 --p_amg_damp 0.668 --p_amg_iter 2 --p_amg_smiter 1 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_coarse 1 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.25 --rey 100 --visc 0 --ang 30 --noel 32

#### Change ang to 67
mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 120 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --f_solver 96 --p_solver 96 --p_amg_coarse 1 --p_amg_str 0.25 --p_amg_sim_smoo 0 --p_amg_damp 0.668 --p_amg_iter 2 --p_amg_smiter 1 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_coarse 1 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.25 --rey 100 --visc 0 --ang 67 --noel 32

mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 120 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --f_solver 96 --p_solver 96 --p_amg_coarse 1 --p_amg_str 0.25 --p_amg_sim_smoo 0 --p_amg_damp 0.668 --p_amg_iter 2 --p_amg_smiter 1 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_coarse 1 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.668 --rey 100 --visc 1 --ang 30 --noel 32

mpirun -np 1 ./sq_temp --dist_prob --prob_id 11  --max_solver_iter 120 --itstimedir res_iterations --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --f_solver 96 --p_solver 96 --p_amg_coarse 1 --p_amg_str 0.25 --p_amg_sim_smoo 0 --p_amg_damp 0.668 --p_amg_iter 2 --p_amg_smiter 1 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_coarse 1 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.668 --rey 100 --visc 1 --ang 67 --noel 32

}




rm -rf res_temp && \
cd ${BASEDIR} && \
./non_interactive_autogen.sh -j 4 && \
cd ${CURRDIR} && \
make ${PROG} && \
runtest > res_temp && \
grep "RAYITS" res_temp



