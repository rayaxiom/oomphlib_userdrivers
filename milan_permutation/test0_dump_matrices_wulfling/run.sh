#!/bin/bash

CURR_DIR=`pwd`
PROG="poisson_3d_no_bpf"

cd ..
make $PROG

cd $CURR_DIR
mv ./../$PROG .



mpirun -np 1 ./poisson_3d_no_bpf --amg_iter 1 --amg_smiter 1 --amg_sim_smoo 0 --amg_strn 0.25 --amg_coarse 0 --amg_damp 0.8 --tetgenfile 11
mpirun -np 2 ./poisson_3d_no_bpf --amg_iter 1 --amg_smiter 1 --amg_sim_smoo 0 --amg_strn 0.25 --amg_coarse 0 --amg_damp 0.8 --tetgenfile 12
mpirun -np 4 ./poisson_3d_no_bpf --amg_iter 1 --amg_smiter 1 --amg_sim_smoo 0 --amg_strn 0.25 --amg_coarse 0 --amg_damp 0.8 --tetgenfile 13
mpirun -np 8 ./poisson_3d_no_bpf --amg_iter 1 --amg_smiter 1 --amg_sim_smoo 0 --amg_strn 0.25 --amg_coarse 0 --amg_damp 0.8 --tetgenfile 14





