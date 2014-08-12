#!/bin/bash


mpirun -np 2 ./poisson_3d --amg_iter 1 --amg_smiter 2 --amg_sim_smoo 0 --amg_damp 0.8 --amg_strn 0.5 --amg_coarse 6 --noel 51





