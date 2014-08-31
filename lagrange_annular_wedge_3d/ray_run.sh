#!/bin/bash

mpirun -np 1 ./annular_wedge_threed --time_type 1 --solver_type 1 --dist_prob --max_solver_iter 1000 --time_start 0.0 --time_end 1.0 --rey 200 --w_solver 0 --ns_solver 1 --f_solver 1 --p_solver 0 --f_amg_iter 2 --f_amg_smiter 2 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.75 --f_amg_coarse 1 --print_f_hypre --prob_id 0 --visc 0 --noel 4



