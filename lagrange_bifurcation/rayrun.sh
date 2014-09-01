#!/bin/bash


mpirun -np 1 ./unstructured_three_d_fluid --time_type 1 --solver_type 1 --dist_prob --max_solver_iter 1000  --time_start 0.0 --time_end 1.0  --itstimedir res_iter_times --visc 0 --rey 100 --prob_id 0 --w_solver 0 --ns_solver 1 --f_solver 1 --p_solver 1 --p_amg_iter 1 --p_amg_smiter 2 --p_amg_sim_smoo 1   --p_amg_str 0.7 --p_amg_coarse 0 --print_p_hypre --f_amg_iter 2 --f_amg_smiter 2 --print_f_hypre --f_amg_sim_smoo 1 --f_amg_damp -1.0 --f_amg_coarse 1 --f_amg_str 0.25 --mesh_area 0.2


