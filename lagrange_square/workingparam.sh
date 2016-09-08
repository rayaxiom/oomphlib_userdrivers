#!/bin/bash

make sq_lgr && \
mpirun -np 1 ./sq_lgr --prob_id 11 --solver_type 1 --w_solver 0 --max_solver_iter 100 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 30 --rey 100 --noel 8
