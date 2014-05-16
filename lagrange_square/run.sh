#!/bin/bash

CURRENTDIR=`pwd`
#cd ../../src/ && make && make install && cd $CURRENTDIR && \
#make sq_lgr && mpirun -np 2 sq_lgr --prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 3

#make sq_lgr && mpirun -np 1 ./sq_lgr --dist_prob --prob_id 11 \
#  --doc_soln tmp_soln_dir --itstimedir tmp_its_dir \
#  --visc 0 --rey 100 \
#  --w_solver 0 --ns_solver 0 --ang 42 --noel 3 --doc_prec tmp_doc_prec_dir \
#  --bdw --sigma 42.001

#make sq_lgr && mpirun -np 1 ./sq_lgr --dist_prob --prob_id 11 \
#  --doc_soln tmp_soln_dir --itstimedir tmp_its_dir \
#  --visc 0 --rey 100 \
#  --w_solver 0 --ns_solver 0 --ang 42 --noel 3 --doc_prec tmp_doc_prec_dir \
#  --bdw --sigma 42.001

#make sq_lgr && mpirun -np 1 ./sq_lgr --dist_prob --prob_id 11 \
#  --itstimedir tmp_its_dir \
#  --visc 1 --rey 100 \
#  --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --ang 42 --noel 8

FAMGPARAM="--f_amg_iter 1 --f_amg_smiter 2 --f_amg_str 0.668 --f_amg_damp -1 --f_amg_coarse 1 --f_amg_sim_smoo 1"
make sq_lgr && mpirun -np 1 ./sq_lgr --prob_id 11 --w_solver 0 \
  --ns_solver 1 --p_solver 1 --f_solver 96 $FAMGPARAM --visc 1 --ang 67 --rey 200 --noel 128 \
	--print_hypre

#  --w_solver 1 --ns_solver 0 --ang 42 --noel 3 --doc_prec tmp_doc_prec_dir

