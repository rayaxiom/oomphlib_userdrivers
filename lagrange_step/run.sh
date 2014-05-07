#!/bin/bash

CURRENTDIR=`pwd`

#make step_po && mpirun -np 1 ./step_po --dist_prob --prob_id 11 \
#  --itstimedir tmp_its_dir \
#  --visc 0 --rey_start 0 --rey_end 200 --rey_incre 50 \
#  --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0 --ang 67 --noel 16
make step_po && mpirun -np 1 ./step_po --dist_prob --prob_id 11 \
  --itstimedir tmp_its_dir --doc_soln tmp_soln_dir \
  --visc 0 --rey_start 0 --rey_end 150 --rey_incre 50 \
  --w_solver 0 --ns_solver 0 --ang 0 --noel 4


#  --w_solver 1 --ns_solver 0 --ang 42 --noel 3 --doc_prec tmp_doc_prec_dir

