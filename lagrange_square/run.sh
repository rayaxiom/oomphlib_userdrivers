#!/bin/bash

CURRENTDIR=`pwd`
#cd ../../src/ && make && make install && cd $CURRENTDIR && \
#make sq_lgr && mpirun -np 2 sq_lgr --prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 3

make sq_lgr && mpirun -np 1 ./sq_lgr --dist_prob --prob_id 11 \
  --doc_soln tmp_soln_dir --itstimedir tmp_its_dir \
  --visc 0 --rey 100 \
  --w_solver 0 --ns_solver 0 --ang 42 --noel 3 --doc_prec tmp_doc_prec_dir \
  --bdw --sigma 42.001

#  --w_solver 1 --ns_solver 0 --ang 42 --noel 3 --doc_prec tmp_doc_prec_dir

