#!/bin/bash

CURRENTDIR=`pwd`

#function make_src{
#cd ../../src && make && make install && cd $CURRENTDIR
#}


cd ./../../src && make && make install && cd $CURRENTDIR && \
make sq_lgr && \
mpirun -np 2 ./sq_lgr --prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 16


