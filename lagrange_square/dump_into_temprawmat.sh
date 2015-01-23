#!/bin/bash
set -x

DUMPFOLDER="temprawmat"

# Run tests for:
# Ang = 30
# Re = 0, 200, 500
# Visc = 0, 1
# noel: 4 8, 16 32 64
rey="500"
vis="1"
CURR="A30V${vis}R${rey}N"
NOEL="8"
rm -rf $DUMPFOLDER
mkdir $DUMPFOLDER
mpirun -np 1 ./sq_lgr --max_solver_iter 110 --dist_prob --solver_type 1 --prob_id 11 --w_solver 0 --ns_solver 0 --rey ${rey} --visc ${vis} --ang 30 --noel $NOEL
FOLDER="${CURR}${NOEL}test"
rm -rf $FOLDER
mv $DUMPFOLDER $FOLDER 




