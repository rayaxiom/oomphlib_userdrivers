#!/bin/bash


PROGRAM="two_d_annular_wedge"

MAX_ITER="--max_solver_iter 110"
DIST_PROB="--dist_prob"
TRILINOS_SOLVER="--trilinos_solver"
PHI_LO="--phi_lo 0.0"
PHI_HI="--phi_hi 90.0"
R_LO="--r_lo 1.0"
R_HI="--r_hi 3.0"
PROB_ID="--prob_id 11"
#ANG="--ang 42"
REY="--rey 200"
VIS="--visc 0"
NOEL="--noel 128"

W_SOLVER="--w_solver 0"
NS_SOLVER="--ns_solver 1"
#F_SOLVER="--f_solver 0"
F_SOLVER="--f_solver 96 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_com_smoo 9 --f_amg_damp -1 --f_amg_str 0.25 --f_amg_coarse 1"

P_SOLVER="--p_solver 1"

PREC="$W_SOLVER $NS_SOLVER $F_SOLVER $P_SOLVER"

PARAM="$MAX_ITER $DIST_PROB $TRILINOS_SOLVER $PROB_ID $PHI_LO $PHI_HI $R_LO $R_HI $REY $VIS $NOEL $PREC"


make $PROGRAM && mpirun -np 1 ./$PROGRAM $PARAM


