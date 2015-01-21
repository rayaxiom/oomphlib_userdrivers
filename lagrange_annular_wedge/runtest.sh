#!/bin/bash

# --f_amg_iter is generally 1, --f_amg_smiter is generally 2
# --f_amg_str is 0.668 for stress divergence or 0.25 for simple
#
# --f_coarse:
# 0 = CLJP (parallel coarsening using independent sets)
# 1 = classical RS with no boundary treatment (not recommended
#     in parallel)
# 3 = modified RS with 3rd pass to add C points on the boundaries
# 6 = Falgout (uses 1 then CLJP using interior coarse points as
#     first independent set)
# 8 = PMIS (parallel coarsening using independent sets - lower
#     complexities than 0, maybe also slower convergence)
# 10= HMIS (one pass RS on each processor then PMIS on interior
#     coarse points as first independent set)
# 11= One pass RS on each processor (not recommended)
#
# --f_amg_sim_smoo:
# 0 - Jacobi (Need to set damping as well) (damping needed)
# 1 - Gauss-Seidel, sequential, very slow in parallel
# 2 - GS - interior parallel, serial on boundary.
# 3 - hybrid GS or SOR, forward solve (damping needed ?)
# 4 - hybrid GS or SOR, backwards solve (damping needed ?)
# 6 - hybrid symmetric GS or SSOR. (damping needed ?)
#
# --f_amg_com_smoo:
# 6 - Schwarz
# 7 - Pilut
# 8 - ParaSails
# 9 - Euclid


# Set from NSPP
DIST_PROB="--dist_prob"
PROB_ID="--prob_id 11" ## check if this is correct
DOC_SOLN="--doc_soln soln_temp"
VISC="--visc 1"
REY="--rey 200"
MAX_SOLVER_ITER="--max_solver_iter 100"
ITSTIMEDIR="--itstimedir res_temp"
SOLVER_TYPE="--solver_type 2"

NSPP="$DIST_PROB $PROB_ID $DOC_SOLN $VISC $REY $MAX_SOLVER_ITER $ITSTIMEDIR $SOLVER_TYPE"

#############################################################################

# Set from LPH
DOC_PREC=""
LSC_ONLY=""
SIGMA=""
W_SOLVER="--w_solver 0"
BDW=""
NS_SOLVER="--ns_solver 1"

#############################################################################
P_SOLVER="--p_solver 1"

P_PREC="$P_SOLVER"

##########################
F_SOLVER="--f_solver 96"
F_AMG_STR="--f_amg_str 0.668"
F_AMG_DAMP=""
F_AMG_COARSE="--f_amg_coarse 1"

#F_AMG_SIM_SMOO="--f_amg_sim_smoo 1" # GS
F_AMG_SIM_SMOO=""
F_AMG_COM_SMOO="--f_amg_com_smoo 9" # Euclid
#F_AMG_COM_SMOO=""

F_AMG_SMOO="$F_AMG_SIM_SMOO $F_AMG_COM_SMOO"

F_AMG_ITER="--f_amg_iter 1"
F_AMG_SMITER="--f_amg_smiter 2"

F_PREC="$F_SOLVER $F_AMG_STR $F_AMG_DAMP $F_AMG_COARSE $F_AMG_SMOO $F_AMG_ITER $F_AMG_SMITER"
##########################

RESULT_PREC="$LSC_ONLY $SIGMA $W_SOLVER $BDW $NS_SOLVER $P_PREC $F_PREC"

LPH="$DOC_PREC $RESULT_PREC"



# Specific problem parameters
NOEL="--noel 4"
PHI_LO="--phi_lo 0.0"
PHI_HI="--phi_hi 90.0"
R_LO="--r_lo 1.0"
R_HI="--r_hi 3.0"
BC="--bc 0"

SPECIFIC_PROB="$NOEL $PHI_LO $PHI_HI $R_LO $R_HI $BC"


PARAM="$NSPP $LPH $SPECIFIC_PROB"


mpirun -np 1 ./two_d_annular_wedge $PARAM






