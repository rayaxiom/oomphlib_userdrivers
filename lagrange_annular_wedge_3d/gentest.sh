#!/bin/bash

OOMPHBASE="/home/ray/oomphlib/mpi_debug_paranoid"

NPROC="1"
RUNCOMMAND="mpirun -np" 

PROGRAM="./annular_wedge_threed"

## Combine MOVED TO FUNCTION
#RUNCOMMAND="$RUNCOMMAND $PROGRAM"

############################################################################

PROBID="--prob_id 0"

# 0 steady, 1 adaptive, 2 fixed.
TIMETYPE="--time_type 1"

# 0 exact, 1, oomph gmres, 2, trilinos gmres
# NOTE: This is important, we ARE using oomph gmres when doing serial, this
# is not the case with parallel, where we use trilinos gmres
SOLVERTYPE="--solver_type 2"

# distribute problem
DISTPROB="--dist_prob"

# an integer
MAXITER="--max_solver_iter 100"

# Only set if doing fixed time stepping???
DT="--dt 0.01 "
TSTART="--time_start 0.0"
TEND="--time_end 0.5"


ITSTIMEDIR="--itstimedir res_iter_times"

VISC="--visc 1"
REY="--rey 100" # CHANGE
NOEL="--noel 4" # CHANGE

## Combine MOVED TO FUNCTION
#PROBPARAM="$PROBID $TIMETYPE $SOLVERTYPE $DISTPROB $MAXITER"
#PROBPARAM="$PROBPARAM $DT $TSTART $TEND"
#PROBPARAM="$PROBPARAM $ITSTIMEDIR"
#PROBPARAM="$PROBPARAM $VISC $REY $NOEL"

############################################################################

#--replace_all_f_blocks
#--replace_modified_blocks_only
REPLACETYPE="--replace_modified_blocks_only"
WSOLVER="--w_solver 0"
NSSOLVER="--ns_solver 1"
FSOLVER="--f_solver 1"
PSOLVER="--p_solver 1"

# Part combine MOVED TO FUNCTION
#LGRPARAM="$REPLACETYPE $WSOLVER $NSSOLVER $FSOLVER $PSOLVER"

######################################

PAMGITER="--p_amg_iter 1"
PAMGSMITER="--p_amg_smiter 2"

# 0 Jac, 1 sGS, 2 pGS, 3 SORf, 4 SORb, 6 SSOR
#PAMGSMOO="--p_amg_sim_smoo 1"
PAMGSMOO="--p_amg_sim_smoo 4"
#PAMGDAMP=""
PAMGDAMP="--p_amg_damp 1"
PAMGSTR="--p_amg_str 0.7"

# 0 CLJP, 1 cRS, 3 mRS, 6 Falgout, 8 PMIS. 10 HMIS, 11 opRS
PAMGCOARSE="--p_amg_coarse 6"
PRINTP="--print_p_hypre"

# Part combine MOVED TO FUNCTION
#PPARAM="$PAMGITER $PAMGSMITER $PAMGSMOO $PAMGSTR $PAMGCOARSE $PRINTP"

######################################

FAMGITER="--f_amg_iter 1"
FAMGSMITER="--f_amg_smiter 2"

# 0 Jac, 1 sGS, 2 pGS, 3 SORf, 4 SORb, 6 SSOR
FAMGSMOO="--f_amg_sim_smoo 4"
#FAMGDAMP=""
FAMGDAMP="--f_amg_damp 1.0"

# 0 CLJP, 1 cRS, 3 mRS, 6 Falgout, 8 PMIS. 10 HMIS, 11 opRS
FAMGCOARSE="--f_amg_coarse 0"
FAMGSTR="--f_amg_str 0.9"
PRINTF="--print_f_hypre"

# Part combine
#FPARAM="$FAMGITER $FAMGSMITER $FAMGSMOO $FAMGDAMP $FAMGCOARSE $FAMGSTR $PRINTF"

######################################

## Combine
#PRECPARAM="$LGRPARAM $PPARAM $FPARAM"

############################################################################

FULLRUNCOMMAND=""

############################################################################

combine_args()
{
## Combine
RUNPROGRAMCOMMAND="$RUNCOMMAND $NPROC $PROGRAM"

## Combine
PROBPARAM="$PROBID $TIMETYPE $SOLVERTYPE $DISTPROB $MAXITER"
PROBPARAM="$PROBPARAM $DT $TSTART $TEND"
PROBPARAM="$PROBPARAM $ITSTIMEDIR"
PROBPARAM="$PROBPARAM $VISC $REY $NOEL"

# Part combine
LGRPARAM="$REPLACETYPE $WSOLVER $NSSOLVER $FSOLVER $PSOLVER" 
# Part combine
PPARAM="$PAMGITER $PAMGSMITER $PAMGSMOO $PAMGDAMP $PAMGCOARSE $PAMGSTR $PRINTP"
# Part combine
FPARAM="$FAMGITER $FAMGSMITER $FAMGSMOO $FAMGDAMP $FAMGCOARSE $FAMGSTR $PRINTF"
## Combine
PRECPARAM="$LGRPARAM $PPARAM $FPARAM"

FULLRUNCOMMAND="$RUNPROGRAMCOMMAND $PROBPARAM $PRECPARAM"
} ## Enf of function combine_args()

makesrc()
{
  CURRENTDIR=`pwd`
  cd $OOMPHBASE
  ./autogen.sh -j 2
  cd $CURRENTDIR
}


############################################################################

#rm -rf ./matrixdata/*
#make $PROGRAM

REPLACETYPE="--replace_modified_blocks_only"
REY="--rey 100"
NOEL="--noel 6"
NPROC="3"

combine_args && $FULLRUNCOMMAND 

REPLACETYPE="--replace_all_f_blocks"
combine_args && $FULLRUNCOMMAND 


#makesrc && make $PROGRAM && combine_args && \
#$FULLRUNCOMMAND





