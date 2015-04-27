#!/bin/bash


PROGRAM="annular_wedge_threed"
OOMPHROOT="/home/ray/oomphlib/mpi_debug_paranoid"
CURRDIR=`pwd`

runprog()
{
############################################################################

# 0 steady, 1 adaptive, 2 fixed.
TIMETYPE="--time_type 1"

# 0 exact, 1, oomph gmres, 2, trilinos gmres
# NOTE: This is important, we ARE using oomph gmres when doing serial, this
# is not the case with parallel, where we use trilinos gmres
SOLVERTYPE="--solver_type 2"

# distribute problem
DISTPROB="--dist_prob"

# an integer
MAXSOLVERITER="--max_solver_iter 100"

# Only set if doing fixed time stepping???
DT="--dt 0.01"

TIMESTART="--time_start 0.0"
TIMEEND="--time_end 0.5"

ITSTIMEDIR="--itstimedir res_iter_times"

VISC="--visc 1"

REY="--rey 500"

PROBID="--prob_id 0"

NOEL="--noel 4"
GENPARAM="$TIMETYPE $SOLVERTYPE $DISTPROB $MAXSOLVERITER $DT $TIMESTART $TIMEEND $ITSTIMEDIR $VISC $REY $PROBID $NOEL"
##########################################

WSOLVER="--w_solver 0"
NSSOLVER="--ns_solver 1"

# 
FSOLVER="--f_solver 1" 
FITER="--f_amg_iter 2"
FSMITER="--f_amg_smiter 2"

# 0 Jac, 1 sGS, 2 pGS, 3 SORf, 4 SORb, 6 SSOR
FSMOO="--f_amg_sim_smoo 1"

FPRINT="--print_f_hypre"
FDAMP="--f_amg_damp 1.0"

# 0 CLJP, 1 cRS, 3 mRS, 6 Falgout, 8 PMIS. 10 HMIS, 11 opRS
FCOARSE="--f_amg_coarse 1"
FSTR="--f_amg_str 0.9"
FPARAM="$FSOLVER $FITER $FSMITER $FSMOO $FDAMP $FPRINT $FCOARSE $FSTR"

PSOLVER="--p_solver 1"
PITER="--p_amg_iter 1"
PSMITER="--p_amg_smiter 2"
PSMOO="--p_amg_sim_smoo 1"
PDAMP="--p_amg_damp 1.0"
PSTR="--p_amg_str 0.7"
PCOARSE="--p_amg_coarse 0"
PPRINT="--print_p_hypre"
PPARAM="$PSOLVER $PITER $PSMITER $PSMOO $PDAMP $PSTR $PCOARSE $PPRINT"

PRECPARAM="$WSOLVER $NSSOLVER $FPARAM $PPARAM"

##########################################

ALLPARAM="$GENPARAM $PRECPARAM"


mpirun -np 2 ./$PROGRAM $ALLPARAM
}

############################################################################
####### Make source ########################################################
makesrc()
{
  TMPDIR=`pwd` &&
  cd $OOMPHROOT/src &&
  make && make install &&
  cd $TMPDIR
}

############################################################################
###### Make program ########################################################
makeprog()
{
  make $PROGRAM
}

makesrc && makeprog && runprog 

# This is with Noel = 4, on two procs, trilinos gmres.
#RAYAVGITS:		24.7(42)
#RAYAVGAVGITS:		24.7(3.0)(14)


#runprog 

# Results should be 
# RAYAVGITS:		19.7(27)
# RAYAVGAVGITS:		19.7(1.9)(14)


# These are working parameters for parallel run:
#mpirun -np 2 ./annular_wedge_threed --time_type 1 --solver_type 2 --dist_prob --max_solver_iter 100 --dt 0.01 --time_start 0.0 --time_end 0.5  --itstimedir res_iter_times --visc 0 --rey 200 --prob_id 0 --noel 6 --w_solver 0 --ns_solver 1 --p_amg_iter 1 --p_amg_smiter 2 --p_amg_sim_smoo 4 --p_amg_damp 1 --p_amg_coarse 6 --p_amg_str 0.7 --print_p_hypre --f_amg_iter 1 --f_amg_smiter 2 --f_amg_sim_smoo 4 --f_amg_damp 1 --f_amg_coarse 0 --f_amg_str 0.75 --print_f_hypre


