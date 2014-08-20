#!/bin/bash

#############################################################################
# We want to include Jacobi smoother in our parameter studies
# But first we want to determine the optimal damping parameter.
#
# For unit square, run tests for:
# Sim, Stress
# Ang 0, 30, 67
# Rey 0 200
# Noel 128 256
# Smoother Jacobi with damping factor:
# 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#
#
#############################################################################
## These are handled by the NavierStokesProblemParameters namespace

RES_DIR="res_its_time"

DIST_PROB="--dist_prob"
PROB_ID="--prob_id 11" # RAYRAY SET
DOC_SOLN="" # RAYRAY SET
VISC="" # RAYRAY SET
REY="" # RAYRAY SET
MAX_SOLVER_ITER="--max_solver_iter 100"
ITSTIMEDIR="--itstimedir $RES_DIR"
SOLVER_TYPE="--solver_type 2"
DT=""
TIME_START=""
TIME_END=""
MESH_TYPE=""

## Concatenate the above.
#NSPP="$DIST_PROB $PROB_ID $DOC_SOLN $VISC $REY $MAX_SOLVER_ITER $ITSTIMEDIR "
#NSPP+="$SOLVER_TYPE $DT $TIME_START $TIME_END $MESH_TYPE"

#############################################################################
## These are handled by LagrangianPreconditionerHelpers namespace
DOC_PREC=""
LSC_ONLY=""
SIGMA=""
PRINT_HYPRE="--print_hypre"

## Settiing the W solver options ############################################
W_SOLVER_EXACT="--w_solver 0"
W_SOLVER_CG="" # Using 4 iteration of CG as an inner conditioner

BDW="" # We do not use this.

# Settings for the Navier Stokes block
NS_SOLVER_EXACT=""
NS_SOLVER_LSC="--ns_solver 1"

P_SOLVER_EXACT="--p_solver 0"
P_SOLVER_AMG2D=""
P_SOLVER_AMG3D=""

F_SOLVER_EXACT=""

F_SOLVER_AMG="--f_solver 96"

F_AMG_STRN_25="--f_amg_str 0.25"
F_AMG_STRN_668="--f_amg_str 0.668"
#F_AMG_STRN_75="--f_amg_str 0.75"


#//========================================================================
#  /// \short Default AMG coarsening strategy. Coarsening types include:
#  ///  0 = CLJP (parallel coarsening using independent sets)
#  ///  1 = classical RS with no boundary treatment (not recommended
#  ///      in parallel)
#  ///  3 = modified RS with 3rd pass to add C points on the boundaries
#  ///  6 = Falgout (uses 1 then CLJP using interior coarse points as
#  ///      first independent set)
#  ///  8 = PMIS (parallel coarsening using independent sets - lower
#  ///      complexities than 0, maybe also slower convergence)
#  ///  10= HMIS (one pass RS on each processor then PMIS on interior
#  ///      coarse points as first independent set)
#  ///  11= One pass RS on each processor (not recommended)
#//========================================================================
F_AMG_COARSE_CLJP="--f_amg_coarse 0"
F_AMG_COARSE_RS="--f_amg_coarse 1"
F_AMG_COARSE_MRS="--f_amg_coarse 3"
F_AMG_COARSE_FALGOUT="--f_amg_coarse 6"
F_AMG_COARSE_PMIS="--f_amg_coarse 8"
F_AMG_COARSE_HMIS="--f_amg_coarse 10"
F_AMG_COARSE_1passRS="--f_amg_coarse 11"


#   /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
#   /// include:
#   ///  0 = Jacobi 
#   ///  1 = Gauss-Seidel, sequential
#   ///      (very slow in parallel!)
#   ///  2 = Gauss-Seidel, interior points in parallel, boundary sequential
#   ///      (slow in parallel!)
#   ///  3 = hybrid Gauss-Seidel or SOR, forward solve
#   ///  4 = hybrid Gauss-Seidel or SOR, backward solve
#   ///  6 = hybrid symmetric Gauss-Seidel or SSOR
#   /// To use these methods set AMG_using_simple_smoothing to true
F_AMG_SIM_SMOO_J="--f_amg_sim_smoo 0"
F_AMG_DAMP="--f_amg_damp -1.0"

F_AMG_SIM_SMOO_GS="--f_amg_sim_smoo 1"
F_AMG_SIM_SMOO_GSiterP="--f_amg_sim_smoo 2"
F_AMG_SIM_SMOO_SOR_forward="--f_amg_sim_smoo 3"
F_AMG_SIM_SMOO_SOR_backward="--f_amg_sim_smoo 4"
F_AMG_SIM_SMOO_SSOR="--f_amg_sim_smoo 6"

#   /// \short Complex smoothing methods used in BoomerAMG. Relaxation types
#   /// are:
#   ///  6 = Schwarz
#   ///  7 = Pilut
#   ///  8 = ParaSails
#   ///  9 = Euclid
#   /// To use these methods set AMG_using_simple_smoothing to false
F_AMG_COM_SMOO_SCHWARZ="--f_amg_com_smoo 6"
F_AMG_COM_SMOO_PILUT="--f_amg_com_smoo 7"
F_AMG_COM_SMOO_PARASAILS="--f_amg_com_smoo 8"
F_AMG_COM_SMOO_EUCLID="--f_amg_com_smoo 9"

#F_AMG_ITER="--f_amg_iter 1"
#F_AMG_SMITER="--f_amg_smiter 2"
F_AMG_CYCLE_1V22="--f_amg_iter 1 --f_amg_smiter 2"
F_AMG_CYCLE_2V22="--f_amg_iter 2 --f_amg_smiter 2"


## Concatenate the AMG parameters...
PREC_COMMON="$PRINT_HYPRE $W_SOLVER_EXACT $NS_SOLVER_LSC $P_SOLVER_EXACT $F_SOLVER_AMG"

F_1V22_RS_J_SIM="$PREC_COMMON $F_AMG_CYCLE_1V22 $F_AMG_COARSE_RS $F_AMG_SIM_SMOO_J $F_AMG_STRN_25"
F_1V22_RS_J_STR="$PREC_COMMON $F_AMG_CYCLE_1V22 $F_AMG_COARSE_RS $F_AMG_SIM_SMOO_J $F_AMG_STRN_668"
F_2V22_RS_J_SIM="$PREC_COMMON $F_AMG_CYCLE_2V22 $F_AMG_COARSE_RS $F_AMG_SIM_SMOO_J $F_AMG_STRN_25"
F_2V22_RS_J_STR="$PREC_COMMON $F_AMG_CYCLE_2V22 $F_AMG_COARSE_RS $F_AMG_SIM_SMOO_J $F_AMG_STRN_668"

# We need to append --f_amg_damp x to this.

#########
## Problem specific parameters
ANG="" # Set below
NOEL="" # Set below

PROGRAM="sq_lgr"

TESTLIST_FILE="testlist.list"
function generate_tests()
{
CYCLELIST="0 1" #0 - 1V22, 1 - 2V22
VISCLIST="0 1"
ANGLIST="0 67"
REYLIST="0 200"
NOELLIST="128 256"
DAMPLIST="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"

for CYCLE in $CYCLELIST
do
  for VISC in $VISCLIST
  do
    for ANG in $ANGLIST
    do
      for REY in $REYLIST
      do
        for NOEL in $NOELLIST
        do
          for DAMP in $DAMPLIST
          do

# Set the preconditioner param
PREC_PARAM=""
if [ "$CYCLE" = "0" ]; then
  if [ "$VISC" = "0" ]; then
    PREC_PARAM="$F_1V22_RS_J_SIM --f_amg_damp $DAMP"
  else
    PREC_PARAM="$F_1V22_RS_J_STR --f_amg_damp $DAMP"
  fi
else
  if [ "$VISC" = "0" ]; then
    PREC_PARAM="$F_2V22_RS_J_SIM --f_amg_damp $DAMP"
  else
    PREC_PARAM="$F_2V22_RS_J_STR --f_amg_damp $DAMP"
  fi
fi

# Set up NSPP
NSPP="$DIST_PROB $PROB_ID $DOC_SOLN $MAX_SOLVER_ITER $ITSTIMEDIR "
NSPP+="$SOLVER_TYPE $DT $TIME_START $TIME_END $MESH_TYPE "
NSPP+="--rey $REY --visc $VISC"

# Set up problem parameters
PROB_PARAM="--ang $ANG --noel $NOEL"

echo "mpirun -np 1 ./$PROGRAM $NSPP $PREC_PARAM $PROB_PARAM" >> $TESTLIST_FILE
          done
        done
      done
    done
  done
done

}

mkdir -p $RES_DIR
rm -rf ./$RES_DIR/*
rm -rf $TESTLIST_FILE
generate_tests

CURRENT_DIR=`pwd`
cd ..
make $PROGRAM
cd $CURRENT_DIR
cp ./../$PROGRAM .


SCRATCH_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/lagrange_square/test_damping/"

rsync -av $PROGRAM $SCRATCH_DIR
rsync -av $TESTLIST_FILE $SCRATCH_DIR
rsync -av testlist.qsub $SCRATCH_DIR
rsync -av gen_tests.sh $SCRATCH_DIR
















