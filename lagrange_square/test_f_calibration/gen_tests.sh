#!/bin/bash

## Get a handle on the directories I may require
ABS_TEST_DIR=`pwd` # This is the current directory

# Go one up to get the program directory
cd ..
ABS_PROGRAM_DIR=`pwd`

# Go one up to get the user drivers directory
cd ..
ABS_USERDRIVERS_DIR=`pwd`

# Go one up to get the oomph root directory
cd ..
ABS_OOMPHROOT_DIR=`pwd`


# Now take note of the user driver revisions and oomph-lib revisions, so
# we can reproduce these results
cd $ABS_OOMPHROOT_DIR
git log -1 > $ABS_TEST_DIR/git_rev_oomphlib

cd $ABS_USERDRIVERS_DIR
git log -1 > $ABS_TEST_DIR/git_rev_user_drivers


cd $ABS_TEST_DIR
### Now continue with the rest of the script



#############################################################################
# We want to include Jacobi smoother in our parameter studies
# But first we want to determine the optimal damping parameter.
#
# For unit square, run tests for:
# Prec: F - 
#   1)1v22 Jacobi RS, 2) 2v22 Jacobi RS, 
#   3)1v22 GS RS, 4) 2v22 GS RS
#   5) 1v22 Euclid RS
#
# Sim, Stress
# Ang 0, 30, 67
# Rey 0 200
# Noel 4, 8, 16, 32, 64, 128 256, 512
#
# NOTE: We are using OOMPH-LIB's GMRES since we are still unsure of the
# iteration counts.
#
#
#############################################################################
## These are handled by the NavierStokesProblemParameters namespace

RES_DIR="res_its_time"
QSUB_DIR="qsub_output"

DIST_PROB="--dist_prob"
PROB_ID="--prob_id 11" # RAYRAY SET
DOC_SOLN="" # RAYRAY SET
VISC="" # RAYRAY SET
REY="" # RAYRAY SET
MAX_SOLVER_ITER="--max_solver_iter 1000"
ITSTIMEDIR="--itstimedir $RES_DIR"
SOLVER_TYPE="--solver_type 1"
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
F_AMG_DAMP="--f_amg_damp 1.0"

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

# PREC 0
F_1V22_RS_J="$PREC_COMMON $F_AMG_CYCLE_1V22 $F_AMG_COARSE_RS $F_AMG_SIM_SMOO_J $F_AMG_DAMP"
# PREC 1
F_2V22_RS_J="$PREC_COMMON $F_AMG_CYCLE_2V22 $F_AMG_COARSE_RS $F_AMG_SIM_SMOO_J $F_AMG_DAMP"

# PREC 2
F_1V22_RS_GS="$PREC_COMMON $F_AMG_CYCLE_1V22 $F_AMG_COARSE_RS $F_AMG_SIM_SMOO_GS"
# PREC 3
F_2V22_RS_GS="$PREC_COMMON $F_AMG_CYCLE_2V22 $F_AMG_COARSE_RS $F_AMG_SIM_SMOO_GS"

# PREC 4
F_1V22_RS_EUCLID="$PREC_COMMON $F_AMG_CYCLE_1V22 $F_AMG_COARSE_RS $F_AMG_COM_SMOO_EUCLID"

# REMEMBER TO ADD STRN PARAMETER TO THESE!

#########
## Problem specific parameters
ANG="" # Set below
NOEL="" # Set below

PROGRAM="sq_lgr"

TESTLIST_FILE="testlist.list"
function generate_tests()
{
PRECLIST="0 1 2 3 4"
PREC_PARAM="" # Additional things needs to be set below.

VISCLIST="0 1"
ANGLIST="0 30 67"
REYLIST="0 200"
NOELLIST="4 8 16 32 64 128 256 512"

for PREC in $PRECLIST
do



  for VISC in $VISCLIST
  do


if [ "$PREC" = "0" ]; then
  PREC_PARAM=$F_1V22_RS_J
elif [ "$PREC" = "1" ]; then
  PREC_PARAM=$F_2V22_RS_J
elif [ "$PREC" = "2" ]; then
  PREC_PARAM=$F_1V22_RS_GS
elif [ "$PREC" = "3" ]; then
  PREC_PARAM=$F_2V22_RS_GS
elif [ "$PREC" = "4" ]; then
  PREC_PARAM=$F_1V22_RS_EUCLID
else
  PREC_PARAM="null"
fi



if [ "$VISC" = "0" ]; then
  PREC_PARAM+=" $F_AMG_STRN_25"
elif [ "$VISC" = "1" ]; then
  PREC_PARAM+=" $F_AMG_STRN_668"
else
  PREC_PARAM="null"
fi

    for ANG in $ANGLIST
    do
      for REY in $REYLIST
      do
        for NOEL in $NOELLIST
        do


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

}

# Remove the test list file, so we do not add new tests with old tests
rm -rf $TESTLIST_FILE

# Generate the tests 
generate_tests

# Make the program and copy it into here.
CURRENT_DIR=`pwd`
cd ..
make $PROGRAM
cd $CURRENT_DIR
cp ./../$PROGRAM .


#################################################################
ABS_SCRATCH_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/lagrange_square/test_f_calibration/"

# Make and empty and scratch directory
mkdir -p $ABS_SCRATCH_DIR
rm -rf $ABS_SCRATCH_DIR/*

# make the results output folders
mkdir -p $ABS_SCRATCH_DIR/$RES_DIR
mkdir -p $ABS_SCRATCH_DIR/$QSUB_DIR

# Copy the revisions into the RES_DIR
rsync -av git_rev_oomphlib $ABS_SCRATCH_DIR/$RES_DIR/

rsync -av $PROGRAM $ABS_SCRATCH_DIR
rsync -av $TESTLIST_FILE $ABS_SCRATCH_DIR
rsync -av testlist.qsub $ABS_SCRATCH_DIR



