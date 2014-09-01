#!/bin/bash

set -u

#############################################################################
#############################################################################
#############################################################################

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
#############################################################################
#############################################################################

PROGRAM="navier_stokes_cube"
#parallel R-S, CLJP, HIMS and PIMS

# Initially an empty array. Fill this up later
Files_to_copy=()

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
#
#   /// \short Complex smoothing methods used in BoomerAMG. Relaxation types
#   /// are:
#   ///  6 = Schwarz
#   ///  7 = Pilut
#   ///  8 = ParaSails
#   ///  9 = Euclid
#   /// To use these methods set AMG_using_simple_smoothing to false

TESTLIST=""
NPROC=""
NDOF=""
RUN_COMMMAND=""


########### New stuff for NS



##########




generate_tests()
{
#NPROC_STR="NP${NPROC}"


AMG_COARSE_LIST="0 3 6 8 10"
#AMG_COARSE_STR=""

OUTFILE=""

  for AMG_COARSE in $AMG_COARSE_LIST
  do

# For NDOF= 50,000 base, we want:
# NPROC=1, NOEL=38
# NPROC=2, NOEL=48
# NPROC=4, NOEL=60
# NPROC=8, NOEL=76
# NPROC=16, NOEL=96

# For NDOF= 200,000 base, we want:
# NPROC=1, NOEL=60
# NPROC=2, NOEL=76
# NPROC=4, NOEL=95
# NPROC=8, NOEL=120
# NPROC=16, NOEL=151

# For NDOF= 400,000 base, we want:
# NPROC=1, NOEL=75
# NPROC=2, NOEL=94
# NPROC=4, NOEL=119
# NPROC=8, NOEL=150
# NPROC=16, NOEL=189

# For NDOF= 800,000 base, we want:
# NPROC=1, NOEL=94
# NPROC=2, NOEL=118
# NPROC=4, NOEL=149
# NPROC=8, NOEL=188
# NPROC=16, NOEL=237
NOEL=""
if [ "$NDOF" = "50000" ]; then
  if [ "$NPROC" = "1" ]; then
    NOEL="38"
  elif [ "$NPROC" = "2" ]; then
    NOEL="48"
  elif [ "$NPROC" = "4" ]; then
    NOEL="60"
  elif [ "$NPROC" = "8" ]; then
    NOEL="76"
  else
    NOEL="96"
  fi
elif [ "$NDOF" = "200000" ]; then
  if [ "$NPROC" = "1" ]; then
    NOEL="20"
  elif [ "$NPROC" = "2" ]; then
    NOEL="25"
  elif [ "$NPROC" = "4" ]; then
    NOEL="31"
  elif [ "$NPROC" = "8" ]; then
    NOEL="39"
  else
    NOEL="NULL"
  fi
elif [ "$NDOF" = "400000" ]; then
  if [ "$NPROC" = "1" ]; then
    NOEL="75"
  elif [ "$NPROC" = "2" ]; then
    NOEL="94"
  elif [ "$NPROC" = "4" ]; then
    NOEL="119"
  elif [ "$NPROC" = "8" ]; then
    NOEL="150"
  else
    NOEL="189"
  fi
elif [ "$NDOF" = "800000" ]; then
  if [ "$NPROC" = "1" ]; then
    NOEL="94"
  elif [ "$NPROC" = "2" ]; then
    NOEL="118"
  elif [ "$NPROC" = "4" ]; then
    NOEL="149"
  elif [ "$NPROC" = "8" ]; then
    NOEL="188"
  else
    NOEL="237"
  fi
else
  NOEL="NoelNULL"
fi



TIME_TYPE="--time_type 1" # Adaptive
SOLVER_TYPE="--solver_type 2" # oomph
DIST_PROB="--dist_prob" # problem will take care of this.
MAX_SOLVER_ITER="--max_solver_iter 100" # Allow for high iteration count.

DT="" # Do not set for adaptive
TIME_START="--time_start 0" # Start from 0
TIME_END="--time_end 1" # End at 1

SOLN_DIR="tmp_soln" # Not used
RES_DIR="res_iter_times"
#DOC_SOLN="--doc_soln $SOLN_DIR"
DOC_SOLN=""
ITSTIMEDIR="--itstimedir $RES_DIR"

rm -rf $SOLN_DIR $RES_DIR
mkdir $SOLN_DIR
mkdir $RES_DIR

##
GenProbHelper="$TIME_TYPE $SOLVER_TYPE $DIST_PROB $MAX_SOLVER_ITER "
GenProbHelper+="$DT $TIME_START $TIME_END "
GenProbHelper+="$DOC_SOLN $ITSTIMEDIR"

#################################

VISC="--visc 0"
REY="--rey 200"
REY_START=""
REY_INCRE=""
REY_END=""

##
NavierStokesHelper="$VISC $REY $REY_START $REY_INCRE $REY_END"
###############################

PROB_ID="--prob_id 0"

##
ProbSpecificParam="--noel $NOEL $PROB_ID"
###############################

SIGMA=""
W_SOLVER=""
NS_SOLVER=""
F_SOLVER="--f_solver 1"
P_SOLVER="--p_solver 1"

LGR_PREC="$SIGMA $W_SOLVER $NS_SOLVER $F_SOLVER $P_SOLVER"



#  void set_defaults_for_navier_stokes_momentum_block(
#   HyprePreconditioner* hypre_preconditioner_pt)
#  {
#   // Set default settings as for 2D Poisson
#   set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
#   
#   // Change smoother type:
#   //           0=Jacobi
#   //           1=Gauss-Seidel
#   hypre_preconditioner_pt->amg_simple_smoother() = 0;
#    
#   // Set smoother damping
#   hypre_preconditioner_pt->amg_damping() = 0.5;
#   
#   // Change strength parameter for amg
#   hypre_preconditioner_pt->amg_strength() = 0.75;
#  }


F_AMG_ITER="--f_amg_iter 1"
F_AMG_SMITER="--f_amg_smiter 2"
F_AMG_SIM_SMOO="--f_amg_sim_smoo 0"
F_AMG_COM_SMOO=""
F_AMG_DAMP="--f_amg_damp 0.5"
F_AMG_STR="--f_amg_str 0.75"

# 0 = CLJP (parallel coarsening using independent sets)
# 1 = classical RS with no boundary treatment (not recommended
#     in parallel)
# 3 = modified RS with 3rd pass to add C points on the boundaries
# 6 = Falgout (uses 1 then CLJP using interior coarse points as
#     first independent set) THIS IS DEFAULT ON DOCUMENTATION
# 8 = PMIS (parallel coarsening using independent sets - lower
#     complexities than 0, maybe also slower convergence)
# 10= HMIS (one pass RS on each processor then PMIS on interior
#     coarse points as first independent set)
# 11= One pass RS on each processor (not recommended)
F_AMG_COARSE="--f_amg_coarse $AMG_COARSE"
F_AMG_PRINT="--print_f_hypre"

##
F_AMG_PREC="$F_AMG_ITER $F_AMG_SMITER "
F_AMG_PREC+="$F_AMG_SIM_SMOO $F_AMG_COM_SMOO $F_AMG_DAMP "
F_AMG_PREC+="$F_AMG_STR $F_AMG_COARSE $F_AMG_PRINT"


###### 2D poisson problem
#   // Set iterations to 1
#   hypre_preconditioner_pt->set_amg_iterations(1);
#
#   // Use simple smoother
#   hypre_preconditioner_pt->amg_using_simple_smoothing();
#   
#   // Smoother types:
#   //           0=Jacobi
#   //           1=Gauss-Seidel
#   hypre_preconditioner_pt->amg_simple_smoother() = 1;
#   
#   // AMG preconditioner
#   hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
#   
#   // Choose strength parameter for amg
#   hypre_preconditioner_pt->amg_strength() = 0.25;
#
#   // Coarsening type
#   hypre_preconditioner_pt->amg_coarsening() = 0; 

#####  3D poisson problem
#   // Set default settings as for 2D Poisson
#   set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
#   
#   // Change strength parameter for amg
#   hypre_preconditioner_pt->amg_strength() = 0.7;


P_AMG_ITER="--p_amg_iter 1"
P_AMG_SMITER="--p_amg_smiter 2"
P_AMG_SIM_SMOO="--p_amg_sim_smoo 0" # 0 Jacobi, 1 GS, 
P_AMG_COM_SMOO=""
P_AMG_DAMP="--p_amg_damp 0.8"
P_AMG_STR="--p_amg_str 0.7"
P_AMG_COARSE="--p_amg_coarse 0"
P_AMG_PRINT="--print_p_hypre"


##
P_AMG_PREC="$P_AMG_ITER $P_AMG_SMITER "
P_AMG_PREC+="$P_AMG_SIM_SMOO $P_AMG_COM_SMOO $P_AMG_DAMP "
P_AMG_PREC+="$P_AMG_STR $P_AMG_COARSE $P_AMG_PRINT"

PrecHelper="$LGR_PREC $F_AMG_PREC $P_AMG_PREC"



PARAM="$GenProbHelper $NavierStokesHelper $ProbSpecificParam $PrecHelper"


echo "$RUN_COMMAND ./$PROGRAM $PARAM" >> $TESTLIST

 done

}

#############################################################################
# First compile the program and move it into this folder.
cd ..
make $PROGRAM
cd $ABS_TEST_DIR
cp ./../$PROGRAM .
Files_to_copy+=($PROGRAM)

#############################################################################

# First generate tests for NDOF=200000#######################################
NDOF="200000"
NPROC="1"
TESTLIST="testlist_ndof${NDOF}_np${NPROC}.list"
Files_to_copy+=($TESTLIST)
rm -rf $TESTLIST
RUN_COMMAND="mpirun -np ${NPROC} taskset -c 0"
generate_tests

NPROC="2"
TESTLIST="testlist_ndof${NDOF}_np${NPROC}.list"
Files_to_copy+=($TESTLIST)
rm -rf $TESTLIST
RUN_COMMAND="mpirun -np ${NPROC} taskset -c 0,1"
generate_tests

NPROC="4"
TESTLIST="testlist_ndof${NDOF}_np${NPROC}.list"
Files_to_copy+=($TESTLIST)
rm -rf $TESTLIST
RUN_COMMAND="mpirun -np ${NPROC} taskset -c 0,1,2,3"
generate_tests

NPROC="8"
TESTLIST="testlist_ndof${NDOF}_np${NPROC}.list"
Files_to_copy+=($TESTLIST)
rm -rf $TESTLIST
RUN_COMMAND="mpirun -np ${NPROC} taskset -c 0,1,2,3,4,5,6,7"
generate_tests

###############################################################################
## Add the qsub files, these are created manually, remember that!

# First generate tests for NDOF=200000
NDOF="200000"
NPROC="1"
TESTLIST="ns_ndof${NDOF}_np${NPROC}.qsub"
Files_to_copy+=($TESTLIST)

NPROC="2"
TESTLIST="ns_ndof${NDOF}_np${NPROC}.qsub"
Files_to_copy+=($TESTLIST)

NPROC="4"
TESTLIST="ns_ndof${NDOF}_np${NPROC}.qsub"
Files_to_copy+=($TESTLIST)

NPROC="8"
TESTLIST="ns_ndof${NDOF}_np${NPROC}.qsub"
Files_to_copy+=($TESTLIST)

##############################################################################
#
#ABS_SCRATCH_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/navier_stokes/test1_ndof200000_parallel/"
#
#rm -rf $ABS_SCRATCH_DIR
#mkdir -p $ABS_SCRATCH_DIR
#
#RES_DIR="qsub_output"
#
#mkdir -p $ABS_SCRATCH_DIR$RES_DIR
#
## The loop below will copy to scratch the following:
##
## The program
## The grepstuff.sh file
## All testlist.list
## All qsubs
##
##
### now loop through the above array
#for i in "${Files_to_copy[@]}"
#do
#  rsync -av $i $ABS_SCRATCH_DIR
#done
#
## Copy the revisions into the RES_DIR
#rsync -av git_rev_oomphlib $ABS_SCRATCH_DIR$RES_DIR/
#rsync -av git_rev_user_drivers $ABS_SCRATCH_DIR$RES_DIR/


#
### Loop through the files to copy
#rsync -av $PROGRAM $SCRATCH_PATH
#rsync -av grepstuff.sh $SCRATCH_PATH
#
#
#rsync -av testlist_np1.list $SCRATCH_PATH
#rsync -av testlist_np2.list $SCRATCH_PATH
#rsync -av testlist_np4.list $SCRATCH_PATH
#rsync -av testlist_np8.list $SCRATCH_PATH
#
# Remember to write these, these are not generated by this script or any
## other script
#rsync -av testlist_np1.qsub $SCRATCH_PATH
#rsync -av testlist_np2.qsub $SCRATCH_PATH
#rsync -av testlist_np3.qsub $SCRATCH_PATH
#rsync -av testlist_np4.qsub $SCRATCH_PATH
#rsync -av testlist_np8.qsub $SCRATCH_PATH





