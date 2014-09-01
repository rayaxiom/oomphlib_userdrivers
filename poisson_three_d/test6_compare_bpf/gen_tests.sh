#!/bin/bash

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

PROGRAM="poisson_3d_no_bpf"
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

USE_BPF=""
generate_tests()
{
#NPROC_STR="NP${NPROC}"

USE_BPF_LIST="0 1"

AMG_STRN_LIST="0.5 0.7"

AMG_COARSE_LIST="0 3 6 8 10"

for BPF in $USE_BPF_LIST
do

for AMG_STRN in $AMG_STRN_LIST
do
  for AMG_COARSE in $AMG_COARSE_LIST
  do

USE_BPF=""
if [ "$BPF" = "1" ]; then
  USE_BPF="--use_bpf"
else
  USE_BPF=""
fi

NOEL=""
if [ "$NDOF" = "200000" ]; then
  if [ "$NPROC" = "1" ]; then
    NOEL="60"
  elif [ "$NPROC" = "2" ]; then
    NOEL="76"
  elif [ "$NPROC" = "4" ]; then
    NOEL="95"
  elif [ "$NPROC" = "8" ]; then
    NOEL="120"
  else
    NOEL="NULL"
  fi
else
  NOEL="NoelNULL"
fi



PARAM="--amg_iter 1 --amg_smiter 2 --amg_sim_smoo 0 --amg_damp 0.8 --amg_strn ${AMG_STRN} --amg_coarse ${AMG_COARSE} --noel $NOEL $USE_BPF"

echo "$RUN_COMMAND ./$PROGRAM $PARAM" >> $TESTLIST

 done
done

done

}

#############################################################################
# First compile the program and move it into this folder.
cd ..
make $PROGRAM
cd $ABS_TEST_DIR
cp ./../$PROGRAM .
cp ./../grepstuff.sh .
Files_to_copy+=($PROGRAM)
Files_to_copy+=("grepstuff.sh")

#############################################################################

# First generate tests for NDOF=200000#######################################
NDOF="200000"
USE_BPF=""

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
TESTLIST="poisson3d_ndof${NDOF}_np${NPROC}.qsub"
Files_to_copy+=($TESTLIST)

NPROC="2"
TESTLIST="poisson3d_ndof${NDOF}_np${NPROC}.qsub"
Files_to_copy+=($TESTLIST)

NPROC="4"
TESTLIST="poisson3d_ndof${NDOF}_np${NPROC}.qsub"
Files_to_copy+=($TESTLIST)

NPROC="8"
TESTLIST="poisson3d_ndof${NDOF}_np${NPROC}.qsub"
Files_to_copy+=($TESTLIST)

#############################################################################

ABS_SCRATCH_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/poisson_three_d/test6_compare_bpf/"

rm -rf $ABS_SCRATCH_DIR
mkdir -p $ABS_SCRATCH_DIR

RES_DIR="qsub_output"

mkdir -p $ABS_SCRATCH_DIR$RES_DIR

# The loop below will copy to scratch the following:
#
# The program
# The grepstuff.sh file
# All testlist.list
# All qsubs
#
#
## now loop through the above array
for i in "${Files_to_copy[@]}"
do
  rsync -av $i $ABS_SCRATCH_DIR
done

# Copy the revisions into the RES_DIR
rsync -av git_rev_oomphlib $ABS_SCRATCH_DIR$RES_DIR/
rsync -av git_rev_user_drivers $ABS_SCRATCH_DIR$RES_DIR/




