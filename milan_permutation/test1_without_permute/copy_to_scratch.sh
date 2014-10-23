#!/bin/bash
set -e

# For this test, all files are hand crafted.

############################################################################
############################################################################
############################################################################

## Get a handle on the directories I may require

# This is the current directory
# I.e. ~/oomph/user_drivers/milan_permutation/test1_without_permute
ABS_TEST_DIR=`pwd`

# Go one up to get the program directory
# I.e. ~/oomph/user_drivers/milan_permutation
cd ..
ABS_PROGRAM_DIR=`pwd`

# Go one up to get the user drivers directory
# I.e. ~/oomph/user_drivers
cd ..
ABS_USERDRIVERS_DIR=`pwd`

# Go one up to get the oomph root directory
# I.e. ~/oomph
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

# A list of files to copy to the scratch directory.
# Initially an empty array. Fill this up later
Files_to_copy=()

### Put in the rev. so we can re-create the results (hopefully!)
Files_to_copy+=("git_rev_oomphlib")
Files_to_copy+=("git_rev_user_drivers")

############################################################################
############################################################################
############################################################################

PROGRAM="poisson_3d_no_bpf"

# First compile the program and move it into the test folder.
cd ..
make $PROGRAM
cd $ABS_TEST_DIR
cp ./../$PROGRAM .
Files_to_copy+=($PROGRAM)


### Next we generate the list of tests.

### Recap on list of coarsening strategies:
#//========================================================================
# Default AMG coarsening strategy. Coarsening types include:
# 0 = CLJP (parallel coarsening using independent sets)
# 1 = classical RS wilibth no boundary treatment (not recommended
#     in parallel)
# 3 = modified RS with 3rd pass to add C points on the boundaries
# 6 = Falgout (uses 1 then CLJP using interior coarse points as
#     first independent set)
# 8 = PMIS (parallel coarsening using independent sets - lower
#     complexities than 0, maybe also slower convergence)
# 10= HMIS (one pass RS on each processor then PMIS on interior
#     coarse points as first independent set)
# 11= One pass RS on each processor (not recommended)
#//========================================================================


### Recap on list of smoothers:
# Simple smoothing methods used in BoomerAMG. Relaxation types
# include:
#  0 = Jacobi 
#  1 = Gauss-Seidel, sequential
#      (very slow in parallel!)
#  2 = Gauss-Seidel, interior points in parallel, boundary sequential
#      (slow in parallel!)
#  3 = hybrid Gauss-Seidel or SOR, forward solve
#  4 = hybrid Gauss-Seidel or SOR, backward solve
#  6 = hybrid symmetric Gauss-Seidel or SSOR
# To use these methods set AMG_using_simple_smoothing to true.
#
#
# Complex smoothing methods used in BoomerAMG. Relaxation types
# are:
# 6 = Schwarz
# 7 = Pilut
# 8 = ParaSails
# 9 = Euclid
# To use these methods set AMG_using_simple_smoothing to false


### But the test list is already written, they are hand crafted.
# They are in:
# testlist_ndof200000_np1.list
# testlist_ndof200000_np2.list
# testlist_ndof200000_np4.list
# testlist_ndof200000_np8.list

# Put these in the Files to copy.

Files_to_copy+=("testlist_ndof200000_np1.list")
Files_to_copy+=("testlist_ndof200000_np2.list")
Files_to_copy+=("testlist_ndof200000_np4.list")
Files_to_copy+=("testlist_ndof200000_np8.list")

# Now add the qsub files to the files to copy.
Files_to_copy+=("poisson3d_ndof200000_np1.qsub")
Files_to_copy+=("poisson3d_ndof200000_np2.qsub")
Files_to_copy+=("poisson3d_ndof200000_np4.qsub")
Files_to_copy+=("poisson3d_ndof200000_np8.qsub")


############################################################################
ABS_SCRATCH_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/milan_permutation/test1_without_permute/"

rm -rf $ABS_SCRATCH_DIR
mkdir -p $ABS_SCRATCH_DIR

QSUB_OUTPUT_DIR="qsub_output"

mkdir -p $ABS_SCRATCH_DIR$QSUB_OUTPUT_DIR

# The loop below will copy to scratch the following:
#
# Revision files:
# git_rev_oomphlib
# git_rev_user_drivers  
# 
# The program:
# poisson_3d_no_bpf
#
# testlist.list:
# testlist_ndof200000_np1.list
# testlist_ndof200000_np2.list
# testlist_ndof200000_np4.list
# testlist_ndof200000_np8.list
# 
# qsubs:
# poisson3d_ndof200000_np1.qsub
# poisson3d_ndof200000_np2.qsub
# poisson3d_ndof200000_np4.qsub
# poisson3d_ndof200000_np8.qsub
# 
## now loop through the above array
for i in "${Files_to_copy[@]}"
do
  rsync -av $i $ABS_SCRATCH_DIR
done


## Important: We now link the tetgen folder.
TETGEN_DIR="tetgen_files_unit_cube"

ABS_TETGEN_DIR="$ABS_PROGRAM_DIR/$TETGEN_DIR"
cd $ABS_SCRATCH_DIR

ln -s $ABS_TETGEN_DIR $TETGEN_DIR


