#!/bin/bash

# Compiled dir:
COMPILEDDIR="/mnt/iusers01/mh01/mbax5ml3/mpi_optimized"

# Scratch dir:
SCRATCHDIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized"

## CHANGE THIS ##############################
# Directory where the code/test is. 
TESTDIR="user_drivers/lagrange_annular_wedge_3d/test3_7_reverted_back_lgrprec_rpl_blocks"


### FIles to copy
# Initially an empty array. Fill this up later
Files_to_copy=()
Files_to_copy+=("annular_wedge_threed")
Files_to_copy+=("testlist_ndof200000_np16.list")
Files_to_copy+=("awlgr_ndof200000_np16.qsub")


## New stuff: We run the test 3 times. Each in test1, test2, test3
## Within the folders:
## test3_7_reverted_back_lgrprec_rpl_blocks/test1
## test3_7_reverted_back_lgrprec_rpl_blocks/test2
## test3_7_reverted_back_lgrprec_rpl_blocks/test3

## Each folder has:
# qsub_output folder
# res_iter_times folder folder
# stuff in Files to copy.



# Scratch CODE directory - this is the sub directory.
SCRATCHTESTDIR="$SCRATCHDIR/$TESTDIR"
SCRATCHTESTDIRT1="$SCRATCHDIR/$TESTDIR/test1"
SCRATCHTESTDIRT2="$SCRATCHDIR/$TESTDIR/test2"
SCRATCHTESTDIRT3="$SCRATCHDIR/$TESTDIR/test3"


# Remove the test folder.
rm -rf $SCRATCHTESTDIR

# Make the test folder.
mkdir -p $SCRATCHTESTDIR

# Make the three test folders.
mkdir -p $SCRATCHTESTDIRT1
mkdir -p $SCRATCHTESTDIRT2
mkdir -p $SCRATCHTESTDIRT3

# Make the two folders qsub_output and res_iter_times
mkdir -p $SCRATCHTESTDIRT1/qsub_output
mkdir -p $SCRATCHTESTDIRT1/res_iter_times

mkdir -p $SCRATCHTESTDIRT2/qsub_output
mkdir -p $SCRATCHTESTDIRT2/res_iter_times

mkdir -p $SCRATCHTESTDIRT3/qsub_output
mkdir -p $SCRATCHTESTDIRT3/res_iter_times


## Copy the three files:
# annular_wedge_threed 
# awlgr_ndof200000_np16.qsub
# testlist_ndof200000_np16.list
for i in "${Files_to_copy[@]}"
do
  rsync -av $i $SCRATCHTESTDIRT1
done

for i in "${Files_to_copy[@]}"
do
  rsync -av $i $SCRATCHTESTDIRT2
done

for i in "${Files_to_copy[@]}"
do
  rsync -av $i $SCRATCHTESTDIRT3
done

