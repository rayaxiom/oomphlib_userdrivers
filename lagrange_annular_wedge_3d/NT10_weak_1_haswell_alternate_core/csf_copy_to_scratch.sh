#!/bin/bash

# Compiled dir:
COMPILEDDIR="/mnt/iusers01/mh01/mbax5ml3/mpi_optimized"

# Scratch dir:
SCRATCHDIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized"

## CHANGE THIS ##############################
# Directory where the code/test is. 
TESTDIR="user_drivers/lagrange_annular_wedge_3d/NT10_weak_1_haswell_alternate_core"


### FIles to copy
# Initially an empty array. Fill this up later
Files_to_copy=()
Files_to_copy+=("annular_wedge_threed")
Files_to_copy+=("testlist_ndof200000_np1.list")
Files_to_copy+=("testlist_ndof200000_np2.list")
Files_to_copy+=("testlist_ndof200000_np4.list")
Files_to_copy+=("testlist_ndof200000_np8.list")
Files_to_copy+=("testlist_ndof200000_np16.list")
Files_to_copy+=("testlist_ndof200000_np24.list")
Files_to_copy+=("awlgr_ndof200000_np1.qsub")
Files_to_copy+=("awlgr_ndof200000_np2.qsub")
Files_to_copy+=("awlgr_ndof200000_np4.qsub")
Files_to_copy+=("awlgr_ndof200000_np8.qsub")
Files_to_copy+=("awlgr_ndof200000_np16.qsub")
Files_to_copy+=("awlgr_ndof200000_np24.qsub")

# Scratch CODE directory - this is the sub directory.
SCRATCHTESTDIR="$SCRATCHDIR/$TESTDIR"
SCRATCHTESTDIR1="$SCRATCHTESTDIR/test1"
SCRATCHTESTDIR2="$SCRATCHTESTDIR/test2"
SCRATCHTESTDIR3="$SCRATCHTESTDIR/test3"


#SCRATCHTESTQSUBDIR="$SCRATCHTESTDIR/qsub_output"
#SCRATCHTESTRESITERTIMESDIR="$SCRATCHTESTDIR/res_iter_times"

# Remove the test dir
rm -rf $SCRATCHTESTDIR

# Create the first directory.
mkdir -p $SCRATCHTESTDIR

# The individual test directories
mkdir -p $SCRATCHTESTDIR1
mkdir -p $SCRATCHTESTDIR2
mkdir -p $SCRATCHTESTDIR3

# In each sub test dir, we create res_iter_times and qsub_output
mkdir -p "$SCRATCHTESTDIR1/res_iter_times"
mkdir -p "$SCRATCHTESTDIR1/qsub_output"
mkdir -p "$SCRATCHTESTDIR2/res_iter_times"
mkdir -p "$SCRATCHTESTDIR2/qsub_output"
mkdir -p "$SCRATCHTESTDIR3/res_iter_times"
mkdir -p "$SCRATCHTESTDIR3/qsub_output"



for i in "${Files_to_copy[@]}"
do
  rsync -av $i $SCRATCHTESTDIR1
done

for i in "${Files_to_copy[@]}"
do
  rsync -av $i $SCRATCHTESTDIR2
done

for i in "${Files_to_copy[@]}"
do
  rsync -av $i $SCRATCHTESTDIR3
done






