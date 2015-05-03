#!/bin/bash

# Compiled dir:
COMPILEDDIR="/mnt/iusers01/mh01/mbax5ml3/mpi_optimized"

# Scratch dir:
SCRATCHDIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized"

## CHANGE THIS ##############################
# Directory where the code/test is. 
TESTDIR="user_drivers/lagrange_annular_wedge_3d/test3_6_fixed_lgr_replacement_blocks"


### FIles to copy
# Initially an empty array. Fill this up later
Files_to_copy=()
Files_to_copy+=("annular_wedge_threed")
Files_to_copy+=("testlist_ndof200000_np1.list")
Files_to_copy+=("testlist_ndof200000_np2.list")
Files_to_copy+=("testlist_ndof200000_np4.list")
Files_to_copy+=("testlist_ndof200000_np8.list")
Files_to_copy+=("testlist_ndof200000_np16.list")
Files_to_copy+=("awlgr_ndof200000_np1.qsub")
Files_to_copy+=("awlgr_ndof200000_np2.qsub")
Files_to_copy+=("awlgr_ndof200000_np4.qsub")
Files_to_copy+=("awlgr_ndof200000_np8.qsub")
Files_to_copy+=("awlgr_ndof200000_np16.qsub")

# Scratch CODE directory - this is the sub directory.
SCRATCHTESTDIR="$SCRATCHDIR/$TESTDIR"
SCRATCHTESTQSUBDIR="$SCRATCHTESTDIR/qsub_output"
SCRATCHTESTRESITERTIMESDIR="$SCRATCHTESTDIR/res_iter_times"

rm -rf $SCRATCHTESTDIR
mkdir -p $SCRATCHTESTDIR
mkdir -p $SCRATCHTESTQSUBDIR
mkdir -p $SCRATCHTESTRESITERTIMESDIR

for i in "${Files_to_copy[@]}"
do
  rsync -av $i $SCRATCHTESTDIR
done


