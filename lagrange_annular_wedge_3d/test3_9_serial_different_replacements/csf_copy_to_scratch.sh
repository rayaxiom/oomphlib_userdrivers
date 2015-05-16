#!/bin/bash

# Compiled dir:
COMPILEDDIR="/mnt/iusers01/mh01/mbax5ml3/mpi_optimized"

# Scratch dir:
SCRATCHDIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized"

## CHANGE THIS ##############################
# Directory where the code/test is. 
TESTDIR="user_drivers/lagrange_annular_wedge_3d/test3_9_serial_different_replacements"


### FIles to copy
# Initially an empty array. Fill this up later
Files_to_copy=()
Files_to_copy+=("annular_wedge_threed")
Files_to_copy+=("testlist.list")
Files_to_copy+=("Aw3D_N36_two_diff_replacements.qsub")


## Each folder has:
# qsub_output folder
# res_iter_times folder folder
# stuff in Files to copy.

# Scratch CODE directory - this is the sub directory.
SCRATCHTESTDIR="$SCRATCHDIR/$TESTDIR"

# Remove the test folder.
rm -rf $SCRATCHTESTDIR

# Make the test folder.
mkdir -p $SCRATCHTESTDIR


# Make the two folders qsub_output and res_iter_times
mkdir -p $SCRATCHTESTDIR/qsub_output
mkdir -p $SCRATCHTESTDIR/res_iter_times


## Copy the three files:
# annular_wedge_threed 
# awlgr_ndof200000_np16.qsub
# testlist_ndof200000_np16.list
for i in "${Files_to_copy[@]}"
do
  rsync -av $i $SCRATCHTESTDIR
done


