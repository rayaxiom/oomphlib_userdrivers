#!/bin/bash
PROGRAM="mat_cat"

SCRATCHDIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/concatenate_matrix_without_communication"

make $PROGRAM && cp $PROGRAM $SCRATCHDIR/






