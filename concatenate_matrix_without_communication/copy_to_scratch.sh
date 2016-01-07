#!/bin/bash
PROGRAM="mat_cat"
QSUBFILE="strong_matcat_np16.qsub"
TESTLIST="strong_testlist_np16.list"

SCRATCHDIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/concatenate_matrix_without_communication"

make $PROGRAM && 
cp $PROGRAM $SCRATCHDIR/ && \
cp $QSUBFILE $SCRATCHDIR/ && \
cp $TESTLIST $SCRATCHDIR






