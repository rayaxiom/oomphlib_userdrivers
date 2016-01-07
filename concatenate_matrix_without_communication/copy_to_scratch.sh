#!/bin/bash
PROGRAM="mat_cat"
STRONGQSUBFILE="strong_matcat_np16.qsub"
STRONGTESTLIST="strong_testlist_np16.list"

WEAKQSUBFILE="weak_matcat_np16.qsub"
WEAKTESTLIST="weak_testlist_np16.list"


SCRATCHDIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/concatenate_matrix_without_communication"

make $PROGRAM && 
cp $PROGRAM $SCRATCHDIR/ && \
cp $STRONGQSUBFILE $SCRATCHDIR/ && \
cp $STRONGTESTLIST $SCRATCHDIR/ && \
cp $WEAKQSUBFILE $SCRATCHDIR/ && \
cp $WEAKTESTLIST $SCRATCHDIR/






