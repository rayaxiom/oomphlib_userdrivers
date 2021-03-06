#!/bin/bash

set -e


CURRENT_DIR=`pwd`

TETGEN_DIR="tetgen_files_unit_cube"
RAWMAT_DIR="raw_lin_system"

# Set by if statement
PROG_DEBUG=""
PROG_OPT=""


CURRHOST="$HOSTNAME"
if [ "$CURRHOST" = "onigiri" ]
then
PROG_DEBUG="/home/ray/oomphlib/mpi_debug_paranoid/user_drivers/milan_permutation"
PROG_OPT="/home/ray/oomphlib/mpi_optimized/user_drivers/milan_permutation"
else # Assume that this is csf TODO
echo "not on csf but got to csf loop"
fi
  
TETGEN_DEBUG="$PROG_DEBUG/$TETGEN_DIR"
TETGEN_OPT="$PROG_OPT/$TETGEN_DIR"

RAWMAT_DEBUG="$PROG_DEBUG/$RAWMAT_DIR"
RAWMAT_OPT="$PROG_OPT/$RAWMAT_DIR"

cd $PROG_OPT

rm -rf $RAWMAT_DIR $TETGEN_DIR

ln -s $TETGEN_DEBUG $TETGEN_DIR
ln -s $RAWMAT_DEBUG $RAWMAT_DIR





