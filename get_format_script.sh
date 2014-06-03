#!/bin/bash

# This is a template for getting the format scripts.

##########################################################################
# DEFINE THESE, the rest should be automatic!
PROGRAM_DIR="lagrange_step"
FORMAT_SCRIPT="format_test7_StPo_final.sh"
##########################################################################

CURRENT_DIR=`pwd`

OOMPHOPT="/mnt/iusers01/mh01/mbax5ml3/oomphlib_optimized"
OOMPHOPTUDRI="$OOMPHOPT/user_drivers"
ABS_PROGRAM_DIR="$OOMPHOPTUDRI/$PROGRAM_DIR"

pullonly()
{
  TMPDIR=`pwd`
  cd $OOMPHOPT && \
  git fetch --all && \
  git reset --hard origin/master && \
  cd $OOMPHOPTUDRI && \
  git fetch --all && \
  git reset --hard origin/master && \
  cd $TMPDIR
}

# First get the latest from git.
pullonly

rsync -av $ABS_PROGRAM_DIR/$FORMAT_SCRIPT ./format_results.sh




