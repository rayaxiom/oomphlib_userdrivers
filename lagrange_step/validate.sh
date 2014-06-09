#!/bin/bash

###############################################################################
# Current directory
PROGRAM_DIR=`pwd`

# Where the static validata is.
VALIDATA_DIR="validata"
VALIDATA_TAR="StPo_validata.tar.gz"

# This is where the validation is performed. 
# This will be removed at the beginning of every validation.
VALIDATE_DIR="Validate"
TEMPRES_DIR="temp_validata"



# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)
###############################################################################
## EDIT THIS
############

PROGRAM="step_po"

###############################################################################
###############################################################################
###############################################################################
###############################################################################
## Start of the validation process

make clean
make $PROGRAM

touch $VALIDATE_DIR
rm -rf $VALIDATE_DIR
mkdir $VALIDATE_DIR
cd $VALIDATE_DIR
mkdir $TEMPRES_DIR

## NOW WE ARE RUNNING PROGRAM INSIDE Validate, results are doing into 
## temp_validata

MAX_SOLVER_ITER="--max_solver_iter 110"
DIST_PROB="--dist_prob"
SOLVER_TYPE="--trilinos_solver"
PROB_ID="--prob_id 11"
ANG="--ang 42"
REY="--rey_start 0 --rey_end 50 --rey_incre 25"
VIS="--visc 1"
NOEL='--noel 8'
ITSTIMEDIR="--itstimedir $TEMPRES_DIR"


COMMONRUNARGUMENTS="$MAX_SOLVER_ITER $DIST_PROB $SOLVER_TYPE $PROB_ID $ANG $REY $VIS $NOEL $ITSTIMEDIR"

PREC="--w_solver 0 --ns_solver 0"
RUNARGUMENTS="$COMMONRUNARGUMENTS $PREC"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS

PREC="--w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0"
RUNARGUMENTS="$COMMONRUNARGUMENTS $PREC"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS


###############################################################################
###############################################################################
###############################################################################


# Now I need to compare files... grep for RAYITS
# Generate the array of files!
# Parameter list (to be concatenated into file names):
PROBLIST="StPo"
VISLIST="Str"
RELIST="R0 R25 R50"
PRECLIST="WedNe WedNlFePe"
ANGLIST="A42"
NOELLIST="N8"
PROCLIST="NP1R0 NP2R0 NP2R1 NP3R0 NP3R1 NP3R2 NP4R0 NP4R1 NP4R2 NP4R3"

files=()
for PROB in $PROBLIST
do
  for VIS in $VISLIST
  do
    for RE in $RELIST
    do
      for PREC in $PRECLIST
      do
        for ANG in $ANGLIST
        do
          for NOEL in $NOELLIST
          do
            for PROC in $PROCLIST
            do
              RESFILE="$PROB$VIS$RE$PREC$ANG$NOEL$PROC"
              files+=("$RESFILE")
            done
          done
        done
      done
    done
  done
done

###############################################################################
###############################################################################
###############################################################################
# Now do the comparison

# These needs to be defined:
# PROGRAM_DIR
# VALIDATE_DIR
# VALIDATA_DIR
# VALIDATA_TAR
# files
# TEMPRES_DIR
. ./../../validate_common_code.sh

cd $PROGRAM_DIR

