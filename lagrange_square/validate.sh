#!/bin/bash

###############################################################################
# Current directory
PROGRAM_DIR=`pwd`

# Where the static validata is.
VALIDATA_DIR="validata"
VALIDATA_TAR="SqPo_validata.tar.gz" ## RAYRAY EDIT!!!

# This is where the validation is performed. 
# This will be removed at the beginning of every validation.
VALIDATE_DIR="Validate"
TEMPRES_DIR="temp_results"



# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)
###############################################################################
## EDIT THIS
############

PROGRAM="sq_lgr" ## RAYRAY EDIT!!!

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

COMMONRUNARGUMENTS="--max_solver_iter 110 --dist_prob --solver_type 2 --prob_id 11 --ang 42 --rey 100 --itstimedir $TEMPRES_DIR"
COMMONRUNARGUMENT_SIM="$COMMONRUNARGUMENTS --visc 0"
COMMONRUNARGUMENT_STR="$COMMONRUNARGUMENTS --visc 1"

PREC="--w_solver 0 --ns_solver 0"
RUNARGUMENTS="$COMMONRUNARGUMENT_SIM $PREC --noel 16"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="$COMMONRUNARGUMENT_SIM $PREC --noel 32"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS



PREC="--w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0"
RUNARGUMENTS="$COMMONRUNARGUMENT_STR $PREC --noel 16"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="$COMMONRUNARGUMENT_STR $PREC --noel 32"
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
PROBLIST="SqPo"
VISLIST="x"
RELIST="R100"
PRECLIST="WeNe WeNlF_eP_e"
ANGLIST="A42"
NOELLIST="N16 N32"
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
if [ "$PREC" == "WeNe" ]; then
  VIS="Sim"
elif [ "$PREC" == "WeNlF_eP_e" ]; then
  VIS="Str"
else
  echo "Problem with generating list of files."
  echo "The preconditioner does not make sense"
  exit 1
fi
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

