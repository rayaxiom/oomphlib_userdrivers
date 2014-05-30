#!/bin/bash

###############################################################################
# Current directory
CURRENTDIR=`pwd`

# Where the static validata is.
VALIDATADIR="validata"

# This is where the validation is performed. 
# This will be removed at the beginning of every validation.
VALIDATEDIR="Validate"
TEMPVALIDATADIR="temp_validata"



# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)
###############################################################################
## EDIT THIS
############

PROGRAM="step_po"

###############################################################################
#cd $OOMPH_ROOT_DIR && \
#./autogen_wipebuild_noselftest.sh --rebuild --jobs=4 && \
#cd $CURRENTDIR && \
#make $PROGRAM && \
#mpirun -np 2 ./$PROGRAM $RUNARGUMENTS


#cd $OOMPH_ROOT_DIR/src && make && make install && \
#cd $CURRENTDIR && \
#make $PROGRAM && mpirun -np 1 ./$PROGRAM $RUNARGUMENTS

make clean
make $PROGRAM

# Make the Validation folder, go into it and mkdir temp_validata, this is
# where the program's results will be outputted.
touch $VALIDATEDIR
rm -rf $VALIDATEDIR
mkdir $VALIDATEDIR
cd $VALIDATEDIR
mkdir $TEMPVALIDATADIR


MAX_SOLVER_ITER="--max_solver_iter 110"
DIST_PROB="--dist_prob"
SOLVER_TYPE="--trilinos_solver"
PROB_ID="--prob_id 11"
ANG="--ang 42"
REY="--rey_start 0 --rey_end 50 --rey_incre 25"
VIS="--visc 1"
NOEL='--noel 8'
ITSTIMEDIR="--itstimedir $TEMPVALIDATADIR"


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

touch validation.log

for i in "${files[@]}"
do
	grep "RAYITS" $TEMPVALIDATADIR/$i > RAYITS_new
  grep "RAYITS" ./../$VALIDATADIR/$i > RAYITS_old
  
  DIFF=$(diff RAYITS_new RAYITS_old)
  if [ "$DIFF" != "" ]
  then
    echo "File not the same: $i" >> validation.log
  fi
done

cd $CURRENTDIR

#make $PROGRAM && mpirun -np 2 ./$PROGRAM $RUNARGUMENTS

###############################################################################

#cd $OOMPH_ROOT_DIR && \
#./autogen_wipebuild_noselftest.sh --rebuild --jobs=4 && \
#cd $CURRENTDIR && \
#make $PROGRAM && \
#mpirun -np 1 ./$PROGRAM --ns_solver 1 --visc 0 --ang 0 --rey 0 --noel 64 --itstimedir ray_temp

# I should get:
# RAYITS: 0    19 34    26.5(2)



