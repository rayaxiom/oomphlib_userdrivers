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

PROGRAM="sq_lgr"

###############################################################################
#cd $OOMPH_ROOT_DIR && \
#./autogen_wipebuild_noselftest.sh --rebuild --jobs=4 && \
#cd $CURRENTDIR && \
#make $PROGRAM && \
#mpirun -np 2 ./$PROGRAM $RUNARGUMENTS


#cd $OOMPH_ROOT_DIR/src && make && make install && \
#cd $CURRENTDIR && \
#make $PROGRAM && mpirun -np 1 ./$PROGRAM $RUNARGUMENTS

make $PROGRAM

touch $VALIDATEDIR
rm -rf $VALIDATEDIR
mkdir $VALIDATEDIR
cd $VALIDATEDIR
mkdir $TEMPVALIDATADIR

RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 0 --visc 0 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 0 --visc 0 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 0 --visc 0 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS

RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS




RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 0 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 0 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 0 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS


RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
RUNARGUMENTS="--dist_prob --prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
mpirun -np 4 ../$PROGRAM $RUNARGUMENTS

# Now I need to compare files... grep for RAYITS
files=("SqPoSimR100WedNeA42N16NP1R0" \
"SqPoSimR100WedNeA42N16NP2R0" \
"SqPoSimR100WedNeA42N16NP2R1" \
"SqPoSimR100WedNeA42N16NP3R0" \
"SqPoSimR100WedNeA42N16NP3R1" \
"SqPoSimR100WedNeA42N16NP3R2" \
"SqPoSimR100WedNeA42N16NP4R0" \
"SqPoSimR100WedNeA42N16NP4R1" \
"SqPoSimR100WedNeA42N16NP4R2" \
"SqPoSimR100WedNeA42N16NP4R3" \
"SqPoSimR100WedNeA42N32NP1R0" \
"SqPoSimR100WedNeA42N32NP2R0" \
"SqPoSimR100WedNeA42N32NP2R1" \
"SqPoSimR100WedNeA42N32NP3R0" \
"SqPoSimR100WedNeA42N32NP3R1" \
"SqPoSimR100WedNeA42N32NP3R2" \
"SqPoSimR100WedNeA42N32NP4R0" \
"SqPoSimR100WedNeA42N32NP4R1" \
"SqPoSimR100WedNeA42N32NP4R2" \
"SqPoSimR100WedNeA42N32NP4R3" \
"SqPoSimR100WedNeA42N8NP1R0" \
"SqPoSimR100WedNeA42N8NP2R0" \
"SqPoSimR100WedNeA42N8NP2R1" \
"SqPoSimR100WedNeA42N8NP3R0" \
"SqPoSimR100WedNeA42N8NP3R1" \
"SqPoSimR100WedNeA42N8NP3R2" \
"SqPoSimR100WedNeA42N8NP4R0" \
"SqPoSimR100WedNeA42N8NP4R1" \
"SqPoSimR100WedNeA42N8NP4R2" \
"SqPoSimR100WedNeA42N8NP4R3" \
"SqPoSimR100WedNlFePeA42N16NP1R0" \
"SqPoSimR100WedNlFePeA42N16NP2R0" \
"SqPoSimR100WedNlFePeA42N16NP2R1" \
"SqPoSimR100WedNlFePeA42N16NP3R0" \
"SqPoSimR100WedNlFePeA42N16NP3R1" \
"SqPoSimR100WedNlFePeA42N16NP3R2" \
"SqPoSimR100WedNlFePeA42N16NP4R0" \
"SqPoSimR100WedNlFePeA42N16NP4R1" \
"SqPoSimR100WedNlFePeA42N16NP4R2" \
"SqPoSimR100WedNlFePeA42N16NP4R3" \
"SqPoSimR100WedNlFePeA42N32NP1R0" \
"SqPoSimR100WedNlFePeA42N32NP2R0" \
"SqPoSimR100WedNlFePeA42N32NP2R1" \
"SqPoSimR100WedNlFePeA42N32NP3R0" \
"SqPoSimR100WedNlFePeA42N32NP3R1" \
"SqPoSimR100WedNlFePeA42N32NP3R2" \
"SqPoSimR100WedNlFePeA42N32NP4R0" \
"SqPoSimR100WedNlFePeA42N32NP4R1" \
"SqPoSimR100WedNlFePeA42N32NP4R2" \
"SqPoSimR100WedNlFePeA42N32NP4R3" \
"SqPoSimR100WedNlFePeA42N8NP1R0" \
"SqPoSimR100WedNlFePeA42N8NP2R0" \
"SqPoSimR100WedNlFePeA42N8NP2R1" \
"SqPoSimR100WedNlFePeA42N8NP3R0" \
"SqPoSimR100WedNlFePeA42N8NP3R1" \
"SqPoSimR100WedNlFePeA42N8NP3R2" \
"SqPoSimR100WedNlFePeA42N8NP4R0" \
"SqPoSimR100WedNlFePeA42N8NP4R1" \
"SqPoSimR100WedNlFePeA42N8NP4R2" \
"SqPoSimR100WedNlFePeA42N8NP4R3" \
"SqPoStrR100WedNeA42N16NP1R0" \
"SqPoStrR100WedNeA42N16NP2R0" \
"SqPoStrR100WedNeA42N16NP2R1" \
"SqPoStrR100WedNeA42N16NP3R0" \
"SqPoStrR100WedNeA42N16NP3R1" \
"SqPoStrR100WedNeA42N16NP3R2" \
"SqPoStrR100WedNeA42N16NP4R0" \
"SqPoStrR100WedNeA42N16NP4R1" \
"SqPoStrR100WedNeA42N16NP4R2" \
"SqPoStrR100WedNeA42N16NP4R3" \
"SqPoStrR100WedNeA42N32NP1R0" \
"SqPoStrR100WedNeA42N32NP2R0" \
"SqPoStrR100WedNeA42N32NP2R1" \
"SqPoStrR100WedNeA42N32NP3R0" \
"SqPoStrR100WedNeA42N32NP3R1" \
"SqPoStrR100WedNeA42N32NP3R2" \
"SqPoStrR100WedNeA42N32NP4R0" \
"SqPoStrR100WedNeA42N32NP4R1" \
"SqPoStrR100WedNeA42N32NP4R2" \
"SqPoStrR100WedNeA42N32NP4R3" \
"SqPoStrR100WedNeA42N8NP1R0" \
"SqPoStrR100WedNeA42N8NP2R0" \
"SqPoStrR100WedNeA42N8NP2R1" \
"SqPoStrR100WedNeA42N8NP3R0" \
"SqPoStrR100WedNeA42N8NP3R1" \
"SqPoStrR100WedNeA42N8NP3R2" \
"SqPoStrR100WedNeA42N8NP4R0" \
"SqPoStrR100WedNeA42N8NP4R1" \
"SqPoStrR100WedNeA42N8NP4R2" \
"SqPoStrR100WedNeA42N8NP4R3" \
"SqPoStrR100WedNlFePeA42N16NP1R0" \
"SqPoStrR100WedNlFePeA42N16NP2R0" \
"SqPoStrR100WedNlFePeA42N16NP2R1" \
"SqPoStrR100WedNlFePeA42N16NP3R0" \
"SqPoStrR100WedNlFePeA42N16NP3R1" \
"SqPoStrR100WedNlFePeA42N16NP3R2" \
"SqPoStrR100WedNlFePeA42N16NP4R0" \
"SqPoStrR100WedNlFePeA42N16NP4R1" \
"SqPoStrR100WedNlFePeA42N16NP4R2" \
"SqPoStrR100WedNlFePeA42N16NP4R3" \
"SqPoStrR100WedNlFePeA42N32NP1R0" \
"SqPoStrR100WedNlFePeA42N32NP2R0" \
"SqPoStrR100WedNlFePeA42N32NP2R1" \
"SqPoStrR100WedNlFePeA42N32NP3R0" \
"SqPoStrR100WedNlFePeA42N32NP3R1" \
"SqPoStrR100WedNlFePeA42N32NP3R2" \
"SqPoStrR100WedNlFePeA42N32NP4R0" \
"SqPoStrR100WedNlFePeA42N32NP4R1" \
"SqPoStrR100WedNlFePeA42N32NP4R2" \
"SqPoStrR100WedNlFePeA42N32NP4R3" \
"SqPoStrR100WedNlFePeA42N8NP1R0" \
"SqPoStrR100WedNlFePeA42N8NP2R0" \
"SqPoStrR100WedNlFePeA42N8NP2R1" \
"SqPoStrR100WedNlFePeA42N8NP3R0" \
"SqPoStrR100WedNlFePeA42N8NP3R1" \
"SqPoStrR100WedNlFePeA42N8NP3R2" \
"SqPoStrR100WedNlFePeA42N8NP4R0" \
"SqPoStrR100WedNlFePeA42N8NP4R1" \
"SqPoStrR100WedNlFePeA42N8NP4R2" \
"SqPoStrR100WedNlFePeA42N8NP4R3")

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



