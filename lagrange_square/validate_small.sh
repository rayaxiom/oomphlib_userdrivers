#!/bin/bash

###############################################################################
# Current directory
CURRENTDIR=`pwd`

# Where the static validata is.
VALIDATADIR="validata_small"

# This is where the validation is performed. 
# This will be removed at the beginning of every validation.
VALIDATEDIR="Validate_small"
TEMPVALIDATADIR="temp_validata"



# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)
###############################################################################
## EDIT THIS
############

PROGRAM="sq_lgr"

#### Helper functions ####

function create_folders()
{
  touch $VALIDATEDIR && rm -rf $VALIDATEDIR && mkdir $VALIDATEDIR && \
  cd $VALIDATEDIR && mkdir $TEMPVALIDATADIR && cd $CURRENTDIR
}

## make the program only.
function make_prog()
{
  make $PROGRAM
}

## make src, then make the program.
function make_src_prog()
{
  cd $OOMPH_ROOT_DIR/src && make && make install && \
  cd $CURRENTDIR && make_prog
}

#function cat_test11_result()
#{
#  
#}

function test11()
{
  if [ -d "$TEMPVALIDATADIR" ] 
  then
    RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
    mpirun -np 1 ./$PROGRAM $RUNARGUMENTS
  else
    echo "================================================================"
    echo "The directory $TEMPVALIDATADIR does not exist.                  "
    echo "Please create it in `pwd`"
    echo "================================================================"
  fi
}
function test12()
{ 
  ## RAYRAY check the program folder!!!
  RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
  mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
}
function test21()
{
  if [ -d "$TEMPVALIDATADIR" ] 
  then
    RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 3 --itstimedir $TEMPVALIDATADIR"
    mpirun -np 2 ./$PROGRAM $RUNARGUMENTS
  else
    echo "================================================================"
    echo "The directory $TEMPVALIDATADIR does not exist.                  "
    echo "Please create it in `pwd`"
    echo "================================================================"
  fi
}

function test22()
{
  ## RAYRAY check the program folder!!!
  RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
  mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
}

function validate_test22()
{
# These are the files outputted by test22
files=("SqPoWeNlFePeStrA42R100N16WdNP2R0" \
"SqPoWeNlFePeStrA42R100N16WdNP2R0")

# run test22
cd $VALIDATEDIR && test22 > /dev/null && \
touch validate.log

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
}

function temp_cat_results()
{
  echo " "
  echo " "
  echo "================================================================"
  echo "test11 results: "
  cat $TEMPVALIDATADIR/SqPoWeNeStrA42R100N16WdNP1R0
  echo "================================================================"
  echo "test21 results: "
  cat $TEMPVALIDATADIR/SqPoWeNlFePeStrA42R100N16WdNP1R0
  echo "================================================================"

  ## Check if they are both correct.
files=("SqPoWeNeStrA42R100N16WdNP1R0" \
"SqPoWeNlFePeStrA42R100N16WdNP1R0")

touch validation.log && rm -rf validation.log

for i in "${files[@]}"
do
	grep "RAYITS" ./$TEMPVALIDATADIR/$i > RAYITS_new
  grep "RAYITS" ./$VALIDATADIR/$i > RAYITS_old
  
  DIFF=$(colordiff RAYITS_new RAYITS_old)
  if [ "$DIFF" != "" ]
  then
    echo "$i failed." >> validation.log
    echo $DIFF >> validation.log
  else
    echo "$i passed." >> validation.log
  fi

  rm -rf RAYITS_new RAYITS_old
done

echo " "
echo " "
cat validation.log 
}

function temp_create_temp_validata_dir()
{
  touch $TEMPVALIDATADIR && rm -rf $TEMPVALIDATADIR && mkdir $TEMPVALIDATADIR
}

#make_src_prog && create_folders && \
#test22

#make_src_prog && temp_create_temp_validata_dir && \
#test11 && test21 && \
#temp_cat_results

#make_src_prog && temp_create_temp_validata_dir && \
#test11

make_src_prog && temp_create_temp_validata_dir && \
test21


###############################################################################
#cd $OOMPH_ROOT_DIR && \
#./autogen_wipebuild_noselftest.sh --rebuild --jobs=4 && \
#cd $CURRENTDIR && \
#make $PROGRAM && \
#mpirun -np 2 ./$PROGRAM $RUNARGUMENTS


#cd $OOMPH_ROOT_DIR/src && make && make install && \
#cd $CURRENTDIR && make $PROGRAM
#make $PROGRAM && mpirun -np 1 ./$PROGRAM $RUNARGUMENTS

#make $PROGRAM

#touch $VALIDATEDIR
#rm -rf $VALIDATEDIR
#mkdir $VALIDATEDIR
#cd $VALIDATEDIR
#mkdir $TEMPVALIDATADIR

#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 0 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 0 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 0 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 0 --visc 1 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#
#
#
#
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 0 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 0 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 0 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#
#
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 8 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 16 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#RUNARGUMENTS="--prob_id 11 --w_solver 0 --ns_solver 1 --f_solver 0 --p_solver 0 --visc 1 --ang 42 --rey 100 --noel 32 --itstimedir $TEMPVALIDATADIR"
#mpirun -np 1 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 2 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 3 ../$PROGRAM $RUNARGUMENTS
#mpirun -np 4 ../$PROGRAM $RUNARGUMENTS
#
## Now I need to compare files... grep for RAYITS
#files=("SqPoWeNeFePeSimA42R100N16WdNP1R0" \
#"SqPoWeNeFePeSimA42R100N16WdNP2R0" \
#"SqPoWeNeFePeSimA42R100N16WdNP2R1" \
#"SqPoWeNeFePeSimA42R100N16WdNP3R0" \
#"SqPoWeNeFePeSimA42R100N16WdNP3R1" \
#"SqPoWeNeFePeSimA42R100N16WdNP3R2" \
#"SqPoWeNeFePeSimA42R100N16WdNP4R0" \
#"SqPoWeNeFePeSimA42R100N16WdNP4R1" \
#"SqPoWeNeFePeSimA42R100N16WdNP4R2" \
#"SqPoWeNeFePeSimA42R100N16WdNP4R3" \
#"SqPoWeNeFePeSimA42R100N32WdNP1R0" \
#"SqPoWeNeFePeSimA42R100N32WdNP2R0" \
#"SqPoWeNeFePeSimA42R100N32WdNP2R1" \
#"SqPoWeNeFePeSimA42R100N32WdNP3R0" \
#"SqPoWeNeFePeSimA42R100N32WdNP3R1" \
#"SqPoWeNeFePeSimA42R100N32WdNP3R2" 
#"SqPoWeNeFePeSimA42R100N32WdNP4R0" \
#"SqPoWeNeFePeSimA42R100N32WdNP4R1" \
#"SqPoWeNeFePeSimA42R100N32WdNP4R2" \
#"SqPoWeNeFePeSimA42R100N32WdNP4R3" \
#"SqPoWeNeFePeSimA42R100N8WdNP1R0" \
#"SqPoWeNeFePeSimA42R100N8WdNP2R0" \
#"SqPoWeNeFePeSimA42R100N8WdNP2R1" \
#"SqPoWeNeFePeSimA42R100N8WdNP3R0" \
#"SqPoWeNeFePeSimA42R100N8WdNP3R1" \
#"SqPoWeNeFePeSimA42R100N8WdNP3R2" \
#"SqPoWeNeFePeSimA42R100N8WdNP4R0" \
#"SqPoWeNeFePeSimA42R100N8WdNP4R1" \
#"SqPoWeNeFePeSimA42R100N8WdNP4R2" \
#"SqPoWeNeFePeSimA42R100N8WdNP4R3" \
#"SqPoWeNeFePeStrA42R100N16WdNP1R0" \
#"SqPoWeNeFePeStrA42R100N16WdNP2R0" \
#"SqPoWeNeFePeStrA42R100N16WdNP2R1" \
#"SqPoWeNeFePeStrA42R100N16WdNP3R0" \
#"SqPoWeNeFePeStrA42R100N16WdNP3R1" \
#"SqPoWeNeFePeStrA42R100N16WdNP3R2" \
#"SqPoWeNeFePeStrA42R100N16WdNP4R0" \
#"SqPoWeNeFePeStrA42R100N16WdNP4R1" \
#"SqPoWeNeFePeStrA42R100N16WdNP4R2" \
#"SqPoWeNeFePeStrA42R100N16WdNP4R3" \
#"SqPoWeNeFePeStrA42R100N32WdNP1R0" \
#"SqPoWeNeFePeStrA42R100N32WdNP2R0" \
#"SqPoWeNeFePeStrA42R100N32WdNP2R1" \
#"SqPoWeNeFePeStrA42R100N32WdNP3R0" \
#"SqPoWeNeFePeStrA42R100N32WdNP3R1" \
#"SqPoWeNeFePeStrA42R100N32WdNP3R2" \
#"SqPoWeNeFePeStrA42R100N32WdNP4R0" \
#"SqPoWeNeFePeStrA42R100N32WdNP4R1" \
#"SqPoWeNeFePeStrA42R100N32WdNP4R2" \
#"SqPoWeNeFePeStrA42R100N32WdNP4R3" \
#"SqPoWeNeFePeStrA42R100N8WdNP1R0" \
#"SqPoWeNeFePeStrA42R100N8WdNP2R0" \
#"SqPoWeNeFePeStrA42R100N8WdNP2R1" \
#"SqPoWeNeFePeStrA42R100N8WdNP3R0" \
#"SqPoWeNeFePeStrA42R100N8WdNP3R1" \
#"SqPoWeNeFePeStrA42R100N8WdNP3R2" \
#"SqPoWeNeFePeStrA42R100N8WdNP4R0" \
#"SqPoWeNeFePeStrA42R100N8WdNP4R1" \
#"SqPoWeNeFePeStrA42R100N8WdNP4R2" \
#"SqPoWeNeFePeStrA42R100N8WdNP4R3" \
#"SqPoWeNlFePeSimA42R100N16WdNP1R0" \
#"SqPoWeNlFePeSimA42R100N16WdNP2R0" \
#"SqPoWeNlFePeSimA42R100N16WdNP2R1" \
#"SqPoWeNlFePeSimA42R100N16WdNP3R0" \
#"SqPoWeNlFePeSimA42R100N16WdNP3R1" \
#"SqPoWeNlFePeSimA42R100N16WdNP3R2" \
#"SqPoWeNlFePeSimA42R100N16WdNP4R0" \
#"SqPoWeNlFePeSimA42R100N16WdNP4R1" \
#"SqPoWeNlFePeSimA42R100N16WdNP4R2" \
#"SqPoWeNlFePeSimA42R100N16WdNP4R3" \
#"SqPoWeNlFePeSimA42R100N32WdNP1R0" \
#"SqPoWeNlFePeSimA42R100N32WdNP2R0" \
#"SqPoWeNlFePeSimA42R100N32WdNP2R1" \
#"SqPoWeNlFePeSimA42R100N32WdNP3R0" \
#"SqPoWeNlFePeSimA42R100N32WdNP3R1" \
#"SqPoWeNlFePeSimA42R100N32WdNP3R2" \
#"SqPoWeNlFePeSimA42R100N32WdNP4R0" \
#"SqPoWeNlFePeSimA42R100N32WdNP4R1" \
#"SqPoWeNlFePeSimA42R100N32WdNP4R2" \
#"SqPoWeNlFePeSimA42R100N32WdNP4R3" \
#"SqPoWeNlFePeSimA42R100N8WdNP1R0" \
#"SqPoWeNlFePeSimA42R100N8WdNP2R0" \
#"SqPoWeNlFePeSimA42R100N8WdNP2R1" \
#"SqPoWeNlFePeSimA42R100N8WdNP3R0" \
#"SqPoWeNlFePeSimA42R100N8WdNP3R1" \
#"SqPoWeNlFePeSimA42R100N8WdNP3R2" \
#"SqPoWeNlFePeSimA42R100N8WdNP4R0" \
#"SqPoWeNlFePeSimA42R100N8WdNP4R1" \
#"SqPoWeNlFePeSimA42R100N8WdNP4R2" \
#"SqPoWeNlFePeSimA42R100N8WdNP4R3" \
#"SqPoWeNlFePeStrA42R100N16WdNP1R0" \
#"SqPoWeNlFePeStrA42R100N16WdNP2R0" \
#"SqPoWeNlFePeStrA42R100N16WdNP2R1" \
#"SqPoWeNlFePeStrA42R100N16WdNP3R0" \
#"SqPoWeNlFePeStrA42R100N16WdNP3R1" \
#"SqPoWeNlFePeStrA42R100N16WdNP3R2" \
#"SqPoWeNlFePeStrA42R100N16WdNP4R0" \
#"SqPoWeNlFePeStrA42R100N16WdNP4R1" \
#"SqPoWeNlFePeStrA42R100N16WdNP4R2" \
#"SqPoWeNlFePeStrA42R100N16WdNP4R3" \
#"SqPoWeNlFePeStrA42R100N32WdNP1R0" \
#"SqPoWeNlFePeStrA42R100N32WdNP2R0" \
#"SqPoWeNlFePeStrA42R100N32WdNP2R1" \
#"SqPoWeNlFePeStrA42R100N32WdNP3R0" \
#"SqPoWeNlFePeStrA42R100N32WdNP3R1" \
#"SqPoWeNlFePeStrA42R100N32WdNP3R2" \
#"SqPoWeNlFePeStrA42R100N32WdNP4R0" \
#"SqPoWeNlFePeStrA42R100N32WdNP4R1" \
#"SqPoWeNlFePeStrA42R100N32WdNP4R2" \
#"SqPoWeNlFePeStrA42R100N32WdNP4R3" \
#"SqPoWeNlFePeStrA42R100N8WdNP1R0" \
#"SqPoWeNlFePeStrA42R100N8WdNP2R0" \
#"SqPoWeNlFePeStrA42R100N8WdNP2R1" \
#"SqPoWeNlFePeStrA42R100N8WdNP3R0" \
#"SqPoWeNlFePeStrA42R100N8WdNP3R1" \
#"SqPoWeNlFePeStrA42R100N8WdNP3R2" \
#"SqPoWeNlFePeStrA42R100N8WdNP4R0" \
#"SqPoWeNlFePeStrA42R100N8WdNP4R1" \
#"SqPoWeNlFePeStrA42R100N8WdNP4R2" \
#"SqPoWeNlFePeStrA42R100N8WdNP4R3" )
#
#touch validation.log
#
#for i in "${files[@]}"
#do
#	grep "RAYITS" $TEMPVALIDATADIR/$i > RAYITS_new
#  grep "RAYITS" ./../$VALIDATADIR/$i > RAYITS_old
#  
#  DIFF=$(diff RAYITS_new RAYITS_old)
#  if [ "$DIFF" != "" ]
#  then
#    echo "File not the same: $i" >> validation.log
#  fi
#done
#
#cd $CURRENTDIR

#make $PROGRAM && mpirun -np 2 ./$PROGRAM $RUNARGUMENTS



