#!/bin/bash

set -e

# The name of the program
PROGRAM="step_po"

TEST_DIR=""
PROGRAM_DIR=""
OOMPHROOT_DIR=""

. ./../test_common_code.sh

setup_initial

# folder of where the iteration counts will be.
RESITS_DIR="res_iterations"


cd $TEST_DIR
###############################################################################
####### WE ARE NOW INSIDE THE TEST DIRECTORY ##################################
###############################################################################
mkdir $RESITS_DIR

# There may be two test lists, so we declare the strings up here.
# Then we change this right before we call the gen_testxy functions.
TESTLIST_FILEBASE=""
TEST_LIST=""

# NOTE: The F amg settings are the same as --f_solver 69

## Preconditioner parameters.
LPREC0_LU_LU=""
LPREC1_LU_LSClulu=""
LPREC2_LU_LSCamglu=""
LPREC3_LU_LSCluamgsim=""
LPREC3_LU_LSCluamgstr=""
LPREC4_LU_LSCamgamgsim=""
LPREC4_LU_LSCamgamgstr=""

###### Vanilla LSC
VPREC0_lulu=""
VPREC1_amglu=""
VPREC2_luamg_sim=""
VPREC2_luamg_str=""
VPREC3_amgamg_sim=""
VPREC3_amgamg_str=""

load_prec_param


## Creates test lists for all prec combinations for noel = 4 to 128
function gen_tests_StPo_exacts()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - W SuperLU, NS SuperLU
# 1 - W SuperLU, NS LSC: P SuperLU, F SuperLU
# 2 - W SuperLU, NS LSC: P AMG, F LU
# 3 - W SuperLU, NS LSC: P Lu, F AMG
# 4 - W Super LU, NS LSC: P AMG F AMG

PRECLIST="0 1 2 3"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0 1"
ANGLIST="0 30 67"
RELIST="0"
REPARAM="--rey_start 0 --rey_end 200 --rey_incre 25"
NOELLIST="2 4 8 16 32 64"

for PREC  in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do
case "$PREC" in
  0)
    PRECPARAM="$LPREC0_LU_LU"
    ;;
  1)
    PRECPARAM="$LPREC1_LU_LSClulu"
    ;;
  2)
    PRECPARAM="$LPREC2_LU_LSCamglu"
    ;;
  3)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$LPREC3_LU_LSCluamgsim"
    else
      PRECPARAM="$LPREC3_LU_LSCluamgstr"
    fi
    ;;
  4)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$LPREC3_LU_LSCluamgsim"
    else
      PRECPARAM="$LPREC3_LU_LSCluamgstr"
    fi
    ;;
esac
# Note: I took out --dist_prob and --trilinos_solver because we ARE using OOMPHLIB's GMRES, not trilinos
echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 1000 --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG $REPARAM --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST
        done
      done
    done
  done
done
} # gen_tests function

function gen_tests_StPo_amg()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - W SuperLU, NS SuperLU
# 1 - W SuperLU, NS LSC: P SuperLU, F SuperLU
# 2 - W SuperLU, NS LSC: P AMG, F LU
# 3 - W SuperLU, NS LSC: P Lu, F AMG
# 4 - W Super LU, NS LSC: P AMG F AMG

PRECLIST="4"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0 1"
ANGLIST="0 30 67"
RELIST="0"
REPARAM=""
NOELLIST="2 4 8 16 32 64 128"

for PREC  in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do

if [ "$ANG" -eq "0" ]
then
  REPARAM="--rey_start 0 --rey_end 500 --rey_incre 25"
else
  REPARAM="--rey_start 0 --rey_end 200 --rey_incre 25"
fi

        for NOEL in $NOELLIST
        do
case "$PREC" in
  0)
    PRECPARAM="$LPREC0_LU_LU"
    ;;
  1)
    PRECPARAM="$LPREC1_LU_LSClulu"
    ;;
  2)
    PRECPARAM="$LPREC2_LU_LSCamglu"
    ;;
  3)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$LPREC3_LU_LSCluamgsim"
    else
      PRECPARAM="$LPREC3_LU_LSCluamgstr"
    fi
    ;;
  4)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$LPREC4_LU_LSCamgamgsim"
    else
      PRECPARAM="$LPREC4_LU_LSCamgamgstr"
    fi
    ;;
esac
# Note: I took out --dist_prob and --trilinos_solver because we ARE using OOMPHLIB's GMRES, not trilinos
echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 1000 --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG $REPARAM --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST
        done
      done
    done
  done
done
} # gen_tests function

## Creates test lists for all prec combinations for noel = 4 to 128
function gen_tests_StVa_exacts()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - W SuperLU, NS SuperLU
# 1 - W SuperLU, NS LSC: P SuperLU, F SuperLU
# 2 - W SuperLU, NS LSC: P AMG, F LU
# 3 - W SuperLU, NS LSC: P Lu, F AMG
# 4 - W Super LU, NS LSC: P AMG F AMG

PRECLIST="0 1 2"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0 1"
ANGLIST="0"
RELIST="0"
REPARAM="--rey_start 0 --rey_end 200 --rey_incre 25"
NOELLIST="2 4 8 16 32 64"

for PREC  in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do
case "$PREC" in
  0)
    PRECPARAM="$VPREC0_lulu"
    ;;
  1)
    PRECPARAM="$VPREC1_amglu"
    ;;
  2)
    if [ "$VIS" -eq "0"  ]; then
      PRECPARAM="$VPREC2_luamg_sim"
    else
      PRECPARAM="$VPREC2_luamg_str"
    fi
    ;;
  3)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$VPREC3_amgamg_sim"
    else
      PRECPARAM="$VPREC3_amgamg_str"
    fi
    ;;
esac

# Note: I took out --dist_prob and --trilinos_solver because we ARE using OOMPHLIB's GMRES, not trilinos
echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 1000 --prob_id 88 $PRECPARAM --visc $VIS $REPARAM --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST
        done
      done
    done
  done
done
} # gen_tests function

function gen_tests_StVa_amg()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - W SuperLU, NS SuperLU
# 1 - W SuperLU, NS LSC: P SuperLU, F SuperLU
# 2 - W SuperLU, NS LSC: P AMG, F LU
# 3 - W SuperLU, NS LSC: P Lu, F AMG
# 4 - W Super LU, NS LSC: P AMG F AMG

PRECLIST="3"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0 1"
ANGLIST="0"
RELIST="0"
REPARAM="--rey_start 0 --rey_end 500 --rey_incre 25"
NOELLIST="2 4 8 16 32 64 128"

for PREC  in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do

case "$PREC" in
  0)
    PRECPARAM="$VPREC0_lulu"
    ;;
  1)
    PRECPARAM="$VPREC1_amglu"
    ;;
  2)
    if [ "$VIS" -eq "0"  ]; then
      PRECPARAM="$VPREC2_luamg_sim"
    else
      PRECPARAM="$VPREC2_luamg_str"
    fi
    ;;
  3)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$VPREC3_amgamg_sim"
    else
      PRECPARAM="$VPREC3_amgamg_str"
    fi
    ;;
esac

# Note: I took out --dist_prob and --trilinos_solver because we ARE using OOMPHLIB's GMRES, not trilinos
echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 1000 --prob_id 88 $PRECPARAM --visc $VIS $REPARAM --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST
        done
      done
    done
  done
done
} # gen_tests function


TESTLIST_FILEBASE="tests_StPo_Va"
TEST_LIST="$TESTLIST_FILEBASE.list"

gen_tests_StPo_exacts
gen_tests_StPo_amg
gen_tests_StVa_exacts
gen_tests_StVa_amg

. $PROGRAM_DIR/../generate_qsub_script.sh
# Do this bit in a sub shell
(generate_qsub_script $TEST_LIST)

TEST_RUN="$TESTLIST_FILEBASE.sh"
echo "#!/bin/bash" >> $TEST_RUN
cat $TEST_LIST >> $TEST_RUN

cp ./../$0 .

QSUBFILE="$TESTLIST_FILEBASE.qsub"


###############################################################################
###############################################################################
###############################################################################
###############################################################################
################### Now check if I'm on csf, if so, delete the related scratch
# and copy the current stuff there.
echo -e "\n"
OOMPH_ABS_ROOT_DIR=""
SCRATCH_ABS_ROOT_DIR=""
COPY_TO_SCRATCH=0

if [[ $HOME == *mbax5ml3* ]]
then
  OOMPH_ABS_ROOT_DIR="/mnt/iusers01/mh01/mbax5ml3/oomphlib_optimized"
  SCRATCH_ABS_ROOT_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/oomphlib_optimized"
  echo "This is csf, copying files into:"
  echo "$SCRATCH_ABS_ROOT_DIR"
  COPY_TO_SCRATCH=1
elif [[ $HOSTNAME == *onigiri* ]]
then
  OOMPH_ABS_ROOT_DIR="/home/ray/oomphlib/matthias_bpf_rewrite"
  SCRATCH_ABS_ROOT_DIR="/home/ray/oomphlib/scratch/oomphlib_optimized"
  echo "This is onigiri, testing of files into:"
  echo "$SCRATCH_ABS_ROOT_DIR"
  COPY_TO_SCRATCH=1
else
  echo "Unrecognised machine, not doing copy into scratch."
  COPY_TO_SCRATCH=0
fi

echo -e "\n"

# NOTE, we are still in TEST_DIR
if [ "$COPY_TO_SCRATCH" -eq "1" ]
then
  echo "Doing copy to scratch!"

  cd $PROGRAM_DIR
  LOCAL_PROGRAM_DIR=${PWD##*/}
  cd $TEST_DIR

  OOMPH_PROGRAM_DIR="$OOMPH_ABS_ROOT_DIR/user_drivers/$LOCAL_PROGRAM_DIR"
  SCRATCH_PROGRAM_DIR="$SCRATCH_ABS_ROOT_DIR/user_drivers/$LOCAL_PROGRAM_DIR"

  OOMPH_TEST_DIR="$OOMPH_PROGRAM_DIR/$TEST_DIR"
  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$TEST_DIR"

  echo -e "\n"
  echo "OOMPH_TEST_DIR: $OOMPH_TEST_DIR"
  echo "SCRATCH_TEST_DIR: $SCRATCH_TEST_DIR"
  echo "Removing SCRATCH_TEST_DIR and recopying"

  # Remove the scratch stuff.
  rm -rf $SCRATCH_TEST_DIR
  mkdir -p $SCRATCH_TEST_DIR

  rsync -av $OOMPH_TEST_DIR/$PROGRAM $SCRATCH_TEST_DIR/
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR/
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR/

  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  QSUBOUTPUT_DIR="qsub_output_$TESTLIST_FILEBASE"
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR
  echo -e "\n"
  echo "I have moved the files:"
  echo "$PROGRAM"
  echo "$QSUBFILE"
  echo "$TEST_LIST"
  echo "and created directories:"
  echo "$RESITS_DIR"
  echo "$QSUBOUTPUT_DIR"
fi








