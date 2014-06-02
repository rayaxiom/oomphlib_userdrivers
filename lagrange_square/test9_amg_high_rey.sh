#!/bin/bash

set -e

# The name of the program
PROGRAM="sq_lgr"

# Create the FILEBASE, this is the folder where the testing will be done.
THISFILE=$0 # This contains "./", which we do not want.
THISFILE=${THISFILE:2} # Gets rid of "./"
FILEBASE=${THISFILE%%.*} # Get rid of the extension (in this case, ".sh")

## The test folder is the same as the file base.
TEST_DIR=$FILEBASE


# Create the new folder (remove old one)
touch $TEST_DIR
rm -rf $TEST_DIR
mkdir $TEST_DIR

# Get the current directory and the oomph-base
PROGRAM_DIR=`pwd`
OOMPHROOT_DIR=$(make -s --no-print-directory print-top_builddir)

# folder of where the iteration counts will be.
RESITS_DIR="res_iterations"

# Get version of oomph-lib
# The OOMPHROOT_DIR is relative to the PROGRAM_DIR
# So we have to go into PROGRAM_DIR first (to be safe).
cd $PROGRAM_DIR
cd $OOMPHROOT_DIR
git log -1 > $PROGRAM_DIR/$TEST_DIR/oomphlib_revision
cd $PROGRAM_DIR

# Get version of user drivers
cd $PROGRAM_DIR/..
git log -1 > $PROGRAM_DIR/$TEST_DIR/user_driver_revision
cd $PROGRAM_DIR

# make the program and move it into the test folder.
make $PROGRAM
mv $PROGRAM ./$TEST_DIR
cd $TEST_DIR


###############################################################################
####### WE ARE NOW INSIDE THE TEST DIRECTORY ##################################
###############################################################################

# Make the results directory. I have it in an if statement...
# Because some times when it exists, we wish to reuse it... of course not this 
# time...
if [ ! -d "$RESITS_DIR" ]; then
  mkdir $RESITS_DIR
fi

# There may be two test lists, so we declare the strings up here.
# Then we change this right before we call the gen_testxy functions.
TEST_FILEBASE=""
TEST_LIST=""

# NOTE: The F amg settings are the same as --f_solver 69

## Preconditioner parameters.
Famg_BASE="--f_solver 96"
Famg_ITER="--f_amg_iter 1"
Famg_SMITER="--f_amg_smiter 2"
Famg_SMOOTHER="--f_amg_com_smoo 9" #Euclid
Famg_DAMP="--f_amg_damp -1"
Famg_STRN_SIM="--f_amg_str 0.25"
Famg_STRN_STR="--f_amg_str 0.668"
Famg_COARSE="--f_amg_coarse 1" #RS

Prec_WLu_NSLu="--w_solver 0 --ns_solver 0"
Prec_WLu_NSLscLu="--w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0"
Prec_WLu_NSLscPamgFlu="--w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 0"

Famg_sim="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
Famg_str="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_STR $Famg_COARSE"
Prec_WLu_NSLscPLuFamgsim="--w_solver 0 --ns_solver 1 --p_solver 0 $Famg_sim"
Prec_WLu_NSLscPLuFamgstr="--w_solver 0 --ns_solver 1 --p_solver 0 $Famg_str"
Prec_WLu_NSLscPamgFamgsim="--w_solver 0 --ns_solver 1 --p_solver 1 $Famg_sim"
Prec_WLu_NSLscPamgFamgstr="--w_solver 0 --ns_solver 1 --p_solver 1 $Famg_str"


## Creates test lists for all prec combinations for noel = 4 to 128
function gen_po_tests()
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
ANGLIST="0"
RELIST="0 100 200 500 1000"
NOELLIST="4 8 16 32 64 128 256 512"

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
    PRECPARAM="$Prec_WLu_NSLu"
    ;;
  1)
    PRECPARAM="$Prec_WLu_NSLscLu"
    ;;
  2)
    PRECPARAM="$Prec_WLu_NSLscPamgFlu"
    ;;
  3)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$Prec_WLu_NSLscPLuFamgsim"
    else
      PRECPARAM="$Prec_WLu_NSLscPLuFamgstr"
    fi
    ;;
  4)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$Prec_WLu_NSLscPamgFamgsim"
    else
      PRECPARAM="$Prec_WLu_NSLscPamgFamgstr"
    fi
    ;;
esac

echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 110 --dist_prob --trilinos_solver --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST

        done
      done
    done
  done
done
} # gen_po_tests function


## Prec 0, LU for both P and F solve.
VaPREC0_Plu_Flu="--lsc_only --p_solver 0 --f_solver 0"
## Prec 1 P - AMG, F - LU
VaPREC1_Pamg_Flu="--lsc_only --p_solver 1 --f_solver 0"
## Prec 2 P - LU F - AMG
VaPREC2_Plu_Famgsim="--lsc_only --p_solver 0 $Famg_sim"
VaPREC2_Plu_Famgstr="--lsc_only --p_solver 0 $Famg_str"
## PREC 3 full AMG
VaPREC3_Pamg_Famgsim="--lsc_only --p_solver 1 $Famg_sim"
VaPREC3_Pamg_Famgstr="--lsc_only --p_solver 1 $Famg_str"

## Creates test lists for all prec combinations for noel = 4 to 128
function gen_va_tests()
{
PRECLIST="3"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0 1"
RELIST="0 100 200 500 1000"
NOELLIST="4 8 16 32 64 128 256 512"

for PREC  in $PRECLIST
do
  for VIS in $VISLIST
  do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do
case "$PREC" in
  0)
    PRECPARAM="$VaPREC0_Plu_Flu"
    ;;
  1)
    PRECPARAM="$VaPREC1_Pamg_Flu"
    ;;
  2)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$VaPREC2_Plu_Famgsim"
    else
      PRECPARAM="$VaPREC2_Plu_Famgstr"
    fi

    ;;
  3)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$VaPREC3_Pamg_Famgsim"
    else
      PRECPARAM="$VaPREC3_Pamg_Famgstr"
    fi
    ;;
esac

echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 110 --dist_prob --trilinos_solver --prob_id 88 $PRECPARAM --visc $VIS --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST

        done
      done
  done
done
} # gen_tests function



TEST_FILEBASE="tests_amg_PoVa_high_rey"
TEST_LIST="$TEST_FILEBASE.list"

gen_po_tests
gen_va_tests

echo -e "\n"

. $PROGRAM_DIR/../generate_qsub_script.sh
# Do this bit in a sub shell
(generate_qsub_script $TEST_LIST)

TEST_RUN="$TEST_FILEBASE.sh"
echo "#!/bin/bash" >> $TEST_RUN
cat $TEST_LIST >> $TEST_RUN

cp ./../$0 .

QSUBFILE="$TEST_FILEBASE.qsub"


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

  OOMPH_TEST_DIR="$OOMPH_PROGRAM_DIR/$FILEBASE"
  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE"

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
  QSUBOUTPUT_DIR="qsub_output_$TEST_FILEBASE"
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



