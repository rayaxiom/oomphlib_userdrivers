#!/bin/bash

set -e

# The name of the program
PROGRAM="two_d_annular_wedge"

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
function gen_tests_AwPo_exacts()
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
RELIST="0 100 200 500 1000"
NOELLIST="4 8 16 32 64 128 256"

COMMONPARAM="--max_solver_iter 110 --dist_prob --trilinos_solver --phi_lo 0.0 --phi_hi 90.0 --r_lo 1.0 --r_hi 3.0 --prob_id 11"

for PREC  in $PRECLIST
do
  for VIS in $VISLIST
  do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do
load_LPREC_case
# Note: I took out --dist_prob and --trilinos_solver because we ARE using OOMPHLIB's GMRES, not trilinos
echo "mpirun -np 1 ./$PROGRAM $COMMONPARAM $PRECPARAM --visc $VIS --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST
        done
      done
  done
done
} # gen_tests function
## Creates test lists for all prec combinations for noel = 4 to 128
function gen_tests_AwPo_amg()
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
RELIST="0 100 200 500 1000"
NOELLIST="4 8 16 32 64 128 256 512"

COMMONPARAM="--max_solver_iter 110 --dist_prob --trilinos_solver --phi_lo 0.0 --phi_hi 90.0 --r_lo 1.0 --r_hi 3.0 --prob_id 11"

for PREC  in $PRECLIST
do
  for VIS in $VISLIST
  do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do
load_LPREC_case
# Note: I took out --dist_prob and --trilinos_solver because we ARE using OOMPHLIB's GMRES, not trilinos
echo "mpirun -np 1 ./$PROGRAM $COMMONPARAM $PRECPARAM --visc $VIS --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST
        done
      done
  done
done
} # gen_tests function

TESTLIST_FILEBASE="tests_AwPo_initial"
TEST_LIST="$TESTLIST_FILEBASE.list"

gen_tests_AwPo_exacts
gen_tests_AwPo_amg

. $PROGRAM_DIR/../generate_qsub_script.sh
# Do this bit in a sub shell
(generate_qsub_script $TEST_LIST)

TEST_RUN="$TESTLIST_FILEBASE.sh"
echo "#!/bin/bash" >> $TEST_RUN
cat $TEST_LIST >> $TEST_RUN

cp ./../$0 .

QSUBFILE="$TESTLIST_FILEBASE.qsub"


###############################################################################
## Produce the script to format results. 
function local_format_results()
{
CURRENT_DIR=`pwd`

function format_AwPo_exacts()
{
ITSTIMEDIR="res_iterations"

cd $ITSTIMEDIR

TAG="RAYITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

LINE=""

Prob_str="AwPo"
PRECLIST="WedNe WedNlFePe WedNlFrayPe WedNlFePa"
VISLIST="Sim Str"
ANGLIST="A_"
RELIST="R0 R100 R200 R500 R1000"
NOELLIST="N2 N4 N8 N16 N32 N64 N128 N256"

for PREC in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        LINE=""
        for NOEL in $NOELLIST
        do

          RESFILE="${Prob_str}${VIS}${RE}${PREC}${ANG}${NOEL}NP1R0"
          TOKEN=""

          if [ -f $RESFILE ]
          then
            TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
          else
            TOKEN="x"
          fi
          
          LINE="$LINE $TOKEN"
        done
        echo $LINE
      done
    done
  done
done
}

function format_AwPo_amg()
{
ITSTIMEDIR="res_iterations"

cd $ITSTIMEDIR

TAG="RAYITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

LINE=""

Prob_str="AwPo"
PRECLIST="WedNlFrayPa"
VISLIST="Sim Str"
ANGLIST="A_"
RELIST="R0 R100 R200 R500 R1000"
NOELLIST="N2 N4 N8 N16 N32 N64 N128 N256 N512"

for PREC in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        LINE=""
        for NOEL in $NOELLIST
        do

          RESFILE="${Prob_str}${VIS}${RE}${PREC}${ANG}${NOEL}NP1R0"
          TOKEN=""

          if [ -f $RESFILE ]
          then
            TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
          else
            TOKEN="x"
          fi
          
          LINE="$LINE $TOKEN"
        done
        echo $LINE
      done
    done
  done
done
}


cd $CURRENT_DIR
echo -e "\n"
echo -e "Formatting EXACTS results"
format_AwPo_exacts

cd $CURRENT_DIR
echo -e "\n"
echo -e "Formatting pure AMG results"
format_AwPo_amg
}

FORMAT_RESULTS_SCRIPT="format_results.sh"
echo '#!/bin/bash' >> $FORMAT_RESULTS_SCRIPT
declare -f local_format_results >> $FORMAT_RESULTS_SCRIPT
echo 'local_format_results' >> $FORMAT_RESULTS_SCRIPT
echo -e "\n I have created the format results script: $FORMAT_RESULTS_SCRIPT"
echo -e "\n"

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
  rsync -av $OOMPH_TEST_DIR/$FORMAT_RESULTS_SCRIPT $SCRATCH_TEST_DIR/

  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  QSUBOUTPUT_DIR="qsub_output_$TESTLIST_FILEBASE"
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR
  echo -e "\n"
  echo "I have moved the files:"
  echo "$PROGRAM"
  echo "$QSUBFILE"
  echo "$TEST_LIST"
  echo "$FORMAT_RESULTS_SCRIPT"
  echo "and created directories:"
  echo "$RESITS_DIR"
  echo "$QSUBOUTPUT_DIR"
fi




