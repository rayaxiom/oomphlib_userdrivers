#!/bin/bash

set -e

# The name of the program
PROGRAM="cube"

TEST_DIR=""
PROGRAM_DIR=""
OOMPHROOT_DIR=""

. ./../test_common_code.sh

setup_initial

# folder of where the iteration counts will be.
RESITS_DIR="res_its_RS_GS"


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
function gen_tests_cube_amg()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - W SuperLU, NS SuperLU
# 1 - W SuperLU, NS LSC: P SuperLU, F SuperLU
# 2 - W SuperLU, NS LSC: P AMG, F LU
# 3 - W SuperLU, NS LSC: P Lu, F AMG
# 4 - W Super LU, NS LSC: P AMG F AMG

PRECLIST="4"
# The precs are set according to the PRECLIST above.
PREC_PARAM=""

VISLIST="0 1"
ANGLIST="0 30"
RELIST="0 100 200 500 1000"
NOELLIST="4 6 8 10 12 14 17 27 35 39"

COMMON_PARAM="--dist_prob --print_hypre --prob_id 21 --max_solver_iter 100 --solver_type 2 --dt 0.1"

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
    PREC_PARAM="$LPREC0_LU_LU"
    ;;
  1)
    PREC_PARAM="$LPREC1_LU_LSCLuLu"
    ;;
  2)
    PREC_PARAM="$LPREC2_LU_LSCAmgLu"
    ;;
  3)
    if [ "$VIS" -eq "0" ]; then
      PREC_PARAM="$LPREC3_LU_LSCLu_C1v22_Eu_RS_sim"
    else
      PREC_PARAM="$LPREC3_LU_LSCLu_C1v22_Eu_RS_str"
    fi
    ;;
  4)
    if [ "$VIS" -eq "0" ]; then
      PREC_PARAM="--w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 96 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.25 --f_amg_coarse 1"
    else
      PREC_PARAM="--w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 96 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.668 --f_amg_coarse 1"
    fi
    ;;
esac
# Note: I took out --dist_prob and --trilinos_solver because we ARE using OOMPHLIB's GMRES, not trilinos
echo "mpirun -np 1 ./$PROGRAM $COMMON_PARAM $PREC_PARAM --visc $VIS --ang $ANG --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST
        done
      done
    done
  done
done
} # gen_tests function

TESTLIST_FILEBASE="tests_cube_RS_GS"
TEST_LIST="$TESTLIST_FILEBASE.list"

gen_tests_cube_amg


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

function format_cube_amg()
{
ITSTIMEDIR="res_its_RS_GS"

cd $ITSTIMEDIR

TAG="RAYAVGITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

LINE=""

Prob_str="CuPoQ"
VISLIST="Sim Str"
RELIST="R0 R100 R200 R500 R1000"
PRECLIST="WedNe WedNlFrayPa"
ANGLIST="Ax0y0z0 Ax30y30z30"
NOELLIST="N4 N6 N8 N10 N12 N14 N17 N27 N35 N39"

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
echo -e "Formatting pure AMG results"
format_cube_amg
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




