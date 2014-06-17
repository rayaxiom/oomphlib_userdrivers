#!/bin/bash

set -e

# --f_amg_iter is generally 1, --f_amg_smiter is generally 2
# --f_amg_str is 0.668 for stress divergence or 0.25 for simple
#
# --f_coarse:
# 0 = CLJP (parallel coarsening using independent sets)
# 1 = classical RS with no boundary treatment (not recommended
#     in parallel)
# 3 = modified RS with 3rd pass to add C points on the boundaries
# 6 = Falgout (uses 1 then CLJP using interior coarse points as
#     first independent set)
# 8 = PMIS (parallel coarsening using independent sets - lower
#     complexities than 0, maybe also slower convergence)
# 10= HMIS (one pass RS on each processor then PMIS on interior
#     coarse points as first independent set)
# 11= One pass RS on each processor (not recommended)
#
# --f_amg_sim_smoo:
# 0 - Jacobi (Need to set damping as well) (damping needed)
# 1 - Gauss-Seidel, sequential, very slow in parallel
# 2 - GS - interior parallel, serial on boundary.
# 3 - hybrid GS or SOR, forward solve (damping needed ?)
# 4 - hybrid GS or SOR, backwards solve (damping needed ?)
# 6 - hybrid symmetric GS or SSOR. (damping needed ?)
#
# --f_amg_com_smoo:
# 6 - Schwarz
# 7 - Pilut
# 8 - ParaSails
# 9 - Euclid
#
# --f_amg_damp should only be positive for Jacobi smoother, (f_amg_sim_smoo 0)
# if we are not using Jacobi for smoothing, then the damping value is ignored.
#
#
# --f_solver 69 is equivalent to these settings:
#FAMG_ITER="--f_amg_iter 1"
#FAMG_SMITER="--f_amg_smiter 2"
#FAMG_SSMOOTHER="--f_amg_sim_smoo 1"
#FAMG_CSMOOTHER=""
#FAMG_DAMP="--f_amg_damp -1"
#FAMG_STRN="--f_amg_str 0.668"
#FAMG_COARSE="--f_amg_coarse 1"

# The name of the program
PROGRAM="sq_lgr"

TEST_DIR=""
PROGRAM_DIR=""
OOMPHROOT_DIR=""

. ./../test_common_code.sh

setup_initial

## Load the preconditioners... preconditioner parameters:

## All exact
LPREC0_LU_LU=""

## Exact LSC
LPREC1_LU_LSCLuLu=""

## AMG for pressure block.
LPREC2_LU_LSCAmgLu=""

### AMG for LSC block.
LPREC3_LU_LSCLu_C1v22_GS_RS_sim=""
LPREC3_LU_LSCLu_C1v22_GS_RS_str=""

LPREC3_LU_LSCLu_C2v22_GS_RS_sim=""
LPREC3_LU_LSCLu_C2v22_GS_RS_str=""

LPREC3_LU_LSCLu_C1v22_Eu_RS_sim=""
LPREC3_LU_LSCLu_C1v22_Eu_RS_str=""

### Full AMG
LPREC4_LU_LSCAmg_C1v22_GS_RS_sim=""
LPREC4_LU_LSCAmg_C1v22_GS_RS_str=""

LPREC4_LU_LSCAmg_C2v22_GS_RS_sim=""
LPREC4_LU_LSCAmg_C2v22_GS_RS_str=""

LPREC4_LU_LSCAmg_C1v22_Eu_RS_sim=""
LPREC4_LU_LSCAmg_C1v22_Eu_RS_str=""

################################### Vanilla prec
## Exact for both P and F solve
VPREC0_Lu_Lu=""

## P - AMG, F exact
VPREC1_Amg_Lu=""

## P - Exact, F - AMG
VPREC2_Lu_C1v22_GS_RS_sim=""
VPREC2_Lu_C1v22_GS_RS_str=""

VPREC2_Lu_C2v22_GS_RS_sim=""
VPREC2_Lu_C2v22_GS_RS_str=""

VPREC2_Lu_C1v22_Eu_RS_sim=""
VPREC2_Lu_C1v22_Eu_RS_str=""

## P - Amg, F - AMG
VPREC3_Amg_C1v22_GS_RS_sim=""
VPREC3_Amg_C1v22_GS_RS_str=""

VPREC3_Amg_C2v22_GS_RS_sim=""
VPREC3_Amg_C2v22_GS_RS_str=""

VPREC3_Amg_C1v22_Eu_RS_sim=""
VPREC3_Amg_C1v22_Eu_RS_str=""

load_prec_param



# folder of where the iteration counts will be.
RESITS_DIR=""


cd $TEST_DIR
###############################################################################
####### WE ARE NOW INSIDE THE TEST DIRECTORY ##################################
###############################################################################
## REMEMBER: For this case, this will be made later.
##mkdir $RESITS_DIR

# There may be two test lists, so we declare the strings up here.
# Then we change this right before we call the gen_testxy functions.
TESTLIST_FILEBASE=""
TEST_LIST=""

## Common problem parameters ######################

## The lists (param to loop through)
PRECLIST="4" # 4 is full AMG, check out test_common_code.sh
VISLIST="0 1"
ANGLIST="0 30 67"
RELIST="200"
NOELLIST="4 8 16 32 64 128 256 512"

# Other problem parameters:
PROB_ID="--prob_id 11"
MAX_SOLVER_ITER="--max_solver_iter 110"
DIST_PROB="--dist_prob"
ITERATIVE_SOLVER="--trilinos_solver"
PRINT_HYPRE="--print_hypre"

COMMON_PARAM="$PROB_ID $MAX_SOLVER_ITER $DIST_PROB $ITERATIVE_SOLVER $PRINT_HYPRE"


#Famg_BASE="--f_solver 96"
#Famg_STRN_SIM="--f_amg_str 0.25"
#Famg_STRN_STR="--f_amg_str 0.668"
#Famg_COARSE="--f_amg_coarse 1" #RS



###############################################################################
###############################################################################
## Test batch 1: Cycle: 1V22, Smoothing: GS, Coarsening: RS
###############################################################################
###############################################################################
function gen_tests_1v22_GS_RS()
{
PREC_PARAM=""

for PREC in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do

load_LPREC_C1v22_GS_RS_case

echo "mpirun -np 1 ./$PROGRAM $COMMON_PARAM $PREC_PARAM --visc $VIS --ang $ANG --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST
        done
      done
    done
  done
done
}

###############################################################################
###############################################################################
## Test batch 2: Cycle: 2V22, Smoothing: GS, Coarsening: RS
###############################################################################
###############################################################################
function gen_tests_2v22_GS_RS()
{
PREC_PARAM=""
for PREC in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do


load_LPREC_C2v22_GS_RS_case

echo "mpirun -np 1 ./$PROGRAM $COMMON_PARAM $PREC_PARAM --visc $VIS --ang $ANG --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST
        done
      done
    done
  done
done
}

###############################################################################
###############################################################################
## Test batch 3: Cycle: 1V22, Smoothing: Euclid, Coarsening: RS
###############################################################################
###############################################################################
function gen_tests_1v22_EUCLID_RS()
{
PREC_PARAM=""
for PREC in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do


load_LPREC_C1v22_Eu_RS_case

echo "mpirun -np 1 ./$PROGRAM $COMMON_PARAM $PREC_PARAM --visc $VIS --ang $ANG --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST
        done
      done
    done
  done
done
}

TESTLIST_FILEBASE="tests_varyAmg"
TEST_LIST="$TESTLIST_FILEBASE.list"

RESITS_DIR="res_its_1v22_GS_RS_T1"
gen_tests_1v22_GS_RS
RESITS_DIR="res_its_1v22_GS_RS_T2"
gen_tests_1v22_GS_RS
RESITS_DIR="res_its_1v22_GS_RS_T3"
gen_tests_1v22_GS_RS

RESITS_DIR="res_its_2v22_GS_RS_T1"
gen_tests_2v22_GS_RS
RESITS_DIR="res_its_2v22_GS_RS_T2"
gen_tests_2v22_GS_RS
RESITS_DIR="res_its_2v22_GS_RS_T3"
gen_tests_2v22_GS_RS

RESITS_DIR="res_its_1v22_EUCLID_RS_T1"
gen_tests_1v22_EUCLID_RS
RESITS_DIR="res_its_1v22_EUCLID_RS_T2"
gen_tests_1v22_EUCLID_RS
RESITS_DIR="res_its_1v22_EUCLID_RS_T3"
gen_tests_1v22_EUCLID_RS

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

ITSTIMEDIR_LIST=""

TAG=""
# Set these below:
#TAG="RAYITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

function format_res_all_params()
{
  LINE=""

  Prob_str="SqPo"
  PRECLIST="WedNlFrayPa"
  VISLIST="Sim Str"
  ANGLIST="A0 A30 A67"
  RELIST="R200"
  NOELLIST="N4 N8 N16 N32 N64 N128 N256 N512"

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


            ## Set the first min.
            MIN_TOKEN=""
            MIN_NUM="-1.0"

            ITSTIMEDIR=${ITSTIMEDIR_LIST[0]}
            RESFILE="${ITSTIMEDIR}/${Prob_str}${VIS}${RE}${PREC}${ANG}${NOEL}NP1R0"

            if [ -f $RESFILE ]; then

              MIN_TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
              MIN_NUM=${MIN_TOKEN%(*}

              for ITSTIMEDIR in "${ITSTIMEDIR_LIST[@]}"
              do
                RESFILE="${ITSTIMEDIR}/${Prob_str}${VIS}${RE}${PREC}${ANG}${NOEL}NP1R0"

                CURRENT_TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
                CURRENT_NUM=${CURRENT_TOKEN%(*}

                # Now check to see which is bigger...
                bools=$(echo "$CURRENT_NUM < $MIN_NUM" | bc)

                if [ "$bools" == "1" ]; then
                  MIN_TOKEN=$CURRENT_TOKEN
                  MIN_NUM=$CURRENT_NUM
                fi

              done

              LINE="$LINE $MIN_TOKEN"
            else
              LINE="$LINE x"
            fi

          done
          echo $LINE
        done
      done
    done
  done
}

TAG="RAYITS"

cd $CURRENT_DIR
echo -e "\n"
ITSTIMEDIR_BASE="res_its_1v22_GS_RS"
echo -e "Formatting $ITSTIMEDIR_BASE"
ITSTIMEDIR_LIST=("${ITSTIMEDIR_BASE}_T1" "${ITSTIMEDIR_BASE}_T2" "${ITSTIMEDIR_BASE}_T3" )
format_res_all_params

cd $CURRENT_DIR
echo -e "\n"
ITSTIMEDIR_BASE="res_its_2v22_GS_RS"
echo -e "Formatting $ITSTIMEDIR_BASE"
ITSTIMEDIR_LIST=("${ITSTIMEDIR_BASE}_T1" "${ITSTIMEDIR_BASE}_T2" "${ITSTIMEDIR_BASE}_T3" )
format_res_all_params

cd $CURRENT_DIR
echo -e "\n"
ITSTIMEDIR_BASE="res_its_1v22_EUCLID_RS"
echo -e "Formatting $ITSTIMEDIR_BASE"
ITSTIMEDIR_LIST=("${ITSTIMEDIR_BASE}_T1" "${ITSTIMEDIR_BASE}_T2" "${ITSTIMEDIR_BASE}_T3" )
format_res_all_params
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
  DIR_LIST=("res_its_1v22_GS_RS_T1" \
            "res_its_1v22_GS_RS_T2" \
            "res_its_1v22_GS_RS_T3" \
            "res_its_2v22_GS_RS_T1" \
            "res_its_2v22_GS_RS_T2" \
            "res_its_2v22_GS_RS_T3" \
            "res_its_1v22_EUCLID_RS_T1" \
            "res_its_1v22_EUCLID_RS_T2" \
            "res_its_1v22_EUCLID_RS_T3" \
            "qsub_output_$TESTLIST_FILEBASE")

  for DIR in "${DIR_LIST[@]}"
  do
    mkdir -p $SCRATCH_TEST_DIR/$DIR
  done

  ## Tell the user that everything is done.
  echo -e "\n"
  echo "I have moved the files:"
  echo "  $PROGRAM"
  echo "  $QSUBFILE"
  echo "  $TEST_LIST"
  echo "  $FORMAT_RESULTS_SCRIPT"
  echo -e "\nand created directories:"
  printf -- '  %s\n' "${DIR_LIST[@]}"
fi




