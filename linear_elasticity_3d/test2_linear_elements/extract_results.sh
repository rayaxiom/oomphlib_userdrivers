#!/bin/bash
set -u

PROC_NUM="0"

# Please only put: poisson3d_ndof200000_np1
FILEBASE="poisson3d_ndof50000_np"
NPROC="1"
FILE1_CLJP_05="${FILEBASE}${NPROC}.qsub.o*.1" # Str 0.5
#FILE2_mRS_05="$FILEBASE${NPROC}.qsub.o*.2" # Str 0.5
#FILE3_Falgout_05="$FILEBASE${NPROC}.qsub.o*.3" # Str 0.5
#FILE4_PMIS_05="$FILEBASE${NPROC}.qsub.o*.4" # Str 0.5
#FILE5_HMIS_05="$FILEBASE${NPROC}.qsub.o*.5" # Str 0.5
#FILE6_CLJP_07="$FILEBASE${NPROC}.qsub.o*.6" # Str 0.7
#FILE7_mRS_07="$FILEBASE${NPROC}.qsub.o*.7" # Str 0.7
#FILE8_Falgout_07="$FILEBASE${NPROC}.qsub.o*.8" # Str 0.7
#FILE9_PMIS_07="$FILEBASE${NPROC}.qsub.o*.9" # Str 0.7
#FILE10_HMIS_07="$FILEBASE${NPROC}.qsub.o*.10" # Str 0.7


GenJac="Time to generate Jacobian"
GenBoomerAMG="Setting up BoomerAMG"
PrecSetup="Time for preconditioner setup"
GenTrilinosMat="Time to generate Trilinos matrix"
LinSolIter="Linear solver iterations"
TrilinosSolTime="Time for trilinos solve itself"
CompleteTrilinosSolTime="Time for complete trilinos solve"
LinSolTime="Time for linear solver"
TotalLinSolTime="Total time for linear solver"
TotalNewtonSolTime="Total time for Newton solver"
TimeOutsideLinSol="Time outside linear solver"

GrepString="$GenJac\|$GenBoomerAMG\|$PrecSetup\|$GenTrilinosMat\|$LinSolIter\|$TrilinosSolTime\|$CompleteTrilinosSolTime\|$LinSolTime\|$TotalLinSolTime\|$TotalNewtonSolTime\|$TimeOutsideLinSol"

## Get GenJac time for File1, nproc = 1, 2, 4, 8



TAG=""
FILEBASE=""
AMG=""
function output_per_file_vary_np()
{
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    RESFILE="${FILEBASE}${NPROC}.qsub.o*.${AMG}"
    TOKEN=$(grep "Processor $PROC_NUM:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    LINE="$LINE $TOKEN"
  done
  echo $LINE
}


# Note: this is specially crafted for the 3D linear elasticity problem.
# and used to grep boomer AMG. The problem is, BoomerAMG is outputted three
# times, as a result, the TOKEN is 1.1 1.2 1.3 (for example).
# But we want them going down:
# 1.1
# 1.2
# 1.3
#
# Furthermore, we want to want to loop through four times (for each nproc),
# nproc is 1, 2, 4 and 8.
#
# Let us say that the token for file nproc 1, 2, 4, 8 gives
# TOKEN1 = 1.1 1.2 1.3
# TOKEN2 = 2.1 2.2 2.3
# TOKEN4 = 4.1 4.2 4.3
# TOKEN8 = 8.1 8.2 8.3
#
# So ultimately, we want:
# 1.1 2.1 4.1 8.1
# 1.2 2.2 4.2 8.2
# 1.3 2.3 4.3 8.3
#
# So we will store the tokens per file. After this, we will form a line
# for the each token.
# 
# i.e.
# LINE= TOKEN1[0], TOKEN2[0] TOKEN4[0] TOKEN8[0]
# LINE= TOKEN1[1], TOKEN2[1] TOKEN4[1] TOKEN8[1]
# LINE= TOKEN1[2], TOKEN2[2] TOKEN4[2] TOKEN8[2]
#
function output_per_file_vary_np_linelas_boomeramg()
{
  declare -a TOKEN1
  declare -a TOKEN2
  declare -a TOKEN4
  declare -a TOKEN8

  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.o*.${AMG}"

    # Use awk and grep magic
    TOKEN=$(grep "Processor $PROC_NUM:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      TOKEN1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      TOKEN2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      TOKEN4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      TOKEN8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done

  # Now loop through the number of elements in the array.
  INDEXLIST="0 1 2"
  for INDEX in $INDEXLIST
  do
    LINE=""
    LINE="${TOKEN1[$INDEX]} ${TOKEN2[$INDEX]} ${TOKEN4[$INDEX]} ${TOKEN8[$INDEX]}"
    echo $LINE
  done
}

function output_per_file_vary_np_m1()
{
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    RESFILE="${FILEBASE}${NPROC}.qsub.o*.${AMG}"
    TOKEN=$(grep "Processor $PROC_NUM:" $RESFILE | grep "$TAG" | awk '{print $(NF-1)}')
    LINE="$LINE $TOKEN"
  done
  echo $LINE
}


FILEBASE="linq_ndof200000_np"
function damp05()
{
echo "Doing 0.5"
AMG="1"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np_linelas_boomeramg

TAG="$PrecSetup"
output_per_file_vary_np

TAG="$GenTrilinosMat"
output_per_file_vary_np

TAG="$LinSolIter"
output_per_file_vary_np

TAG="$TrilinosSolTime"
output_per_file_vary_np

TAG="$CompleteTrilinosSolTime"
output_per_file_vary_np

TAG="$LinSolTime"
output_per_file_vary_np

TAG="$TotalLinSolTime"
output_per_file_vary_np

TAG="$TotalNewtonSolTime"
output_per_file_vary_np

TAG="$TimeOutsideLinSol"
output_per_file_vary_np_m1

echo ""

###############################
AMG="2"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np_linelas_boomeramg

TAG="$PrecSetup"
output_per_file_vary_np

TAG="$GenTrilinosMat"
output_per_file_vary_np

TAG="$LinSolIter"
output_per_file_vary_np

TAG="$TrilinosSolTime"
output_per_file_vary_np

AG="$CompleteTrilinosSolTime"
output_per_file_vary_np

TAG="$LinSolTime"
output_per_file_vary_np

TAG="$TotalLinSolTime"
output_per_file_vary_np

TAG="$TotalNewtonSolTime"
output_per_file_vary_np

TAG="$TimeOutsideLinSol"
output_per_file_vary_np_m1

echo ""
###############################
AMG="3"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np_linelas_boomeramg

TAG="$PrecSetup"
output_per_file_vary_np

TAG="$GenTrilinosMat"
output_per_file_vary_np

TAG="$LinSolIter"
output_per_file_vary_np

TAG="$TrilinosSolTime"
output_per_file_vary_np

TAG="$CompleteTrilinosSolTime"
output_per_file_vary_np

TAG="$LinSolTime"
output_per_file_vary_np

TAG="$TotalLinSolTime"
output_per_file_vary_np

TAG="$TotalNewtonSolTime"
output_per_file_vary_np

TAG="$TimeOutsideLinSol"
output_per_file_vary_np_m1
echo ""
###############################
AMG="4"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np_linelas_boomeramg

TAG="$PrecSetup"
output_per_file_vary_np

TAG="$GenTrilinosMat"
output_per_file_vary_np

AG="$LinSolIter"
output_per_file_vary_np

TAG="$TrilinosSolTime"
output_per_file_vary_np

TAG="$CompleteTrilinosSolTime"
output_per_file_vary_np

TAG="$LinSolTime"
output_per_file_vary_np

TAG="$TotalLinSolTime"
output_per_file_vary_np

TAG="$TotalNewtonSolTime"
output_per_file_vary_np

TAG="$TimeOutsideLinSol"
output_per_file_vary_np_m1
echo ""
###############################
AMG="5"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np_linelas_boomeramg

TAG="$PrecSetup"
output_per_file_vary_np

TAG="$GenTrilinosMat"
output_per_file_vary_np

TAG="$LinSolIter"
output_per_file_vary_np

TAG="$TrilinosSolTime"
output_per_file_vary_np

TAG="$CompleteTrilinosSolTime"
output_per_file_vary_np

TAG="$LinSolTime"
output_per_file_vary_np

TAG="$TotalLinSolTime"
output_per_file_vary_np

TAG="$TotalNewtonSolTime"
output_per_file_vary_np

TAG="$TimeOutsideLinSol"
output_per_file_vary_np_m1
}

function damp07()
{
echo "Doing 0.7"
AMG="6"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np_linelas_boomeramg

TAG="$PrecSetup"
output_per_file_vary_np

TAG="$GenTrilinosMat"
output_per_file_vary_np

TAG="$LinSolIter"
output_per_file_vary_np

TAG="$TrilinosSolTime"
output_per_file_vary_np

TAG="$CompleteTrilinosSolTime"
output_per_file_vary_np

TAG="$LinSolTime"
output_per_file_vary_np

TAG="$TotalLinSolTime"
output_per_file_vary_np

TAG="$TotalNewtonSolTime"
output_per_file_vary_np

TAG="$TimeOutsideLinSol"
output_per_file_vary_np_m1

echo ""

###############################
AMG="7"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np_linelas_boomeramg

TAG="$PrecSetup"
output_per_file_vary_np

TAG="$GenTrilinosMat"
output_per_file_vary_np

TAG="$LinSolIter"
output_per_file_vary_np

TAG="$TrilinosSolTime"
output_per_file_vary_np

TAG="$CompleteTrilinosSolTime"
output_per_file_vary_np

TAG="$LinSolTime"
output_per_file_vary_np

TAG="$TotalLinSolTime"
output_per_file_vary_np

TAG="$TotalNewtonSolTime"
output_per_file_vary_np

TAG="$TimeOutsideLinSol"
output_per_file_vary_np_m1

echo ""
###############################
AMG="8"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np_linelas_boomeramg

TAG="$PrecSetup"
output_per_file_vary_np

TAG="$GenTrilinosMat"
output_per_file_vary_np

TAG="$LinSolIter"
output_per_file_vary_np

TAG="$TrilinosSolTime"
output_per_file_vary_np

TAG="$CompleteTrilinosSolTime"
output_per_file_vary_np

TAG="$LinSolTime"
output_per_file_vary_np

TAG="$TotalLinSolTime"
output_per_file_vary_np

TAG="$TotalNewtonSolTime"
output_per_file_vary_np

TAG="$TimeOutsideLinSol"
output_per_file_vary_np_m1
echo ""
###############################
AMG="9"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np_linelas_boomeramg

TAG="$PrecSetup"
output_per_file_vary_np

TAG="$GenTrilinosMat"
output_per_file_vary_np

TAG="$LinSolIter"
output_per_file_vary_np

TAG="$TrilinosSolTime"
output_per_file_vary_np

TAG="$CompleteTrilinosSolTime"
output_per_file_vary_np

TAG="$LinSolTime"
output_per_file_vary_np

TAG="$TotalLinSolTime"
output_per_file_vary_np

TAG="$TotalNewtonSolTime"
output_per_file_vary_np

TAG="$TimeOutsideLinSol"
output_per_file_vary_np_m1
echo ""
###############################
AMG="10"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np_linelas_boomeramg

TAG="$PrecSetup"
output_per_file_vary_np

TAG="$GenTrilinosMat"
output_per_file_vary_np

TAG="$LinSolIter"
output_per_file_vary_np

TAG="$TrilinosSolTime"
output_per_file_vary_np

TAG="$CompleteTrilinosSolTime"
output_per_file_vary_np

TAG="$LinSolTime"
output_per_file_vary_np

TAG="$TotalLinSolTime"
output_per_file_vary_np

TAG="$TotalNewtonSolTime"
output_per_file_vary_np

TAG="$TimeOutsideLinSol"
output_per_file_vary_np_m1
}

damp05
echo ""
echo ""
damp07





