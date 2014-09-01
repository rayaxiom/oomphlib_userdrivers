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


FILEBASE="poisson3d_ndof400000_np"
function damp05()
{
echo "Doing 0.5"
AMG="1"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np

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
output_per_file_vary_np

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
AMG="3"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np

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
output_per_file_vary_np

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
AMG="5"
TAG="$GenJac"
output_per_file_vary_np

TAG="$GenBoomerAMG"
output_per_file_vary_np

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
output_per_file_vary_np

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
output_per_file_vary_np

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
output_per_file_vary_np

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
output_per_file_vary_np

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
output_per_file_vary_np

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





