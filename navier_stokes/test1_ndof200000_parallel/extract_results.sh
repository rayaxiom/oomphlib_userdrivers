#!/bin/bash
set -u

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

declare -a GenJacArray1
declare -a GenBoomerAMGArray1
declare -a PrecSetupArray1
declare -a GenTrilinosMatArray1
declare -a LinSolIterArray1
declare -a TrilinosSolTimeArray1
declare -a CompleteTrilinosSolTimeArray1
declare -a LinSolTimeArray1
declare -a TotalLinSolTimeArray1
declare -a TotalNewtonSolTimeArray1
declare -a TimeOutsideLinSolArray1

declare -a GenJacArray2
declare -a GenBoomerAMGArray2
declare -a PrecSetupArray2
declare -a GenTrilinosMatArray2
declare -a LinSolIterArray2
declare -a TrilinosSolTimeArray2
declare -a CompleteTrilinosSolTimeArray2
declare -a LinSolTimeArray2
declare -a TotalLinSolTimeArray2
declare -a TotalNewtonSolTimeArray2
declare -a TimeOutsideLinSolArray2

declare -a GenJacArray4
declare -a GenBoomerAMGArray4
declare -a PrecSetupArray4
declare -a GenTrilinosMatArray4
declare -a LinSolIterArray4
declare -a TrilinosSolTimeArray4
declare -a CompleteTrilinosSolTimeArray4
declare -a LinSolTimeArray4
declare -a TotalLinSolTimeArray4
declare -a TotalNewtonSolTimeArray4
declare -a TimeOutsideLinSolArray4

declare -a GenJacArray8
declare -a GenBoomerAMGArray8
declare -a PrecSetupArray8
declare -a GenTrilinosMatArray8
declare -a LinSolIterArray8
declare -a TrilinosSolTimeArray8
declare -a CompleteTrilinosSolTimeArray8
declare -a LinSolTimeArray8
declare -a TotalLinSolTimeArray8
declare -a TotalNewtonSolTimeArray8
declare -a TimeOutsideLinSolArray8

TAG=""
FILEBASE=""
AMG=""
function output_per_file_vary_np()
{
  NPROCLIST="1"
  LINE=""
  for NPROC in $NPROCLIST
  do
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"
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
function store_genjac()
{
  unset GenJacArray1
  unset GenJacArray2
  unset GenJacArray4
  unset GenJacArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      GenJacArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      GenJacArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      GenJacArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      GenJacArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}

function store_genboomeramg()
{
  unset GenBoomerAMGArray1
  unset GenBoomerAMGArray2
  unset GenBoomerAMGArray4
  unset GenBoomerAMGArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      GenBoomerAMGArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      GenBoomerAMGArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      GenBoomerAMGArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      GenBoomerAMGArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}



function store_precsetup()
{
  unset PrecSetupArray1
  unset PrecSetupArray2
  unset PrecSetupArray4
  unset PrecSetupArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      PrecSetupArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      PrecSetupArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      PrecSetupArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      PrecSetupArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}

function store_gentrilinosmat()
{
  unset GenTrilinosMatArray1
  unset GenTrilinosMatArray2
  unset GenTrilinosMatArray4
  unset GenTrilinosMatArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      GenTrilinosMatArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      GenTrilinosMatArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      GenTrilinosMatArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      GenTrilinosMatArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}

function store_linsoliter()
{
  unset LinSolIterArray1
  unset LinSolIterArray2
  unset LinSolIterArray4
  unset LinSolIterArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      LinSolIterArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      LinSolIterArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      LinSolIterArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      LinSolIterArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}

function store_trilinossoltime()
{
  unset TrilinosSolTimeArray1
  unset TrilinosSolTimeArray2
  unset TrilinosSolTimeArray4
  unset TrilinosSolTimeArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      TrilinosSolTimeArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      TrilinosSolTimeArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      TrilinosSolTimeArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      TrilinosSolTimeArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}

function store_completetrilinossoltime()
{
  unset CompleteTrilinosSolTimeArray1
  unset CompleteTrilinosSolTimeArray2
  unset CompleteTrilinosSolTimeArray4
  unset CompleteTrilinosSolTimeArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      CompleteTrilinosSolTimeArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      CompleteTrilinosSolTimeArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      CompleteTrilinosSolTimeArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      CompleteTrilinosSolTimeArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}

function store_linsoltime()
{
  unset LinSolTimeArray1
  unset LinSolTimeArray2
  unset LinSolTimeArray4
  unset LinSolTimeArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      LinSolTimeArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      LinSolTimeArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      LinSolTimeArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      LinSolTimeArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}

function store_totallinsoltime()
{
  unset TotalLinSolTimeArray1
  unset TotalLinSolTimeArray2
  unset TotalLinSolTimeArray4
  unset TotalLinSolTimeArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      TotalLinSolTimeArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      TotalLinSolTimeArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      TotalLinSolTimeArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      TotalLinSolTimeArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}

function store_totalnewtonsoltime()
{
  unset TotalNewtonSolTimeArray1
  unset TotalNewtonSolTimeArray2
  unset TotalNewtonSolTimeArray4
  unset TotalNewtonSolTimeArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $NF}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      TotalNewtonSolTimeArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      TotalNewtonSolTimeArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      TotalNewtonSolTimeArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      TotalNewtonSolTimeArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}

function store_timeoutsidelinsol()
{
  unset TimeOutsideLinSolArray1
  unset TimeOutsideLinSolArray2
  unset TimeOutsideLinSolArray4
  unset TimeOutsideLinSolArray8
  NPROCLIST="1 2 4 8"
  LINE=""
  for NPROC in $NPROCLIST
  do
    # Form the results file, which we grep from
    RESFILE="${FILEBASE}${NPROC}.qsub.${AMG}_proc0"

    # Use awk and grep magic
    TOKEN=$(grep "Processor 0:" $RESFILE | grep "$TAG" | awk '{print $(NF-1)}')
    
    # Split the token into the correct variable
    if [ "$NPROC" = "1" ]; then
      TimeOutsideLinSolArray1=(${TOKEN})
    elif [ "$NPROC" = "2" ]; then
      TimeOutsideLinSolArray2=(${TOKEN})
    elif [ "$NPROC" = "4" ]; then
      TimeOutsideLinSolArray4=(${TOKEN})
    elif [ "$NPROC" = "8" ]; then
      TimeOutsideLinSolArray8=(${TOKEN})
    else
      echo "Each end of trying to split stuff."
    fi
  done
}


function output_stored_results()
{
NTIMESTEP="0 1 2 3"
for TIMESTEP in $NTIMESTEP
do
echo "Time step: $TIMESTEP"
LINE="${GenJacArray1[$TIMESTEP]} ${GenJacArray2[$TIMESTEP]} ${GenJacArray4[$TIMESTEP]} ${GenJacArray8[$TIMESTEP]}"
echo $LINE

boom1=$((TIMESTEP*2))
boom2=$((boom1+1))
LINE="${GenBoomerAMGArray1[$boom1]} ${GenBoomerAMGArray2[$boom1]} ${GenBoomerAMGArray4[$boom1]} ${GenBoomerAMGArray8[$boom1]}"
echo $LINE

LINE="${GenBoomerAMGArray1[$boom2]} ${GenBoomerAMGArray2[$boom2]} ${GenBoomerAMGArray4[$boom2]} ${GenBoomerAMGArray8[$boom2]}"
echo $LINE


LINE="${PrecSetupArray1[$TIMESTEP]} ${PrecSetupArray2[$TIMESTEP]} ${PrecSetupArray4[$TIMESTEP]} ${PrecSetupArray8[$TIMESTEP]}"
echo $LINE

LINE="${GenTrilinosMatArray1[$TIMESTEP]} ${GenTrilinosMatArray2[$TIMESTEP]} ${GenTrilinosMatArray4[$TIMESTEP]} ${GenTrilinosMatArray8[$TIMESTEP]}"
echo $LINE

LINE="${LinSolIterArray1[$TIMESTEP]} ${LinSolIterArray2[$TIMESTEP]} ${LinSolIterArray4[$TIMESTEP]} ${LinSolIterArray8[$TIMESTEP]}"
echo $LINE

LINE="${TrilinosSolTimeArray1[$TIMESTEP]} ${TrilinosSolTimeArray2[$TIMESTEP]} ${TrilinosSolTimeArray4[$TIMESTEP]} ${TrilinosSolTimeArray8[$TIMESTEP]}"
echo $LINE

LINE="${CompleteTrilinosSolTimeArray1[$TIMESTEP]} ${CompleteTrilinosSolTimeArray2[$TIMESTEP]} ${CompleteTrilinosSolTimeArray4[$TIMESTEP]} ${CompleteTrilinosSolTimeArray8[$TIMESTEP]}"
echo $LINE

LINE="${LinSolTimeArray1[$TIMESTEP]} ${LinSolTimeArray2[$TIMESTEP]} ${LinSolTimeArray4[$TIMESTEP]} ${LinSolTimeArray8[$TIMESTEP]}"
echo $LINE

LINE="${TotalLinSolTimeArray1[$TIMESTEP]} ${TotalLinSolTimeArray2[$TIMESTEP]} ${TotalLinSolTimeArray4[$TIMESTEP]} ${TotalLinSolTimeArray8[$TIMESTEP]}"
echo $LINE

LINE="${TotalNewtonSolTimeArray1[$TIMESTEP]} ${TotalNewtonSolTimeArray2[$TIMESTEP]} ${TotalNewtonSolTimeArray4[$TIMESTEP]} ${TotalNewtonSolTimeArray8[$TIMESTEP]}"
echo $LINE

LINE="${TimeOutsideLinSolArray1[$TIMESTEP]} ${TimeOutsideLinSolArray2[$TIMESTEP]} ${TimeOutsideLinSolArray4[$TIMESTEP]} ${TimeOutsideLinSolArray8[$TIMESTEP]}"
echo $LINE
echo ""
done
}


FILEBASE="ns_ndof200000_np"
function damp05()
{
echo "Doing CLJP:"
AMG="1"
TAG="$GenJac"
store_genjac

TAG="$GenBoomerAMG"
store_genboomeramg

TAG="$PrecSetup"
store_precsetup

TAG="$GenTrilinosMat"
store_gentrilinosmat

TAG="$LinSolIter"
store_linsoliter

TAG="$TrilinosSolTime"
store_trilinossoltime

TAG="$CompleteTrilinosSolTime"
store_completetrilinossoltime

TAG="$LinSolTime"
store_linsoltime

TAG="$TotalLinSolTime"
store_totallinsoltime

TAG="$TotalNewtonSolTime"
store_totalnewtonsoltime

TAG="$TimeOutsideLinSol"
store_timeoutsidelinsol

output_stored_results
echo ""
echo ""

###############################
echo "Doing Modified RS"
AMG="2"

TAG="$GenJac"
store_genjac

TAG="$GenBoomerAMG"
store_genboomeramg

TAG="$PrecSetup"
store_precsetup

TAG="$GenTrilinosMat"
store_gentrilinosmat

TAG="$LinSolIter"
store_linsoliter

TAG="$TrilinosSolTime"
store_trilinossoltime

TAG="$CompleteTrilinosSolTime"
store_completetrilinossoltime

TAG="$LinSolTime"
store_linsoltime

TAG="$TotalLinSolTime"
store_totallinsoltime

TAG="$TotalNewtonSolTime"
store_totalnewtonsoltime

TAG="$TimeOutsideLinSol"
store_timeoutsidelinsol

output_stored_results
echo ""
echo ""
##############################
AMG="3"
echo "Doing Falgout"

TAG="$GenJac"
store_genjac

TAG="$GenBoomerAMG"
store_genboomeramg

TAG="$PrecSetup"
store_precsetup

TAG="$GenTrilinosMat"
store_gentrilinosmat

TAG="$LinSolIter"
store_linsoliter

TAG="$TrilinosSolTime"
store_trilinossoltime

TAG="$CompleteTrilinosSolTime"
store_completetrilinossoltime

TAG="$LinSolTime"
store_linsoltime

TAG="$TotalLinSolTime"
store_totallinsoltime

TAG="$TotalNewtonSolTime"
store_totalnewtonsoltime

TAG="$TimeOutsideLinSol"
store_timeoutsidelinsol

output_stored_results
echo ""
echo ""
##############################
AMG="4"
echo "Doing PMIS"

TAG="$GenJac"
store_genjac

TAG="$GenBoomerAMG"
store_genboomeramg

TAG="$PrecSetup"
store_precsetup

TAG="$GenTrilinosMat"
store_gentrilinosmat

TAG="$LinSolIter"
store_linsoliter

TAG="$TrilinosSolTime"
store_trilinossoltime

TAG="$CompleteTrilinosSolTime"
store_completetrilinossoltime

TAG="$LinSolTime"
store_linsoltime

TAG="$TotalLinSolTime"
store_totallinsoltime

TAG="$TotalNewtonSolTime"
store_totalnewtonsoltime

TAG="$TimeOutsideLinSol"
store_timeoutsidelinsol

output_stored_results
echo ""
echo
###############################
AMG="5"
echo "Doing HMIS"
TAG="$GenJac"
store_genjac

TAG="$GenBoomerAMG"
store_genboomeramg

TAG="$PrecSetup"
store_precsetup

TAG="$GenTrilinosMat"
store_gentrilinosmat

TAG="$LinSolIter"
store_linsoliter

TAG="$TrilinosSolTime"
store_trilinossoltime

TAG="$CompleteTrilinosSolTime"
store_completetrilinossoltime

TAG="$LinSolTime"
store_linsoltime

TAG="$TotalLinSolTime"
store_totallinsoltime

TAG="$TotalNewtonSolTime"
store_totalnewtonsoltime

TAG="$TimeOutsideLinSol"
store_timeoutsidelinsol

output_stored_results
echo ""
echo ""
}

damp05
#echo ""
#echo ""
#damp07





