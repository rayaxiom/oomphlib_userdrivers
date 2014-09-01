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


grep "Processor $PROC_NUM:" $FILE1_CLJP_05 | grep "$GrepString"
echo ""

#grep "Processor $PROC_NUM:" $FILE2 | grep "$GrepString"
#echo ""
#
#grep "Processor $PROC_NUM:" $FILE3 | grep "$GrepString"
#echo ""
#
#grep "Processor $PROC_NUM:" $FILE4 | grep "$GrepString"
#echo ""
#
#grep "Processor $PROC_NUM:" $FILE5 | grep "$GrepString"
#echo ""
#
#grep "Processor $PROC_NUM:" $FILE6 | grep "$GrepString"
#echo ""
#
#grep "Processor $PROC_NUM:" $FILE7 | grep "$GrepString"
#echo ""
#
#grep "Processor $PROC_NUM:" $FILE8 | grep "$GrepString"
#echo ""
#
#grep "Processor $PROC_NUM:" $FILE9 | grep "$GrepString"
#echo ""
#
#grep "Processor $PROC_NUM:" $FILE10 | grep "$GrepString"
#echo ""





