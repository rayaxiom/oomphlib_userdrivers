#!/bin/bash

# Create the THISFILEBASE, this is the folder where the testing will be done.
THISFILE=$0 # This contains "./", which we do not want.
THISFILE=${THISFILE:2} # Gets rid of "./"
THISFILEBASE=${THISFILE%%.*} # Get rid of the extension (in this case, ".sh")

TLIST_FILE_LIST="${THISFILEBASE}.list"
rm -rf $TLIST_FILE_LIST


PROGRAM="annular_wedge_threed"

#SOLN_DIR="RESLT_SOLN_tmp"
TIME_DIR="RESLT_ITSTIME"
rm -rf $TIME_DIR
mkdir $TIME_DIR


RUN_COMMAND="mpirun -np 1"


gen_exact_tests()
{
DIST_PROB="--dist_prob"
PROB_ID="--prob_id 20"
#DOC_SOLN="--doc_soln $SOLN_DIR"
DOC_SOLN=""

#VISC="--visc 1" ################################################
VISCLIST="0 1"

#REY="--rey 100" ################################################
REYLIST="0 100 200 500 1000"


MAX_ITER="--max_solver_iter 100" # TO CHECK
ITSTIMEDIR="--itstimedir $TIME_DIR" # TO CHECK
#ITSTIMEDIR="" # TO CHECK
SOLVER_TYPE="--solver_type 2" # TO CHECK
DT="--dt -1.0"
TIME_START="--time_start 0.0"
TIME_END="--time_end 1.0"


###############################################################################
PRINT_HYPRE="--print_hypre"


W_SOLVER_E="--w_solver 0"

NS_SOLVER_E="--ns_solver 0"
NS_SOLVER_LSC="--ns_solver 1"

P_SOLVER_E="--p_solver 0"
P_SOLVER_A="--p_solver 1"


FAMG_ITER="--f_amg_iter 2"
FAMG_SMITER="--f_amg_smiter 2"
FAMG_SSMOOTHER="--f_amg_sim_smoo 1" #GS
#FAMG_CSMOOTHER="--f_amg_com_smoo 9" # Euclid
FAMG_DAMP="--f_amg_damp -1"

FAMG_STRN_SIM="--f_amg_str 0.25"
FAMG_STRN_STR="--f_amg_str 0.668"

FAMG_COARSE="--f_amg_coarse 1" #RS - 1, Falgout - 6

F_SOLVER_E="--f_solver 0"
F_SOLVER_A_SIM="--f_solver 96 $FAMG_ITER $FAMG_SMITER $FAMG_SSMOOTHER $FAMG_CSMOOTHER $FAMG_DAMP $FAMG_STRN_SIM $FAMG_COARSE"
F_SOLVER_A_STR="--f_solver 96 $FAMG_ITER $FAMG_SMITER $FAMG_SSMOOTHER $FAMG_CSMOOTHER $FAMG_DAMP $FAMG_STRN_STR $FAMG_COARSE"

# USE THESE
PREC_EXACT="$PRINT_HYPRE $W_SOLVER_E $NS_SOLVER_E"
PREC_LSC_EXACT="$PRINT_HYPRE $W_SOLVER_E $NS_SOLVER_LSC $P_SOLVER_E $F_SOLVER_E"
PREC_LSC_AMG_SIM="$PRINT_HYPRE $W_SOLVER_E $NS_SOLVER_LSC $P_SOLVER_A $F_SOLVER_A_SIM"
PREC_LSC_AMG_STR="$PRINT_HYPRE $W_SOLVER_E $NS_SOLVER_LSC $P_SOLVER_A $F_SOLVER_A_STR"

PRECLIST="0 1"

# The precs are set according to the PRECLIST above.
PRECPARAM=""

###############################################################################
NOELLIST="4 6 8 10 13"

###############################################################################
###############################################################################
## The main for loops
for PREC in $PRECLIST
do
  for VISC in $VISCLIST
  do
    for REY in $REYLIST
    do
      for NOEL in $NOELLIST
      do

## Parameters set by NSPP:
NSPP="$DIST_PROB $PROB_ID $DOC_SOLN --visc $VISC --rey $REY $MAX_ITER $ITSTIMEDIR $SOLVER_TYPE $DT $TIME_START $TIME_END"


## Parameters get by LPH:
case "$PREC" in
  0)
    PRECPARAM="$PREC_EXACT"
    ;;
  1)
    PRECPARAM="$PREC_LSC_EXACT"
    ;;
  2)
    if [ "$VISC" -eq "0" ]; then
      PRECPARAM="$PREC_LSC_AMG_SIM"
    else
      PRECPARAM="$PREC_LSC_AMG_STR"
    fi
    ;;
esac

LPH="$PRECPARAM"

# Parameters set by the specific problem
PROB_SPECIFIC_PARAM="--noel $NOEL"


# Put them all together:
ALL_PARAM="$NSPP $LPH $PROB_SPECIFIC_PARAM"

##### ECHO
echo "$RUN_COMMAND ./$PROGRAM $ALL_PARAM" >> $TLIST_FILE_LIST
      done
    done
  done
done

} # gen_tests function


gen_small_amg_tests()
{
DIST_PROB="--dist_prob"
PROB_ID="--prob_id 20"
#DOC_SOLN="--doc_soln $SOLN_DIR"
DOC_SOLN=""

#VISC="--visc 1" ################################################
VISCLIST="0 1"

#REY="--rey 100" ################################################
REYLIST="0 100 200 500 1000"


MAX_ITER="--max_solver_iter 100" # TO CHECK
ITSTIMEDIR="--itstimedir $TIME_DIR" # TO CHECK
#ITSTIMEDIR="" # TO CHECK
SOLVER_TYPE="--solver_type 2" # TO CHECK
DT="--dt -1.0"
TIME_START="--time_start 0.0"
TIME_END="--time_end 1.0"


###############################################################################
PRINT_HYPRE="--print_hypre"


W_SOLVER_E="--w_solver 0"

NS_SOLVER_E="--ns_solver 0"
NS_SOLVER_LSC="--ns_solver 1"

P_SOLVER_E="--p_solver 0"
P_SOLVER_A="--p_solver 1"


FAMG_ITER="--f_amg_iter 2"
FAMG_SMITER="--f_amg_smiter 2"
FAMG_SSMOOTHER="--f_amg_sim_smoo 1" #GS
#FAMG_CSMOOTHER="--f_amg_com_smoo 9" # Euclid
FAMG_DAMP="--f_amg_damp -1"

FAMG_STRN_SIM="--f_amg_str 0.25"
FAMG_STRN_STR="--f_amg_str 0.668"

FAMG_COARSE="--f_amg_coarse 1" #RS - 1, Falgout - 6

F_SOLVER_E="--f_solver 0"
F_SOLVER_A_SIM="--f_solver 96 $FAMG_ITER $FAMG_SMITER $FAMG_SSMOOTHER $FAMG_CSMOOTHER $FAMG_DAMP $FAMG_STRN_SIM $FAMG_COARSE"
F_SOLVER_A_STR="--f_solver 96 $FAMG_ITER $FAMG_SMITER $FAMG_SSMOOTHER $FAMG_CSMOOTHER $FAMG_DAMP $FAMG_STRN_STR $FAMG_COARSE"

# USE THESE
PREC_EXACT="$PRINT_HYPRE $W_SOLVER_E $NS_SOLVER_E"
PREC_LSC_EXACT="$PRINT_HYPRE $W_SOLVER_E $NS_SOLVER_LSC $P_SOLVER_E $F_SOLVER_E"
PREC_LSC_AMG_SIM="$PRINT_HYPRE $W_SOLVER_E $NS_SOLVER_LSC $P_SOLVER_A $F_SOLVER_A_SIM"
PREC_LSC_AMG_STR="$PRINT_HYPRE $W_SOLVER_E $NS_SOLVER_LSC $P_SOLVER_A $F_SOLVER_A_STR"

PRECLIST="2"

# The precs are set according to the PRECLIST above.
PRECPARAM=""

###############################################################################
NOELLIST="4 6 8 10 13 17"

###############################################################################
###############################################################################
## The main for loops
for PREC in $PRECLIST
do
  for VISC in $VISCLIST
  do
    for REY in $REYLIST
    do
      for NOEL in $NOELLIST
      do

## Parameters set by NSPP:
NSPP="$DIST_PROB $PROB_ID $DOC_SOLN --visc $VISC --rey $REY $MAX_ITER $ITSTIMEDIR $SOLVER_TYPE $DT $TIME_START $TIME_END"


## Parameters get by LPH:
case "$PREC" in
  0)
    PRECPARAM="$PREC_EXACT"
    ;;
  1)
    PRECPARAM="$PREC_LSC_EXACT"
    ;;
  2)
    if [ "$VISC" -eq "0" ]; then
      PRECPARAM="$PREC_LSC_AMG_SIM"
    else
      PRECPARAM="$PREC_LSC_AMG_STR"
    fi
    ;;
esac

LPH="$PRECPARAM"

# Parameters set by the specific problem
PROB_SPECIFIC_PARAM="--noel $NOEL"


# Put them all together:
ALL_PARAM="$NSPP $LPH $PROB_SPECIFIC_PARAM"

##### ECHO
echo "$RUN_COMMAND ./$PROGRAM $ALL_PARAM" >> $TLIST_FILE_LIST
      done
    done
  done
done

} # gen_tests function

gen_exact_tests
gen_small_amg_tests



