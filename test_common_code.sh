#!/bin/bash

function setup_initial()
{
# Create the THISFILEBASE, this is the folder where the testing will be done.
THISFILE=$0 # This contains "./", which we do not want.
THISFILE=${THISFILE:2} # Gets rid of "./"
THISFILEBASE=${THISFILE%%.*} # Get rid of the extension (in this case, ".sh")

## The test folder is the same as the file base.
TEST_DIR=$THISFILEBASE

# Create the new folder (remove old one)
touch $TEST_DIR
rm -rf $TEST_DIR
mkdir $TEST_DIR


# Get the current directory and the oomph-base
PROGRAM_DIR=`pwd`
OOMPHROOT_DIR=$(make -s --no-print-directory print-top_builddir)

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
}

function load_prec_param()
{
## Preconditioner parameters.
Famg_BASE="--f_solver 96"
Famg_STRN_SIM="--f_amg_str 0.25"
Famg_STRN_STR="--f_amg_str 0.668"

# The cycle
Famg_ITER_1="--f_amg_iter 1"
Famg_ITER_2="--f_amg_iter 2"
Famg_SMITER_2="--f_amg_smiter 2"

# Smoother
Famg_SMOOTHER_GS="--f_amg_sim_smoo 1" # GS
Famg_SMOOTHER_Eu="--f_amg_com_smoo 9" #Euclid
Famg_DAMP="--f_amg_damp -1"

# Coarsener
Famg_COARSE="--f_amg_coarse 1" #RS

#### 1V22 GS RS
Famg_C1v22_GS_RS_sim="$Famg_BASE $Famg_ITER_1 $Famg_SMITER_2 $Famg_SMOOTHER_GS $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
Famg_C1v22_GS_RS_str="$Famg_BASE $Famg_ITER_1 $Famg_SMITER_2 $Famg_SMOOTHER_GS $Famg_DAMP $Famg_STRN_STR $Famg_COARSE"

#### 2V22 GS RS 
Famg_C2v22_GS_RS_sim="$Famg_BASE $Famg_ITER_2 $Famg_SMITER_2 $Famg_SMOOTHER_GS $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
Famg_C2v22_GS_RS_str="$Famg_BASE $Famg_ITER_2 $Famg_SMITER_2 $Famg_SMOOTHER_GS $Famg_DAMP $Famg_STRN_STR $Famg_COARSE"

#### 1V22 Euclid RS
Famg_C1v22_Eu_RS_sim="$Famg_BASE $Famg_ITER_1 $Famg_SMITER_2 $Famg_SMOOTHER_Eu $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
Famg_C1v22_Eu_RS_str="$Famg_BASE $Famg_ITER_1 $Famg_SMITER_2 $Famg_SMOOTHER_Eu $Famg_DAMP $Famg_STRN_STR $Famg_COARSE"




#Famg_sim="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
#Famg_str="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_STR $Famg_COARSE"



Flu="--f_solver 0"

Plu="--p_solver 0"
Pamg="--p_solver 1"

Wlu="--w_solver 0"

NSlu="--ns_solver 0"
NSlsc="--ns_solver 1"

## All exact
LPREC0_LU_LU="$Wlu $NSlu"

## Exact LSC
LPREC1_LU_LSCLuLu="$Wlu $NSlsc $Plu $Flu"

## AMG for pressure block.
LPREC2_LU_LSCAmgLu="$Wlu $NSlsc $Pamg $Flu"

### AMG for LSC block.
LPREC3_LU_LSCLu_C1v22_GS_RS_sim="$Wlu $NSlsc $Plu $Famg_C1v22_GS_RS_sim"
LPREC3_LU_LSCLu_C1v22_GS_RS_str="$Wlu $NSlsc $Plu $Famg_C1v22_GS_RS_str"

LPREC3_LU_LSCLu_C2v22_GS_RS_sim="$Wlu $NSlsc $Plu $Famg_C2v22_GS_RS_sim"
LPREC3_LU_LSCLu_C2v22_GS_RS_str="$Wlu $NSlsc $Plu $Famg_C2v22_GS_RS_str"

LPREC3_LU_LSCLu_C1v22_Eu_RS_sim="$Wlu $NSlsc $Plu $Famg_C1v22_Eu_RS_sim"
LPREC3_LU_LSCLu_C1v22_Eu_RS_str="$Wlu $NSlsc $Plu $Famg_C1v22_Eu_RS_str"

### Full AMG
LPREC4_LU_LSCAmg_C1v22_GS_RS_sim="$Wlu $NSlsc $Pamg $Famg_C1v22_GS_RS_sim"
LPREC4_LU_LSCAmg_C1v22_GS_RS_str="$Wlu $NSlsc $Pamg $Famg_C1v22_GS_RS_str"

LPREC4_LU_LSCAmg_C2v22_GS_RS_sim="$Wlu $NSlsc $Pamg $Famg_C2v22_GS_RS_sim"
LPREC4_LU_LSCAmg_C2v22_GS_RS_str="$Wlu $NSlsc $Pamg $Famg_C2v22_GS_RS_str"

LPREC4_LU_LSCAmg_C1v22_Eu_RS_sim="$Wlu $NSlsc $Pamg $Famg_C1v22_Eu_RS_sim"
LPREC4_LU_LSCAmg_C1v22_Eu_RS_str="$Wlu $NSlsc $Pamg $Famg_C1v22_Eu_RS_str"


###### Vanilla LSC
LSConly="--lsc_only"

## Exact for both P and F solve
VPREC0_Lu_Lu="$LSConly $Plu $Flu"

## P - AMG, F exact
VPREC1_Amg_Lu="$LSConly $Pamg $Flu"

## P - Exact, F - AMG
VPREC2_Lu_C1v22_GS_RS_sim="$LSConly $Plu $Famg_C1v22_GS_RS_sim"
VPREC2_Lu_C1v22_GS_RS_str="$LSConly $Plu $Famg_C1v22_GS_RS_str"

VPREC2_Lu_C2v22_GS_RS_sim="$LSConly $Plu $Famg_C2v22_GS_RS_sim"
VPREC2_Lu_C2v22_GS_RS_str="$LSConly $Plu $Famg_C2v22_GS_RS_str"

VPREC2_Lu_C1v22_Eu_RS_sim="$LSConly $Plu $Famg_C1v22_Eu_RS_sim"
VPREC2_Lu_C1v22_Eu_RS_str="$LSConly $Plu $Famg_C1v22_Eu_RS_str"

## P - Amg, F - AMG
VPREC3_Amg_C1v22_GS_RS_sim="$LSConly $Pamg $Famg_C1v22_GS_RS_sim"
VPREC3_Amg_C1v22_GS_RS_str="$LSConly $Pamg $Famg_C1v22_GS_RS_str"

VPREC3_Amg_C2v22_GS_RS_sim="$LSConly $Pamg $Famg_C2v22_GS_RS_sim"
VPREC3_Amg_C2v22_GS_RS_str="$LSConly $Pamg $Famg_C2v22_GS_RS_str"

VPREC3_Amg_C1v22_Eu_RS_sim="$LSConly $Pamg $Famg_C1v22_Eu_RS_sim"
VPREC3_Amg_C1v22_Eu_RS_str="$LSConly $Pamg $Famg_C1v22_Eu_RS_str"
}

function load_LPREC_C1v22_GS_RS_case()
{
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
      PREC_PARAM="$LPREC3_LU_LSCLu_C1v22_GS_RS_sim"
    else
      PREC_PARAM="$LPREC3_LU_LSCLu_C1v22_GS_RS_str"
    fi
    ;;
  4)
    if [ "$VIS" -eq "0" ]; then
      PREC_PARAM="$LPREC4_LU_LSCAmg_C1v22_GS_RS_sim"
    else
      PREC_PARAM="$LPREC4_LU_LSCAmg_C1v22_GS_RS_str"
    fi
    ;;
esac
}

function load_LPREC_C2v22_GS_RS_case()
{
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
      PREC_PARAM="$LPREC3_LU_LSCLu_C2v22_GS_RS_sim"
    else
      PREC_PARAM="$LPREC3_LU_LSCLu_C2v22_GS_RS_str"
    fi
    ;;
  4)
    if [ "$VIS" -eq "0" ]; then
      PREC_PARAM="$LPREC4_LU_LSCAmg_C2v22_GS_RS_sim"
    else
      PREC_PARAM="$LPREC4_LU_LSCAmg_C2v22_GS_RS_str"
    fi
    ;;
esac
}

function load_LPREC_C1v22_Eu_RS_case()
{
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
      PREC_PARAM="$LPREC4_LU_LSCAmg_C1v22_Eu_RS_sim"
    else
      PREC_PARAM="$LPREC4_LU_LSCAmg_C1v22_Eu_RS_str"
    fi
    ;;
esac
}


