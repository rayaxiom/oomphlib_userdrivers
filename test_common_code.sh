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
Famg_ITER="--f_amg_iter 1"
Famg_SMITER="--f_amg_smiter 2"
Famg_SMOOTHER="--f_amg_com_smoo 9" #Euclid
Famg_DAMP="--f_amg_damp -1"
Famg_STRN_SIM="--f_amg_str 0.25"
Famg_STRN_STR="--f_amg_str 0.668"
Famg_COARSE="--f_amg_coarse 1" #RS

Famg_sim="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
Famg_str="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_STR $Famg_COARSE"
Flu="--f_solver 0"

Plu="--p_solver 0"
Pamg="--p_solver 1"

Wlu="--w_solver 0"

NSlu="--ns_solver 0"
NSlsc="--ns_solver 1"

LPREC0_LU_LU="$Wlu $NSlu"

LPREC1_LU_LSClulu="$Wlu $NSlsc $Plu $Flu"

LPREC2_LU_LSCamglu="$Wlu $NSlsc $Pamg $Flu"

LPREC3_LU_LSCluamgsim="$Wlu $NSlsc $Plu $Famg_sim"
LPREC3_LU_LSCluamgstr="$Wlu $NSlsc $Plu $Famg_str"

LPREC4_LU_LSCamgamgsim="$Wlu $NSlsc $Pamg $Famg_sim"
LPREC4_LU_LSCamgamgstr="$Wlu $NSlsc $Pamg $Famg_str"

###### Vanilla LSC
LSConly="--lsc_only"
VPREC0_lulu="$LSConly $Plu $Flu"
VPREC1_amglu="$LSConly $Pamg $Flu"

VPREC2_luamg_sim="$LSConly $Plu $Famg_sim"
VPREC2_luamg_str="$LSConly $Plu $Famg_str"

VPREC3_amgamg_sim="$LSConly $Pamg $Famg_sim"
VPREC3_amgamg_str="$LSConly $Pamg $Famg_str"
}




