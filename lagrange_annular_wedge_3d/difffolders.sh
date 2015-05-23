#!/bin/bash

MATFOLDER1="matrixdata1"
MATFOLDER2="matrixdata2"
MATFOLDER3="matrixdata3"

CURRDIR=`pwd`

###########################################################################


#CNPR="C1NP3R0_"

ModBlocksOnlyLgrIbBefore="RMBO_lgr_ibbefore_"
ModBlocksOnlyLgrIbAfter="RMBO_lgr_ibafter_"

ModBlocksOnlyLscBlock="RMBO_lsc_b_"
ModBlocksOnlyLscDof="RMBO_lBO_lsc_b_sc_d_"

ModBlocksOnlyLscMasterIblock="RMBO_master_ib_"
ModBlocksOnlyLscMasterBlock="RMBO_master_b_"
ModBlocksOnlyLscMasterDof="RMBO_master_d_"

#########

AllFBlockLgrIbBefore="RAFB_lgr_ibbefore_"
AllFBlockLgrIbAfter="RAFB_lgr_ibafter_"

AllFBlockLscBlock="RAFB_lsc_b_"
AllFBlockLscDof="RAFB_lsc_d_"

AllFBlockLscMasterIblock="RAFB_master_ib_"
AllFBlockLscMasterBlock="RAFB_master_b_"
AllFBlockLscMasterDof="RAFB_master_d_"

#ls -l ${ModBlocksOnlyLgrIbBefore}${CNPR}*
#ls -l ${ModBlocksOnlyLgrIbAfter}${CNPR}*

#ls -l ${ModBlocksOnlyLscBlock}${CNPR}*
#ls -l ${ModBlocksOnlyLscDof}${CNPR}*

#ls -l ${ModBlocksOnlyLscMasterIblock}${CNPR}*
#ls -l ${ModBlocksOnlyLscMasterBlock}${CNPR}*
#ls -l ${ModBlocksOnlyLscMasterDof}${CNPR}*

#########################################################

#ls -l ${AllFBlockLgrIbBefore}${CNPR}*
#ls -l ${AllFBlockLgrIbAfter}${CNPR}*

#ls -l ${AllFBlockLscBlock}${CNPR}*
#ls -l ${AllFBlockLscDof}${CNPR}*

#ls -l ${AllFBlockLscMasterIblock}${CNPR}*
#ls -l ${AllFBlockLscMasterBlock}${CNPR}*
#ls -l ${AllFBlockLscMasterDof}${CNPR}*




############################################################################
F1=""
F2=""
difffiles()
{

  DIFF=$(diff ${F1} ${F2} )
  if [ "$DIFF" != "" ]
  then
    echo "diff ${F1} ${F2}"
  fi
}
############################################################################

## First determine which Counter does the extracted matrices differ.

#diff RMBO_lgr_ibbefore_C0NP3R0_0000 RAFB_lgr_ibbefore_C0NP3R0_0000
for i in {1..30}
do
F1="${MATFOLDER1}/${ModBlocksOnlyLgrIbBefore}C${i}NP3R0_0000"
F2="${MATFOLDER2}/${ModBlocksOnlyLgrIbBefore}C${i}NP3R0_0000"

difffiles
done








cd $CURRDIR






