#!/bin/bash

set -e

RAYNPROC="$1"
FILENUMBER="$2"


TETGEN_DIR="tetgen_files_unit_cube"
TETGEN_FILEBASE="cube"
TETGEN_TARGZ="cube_1_16.tar.gz"

PROGRAM="poisson_3d_no_bpf"

CURRENT_TETGEN_FILE=""

function mprog()
{
  make $PROGRAM
}

function extractmesh()
{
  CURR_DIR=`pwd`
  FILENUM="$1"

  cd $TETGEN_DIR

  TETGEN_FILENUM="${TETGEN_FILEBASE}.${FILENUM}"
  TETGEN_ELE="${TETGEN_FILENUM}.ele"
  TETGEN_FACE="${TETGEN_FILENUM}.face"
  TETGEN_NODE="${TETGEN_FILENUM}.node"

  CURRENT_TETGEN_FILE="${TETGEN_FILENUM}"

  tar -xzf $TETGEN_TARGZ $TETGEN_ELE $TETGEN_FACE $TETGEN_NODE

  cd $CURR_DIR
}

function delmesh()
{
  CURR_DIR=`pwd`
  cd $TETGEN_DIR
  
  rm -rf "${TETGEN_FILEBASE}.*.ele"
  rm -rf "${TETGEN_FILEBASE}.*.face"
  rm -rf "${TETGEN_FILEBASE}.*.node"
  
  cd $CURR_DIR
}

function rprog()
{
NPROC="$1"
MPIRUNCOMMAND="mpirun -np $NPROC"

AMG_ITER="--amg_iter 1"
AMG_SMITER="--amg_smiter 1"
AMG_SIM_SMOO="--amg_sim_smoo 0" # 0 - jacobi, 1 - GS
AMG_DAMP="--amg_damp 0.8"
AMG_STRN="--amg_strn 0.25"
AMG_COARSE="--amg_coarse 6" # 1 - RS. 6 - Falgout

AMG_PARAM="$AMG_ITER $AMG_SMITER $AMG_SIM_SMOO $AMG_STRN $AMG_COARSE $AMG_DAMP"


CURRENT_TETGEN_FILE="$2"
TETGEN_FILE="--tetgenfile ${CURRENT_TETGEN_FILE}"

FULLPARAM="${AMG_PARAM} ${TETGEN_FILE}"

echo $FULLPARAM

${MPIRUNCOMMAND} ./${PROGRAM} ${FULLPARAM}
}


#mprog && extractmesh 1 && rprog 1 && delmesh

TTTFILE="$FILENUMBER"

mprog && rprog $RAYNPROC $TTTFILE



