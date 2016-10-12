#!/bin/bash

CURRENT_DIR=`pwd`

#TETGEN_DIR="tetgenmesh"
TETGEN_FILEBASE="fsi_bifurcation_fluid"

# Create the folder
#touch ${TETGEN_DIR} && rm -rf ${TETGEN_DIR} && mkdir -p $TETGEN_DIR

# copy the poly file.
#cp ${TETGEN_FILEBASE}.poly ./$TETGEN_DIR/

#cd $TETGEN_DIR

# Use default quality -pq
# Which is -pq2.0/0.
#   the first constraint is the maximum allowable radius-edge ratio, 
#     default is 2.0
#   the second constraint is the minimum allowable dihedral angle, 
#     default is 0 (degree)

tetgen -pqa0.8 ${TETGEN_FILEBASE}.poly
tetgen -rqa0.4 ${TETGEN_FILEBASE}.1
tetgen -rqa0.2 ${TETGEN_FILEBASE}.2
tetgen -rqa0.1 ${TETGEN_FILEBASE}.3
tetgen -rqa0.05 ${TETGEN_FILEBASE}.4
tetgen -rqa0.025 ${TETGEN_FILEBASE}.5
tetgen -rqa0.0125 ${TETGEN_FILEBASE}.6
tetgen -rqa0.00625 ${TETGEN_FILEBASE}.7
tetgen -rqa0.003125 ${TETGEN_FILEBASE}.8
tetgen -rqa0.0015625 ${TETGEN_FILEBASE}.9
tetgen -rqa0.00078125 ${TETGEN_FILEBASE}.10
tetgen -rqa0.000390625 ${TETGEN_FILEBASE}.11
tetgen -rqa0.0001953125 ${TETGEN_FILEBASE}.12
#tetgen -rqa0.00009765625 ${TETGEN_FILEBASE}.13
#tetgen -rqa0.000048828125 ${TETGEN_FILEBASE}.15
#tetgen -rqa0.0000244140625 ${TETGEN_FILEBASE}.16
#tetgen -rqa0.00001220703125 ${TETGEN_FILEBASE}.17
#tetgen -rqa0.000006103515625 ${TETGEN_FILEBASE}.18
#tetgen -rqa0.0000030517578125 ${TETGEN_FILEBASE}.19
#tetgen -rqa0.00000152587890625 ${TETGEN_FILEBASE}.20

cd $CURRENT_DIR

