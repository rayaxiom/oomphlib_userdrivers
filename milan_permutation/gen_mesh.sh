#!/bin/bash

# Delete previous files
rm -rf cube.*.edge cube.*.ele cube.*.face cube.*.node

# Use default quality -pq
# Which is -pq2.0/0.
#   the first constraint is the maximum allowable radius-edge ratio, 
#     default is 2.0
#   the second constraint is the minimum allowable dihedral angle, 
#     default is 0 (degree)

tetgen -pqa0.005 cube.poly
tetgen -rqa0.0025 cube.1
tetgen -rqa0.00125 cube.2
tetgen -rqa0.000625 cube.3
tetgen -rqa0.0003125 cube.4
tetgen -rqa0.00015625 cube.5
tetgen -rqa0.000078125 cube.6
tetgen -rqa0.0000390625 cube.7
tetgen -rqa0.00001953125 cube.8
tetgen -rqa0.000009765625 cube.9
tetgen -rqa0.0000048828125 cube.10
tetgen -rqa0.00000244140625 cube.11
tetgen -rqa0.000001220703125 cube.12
tetgen -rqa0.0000006103515625 cube.13
tetgen -rqa0.00000030517578125 cube.14
tetgen -rqa0.000000152587890625 cube.15
tetgen -rqa0.0000000762939453125 cube.16

