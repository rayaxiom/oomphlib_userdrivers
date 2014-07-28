#!/bin/bash


FILE="cube"

# To generate the mesh, do
# tetgen -a0.4 $FILE.poly
#
# This generates three files:
# $FILE.1.ele
# $FILE.1.face
# $FILE.1.node
#



# The folders are arranged such that 
# -a0.4 is in 0d4
# -a0.2 is in 0d2
# So d stands for decimal!

#HEX = 22,763
#TET = 617
AREA="0.4"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX = 35,797
#TET = 1,023
AREA="0.2"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX = 63,367
#TET = 1,901
AREA="0.1"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX = 107,459
#TET = 3,277
AREA="0.05"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX = 237,157
#TET = 7,704
AREA="0.025"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX ~ 474,314
#TET = 16,411
AREA="0.0125"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX ~ 948,628
#TET ~ 32,822
AREA="0.00625"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX ~ 1,897,256
#TET ~ 65,644
AREA="0.003125"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX ~ 3,794,512
#TET ~ 131,288
AREA="0.0015625"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX ~ 7,589,024
#TET ~ 262,576
AREA="0.00078125"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX ~ 15,178,048
#TET ~ 525,152
AREA="0.000390625"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX ~ 30,356,096
#TET ~ 1,050,304
AREA="0.0001953125"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER

#HEX ~ 60,712,192
#TET ~ 2,100,608
AREA="0.00009765625"
FOLDER=${AREA//./d}
touch $FOLDER && rm -rf $FOLDER && mkdir $FOLDER
tetgen -a$AREA $FILE.poly
mv $FILE.1.ele $FOLDER
mv $FILE.1.face $FOLDER
mv $FILE.1.node $FOLDER
