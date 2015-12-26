#!/bin/bash

FILENAME="newtest1_famg_blockdiag.qsub.o923534."
#RAYITS
#RAYPRECSETUP
#RAYLINSOLVER
TAG="RAYITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

# Sim diag
NUMLIST="1 2 3 4 5 6 7"
LINE=""
for NUM in $NUMLIST
do
RESFILE="${FILENAME}${NUM}"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"
done
echo $LINE

# Sim upper
NUMLIST="8 9 10 11 12 13 14"
LINE=""
for NUM in $NUMLIST
do
RESFILE="${FILENAME}${NUM}"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"
done
echo $LINE

# Sim lower
NUMLIST="15 16 17 18 19 20 21"
LINE=""
for NUM in $NUMLIST
do
RESFILE="${FILENAME}${NUM}"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"
done
echo $LINE

# Sim full
NUMLIST="22 23 24 25 26 27 28"
LINE=""
for NUM in $NUMLIST
do
RESFILE="${FILENAME}${NUM}"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"
done
echo $LINE


#################
# Str diag
NUMLIST="29 30 31 32 33 34 35"
LINE=""
for NUM in $NUMLIST
do
RESFILE="${FILENAME}${NUM}"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"
done
echo $LINE

# Str upper
NUMLIST="36 37 38 39 40 41 42"
LINE=""
for NUM in $NUMLIST
do
RESFILE="${FILENAME}${NUM}"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"
done
echo $LINE

# Str lower
NUMLIST="43 44 45 46 47 48 49"
LINE=""
for NUM in $NUMLIST
do
RESFILE="${FILENAME}${NUM}"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"
done
echo $LINE

# Str full
NUMLIST="50 51 52 53 54 55 56"
LINE=""
for NUM in $NUMLIST
do
RESFILE="${FILENAME}${NUM}"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"
done
echo $LINE





