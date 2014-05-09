#!/bin/bash

TESTFOLDER="params_oomphgmres"
ITSTIMEDIR="itstimedir"

cd $TESTFOLDER
cd $ITSTIMEDIR

TAG="RAYITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

LINE=""

Prob_str="SqPo"
PRECLIST="WedNe WedNlFePe WedNlFaPa"
VISLIST="Sim Str"
ANGLIST="A0 A30 A67"
RELIST="R0 R100 R200"
NOELLIST="N4 N8 N16 N32 N64 N128 N256 N512"

for PREC in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        LINE=""
        for NOEL in $NOELLIST
        do

          RESFILE="${Prob_str}${VIS}${RE}${PREC}${ANG}${NOEL}NP1R0"
          TOKEN=""

          if [ -f $RESFILE ]
          then
            TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
          else
            TOKEN="x"
          fi
          
          LINE="$LINE $TOKEN"
        done
        echo $LINE
      done
    done
  done
done
