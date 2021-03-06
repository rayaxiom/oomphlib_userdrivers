#!/bin/bash

RES_TAR_FILE="from_csf/res_its_time.tar.gz"
RES_DIR="res_its_time"
CURRENT_DIR=`pwd`

## We assume that the results are tar and gunzipped up some where.
## In this particular instance, they are in ./from_csf/
rm -rf $RES_DIR
tar -xzf $RES_TAR_FILE
cd $RES_DIR

#TAG="RAYITS"
#TAG="RAYPRECSETUP"
TAG="RAYLINSOLVER"

LINE=""

Prob_str="SqPo"
CYCLELIST="WeNlF_1v2 WeNlF_2v2"
STRN="" # - depends on viscuous term
DAMPLIST="cRSJac0.1P_e cRSJac0.2P_e cRSJac0.3P_e cRSJac0.4P_e cRSJac0.5P_e cRSJac0.6P_e cRSJac0.7P_e cRSJac0.8P_e cRSJac0.9P_e cRSJac1P_e"

VISLIST="Sim Str"
ANGLIST="A0 A67"
RELIST="R0 R200"
NOELLIST="N128 N256"

for CYCLE in $CYCLELIST
do
  for VIS in $VISLIST
  do

if [ "$VIS" = "Sim" ]; then
  STRN="Strn0.25"
elif [ "$VIS" = "Str" ]; then
  STRN="Strn0.668"
else
  STRN="STRN-NULL"
fi

    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do
        LINE=""
          for DAMP in $DAMPLIST
          do

            # Form the preconditioner string.
          PREC="${CYCLE}${STRN}${DAMP}"

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
done

cd $CURRENT_DIR
rm -rf $RES_DIR

