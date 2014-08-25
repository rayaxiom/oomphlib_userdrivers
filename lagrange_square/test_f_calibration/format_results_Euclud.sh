#!/bin/bash

RES_TAR_FILE="from_csf/res_its_time.tar.gz"
RES_DIR="res_its_time"
CURRENT_DIR=`pwd`

## We assume that the results are tar and gunzipped up some where.
## In this particular instance, they are in ./from_csf/
rm -rf $RES_DIR
tar -xzf $RES_TAR_FILE
cd $RES_DIR

TAG="RAYITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

LINE=""

Prob_str="SqPo"
CYCLELIST="WeNlF_1v2"
STRN="" # - depends on viscuous term
PRECLIST="cRSEuclidP_e"

VISLIST="Sim Str"
ANGLIST="A0 A30 A67"
RELIST="R0 R200"
NOELLIST="N4 N8 N16 N32 N64 N128 N256 N512"

for PREC in $PRECLIST
do
  for CYCLE in $CYCLELIST
  do
  for VIS in $VISLIST
  do

if [ "$VIS" = "Sim" ]; then
  STRN="Strn0-25"
elif [ "$VIS" = "Str" ]; then
  STRN="Strn0-668"
else
  STRN="STRN-NULL"
fi

    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        LINE=""
          for NOEL in $NOELLIST
          do

            # Form the preconditioner string.
          REALPREC="${CYCLE}${STRN}${PREC}"

          RESFILE="${Prob_str}${VIS}${RE}${REALPREC}${ANG}${NOEL}NP1R0"
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

