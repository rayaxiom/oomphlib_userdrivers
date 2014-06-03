#!/bin/bash

function local_format_results()
{
CURRENT_DIR=`pwd`

function format_StPo_exacts()
{
ITSTIMEDIR="res_iterations"

cd $ITSTIMEDIR

TAG="RAYITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

LINE=""

Prob_str="StPo"
PRECLIST="WedNe WedNlFePe WedNlFrayPe WedNlFePa"
VISLIST="Sim Str"
ANGLIST="A0 A30 A67"
RELIST="R0 R25 R50 R75 R100 R125 R150 R175 R200"
NOELLIST="N2 N4 N8 N16 N32 N64"

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
}

function format_StPo_amg()
{
ITSTIMEDIR="res_iterations"

cd $ITSTIMEDIR

TAG="RAYITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

LINE=""

Prob_str="StPo"
PRECLIST="WedNlFrayPa"
VISLIST="Sim Str"
ANGLIST="A0 A30 A67"
RELIST=""
NOELLIST="N2 N4 N8 N16 N32 N64 N128"

for PREC in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
if [ "$ANG" == "A0" ]; then
RELIST="R0 R25 R50 R75 R100 R125 R150 R175 R200 R225 R250 R275 R300 R325 R350 R375 R400 R425 R450 R475 R500"
else
RELIST="R0 R25 R50 R75 R100 R125 R150 R175 R200"
fi
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
}

function format_StVa_exacts()
{
ITSTIMEDIR="res_iterations"

cd $ITSTIMEDIR

TAG="RAYITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

LINE=""

Prob_str="StVa"
PRECLIST="FePe FrayPe FePa"
VISLIST="Sim Str"
ANGLIST="A_"
RELIST="R0 R50 R100 R150 R200"
NOELLIST="N2 N4 N8 N16 N32 N64"

for PREC in $PRECLIST
do
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
#if [ "$ANG" == "A0" ]; then
#RELIST="R0 R25 R50 R75 R100 R125 R150 R175 R200 R225 R250 R275 R300 R325 R350 R375 R400 R425 R450 R475 R500"
#else
#RELIST="R0 R25 R50 R75 R100 R125 R150 R175 R200"
#fi
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
}

cd $CURRENT_DIR
echo -e "\n"
echo -e "Formatting StPo EXACTS results"
format_StPo_exacts

cd $CURRENT_DIR
echo -e "\n"
echo -e "Formatting StPo pure AMG results"
format_StPo_amg

cd $CURRENT_DIR
echo -e "\n"
echo -e "Formatting StVa EXACTS results"
format_StVa_exacts



}

local_format_results

