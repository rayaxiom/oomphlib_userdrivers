#!/bin/bash

FILENAME="newtest07_smoothers.qsub.o*."
#RAYITS
#RAYPRECSETUP
#RAYLINSOLVER
TAG="RAYITS"
#TAG="RAYPRECSETUP"
#TAG="RAYLINSOLVER"

function local_format_results()
{
# We loop through the files via artificial loops
# Sim, str
Nvisc="0 1"
# 0 100, 200, 
Nrey="1 2 3"
# J1v22, J2v22, GS1v22, GS2v22, Euclid
Nprec="1 2 3 4 5"
# noel 4 8 16 32 64 128 256 512
Nnoel="1 2 3 4 5 6 7 8"

COUNTER="1"

#RESFILE="${RESFOLDER_LIST[0]}/${FILENAME}${COUNTER}"
#echo "$RESFILE"

for visc in $Nvisc
do
  for rey in $Nrey
  do
    for prec in $Nprec
    do
LINE=""
      for noel in $Nnoel
      do


## Set the first min token
#MIN_TOKEN=""
#MIN_NUM="-1.0"

# Get the first results file and see if it exists
RESFILE="${RESFOLDER_LIST[0]}/${FILENAME}${COUNTER}"

if [ -f $RESFILE ]; then
  MIN_TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
  MIN_NUM=${MIN_TOKEN%(*}
  for RESFOLDER in "${RESFOLDER_LIST[@]}"
  do
    RESFILE="${RESFOLDER}/${FILENAME}${COUNTER}"

    CURRENT_TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
    CURRENT_NUM=${CURRENT_TOKEN%(*}

    # now check to see which is bigger...
    bools=$(echo "$CURRENT_NUM < $MIN_NUM" | bc)
    
    if [ "$bools" == "1" ]; then
      MIN_TOKEN=$CURRENT_TOKEN
      MIN_NUM=$CURRENT_NUM
    fi
  done

  LINE="$LINE $MIN_TOKEN"
else
  echo "FILE NOT FOUND: $RESFILE"
fi

#RESFILE="${RESFOLDER_LIST[0]}/${FILENAME}${COUNTER}"
#TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
#LINE="$LINE $TOKEN"

# Increment counter
COUNTER=$((COUNTER+1))
      done
echo $LINE
    done
  done
done

}


#RESFOLDER_LIST=("test1/qsub_output" "test2/qsub_output" "test3/qsub_output" )
RESFOLDER_LIST=("test1/qsub_output" )
local_format_results

