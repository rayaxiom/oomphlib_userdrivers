#!/bin/bash

NEWBASE="/home/ray/oomphlib/mpi_debug_paranoid/user_drivers"
OLDBASE="/media/ray/PENDRIVE/mpi_debug_paranoid/user_drivers_temp"

CURRENTFOLDER="lagrange_cube/tetgen_files"

declare -a folderarray=(
"0d00009765625"
"0d0001953125"
"0d000390625"
"0d00078125"
"0d0015625"
"0d003125"
"0d00625"
"0d0125"
"0d025"
"0d05"
"0d1"
"0d2"
"0d4")


THISDIR="$NEWBASE/$CURRENTFOLDER"
OLDDIR="$OLDBASE/$CURRENTFOLDER"

for foldern in "${folderarray[@]}"
do

  NEWFOLDER="$THISDIR/$foldern"
  OLDFOLDER="$OLDDIR/$foldern"

  declare -a filerarray=(
  "cube.1.ele"
  "cube.1.face"
  "cube.1.node")

  echo "$foldern"
  for filen in "${filerarray[@]}"
  do
#    echo "$filen"
    NEWFULLFILE="$NEWFOLDER/$filen"
    OLDFULLFILE="$OLDFOLDER/$filen"
   
    if [ ! -f $NEWFULLFILE ]; then
      echo "   NOT ON CURRENT: $filen"
    elif [ ! -f $OLDFULLFILE ]; then
      echo "   NOT ON OLD: $filen"
    else
      DIFF=$(diff $NEWFULLFILE $OLDFULLFILE)
      if [ "$DIFF" != "" ]
      then
        echo "   Modified: $filen"
      else
        echo "Delete: $filen"
#        rm $FULLOLDFILE
      fi
    fi
  done # for - files
done # for - folders












