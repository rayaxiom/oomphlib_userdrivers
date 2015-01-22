#!/bin/bash

NEWDIRBASE="/home/ray/oomphlib/mpi_debug_paranoid/user_drivers"
OLDDIRBASE="/media/ray/PENDRIVE/mpi_debug_paranoid/user_drivers_temp"

# CHANGE HERE
DIRNAME="time_stepping"

# CHANGE HERE
declare -a filearray=(
"Makefile.am"
"rayleigh_channel.cc"
"rayleigh_traction_channel.cc"
"sq_lgr.cc"
"three_d_fp.cc"
"two_d_fp.cc"
"unstructured_three_d_fluid.cc")


THISDIR="$NEWDIRBASE/$DIRNAME"
OLDDIR="$OLDDIRBASE/$DIRNAME"

DELFLAG="$1"

for filen in "${filearray[@]}"
do
  FULLTHISFILE="$THISDIR/$filen"
  FULLOLDFILE="$OLDDIR/$filen"

if [ ! -f $FULLTHISFILE ]; then
  NEWFOLDER="$THISDIR/"
  cp $FULLOLDFILE $NEWFOLDER
  echo "   MOVED: $filen to current"
elif [ ! -f $FULLOLDFILE ]; then
  echo "   NOT ON OLD: $filen"
else
  DIFF=$(diff $FULLTHISFILE $FULLOLDFILE)
  if [ "$DIFF" != "" ] 
  then
    echo "   Modified: $filen"
  else
    echo "Delete: $filen"
#    if [ "$DELFLAG" = "del" ]; then
      rm $FULLOLDFILE
#    fi
  fi
fi
done












