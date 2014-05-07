#!/bin/bash

set -e



CURRENTDIR=`pwd`
# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)
PROGRAM="sq_lgr"

cd $OOMPH_ROOT_DIR
./autogen_quick.sh --rebuild --jobs=4
cd $CURRENTDIR
make $PROGRAM


TESTFOLDER="param_test_4_128"
touch $TESTFOLDER
rm -rf $TESTFOLDER
mkdir $TESTFOLDER

cd $OOMPH_ROOT_DIR
git log -1 > $CURRENTDIR/oomphlib_revision
cd $CURRENTDIR

git log -1 > user_driver_revision

mv oomphlib_revision ./$TESTFOLDER/
mv user_driver_revision ./$TESTFOLDER/

cd $TESTFOLDER


ITSTIMEDIR="itstimedir"
touch $ITSTIMEDIR
rm -rf $ITSTIMEDIR
mkdir $ITSTIMEDIR

TESTLIST_FILE="testlist.list"
TESTOUT_FILE="test_output"
TESTLISTFIN_FILE="tests_done"




gen_tests()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - W SuperLU, NS SuperLU
# 1 - W SuperLU, NS LSC: P SuperLU, F SuperLU
# 2 - W SuperLU, NS LSC: P AMG, F AMG

PRECLIST="0 1 2"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0 1"
ANGLIST="0 30 67"
RELIST="0 100 200"
NOELLIST="4 8 16 32 64 128 256"

for PREC  in $PRECLIST
do
  case "$PREC" in
    0)
      PRECPARAM="--w_solver 0 --ns_solver 0"
      ;;
    1)
      PRECPARAM="--w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0"
      ;;
    2)
      PRECPARAM="--w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 69"
      ;;
  esac
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do
          echo "mpirun -np 1 ./$PROGRAM --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG --rey $RE --noel $NOEL --itstimedir $ITSTIMEDIR" >> $TESTLIST_FILE
        done
      done
    done
  done
done
} # gen_tests function

run_tests()
{
  # State that this is the initial test run (as opposed to run_tests_con)
  # Total number of tests
  NUMTESTS=$(wc $TESTLIST_FILE | awk '{print $1;}')
  echo "RAYRAY: initial test run, total tests: $NUMTESTS"
  
#  # Check if SHELLOUTPUT_FILE contains previous runs.
#  LINESINSHELLOUT=$(grep "RAYRAY: initial" $SHELLOUTPUT_FILE | wc | awk '{print $1;}')
#
#  if [ -f $SHELLOUTPUT_FILE ]; then
#    if [ "$LINESINSHELLOUT" -ne "1" ]; then
#      echo "This is the first run but you have not removed the previous $SHELLOUTPUT_FILE file."
#      echo "Please delete both the file $SHELLOUTPUT_FILE and the directory $ITSTIMEDIR"
#      echo "Then re-run this script."
#      exit 1
#    fi
#  fi


  # Current test number.
  TESTNUM=1
  
  while read RUNCOMMAND; do
    echo "Test $TESTNUM/$NUMTESTS: \"$RUNCOMMAND\" on $(date)"
    bash -c "$RUNCOMMAND" </dev/null >> $TESTOUT_FILE

    # For some reason I have to echo to stdout, otherwise it will execute
    # the rest of the loop, updating TESTLISTFIN_FILE even when I ctrl+c
    echo ""
    TESTNUM=$[$TESTNUM+1]
    echo "$RUNCOMMAND" >> $TESTLISTFIN_FILE
  done < $TESTLIST_FILE
}

gen_tests
mv ./../$PROGRAM .
run_tests







