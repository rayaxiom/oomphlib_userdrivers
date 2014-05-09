#!/bin/bash

set -e

# CHANGE THIS
TESTFOLDER="params_oomphgmres"
PROGRAM="step_po"


CURRENTDIR=`pwd`
# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

# Rebuild oomph-lib.
#cd $OOMPH_ROOT_DIR
#./autogen_quick.sh --rebuild --jobs=4
#cd $CURRENTDIR
#make $PROGRAM


touch $TESTFOLDER
rm -rf $TESTFOLDER
mkdir $TESTFOLDER

# Get version of oomph-lib
cd $OOMPH_ROOT_DIR
git log -1 > $CURRENTDIR/$TESTFOLDER/oomphlib_revision

# Get version of user drivers
cd $CURRENTDIR
git log -1 > $CURRENTDIR/$TESTFOLDER/user_driver_revision

# make the program and move it into the test folder.
make $PROGRAM
mv $PROGRAM ./$TESTFOLDER

cd $TESTFOLDER


ITSTIMEDIR="itstimedir"
touch $ITSTIMEDIR
rm -rf $ITSTIMEDIR
mkdir $ITSTIMEDIR

TEST1_FILEBASE="tests_2_64"
TEST1_FILELIST="$TEST1_FILEBASE.list"
gen_tests1()
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
RELIST="0"
NOELLIST="2 4 8 16 32 64"

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
          case "$RE" in
            0)
              REYPARAM="--rey_start 0 --rey_end 200 --rey_incre 25"
              ;;
          esac
        for NOEL in $NOELLIST
        do
          echo "mpirun -np 1 ./$PROGRAM --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG $REYPARAM --noel $NOEL --itstimedir $ITSTIMEDIR" >> $TEST1_FILELIST
        done
      done
    done
  done
done
} # gen_tests function



TEST2_FILEBASE="tests_128"
TEST2_FILELIST="$TEST2_FILEBASE.list"
gen_tests2()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - W SuperLU, NS SuperLU
# 1 - W SuperLU, NS LSC: P SuperLU, F SuperLU
# 2 - W SuperLU, NS LSC: P AMG, F AMG

PRECLIST="2"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0 1"
ANGLIST="0 30 67"
RELIST="0"
NOELLIST="128"

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
          case "$RE" in
            0)
              REYPARAM="--rey_start 0 --rey_end 200 --rey_incre 25"
              ;;
          esac
        for NOEL in $NOELLIST
        do
          echo "mpirun -np 1 ./$PROGRAM --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG $REYPARAM --noel $NOEL --itstimedir $ITSTIMEDIR" >> $TEST2_FILELIST
        done
      done
    done
  done
done
} # gen_tests function

gen_tests1
TEST1_RUN="$TEST1_FILEBASE.sh"
echo "#!/bin/bash" >> $TEST1_RUN
cat $TEST1_FILELIST >> $TEST1_RUN

gen_tests2
TEST2_RUN="$TEST2_FILEBASE.sh"
echo "#!/bin/bash" >> $TEST2_RUN
cat $TEST2_FILELIST >> $TEST2_RUN


cp ./../$0 .






