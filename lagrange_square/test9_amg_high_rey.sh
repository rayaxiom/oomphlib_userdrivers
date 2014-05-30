#!/bin/bash

set -e

# The name of the program
PROGRAM="sq_lgr"

# Create the FILEBASE, this is the folder where the testing will be done.
THISFILE=$0 # This contains "./", which we do not want.
THISFILE=${THISFILE:2} # Gets rid of "./"
FILEBASE=${THISFILE%%.*} # Get rid of the extension (in this case, ".sh")

# Create the new folder (remove old one)
touch $FILEBASE
rm -rf $FILEBASE
mkdir $FILEBASE

# Get the current directory and the oomph-base
CURRENT_DIR=`pwd`
OOMPHROOT_DIR=$(make -s --no-print-directory print-top_builddir)

# folder of where the iteration counts will be.
RESITS_DIR="res_iterations"

# Get version of oomph-lib
cd $CURRENT_DIR
cd $OOMPHROOT_DIR
git log -1 > $CURRENT_DIR/$FILEBASE/oomphlib_revision
cd $CURRENT_DIR

# Get version of user drivers
cd $OOMPHROOT_DIR/user_drivers
git log -1 > $CURRENT_DIR/$FILEBASE/user_driver_revision
cd $CURRENT_DIR

# make the program and move it into the test folder.
make $PROGRAM
mv $PROGRAM ./$FILEBASE
cd $FILEBASE

###############################################################################
# Now we are inside the test folder (FILEBASE)

# Make the results directory. I have it in an if statement...
# Because some times when it exists, we wish to reuse it... of course not this 
# time...
if [ ! -d "$RESITS_DIR" ]; then
  mkdir $RESITS_DIR
fi

# There may be two test lists, so we declare the strings up here.
# Then we change this right before we call the gen_testxy functions.
TEST_FILEBASE=""
TEST_LIST=""

# NOTE: The F amg settings are the same as --f_solver 69

## Preconditioner parameters.
Famg_BASE="--f_solver 96"
Famg_ITER="--f_amg_iter 1"
Famg_SMITER="--f_amg_smiter 2"
Famg_SMOOTHER="--f_amg_com_smoo 9" #Euclid
Famg_DAMP="--f_amg_damp -1"
Famg_STRN_SIM="--f_amg_str 0.25"
Famg_STRN_STR="--f_amg_str 0.668"
Famg_COARSE="--f_amg_coarse 1" #RS

Prec_WLu_NSLu="--w_solver 0 --ns_solver 0"
Prec_WLu_NSLscLu="--w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0"
Prec_WLu_NSLscPamgFlu="--w_solver 0 --ns_solver 1 --p_solver 1 --f_solver 0"

Famg_sim="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
Famg_str="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_STR $Famg_COARSE"
Prec_WLu_NSLscPLuFamgsim="--w_solver 0 --ns_solver 1 --p_solver 0 $Famg_sim"
Prec_WLu_NSLscPLuFamgstr="--w_solver 0 --ns_solver 1 --p_solver 0 $Famg_str"
Prec_WLu_NSLscPamgFamgsim="--w_solver 0 --ns_solver 1 --p_solver 1 $Famg_sim"
Prec_WLu_NSLscPamgFamgstr="--w_solver 0 --ns_solver 1 --p_solver 1 $Famg_str"


## Creates test lists for all prec combinations for noel = 4 to 128
function gen_po_tests()
{
#PRECLIST="0 1 2" # Doing either full exact or Exact Navier Stokes
# 0 - W SuperLU, NS SuperLU
# 1 - W SuperLU, NS LSC: P SuperLU, F SuperLU
# 2 - W SuperLU, NS LSC: P AMG, F LU
# 3 - W SuperLU, NS LSC: P Lu, F AMG
# 4 - W Super LU, NS LSC: P AMG F AMG

PRECLIST="4"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0 1"
ANGLIST="0"
RELIST="0 100 200 500 1000"
NOELLIST="4 8 16 32 64 128 256 512"

for PREC  in $PRECLIST
do

  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do
case "$PREC" in
  0)
    PRECPARAM="$Prec_WLu_NSLu"
    ;;
  1)
    PRECPARAM="$Prec_WLu_NSLscLu"
    ;;
  2)
    PRECPARAM="$Prec_WLu_NSLscPamgFlu"
    ;;
  3)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$Prec_WLu_NSLscPLuFamgsim"
    else
      PRECPARAM="$Prec_WLu_NSLscPLuFamgstr"
    fi
    ;;
  4)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$Prec_WLu_NSLscPamgFamgsim"
    else
      PRECPARAM="$Prec_WLu_NSLscPamgFamgstr"
    fi
    ;;
esac

echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 110 --dist_prob --trilinos_solver --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST

        done
      done
    done
  done
done
} # gen_po_tests function


## Prec 0, LU for both P and F solve.
VaPREC0_Plu_Flu="--lsc_only --p_solver 0 --f_solver 0"
## Prec 1 P - AMG, F - LU
VaPREC1_Pamg_Flu="--lsc_only --p_solver 1 --f_solver 0"
## Prec 2 P - LU F - AMG
VaPREC2_Plu_Famgsim="--lsc_only --p_solver 0 $Famg_sim"
VaPREC2_Plu_Famgstr="--lsc_only --p_solver 0 $Famg_str"
## PREC 3 full AMG
VaPREC3_Pamg_Famgsim="--lsc_only --p_solver 1 $Famg_sim"
VaPREC3_Pamg_Famgstr="--lsc_only --p_solver 1 $Famg_str"

## Creates test lists for all prec combinations for noel = 4 to 128
function gen_va_tests()
{
PRECLIST="3"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0 1"
RELIST="0 100 200 500 1000"
NOELLIST="4 8 16 32 64 128 256 512"

for PREC  in $PRECLIST
do
  for VIS in $VISLIST
  do
      for RE in $RELIST
      do
        for NOEL in $NOELLIST
        do
case "$PREC" in
  0)
    PRECPARAM="$PREC0_Plu_Flu"
    ;;
  1)
    PRECPARAM="$PREC1_Pamg_Flu"
    ;;
  2)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$PREC2_Plu_Famgsim"
    else
      PRECPARAM="$PREC2_Plu_Famgstr"
    fi

    ;;
  3)
    if [ "$VIS" -eq "0" ]; then
      PRECPARAM="$PREC3_Pamg_Famgsim"
    else
      PRECPARAM="$PREC3_Pamg_Famgstr"
    fi
    ;;
esac

echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 110 --dist_prob --trilinos_solver --prob_id 88 $PRECPARAM --visc $VIS --rey $RE --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST

        done
      done
  done
done
} # gen_tests function



TEST_FILEBASE="tests_amg_PoVa_high_rey"
TEST_LIST="$TEST_FILEBASE.list"

gen_po_tests
gen_va_tests

genqsub_script $TEST_LIST

TEST_RUN="$TEST_FILEBASE.sh"
echo "#!/bin/bash" >> $TEST_RUN
cat $TEST_LIST >> $TEST_RUN

cp ./../$0 .


################################################################################
################################################################################
################################################################################
################################################################################
#### Now create the qsub file.
#QSUBFILE="$FILEBASE.qsub"
#NUMTESTS=$(cat $TEST_LIST | wc -l)
#echo '#!/bin/bash' >> $QSUBFILE
#echo '#$ -S /bin/bash' >> $QSUBFILE
#echo '#$ -cwd' >> $QSUBFILE
#echo '#$ -V' >> $QSUBFILE
#
#echo -e "\n" >> $QSUBFILE
#
#echo "#$ -t 1-$NUMTESTS" >> $QSUBFILE
#
#echo -e "\n" >> $QSUBFILE
#
#echo "# This should be ran in the scratch file system." >> $QSUBFILE
#echo "# Thus results directory '$RESITS_DIR' may not exist." >> $QSUBFILE
#echo -e "# We create it if it does not exist.\n" >> $QSUBFILE
#echo "if [ ! -d \"$RESITS_DIR\" ]; then" >> $QSUBFILE
#echo "  mkdir $RESITS_DIR" >> $QSUBFILE
#echo "fi" >> $QSUBFILE
#
#echo -e "\n" >> $QSUBFILE
#
## Do the same thing with the output directory for qsub
#QSUBOUTPUT_DIR="qsub_output"
#echo "if [ ! -d \"$QSUBOUTPUT_DIR\" ]; then" >> $QSUBFILE
#echo "  mkdir $QSUBOUTPUT_DIR" >> $QSUBFILE
#echo "fi" >> $QSUBFILE
#
#echo -e "\n" >> $QSUBFILE
#
### Some comments for the script.
#echo "# Task id 1 will read line 1 from $TEST_LIST" >> $QSUBFILE
#echo "# Task id 2 will read line 2 from $TEST_LIST" >> $QSUBFILE
#echo "# and so on..." >> $QSUBFILE
#echo "# Each line contains the run command with a different set of parameters" >> $QSUBFILE
#
#echo -e "\n" >> $QSUBFILE
#
### Get the run command from TEST_LIST
#RUNLINE='FULL_RUNCOMMAND=`awk "NR==$SGE_TASK_ID" '
#RUNLINE+="$TEST_LIST"
#RUNLINE+='`'
#echo $RUNLINE >> $QSUBFILE
#
## Now run the command!
#echo '$FULL_RUNCOMMAND' >> $QSUBFILE
#
#echo -e "\n" >> $QSUBFILE
#
## Clean up, move the qsub output and error files into QSUBOUTPUT_DIR
#CLEANUPLINE="mv $QSUBFILE"
#CLEANUPLINE+='.*.$SGE_TASK_ID '
#CLEANUPLINE+=" ./$QSUBOUTPUT_DIR/"
#echo $CLEANUPLINE >> $QSUBFILE



###############################################################################
###############################################################################
###############################################################################
###############################################################################
################### Now check if I'm on csf, if so, delete the related scratch
# and copy the current stuff there.
if [[ $HOME == *mbax5ml3* ]]
then
  CURRENT_DIR=`pwd`
  cd ..
  PROGRAM_DIR=${PWD##*/}
  cd $CURRENT_DIR
  OOMPH_PROGRAM_DIR="/mnt/iusers01/mh01/mbax5ml3/oomphlib_optimized/user_drivers/$PROGRAM_DIR"
  SCRATCH_PROGRAM_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/oomphlib_optimized/user_drivers/$PROGRAM_DIR"

  OOMPH_TEST_DIR="$OOMPH_PROGRAM_DIR/$FILEBASE"
  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE"

  echo "OOMPH_TEST_DIR: $OOMPH_TEST_DIR"
  echo "SCRATCH_TEST_DIR: $SCRATCH_TEST_DIR"

  # Remove the scratch stuff.
  rm -rf $SCRATCH_TEST_DIR

  rsync -av $OOMPH_TEST_DIR/$PROGRAM $SCRATCH_TEST_DIR/
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR/
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR/

  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR
fi



