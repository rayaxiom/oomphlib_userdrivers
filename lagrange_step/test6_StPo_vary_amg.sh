#!/bin/bash

# The name of the program
PROGRAM="step_po"

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

# Folder of where the iteration counts will be.
# This will be set later since we are doing amg sweeps, the folder name
# will reflect the amg parameters used.
RESITS_DIR=""

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
###############################################################################

# There may be two test lists, so we declare the strings up here.
# Then we change this right before we call the gen_testxy functions.
TEST_FILEBASE=""
TEST_LIST=""

# NOTE: The F amg settings are the same as --f_solver 69

## Preconditioner parameters.
Famg_BASE="--f_solver 96"
Famg_STRN_SIM="--f_amg_str 0.25"
Famg_STRN_STR="--f_amg_str 0.668"
Famg_COARSE="--f_amg_coarse 1" #RS

## Creates test lists for all prec combinations for noel = 2 to 128
function gen_tests_1v22_GS_RS_SIM()
{

RESITS_DIR="It1Sit2_GSs_RSc_Strn0.25"

Famg_ITER="--f_amg_iter 1"
Famg_SMITER="--f_amg_smiter 2"
Famg_SMOOTHER="--f_amg_sim_smoo 1" #GS
Famg_DAMP="--f_amg_damp -1"

Famg_sim="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
Prec_WLu_NSLscPamgFamgsim="--w_solver 0 --ns_solver 1 --p_solver 1 $Famg_sim"

PRECLIST="0"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0"
ANGLIST="0 30 67"
RELIST="0"
REPARAM=""
NOELLIST="2 4 8 16 32 64 128"

for PREC  in $PRECLIST
do
  PRECPARAM=$Prec_WLu_NSLscPamgFamgsim
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        REPARAM="--rey_start 0 --rey_end 200 --rey_incre 25"
        for NOEL in $NOELLIST
        do
echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 100 --dist_prob --trilinos_solver --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG $REPARAM --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST

        done
      done
    done
  done
done
} # gen_tests function

function gen_tests_2v22_GS_RS_SIM()
{

RESITS_DIR="It2Sit2_GSs_RSc_Strn0.25"

Famg_ITER="--f_amg_iter 2"
Famg_SMITER="--f_amg_smiter 2"
Famg_SMOOTHER="--f_amg_sim_smoo 1" #GS
Famg_DAMP="--f_amg_damp -1"

Famg_sim="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
Prec_WLu_NSLscPamgFamgsim="--w_solver 0 --ns_solver 1 --p_solver 1 $Famg_sim"

PRECLIST="0"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0"
ANGLIST="0 30 67"
RELIST="0"
REPARAM=""
NOELLIST="2 4 8 16 32 64 128"

for PREC  in $PRECLIST
do
  PRECPARAM=$Prec_WLu_NSLscPamgFamgsim
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        REPARAM="--rey_start 0 --rey_end 200 --rey_incre 25"
        for NOEL in $NOELLIST
        do
echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 100 --dist_prob --trilinos_solver --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG $REPARAM --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST

        done
      done
    done
  done
done
} # gen_tests function

function gen_tests_1v22_EUCLID_RS_SIM()
{

RESITS_DIR="It1Sit2_EUCLIDs_RSc_Strn0.25"

Famg_ITER="--f_amg_iter 1"
Famg_SMITER="--f_amg_smiter 2"
Famg_SMOOTHER="--f_amg_com_smoo 9" #GS
Famg_DAMP="--f_amg_damp -1"

Famg_sim="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_SIM $Famg_COARSE"
Prec_WLu_NSLscPamgFamgsim="--w_solver 0 --ns_solver 1 --p_solver 1 $Famg_sim"

PRECLIST="0"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="0"
ANGLIST="0 30 67"
RELIST="0"
REPARAM=""
NOELLIST="2 4 8 16 32 64 128"

for PREC  in $PRECLIST
do
  PRECPARAM=$Prec_WLu_NSLscPamgFamgsim
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        REPARAM="--rey_start 0 --rey_end 200 --rey_incre 25"
        for NOEL in $NOELLIST
        do
echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 100 --dist_prob --trilinos_solver --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG $REPARAM --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST

        done
      done
    done
  done
done
} # gen_tests function

function gen_tests_1v22_EUCLID_RS_STR()
{

RESITS_DIR="It1Sit2_EUCLIDs_RSc_Strn0.668"

Famg_ITER="--f_amg_iter 1"
Famg_SMITER="--f_amg_smiter 2"
Famg_SMOOTHER="--f_amg_com_smoo 9" #GS
Famg_DAMP="--f_amg_damp -1"

Famg_str="$Famg_BASE $Famg_ITER $Famg_SMITER $Famg_SMOOTHER $Famg_DAMP $Famg_STRN_STR $Famg_COARSE"
Prec_WLu_NSLscPamgFamgstr="--w_solver 0 --ns_solver 1 --p_solver 1 $Famg_str"

PRECLIST="0"
# The precs are set according to the PRECLIST above.
PRECPARAM=""

VISLIST="1"
ANGLIST="0 30 67"
RELIST="0"
REPARAM=""
NOELLIST="2 4 8 16 32 64 128"

for PREC  in $PRECLIST
do
  PRECPARAM=$Prec_WLu_NSLscPamgFamgstr
  for VIS in $VISLIST
  do
    for ANG in $ANGLIST
    do
      for RE in $RELIST
      do
        REPARAM="--rey_start 0 --rey_end 200 --rey_incre 25"
        for NOEL in $NOELLIST
        do
echo "mpirun -np 1 ./$PROGRAM --max_solver_iter 100 --dist_prob --trilinos_solver --prob_id 11 $PRECPARAM --visc $VIS --ang $ANG $REPARAM --noel $NOEL --itstimedir $RESITS_DIR" >> $TEST_LIST

        done
      done
    done
  done
done
} # gen_tests function

TEST_FILEBASE="tests_StPo_vary_amg"
TEST_LIST="$TEST_FILEBASE.list"
gen_tests_1v22_GS_RS_SIM
gen_tests_2v22_GS_RS_SIM
gen_tests_1v22_EUCLID_RS_SIM
gen_tests_1v22_EUCLID_RS_STR
TEST_RUN="$TEST_FILEBASE.sh"
echo "#!/bin/bash" >> $TEST_RUN
cat $TEST_LIST >> $TEST_RUN

cp ./../$0 .


###############################################################################
###############################################################################
###############################################################################
###############################################################################
### Now create the qsub file.
QSUBFILE="$FILEBASE.qsub"
NUMTESTS=$(cat $TEST_LIST | wc -l)
echo '#!/bin/bash' >> $QSUBFILE
echo '#$ -S /bin/bash' >> $QSUBFILE
echo '#$ -cwd' >> $QSUBFILE
echo '#$ -V' >> $QSUBFILE

echo -e "\n" >> $QSUBFILE

echo "#$ -t 1-$NUMTESTS" >> $QSUBFILE

echo -e "\n" >> $QSUBFILE

# Do the same thing with the output directory for qsub
QSUBOUTPUT_DIR="qsub_output"
echo "if [ ! -d \"$QSUBOUTPUT_DIR\" ]; then" >> $QSUBFILE
echo "  mkdir $QSUBOUTPUT_DIR" >> $QSUBFILE
echo "fi" >> $QSUBFILE

echo -e "\n" >> $QSUBFILE

## Some comments for the script.
echo "# Task id 1 will read line 1 from $TEST_LIST" >> $QSUBFILE
echo "# Task id 2 will read line 2 from $TEST_LIST" >> $QSUBFILE
echo "# and so on..." >> $QSUBFILE
echo "# Each line contains the run command with a different set of parameters" >> $QSUBFILE

echo -e "\n" >> $QSUBFILE

## Get the run command from TEST_LIST
RUNLINE='FULL_RUNCOMMAND=`awk "NR==$SGE_TASK_ID" '
RUNLINE+="$TEST_LIST"
RUNLINE+='`'
echo $RUNLINE >> $QSUBFILE

# Now run the command!
echo '$FULL_RUNCOMMAND' >> $QSUBFILE

echo -e "\n" >> $QSUBFILE

# Clean up, move the qsub output and error files into QSUBOUTPUT_DIR
CLEANUPLINE="mv $QSUBFILE"
CLEANUPLINE+='.*.$SGE_TASK_ID '
CLEANUPLINE+=" ./$QSUBOUTPUT_DIR/"
echo $CLEANUPLINE >> $QSUBFILE



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
  mkdir -p $SCRATCH_TEST_DIR/It1Sit2_GSs_RSc_Strn0.25
  mkdir -p $SCRATCH_TEST_DIR/It2Sit2_GSs_RSc_Strn0.25
  mkdir -p $SCRATCH_TEST_DIR/It1Sit2_EUCLIDs_RSc_Strn0.25
  mkdir -p $SCRATCH_TEST_DIR/It1Sit2_EUCLIDs_RSc_Strn0.668
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR
fi



