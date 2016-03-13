#!/bin/bash

PROGRAM="annular_wedge_threed"


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

## RRNEW we run this test 3 times then collect the lowest execution times
# So we create three folder, test1, test2 and test3
# Within each, we need to copy over the files:
# .qsub, .list and the program file, then in the scratch directory, we create
# the directories qsub_output and res_iteration in each of the three folders.
TESTFOLDER1="test1"
TESTFOLDER2="test2"
TESTFOLDER3="test3"
mkdir $TESTFOLDER1
mkdir $TESTFOLDER2
mkdir $TESTFOLDER3

############################################################################

############################################################################
# Now we are inside the test folder (FILEBASE)

# Make the results directory. I have it in an if statement...
# Because some times when it exists, we wish to reuse it... 
# of course not this time...
if [ ! -d "$RESITS_DIR" ]; then
  mkdir $RESITS_DIR
fi

# There may be two test lists, so we declare the strings up here.
# Then we change this right before we call the gen_testxy functions.
TEST_FILEBASE=""
TEST_LIST=""




PROBPARAM="--time_type 1 --solver_type 2 --dist_prob --max_solver_iter 100 --dt 0.01 --time_start 0.0 --time_end 0.5  --itstimedir $RESITS_DIR --visc 0 --rey 200 --prob_id 0"
# --noel taken out

LGRPREC="--w_solver 0 --ns_solver 1"
PPREC="--p_solver 1 --p_amg_iter 2 --p_amg_smiter 1 --p_amg_sim_smoo 4 --p_amg_damp 0.668 --p_amg_coarse 6 --p_amg_str 0.7 --print_p_hypre"
FPREC="--f_solver 1 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_sim_smoo 4 --f_amg_damp 1 --f_amg_str 0.668 --print_f_hypre"
# f_amg_coarse taken out.

#mpirun -np 2 taskset -c 0,12 ./annular_wedge_threed       
## Now, things that vary are:
# NPROC
# With this change, I assign:
# * mpirun -np NPROC taskset -c 0,12
# * --noel
# Generate a test list name.

# f_amg_coarse

NPROC=""
function gen_tests()
{
#FAMGCOARSELIST="0 3 6 8 10"
FAMGCOARSELIST="6"

MPIRUN=""
NOEL=""
case "$NPROC" in
  1)
    MPIRUN="mpirun -np 1 taskset -c 0"
    NOEL="35"
    ;;
  2)
    MPIRUN="mpirun -np 2 taskset -c 0,12"
    NOEL="35"
    ;;
  4)
    MPIRUN="mpirun -np 4 taskset -c 0-1,12-13"
    NOEL="35"
    ;;
  8)
    MPIRUN="mpirun -np 8 taskset -c 0-3,12-15"
    NOEL="35"
    ;;
  16)
    MPIRUN="mpirun -np 16 taskset -c 0-7,12-19"
    NOEL="35"
    ;;
#  24)
#    MPIRUN="mpirun -np 24 taskset -c 0-23"
#    NOEL="41"
#    ;;
esac

for FAMGCOARSE in $FAMGCOARSELIST
do
echo "${MPIRUN} ./${PROGRAM} ${PROBPARAM} --noel ${NOEL} ${LGRPREC} ${PPREC} ${FPREC} --f_amg_coarse ${FAMGCOARSE}" >> $TEST_LIST
done
}

############################################################################
############################################################################
############################################################################
# Now we generate the test lists.
NPROC="1"
TEST_FILEBASE="testlist_np${NPROC}"
TEST_LIST="${TEST_FILEBASE}.list"
gen_tests

NPROC="2"
TEST_FILEBASE="testlist_np${NPROC}"
TEST_LIST="${TEST_FILEBASE}.list"
gen_tests

NPROC="4"
TEST_FILEBASE="testlist_np${NPROC}"
TEST_LIST="${TEST_FILEBASE}.list"
gen_tests

NPROC="8"
TEST_FILEBASE="testlist_np${NPROC}"
TEST_LIST="${TEST_FILEBASE}.list"
gen_tests

NPROC="16"
TEST_FILEBASE="testlist_np${NPROC}"
TEST_LIST="${TEST_FILEBASE}.list"
gen_tests

#NPROC="24"
#TEST_FILEBASE="testlist_np${NPROC}"
#TEST_LIST="${TEST_FILEBASE}.list"
#gen_tests

############################################################################
############################################################################
############################################################################

# Now copy THIS script into the test directory.
cp ./../$0 .

# Copy the test lists into test1, test2 and test3
# First do test1
NPROC="1"
TEST_LIST="testlist_np${NPROC}.list"
cp $TEST_LIST ./$TESTFOLDER1/
cp $TEST_LIST ./$TESTFOLDER2/
cp $TEST_LIST ./$TESTFOLDER3/

NPROC="2"
TEST_LIST="testlist_np${NPROC}.list"
cp $TEST_LIST ./$TESTFOLDER1/
cp $TEST_LIST ./$TESTFOLDER2/
cp $TEST_LIST ./$TESTFOLDER3/

NPROC="4"
TEST_LIST="testlist_np${NPROC}.list"
cp $TEST_LIST ./$TESTFOLDER1/
cp $TEST_LIST ./$TESTFOLDER2/
cp $TEST_LIST ./$TESTFOLDER3/

NPROC="8"
TEST_LIST="testlist_np${NPROC}.list"
cp $TEST_LIST ./$TESTFOLDER1/
cp $TEST_LIST ./$TESTFOLDER2/
cp $TEST_LIST ./$TESTFOLDER3/

NPROC="16"
TEST_LIST="testlist_np${NPROC}.list"
cp $TEST_LIST ./$TESTFOLDER1/
cp $TEST_LIST ./$TESTFOLDER2/
cp $TEST_LIST ./$TESTFOLDER3/

#NPROC="24"
#TEST_LIST="testlist_np${NPROC}.list"
#cp $TEST_LIST ./$TESTFOLDER1/
#cp $TEST_LIST ./$TESTFOLDER2/
#cp $TEST_LIST ./$TESTFOLDER3/

############################################################################
############################################################################
############################################################################

# Now create the qsub script.
NPROC=""
function gen_qsub_file()
{
QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
TEST_LIST="testlist_np${NPROC}.list"
NUMTESTS=$(cat $TEST_LIST | wc -l)

echo '#!/bin/bash' >> $QSUBFILE
echo '#$ -S /bin/bash' >> $QSUBFILE
echo '#$ -cwd' >> $QSUBFILE
echo '#$ -V' >> $QSUBFILE
echo '#$ -pe smp.pe 24' >> $QSUBFILE
echo '#$ -l haswell' >> $QSUBFILE
#echo '#$ -l timing' >> $QSUBFILE

echo -e "\n" >> $QSUBFILE

echo "#$ -t 1-$NUMTESTS" >> $QSUBFILE

echo -e "\n" >> $QSUBFILE

echo "# This should be ran in the scratch file system." >> $QSUBFILE
echo "# Thus results directory '$RESITS_DIR' may not exist." >> $QSUBFILE
echo -e "# We create it if it does not exist.\n" >> $QSUBFILE
echo "if [ ! -d \"$RESITS_DIR\" ]; then" >> $QSUBFILE
echo "  mkdir $RESITS_DIR" >> $QSUBFILE
echo "fi" >> $QSUBFILE

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

}

############################################################################
############################################################################
############################################################################

NPROC="1"
gen_qsub_file

NPROC="2"
gen_qsub_file

NPROC="4"
gen_qsub_file

NPROC="8"
gen_qsub_file

NPROC="16"
gen_qsub_file

#NPROC="24"
#gen_qsub_file

############################################################################
############################################################################
############################################################################

# Now copy the qsub files into test1, test2 and test3
NPROC="1"
QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
cp $QSUBFILE ./$TESTFOLDER1
cp $QSUBFILE ./$TESTFOLDER2
cp $QSUBFILE ./$TESTFOLDER3

NPROC="2"
QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
cp $QSUBFILE ./$TESTFOLDER1
cp $QSUBFILE ./$TESTFOLDER2
cp $QSUBFILE ./$TESTFOLDER3

NPROC="4"
QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
cp $QSUBFILE ./$TESTFOLDER1
cp $QSUBFILE ./$TESTFOLDER2
cp $QSUBFILE ./$TESTFOLDER3

NPROC="8"
QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
cp $QSUBFILE ./$TESTFOLDER1
cp $QSUBFILE ./$TESTFOLDER2
cp $QSUBFILE ./$TESTFOLDER3

NPROC="16"
QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
cp $QSUBFILE ./$TESTFOLDER1
cp $QSUBFILE ./$TESTFOLDER2
cp $QSUBFILE ./$TESTFOLDER3

#NPROC="24"
#QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
#cp $QSUBFILE ./$TESTFOLDER1
#cp $QSUBFILE ./$TESTFOLDER2
#cp $QSUBFILE ./$TESTFOLDER3

############################################################################
############################################################################
############################################################################
# Now make the qsub_output and reslt folders within each of the three
# test folders.

QSUBOUTPUT_DIR="qsub_output"

mkdir ./$TESTFOLDER1/$QSUBOUTPUT_DIR
mkdir ./$TESTFOLDER2/$QSUBOUTPUT_DIR
mkdir ./$TESTFOLDER3/$QSUBOUTPUT_DIR

mkdir ./$TESTFOLDER1/$RESITS_DIR
mkdir ./$TESTFOLDER2/$RESITS_DIR
mkdir ./$TESTFOLDER3/$RESITS_DIR

############################################################################
############################################################################
############################################################################

################### Now check if I'm on csf, if so, delete the related scratch
# and copy the current stuff there.
if [[ $HOME == *mbax5ml3* ]]
then
  CURRENT_DIR=`pwd`
  cd ..
  PROGRAM_DIR=${PWD##*/}
  cd $CURRENT_DIR
  OOMPH_PROGRAM_DIR="/mnt/iusers01/mh01/mbax5ml3/mpi_optimized/user_drivers/$PROGRAM_DIR"
  SCRATCH_PROGRAM_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/$PROGRAM_DIR"

  OOMPH_TEST_DIR="$OOMPH_PROGRAM_DIR/$FILEBASE"
  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE"

  echo "OOMPH_TEST_DIR: $OOMPH_TEST_DIR"
  echo "SCRATCH_TEST_DIR: $SCRATCH_TEST_DIR"

  # Remove the scratch stuff.
  rm -rf $SCRATCH_TEST_DIR
  mkdir -p $SCRATCH_TEST_DIR

  ##########################################################################
  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE/$TESTFOLDER1"
  mkdir -p $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$PROGRAM $SCRATCH_TEST_DIR
  NPROC="1"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="2"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="4"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="8"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="16"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
#  NPROC="24"
#  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
#  TEST_LIST="testlist_np${NPROC}.list"
#  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
#  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR

  ##########################################################################
  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE/$TESTFOLDER2"
  mkdir -p $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$PROGRAM $SCRATCH_TEST_DIR
  NPROC="1"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="2"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="4"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="8"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="16"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
#  NPROC="24"
#  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
#  TEST_LIST="testlist_np${NPROC}.list"
#  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
#  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR

  ###
  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR

  ##########################################################################
  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE/$TESTFOLDER3"
  mkdir -p $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$PROGRAM $SCRATCH_TEST_DIR
  NPROC="1"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="2"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="4"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="8"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  NPROC="16"
  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
  TEST_LIST="testlist_np${NPROC}.list"
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
#  NPROC="24"
#  QSUBFILE="${FILEBASE}_np${NPROC}.qsub"
#  TEST_LIST="testlist_np${NPROC}.list"
#  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
#  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR

  ###
  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR
fi



















