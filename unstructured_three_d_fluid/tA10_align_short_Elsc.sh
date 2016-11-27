#!/bin/bash

PROGRAM="unstructured_lgr_trac"


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


function gen_tests()
{

# Run command
RUNCOMMAND="mpirun -np 1"

USEBRICK="--use_brick"
LINSOLVER="--use_iterative_lin_solver"
TRILINOSSOLVER="--use_trilinos"
NSSOLVER="--use_lsc"
#FSOLVER="--use_amg_for_f"
#PSOLVER="--use_amg_for_p"
VISCTERM="--use_stress_div"
RE="--re 500.0"
OUTPRES="--set_outpres"
#DIR="RESLT"
#DOCDIR="--doc_dir ${DIR}"


DOCLABEL="--doc_label fluid_soln"
TETGENLABEL="--tetgen_label tetgen_original_short/fsi_bifurcation_fluid"

#### Time stepping stuff
DO_UNSTEADY="--do_unsteady"
TSTART="--tstart 0.0"
TEND="--tend 1.0"
DT="--dt 0.025"
#DO_ADAPTTIME="--do_adapt_time"
#TIMETOL="--time_tol 0.0001"
TIMEPARAM="${DO_UNSTEADY} ${TSTART} ${TEND} ${DT} ${DO_ADAPTTIME} ${TIMETOL}"



# Concatenate the above parameters
PARAM="${USEBRICK} "
PARAM="${PARAM} ${TRILINOSSOLVER} ${LINSOLVER} "
PARAM="${PARAM} ${NSSOLVER} ${FSOLVER} ${PSOLVER} "
PARAM="${PARAM} ${VISCTERM} ${RE} ${OUTPRES} "
PARAM="${PARAM} ${DOCDIR} ${DOCNUM} ${DOCLABEL} "
PARAM="${PARAM} ${TETGENLABEL} "
PARAM="${PARAM} ${TIMEPARAM} "

TETGENNUM="--tetgen_num 1"
DOCNUM="--doc_num 1"
PARAM2="${PARAM} ${TETGENNUM} ${DOCNUM}"
echo "$RUNCOMMAND ./$PROGRAM ${PARAM2}" >> ${TEST_LIST}

TETGENNUM="--tetgen_num 2"
DOCNUM="--doc_num 2"
PARAM2="${PARAM} ${TETGENNUM} ${DOCNUM}"
echo "$RUNCOMMAND ./$PROGRAM ${PARAM2}" >> ${TEST_LIST}

TETGENNUM="--tetgen_num 3"
DOCNUM="--doc_num 3"
PARAM2="${PARAM} ${TETGENNUM} ${DOCNUM}"
echo "$RUNCOMMAND ./$PROGRAM ${PARAM2}" >> ${TEST_LIST}

TETGENNUM="--tetgen_num 4"
DOCNUM="--doc_num 4"
PARAM2="${PARAM} ${TETGENNUM} ${DOCNUM}"
echo "$RUNCOMMAND ./$PROGRAM ${PARAM2}" >> ${TEST_LIST}

#TETGENNUM="--tetgen_num 5"
#DOCNUM="--doc_num 5"
#PARAM2="${PARAM} ${TETGENNUM} ${DOCNUM}"
#echo "$RUNCOMMAND ./$PROGRAM ${PARAM2}" >> ${TEST_LIST}

#TETGENNUM="--tetgen_num 6"
#DOCNUM="--doc_num 6"
#PARAM2="${PARAM} ${TETGENNUM} ${DOCNUM}"
#echo "$RUNCOMMAND ./$PROGRAM ${PARAM2}" >> ${TEST_LIST}
#
#TETGENNUM="--tetgen_num 7"
#DOCNUM="--doc_num 7"
#PARAM2="${PARAM} ${TETGENNUM} ${DOCNUM}"
#echo "$RUNCOMMAND ./$PROGRAM ${PARAM2}" >> ${TEST_LIST}
#
#TETGENNUM="--tetgen_num 8"
#DOCNUM="--doc_num 8"
#PARAM2="${PARAM} ${TETGENNUM} ${DOCNUM}"
#echo "$RUNCOMMAND ./$PROGRAM ${PARAM2}" >> ${TEST_LIST}

} # End of gen_tests

############################################################################
############################################################################
############################################################################

TEST_FILEBASE="test_list"
TEST_LIST="$TEST_FILEBASE.list"
gen_tests

# Create a shell script for some reason...
TEST_RUN="$TEST_FILEBASE.sh"
echo "#!/bin/bash" >> $TEST_RUN
cat $TEST_LIST >> $TEST_RUN

cp ./../$0 .

# RRNEW
# copy over the .qsub, .list and program file into test1, test2 and test3
cp $TEST_LIST ./$TESTFOLDER1/
cp $TEST_LIST ./$TESTFOLDER2/
cp $TEST_LIST ./$TESTFOLDER3/

# The qsub is done later and the program file is only copied into the 
# scratch directories

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
echo '#$ -l vhighmem' >> $QSUBFILE

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

############################################################################
############################################################################
############################################################################
############################################################################
cp $QSUBFILE ./$TESTFOLDER1
cp $QSUBFILE ./$TESTFOLDER2
cp $QSUBFILE ./$TESTFOLDER3

mkdir ./$TESTFOLDER1/$QSUBOUTPUT_DIR
mkdir ./$TESTFOLDER2/$QSUBOUTPUT_DIR
mkdir ./$TESTFOLDER3/$QSUBOUTPUT_DIR

mkdir ./$TESTFOLDER1/$RESITS_DIR
mkdir ./$TESTFOLDER2/$RESITS_DIR
mkdir ./$TESTFOLDER3/$RESITS_DIR

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
  OOMPH_PROGRAM_DIR="/mnt/iusers01/mh01/mbax5ml3/trunk_optimized/user_drivers/$PROGRAM_DIR"
  SCRATCH_PROGRAM_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/trunk_optimized/user_drivers/$PROGRAM_DIR"

  OOMPH_TEST_DIR="$OOMPH_PROGRAM_DIR/$FILEBASE"
  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE"

  echo "OOMPH_TEST_DIR: $OOMPH_TEST_DIR"
  echo "SCRATCH_TEST_DIR: $SCRATCH_TEST_DIR"

  # Remove the scratch stuff.
  rm -rf $SCRATCH_TEST_DIR
  mkdir -p $SCRATCH_TEST_DIR

  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE/$TESTFOLDER1"
  mkdir -p $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$PROGRAM $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR

  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE/$TESTFOLDER2"
  mkdir -p $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$PROGRAM $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR

  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE/$TESTFOLDER3"
  mkdir -p $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$PROGRAM $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR
  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR

  # symlink in folder 1
  SCRATCH_TEST_DIR="$SCRATCH_PROGRAM_DIR/$FILEBASE/$TESTFOLDER1"
  cd ${SCRATCH_TEST_DIR} 
  mkdir RESLT
  ln -s /mnt/iusers01/mh01/mbax5ml3/trunk_optimized/user_drivers/unstructured_three_d_fluid/tetgen_original_short tetgen_original_short

fi




