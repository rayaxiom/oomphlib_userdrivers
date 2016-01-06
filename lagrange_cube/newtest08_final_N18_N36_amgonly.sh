#!/bin/bash

PROGRAM="lgr_cube"


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

# Declare generic params here.

# All the p params are here, I just need to set them
#PPARAM="--p_solver 96 --p_amg_str 0.668 --p_amg_damp double --p_amg_coarse int --p_amg_sim_smoo int --p_amg_com_smoo int --p_amg_iter int --p_amg_smiter int"

# Setting p param to 2D poisson per Richard p91
# --p_amg_coarse 1: RS (0 is CLJP)
# --p_amg_str 0.25
# --p_amg_sim_smoo 0: Jacobi (1 is GS)
# --p_amg_damp 0.668 (2/3)
# --p_amg_iter 2
# --p_amg_smiter 1 2XV(1,1)

PPARAM="--p_solver 1 --p_amg_iter 2 --p_amg_smiter 1 --p_amg_sim_smoo 0 --p_amg_damp 0.668  --p_amg_str 0.7 --p_amg_coarse 1 --print_p_hypre"

# Create preconditioner params for:
# LU LU
# LU ELSC
# LU ALSC SIM
# LU ALSC STR
Prec_WLu_NSLu="--w_solver 0 --ns_solver 0"
Prec_WLu_NSLSCExact="--w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 0"
Prec_WLu_NSLSCAMGSim="--w_solver 0 --ns_solver 1 $PPARAM --f_solver 1 --f_amg_iter 2 --f_amg_smiter 2 --f_amg_coarse 1 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.25 --print_f_hypre"
Prec_WLu_NSLSCAMGStr="--w_solver 0 --ns_solver 1 $PPARAM --f_solver 1 --f_amg_iter 2 --f_amg_smiter 2 --f_amg_coarse 1 --f_amg_sim_smoo 1 --f_amg_damp -1 --f_amg_str 0.75 --print_f_hypre"


PARAM="--time_type 1 --solver_type 2 --dist_prob --time_start 0.0 --time_end 1.0 --prob_id 0  --max_solver_iter 100 --itstimedir $RESITS_DIR --solver_type 2"



function gen_tests()
{
# Loop for:
# --visc 0 1

# 8090 block diagonal
# 8091 upper triangular
# 8092 lower triangular
# 8093 full AMG
# --f_solver 8090, 8091, 8092, 8093,

# --noel 4, 8, 16, 32, 64, 128

# This is set according to the list above.
# Only LU and ELSC for this test, up to Noel = 128
PRECPARAM=""
PRECLIST="3"

# Sim / Str
VISLIST="0 1"

ANGLIST="0 67"

# As per FIS p364,
REYLIST="10 100 200 500"

# Up to 
NOELLIST="18 20 22 24 26 28 30 32 34 36"

for PREC in $PRECLIST
do
for VIS in $VISLIST
do

if [ "$VIS" -eq "0" ]; then
  case "$PREC" in
    1)
      PRECPARAM="$Prec_WLu_NSLu"
      ;;
    2)
      PRECPARAM="$Prec_WLu_NSLSCExact"
      ;;
    3)
      PRECPARAM="$Prec_WLu_NSLSCAMGSim"
      ;;
  esac
else
  case "$PREC" in
    1)
      PRECPARAM="$Prec_WLu_NSLu"
      ;;
    2)
      PRECPARAM="$Prec_WLu_NSLSCExact"
      ;;
    3)
      PRECPARAM="$Prec_WLu_NSLSCAMGStr"
      ;;
  esac
fi
  for ANG in $ANGLIST
  do
  for REY in $REYLIST
  do
    for NOEL in $NOELLIST
    do
echo "mpirun -np 1 ./$PROGRAM $PARAM $PRECPARAM --ang $ANG --rey $REY --visc $VIS --noel $NOEL" >> $TEST_LIST
    done # NOEL
  done # REY
  done #ANG
done # VIS
done # PREC
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
  OOMPH_PROGRAM_DIR="/mnt/iusers01/mh01/mbax5ml3/mpi_optimized/user_drivers/$PROGRAM_DIR"
  SCRATCH_PROGRAM_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/$PROGRAM_DIR"

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
fi




