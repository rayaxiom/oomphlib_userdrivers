#!/bin/bash

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
#PARAM="--dist_prob --prob_id 11  --max_solver_iter 300 --itstimedir $RESITS_DIR --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 96 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_coarse 1 --f_amg_sim_smoo 1 --f_amg_str 0.25 --f_amg_damp 0.1 --ang 30"
PARAM="--dist_prob --prob_id 11  --max_solver_iter 300 --itstimedir $RESITS_DIR --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --p_solver 0 --f_solver 96 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_coarse 1 --f_amg_sim_smoo 1 --ang 30"
#PARAM="--dist_prob --prob_id 11  --max_solver_iter 300 --itstimedir $RESITS_DIR --solver_type 2 --print_hypre --w_solver 0 --ns_solver 1 --p_solver 0 --ang 30"


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
PRECPARAM=""
VISLIST="1"

# As per FIS p364
REYLIST="100 200"
STRNLIST="0.66 0.661 0.662 0.663 0.664 0.665 0.666 0.667 0.668 0.669 0.67 0.671 0.672 0.673 0.674 0.675 0.676 0.677 0.678 0.679 0.68 0.681 0.682 0.683 0.684 0.685 0.686 0.687 0.688 0.689 0.69 0.691 0.692 0.693 0.694 0.695 0.696 0.697 0.698 0.699 0.7 0.701 0.702 0.703 0.704 0.705 0.706 0.707 0.708 0.709 0.71 0.711 0.712 0.713 0.714 0.715 0.716 0.717 0.718 0.719 0.72 0.721 0.722 0.723 0.724 0.725 0.726 0.727 0.728 0.729 0.73 0.731 0.732 0.733 0.734 0.735 0.736 0.737 0.738 0.739 0.74 0.741 0.742 0.743 0.744 0.745 0.746 0.747 0.748 0.749 0.75 0.751 0.752 0.753 0.754 0.755 0.756 0.757 0.758 0.759 0.76 0.761 0.762 0.763 0.764 0.765 0.766 0.767 0.768 0.769 0.77 0.771 0.772 0.773 0.774 0.775 0.776 0.777 0.778 0.779 0.78 0.781 0.782 0.783 0.784 0.785 0.786 0.787 0.788 0.789 0.79 0.791 0.792 0.793 0.794 0.795 0.796 0.797 0.798 0.799 0.8"
NOELLIST="128"

for VIS in $VISLIST
do
  for REY in $REYLIST
  do
    for STRN in $STRNLIST
    do
      for NOEL in $NOELLIST
      do
echo "mpirun -np 1 ./$PROGRAM $PARAM --rey $REY --visc $VIS --f_amg_str $STRN --noel $NOEL" >> $TEST_LIST
      done
    done
  done
done
} # End of gen_tests

TEST_FILEBASE="test_list"
TEST_LIST="$TEST_FILEBASE.list"
gen_tests

# Create a shell script for some reason...
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
echo '#$ -l highmem' >> $QSUBFILE

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

  rsync -av $OOMPH_TEST_DIR/$PROGRAM $SCRATCH_TEST_DIR/
  rsync -av $OOMPH_TEST_DIR/$QSUBFILE $SCRATCH_TEST_DIR/
  rsync -av $OOMPH_TEST_DIR/$TEST_LIST $SCRATCH_TEST_DIR/

  ## Create the res_its and qsub output directories in scratch.
  mkdir -p $SCRATCH_TEST_DIR/$RESITS_DIR
  mkdir -p $SCRATCH_TEST_DIR/$QSUBOUTPUT_DIR
fi




