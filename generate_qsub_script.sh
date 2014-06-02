#!/bin/bash

function generate_qsub_script()
{
  local TEST_LIST=$1

  ## Check if the list exists
  if [ ! -e $TEST_LIST ]
  then
    echo "No such file: $TEST_LIST"
  fi

  local FILEBASE="${TEST_LIST%.*}"

  local QSUBFILE="$FILEBASE.qsub"

  ## Check if the list exists
  if [ -e $QSUBFILE ]
  then
    echo "qsub script already exists: $QSUBFILE"
    echo "Removing it..."
    echo -e "\n"
    rm -rf $QSUBFILE
    touch $QSUBFILE
  fi

  local NUMTESTS=$(cat $TEST_LIST | wc -l)

  echo '#!/bin/bash' >> $QSUBFILE
  echo '#$ -S /bin/bash' >> $QSUBFILE
  echo '#$ -cwd' >> $QSUBFILE
  echo '#$ -V' >> $QSUBFILE

  echo -e "\n" >> $QSUBFILE

  echo "#$ -t 1-$NUMTESTS" >> $QSUBFILE

  echo -e "\n" >> $QSUBFILE

  # Do the same thing with the output directory for qsub
  local QSUBOUTPUT_DIR="qsub_output_$FILEBASE"
  if [ -d "$QSUBOUTPUT_DIR" ]; then
    echo "The directory $QSUBOUTPUT_DIR exists"
    echo "The qsub stdout and stderr files will be moved in there by the qsub script."
    echo "If a file exists, it will be overwritten, so be careful!"
    echo -e "\n"
  else
    mkdir $QSUBOUTPUT_DIR
    echo "Created directory $QSUBOUTPUT_DIR"
    echo "qsub stdout and stderr files will be moved into there by the qsub script"
    echo -e "\n"
  fi

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
  local CLEANUPLINE="mv --force $QSUBFILE"
  CLEANUPLINE+='.*.$SGE_TASK_ID '
  CLEANUPLINE+=" ./$QSUBOUTPUT_DIR/"
  echo $CLEANUPLINE >> $QSUBFILE

  echo "qsub script generated: $QSUBFILE"
  echo "Remember to create any results/solution directory!"
}


