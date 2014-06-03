#!/bin/bash

CURRENT_DIR=`pwd`

OOMPHOPT="/mnt/iusers01/mh01/mbax5ml3/oomphlib_optimized"
SCRATCH="/mnt/iusers01/mh01/mbax5ml3/scratch/oomphlib_optimized"

OOMPHOPTUDRI="$OOMPHOPT/user_drivers"
SCRATCHUDRI="$SCRATCH/user_drivers"

pullonly()
{
  TMPDIR=`pwd`
  cd $OOMPHOPT && \
  git fetch --all && \
  git reset --hard origin/master && \
  cd $OOMPHOPTUDRI && \
  git fetch --all && \
  git reset --hard origin/master && \
  cd $TMPDIR
}

pullonly

# Source files, INCLUDING the actual file (filename)
declare -a OOMPH_SOURCE=(\
"lagrange_step/format_test7_StPo_final.sh")

# Destination folder (without the filename)
# The source files will all be renamed "format_results.sh" in the destination
# folder.
declare -a SCRATCH_DEST=(\
"lagrange_step/test7_StPo_final")

# Check if the number of source and destination are correct (equal)
if [ "${#OOMPH_SOURCE[@]}" -ne ${#SCRATCH_DEST[@]} ]
then
  echo "Number of source and target differ"
  exit 1
fi

for ((i=0;i<${#OOMPH_SOURCE[@]};++i)); do
  cp $OOMPHOPTUDRI/${OOMPH_SOURCE[i]} $SCRATCHUDRI/${SCRATCH_DEST[i]}/format_results.sh
done


