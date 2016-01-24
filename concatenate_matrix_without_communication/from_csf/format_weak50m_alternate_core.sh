#!/bin/bash
# First untar strong50m.tar.gz with:
# Note: This has been tar'd with:
# 23:16:51[mbax5ml3@login1.prv.csf.compute.estate:/dev/pts/30 +1] ~/scratch/mpi_optimized/user_drivers/concatenate_matrix_without_communication 
# $ tar cvpzf weak50m.tar.gz weak50m
# So we untar it with:
# tar xpvzf strong50m.tar.gz

############################################################################
############################################################################
############################################################################

# In fact, I'll just do it for you.
TESTDIR="weak50m_alternate_core"

# Check if the above directory exists.
if [ ! -d "$TESTDIR" ]; then
  # So the above directory doesn't exists, now check if the tar file exists.

  TESTTAR="${TESTDIR}.tar.gz"
  echo "${TESTDIR} does not exist, now trying to untar ${TESTTAR}"
  
  if [ ! -e $TESTTAR ]; then
    echo "${TESTTAR} does not exist, exiting..."
    exit 1
  else
    tar xpzf ${TESTTAR}
  fi
fi

############################################################################
############################################################################
############################################################################
# First get the correct qsub directory.

TESTDIRNUM=""
RESFILENUM=""
function get_average_times()
{
TEST1DIR="test${TESTDIRNUM}"
QSUBDIR="${TESTDIR}/${TEST1DIR}/qsub_output"

# Inside the qsub directory, the files are:
#"strong50m.qsub.o69441.1"
#"strong50m.qsub.o69441.2"
#"strong50m.qsub.o69441.3"
#"strong50m.qsub.o69441.4"
#"strong50m.qsub.o69441.5"
#"strong50m.qsub.o69441.6"

RESFILE="weak50m_alternate_core.qsub.o*.${RESFILENUM}"
ABSRESFILE="${QSUBDIR}/${RESFILE}"

## Now determine which column to average.
TAG="Time to cat"
GREPOUTPUT=$(grep "${TAG}" $ABSRESFILE)
#echo "${GREPOUTPUT}"
#echo "${GREPOUTPUT}" | awk '{ print $7 }'

# From the above, I have determined that 7 is the correct column
bashcol="7"
echo "${GREPOUTPUT}" | awk -v awkcol="$bashcol" '{ sum += $awkcol } END { if (NR > 0) print sum / NR }'
}

function loop_through_testdir()
{
TESTDIRNUM="1"
get_average_times
TESTDIRNUM="2"
get_average_times
TESTDIRNUM="3"
get_average_times
}


RESFILENUM="1"
OUTPUT2=$(loop_through_testdir)
echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'

RESFILENUM="2"
OUTPUT2=$(loop_through_testdir)
echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'

RESFILENUM="3"
OUTPUT2=$(loop_through_testdir)
echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'

RESFILENUM="4"
OUTPUT2=$(loop_through_testdir)
echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'

RESFILENUM="5"
OUTPUT2=$(loop_through_testdir)
echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'

RESFILENUM="6"
OUTPUT2=$(loop_through_testdir)
echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'




############################################################################
############################################################################
############################################################################



############################################################################
############################################################################
############################################################################

