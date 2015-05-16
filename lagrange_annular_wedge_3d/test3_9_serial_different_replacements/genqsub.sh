#!/bin/bash

## Input file and output file respectively
############################################################################
## The testlist file, this should exist already
testlistfile="testlist.list"

## The qsub file, this output file, existing one will be rm.
qsubfile="Aw3D_N36_two_diff_replacements.qsub"

############################################################################

# Remove existing qsubfile.
rm -rf $qsubfile

## Get number of tests (number of lines)
numtests=$(cat ${testlistfile} | wc -l)


### First echo the top matter
TM="#!/bin/bash\n"
TM="${TM}#$ -S /bin/bash\n"
TM="${TM}#$ -cwd # Job runs in current directory (where you run qsub)\n"
TM="${TM}#$ -V # Job inherits environment (settings from loaded modules etc)\n"
TM="${TM}#$ -l vhighmem\n"
TM="${TM}#$ -t 1-${numtests}\n"
TM="${TM}\n"
TM="${TM}FULL_RUNCOMMAND=\`awk \"NR==\$SGE_TASK_ID\" ${testlistfile}\`\n"
TM="${TM}\$FULL_RUNCOMMAND\n"
TM="${TM}\n"
TM="${TM}mv ${qsubfile}.*.\$SGE_TASK_ID ./qsub_output/"
TM="${TM}\n"

echo -e $TM > ${qsubfile}


