#!/bin/bash

set -e


function get_validate_nlines()
{
TETNUM="$1"
NPROC="$2"
RANK="$3"

jac_val_file="jac${TETNUM}_np${NPROC}r${RANK}_val"
jac_col_file="jac${TETNUM}_np${NPROC}r${RANK}_col"
jac_row_file="jac${TETNUM}_np${NPROC}r${RANK}_row"

# Check that files are there
[ ! -f ./$jac_val_file ] && echo "File $jac_val_file not found!"
[ ! -f ./$jac_col_file ] && echo "File $jac_col_file not found!"
[ ! -f ./$jac_row_file ] && echo "File $jac_row_file not found!"

# create the nlines files
jac_val_file_nlines="${jac_val_file}_nlines"
jac_col_file_nlines="${jac_col_file}_nlines"
jac_row_file_nlines="${jac_row_file}_nlines"

wc -l < ${jac_val_file} > ${jac_val_file_nlines}
wc -l < ${jac_col_file} > ${jac_col_file_nlines}
wc -l < ${jac_row_file} > ${jac_row_file_nlines}

DIFF=$(diff $jac_val_file_nlines $jac_col_file_nlines) 
if [ "$DIFF" != "" ] 
then
    echo "$jac_val_file_nlines and $jac_col_file_nlines are not the same"
fi
}

#get_validate_nlines 12 1 0


#get_validate_nlines 1 2 0
#get_validate_nlines 1 2 1

get_validate_nlines 12 1 0

get_validate_nlines 13 2 0
get_validate_nlines 13 2 1

get_validate_nlines 14 4 0
get_validate_nlines 14 4 1
get_validate_nlines 14 4 2
get_validate_nlines 14 4 3


