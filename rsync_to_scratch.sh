#!/bin/bash

UDRI="/mnt/iusers01/mh01/mbax5ml3/oomphlib_optimized/user_drivers"
SCRATCHUDRI="/mnt/iusers01/mh01/mbax5ml3/scratch/oomphlib_optimized/user_drivers"

FILELIST="scratch_file_list.txt"

rsync -av --files-from=$FILELIST $UDRI $SCRATCHUDRI



