#!/bin/bash

PROGRAM="poisson_3d"
CURRENT_DIR=`pwd`
#parallel R-S, CLJP, HIMS and PIMS



#//========================================================================
#  /// \short Default AMG coarsening strategy. Coarsening types include:
#  ///  0 = CLJP (parallel coarsening using independent sets)
#  ///  1 = classical RS with no boundary treatment (not recommended
#  ///      in parallel)
#  ///  3 = modified RS with 3rd pass to add C points on the boundaries
#  ///  6 = Falgout (uses 1 then CLJP using interior coarse points as
#  ///      first independent set)
#  ///  8 = PMIS (parallel coarsening using independent sets - lower
#  ///      complexities than 0, maybe also slower convergence)
#  ///  10= HMIS (one pass RS on each processor then PMIS on interior
#  ///      coarse points as first independent set)
#  ///  11= One pass RS on each processor (not recommended)
#//========================================================================




#   /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
#   /// include:
#   ///  0 = Jacobi 
#   ///  1 = Gauss-Seidel, sequential
#   ///      (very slow in parallel!)
#   ///  2 = Gauss-Seidel, interior points in parallel, boundary sequential
#   ///      (slow in parallel!)
#   ///  3 = hybrid Gauss-Seidel or SOR, forward solve
#   ///  4 = hybrid Gauss-Seidel or SOR, backward solve
#   ///  6 = hybrid symmetric Gauss-Seidel or SSOR
#   /// To use these methods set AMG_using_simple_smoothing to true
#
#   /// \short Complex smoothing methods used in BoomerAMG. Relaxation types
#   /// are:
#   ///  6 = Schwarz
#   ///  7 = Pilut
#   ///  8 = ParaSails
#   ///  9 = Euclid
#   /// To use these methods set AMG_using_simple_smoothing to false


# Loop through (3 times)

TESTLIST=""
RUN_COMMMAND=""

gen_tests()
{
#NPROC_STR="NP${NPROC}"

AMG_STRN_LIST="0.5 0.7"
#AMG_STRN_STR=""

AMG_COARSE_LIST="0 3 6 8 10"
#AMG_COARSE_STR=""


OUTFILE=""

for AMG_STRN in $AMG_STRN_LIST
do
  for AMG_COARSE in $AMG_COARSE_LIST
  do

#if [ "$AMG_STRN" = "0.5" ]; then
#  AMG_STRN_STR="Strn05"
#else
#  AMG_STRN_STR="Strn07"
#fi

#if [ "$AMG_COARSE" = "0" ]; then
#  AMG_COARSE_STR="CLJP"
#elif [ "$AMG_COARSE" = "3" ]; then
#  AMG_COARSE_STR="mRS"
#elif [ "$AMG_COARSE" = "6" ]; then
#  AMG_COARSE_STR="Falgout"
#elif [ "$AMG_COARSE" = "8" ]; then
#  AMG_COARSE_STR="PMIS"
#else
#  AMG_COARSE_STR="HMIS"
#fi

#OUTFILE="three_d_poisson_${AMG_STRN_STR}_${AMG_COARSE_STR}_${NPROC_STR}"

PARAM="--amg_iter 1 --amg_smiter 2 --amg_sim_smoo 0 --amg_damp 0.8 --amg_strn ${AMG_STRN} --amg_coarse ${AMG_COARSE} --noel 101"

echo "$RUN_COMMAND ./$PROGRAM $PARAM" >> $TESTLIST

 done
done

}

###############################################################################
# First compile the program and move it into this folder.
#cd ..
#make $PROGRAM
#cd $CURRENT_DIR
#cp ./../$PROGRAM .

###############################################################################
TESTLIST1="testlist_np1.list"
TESTLIST="$TESTLIST1"
rm -rf $TESTLIST
RUN_COMMAND="mpirun -np 1 taskset -c 0"
gen_tests

TESTLIST2="testlist_np2.list"
TESTLIST="$TESTLIST2"
rm -rf $TESTLIST
RUN_COMMAND="mpirun -np 2 taskset -c 0,8"
gen_tests

TESTLIST3="testlist_np4.list"
TESTLIST="$TESTLIST3"
rm -rf $TESTLIST
RUN_COMMAND="mpirun -np 4 taskset -c 0,4,8,12"
gen_tests

TESTLIST4="testlist_np8.list"
TESTLIST="$TESTLIST4"
rm -rf $TESTLIST
RUN_COMMAND="mpirun -np 8 taskset -c 0,2,4,6,8,10,12,14"
gen_tests


###############################################################################

### Copy to scratch
#SCRATCH_PATH="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/poisson_three_d/test_2_using_single_node"
#
#rsync -av $PROGRAM $SCRATCH_PATH
#
#rsync -av $TESTLIST1 $SCRATCH_PATH
#rsync -av $TESTLIST2 $SCRATCH_PATH
#rsync -av $TESTLIST3 $SCRATCH_PATH
#rsync -av $TESTLIST4 $SCRATCH_PATH
#
#rsync -av testlist_np1.qsub $SCRATCH_PATH
#rsync -av testlist_np2.qsub $SCRATCH_PATH
#rsync -av testlist_np3.qsub $SCRATCH_PATH
#rsync -av testlist_np4.qsub $SCRATCH_PATH





