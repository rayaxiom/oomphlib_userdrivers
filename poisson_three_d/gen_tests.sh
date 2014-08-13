#!/bin/bash

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

gen_tests()
{
NPROC=$1
NPROC_STR="NP${NPROC}"

PROGRAM="poisson_3d"
AMG_STRN_LIST="0.5 0.7"
AMG_STRN_STR=""

AMG_COARSE_LIST="0 3 6 8 10"
AMG_COARSE_STR=""

OUTFILE=""

for AMG_STRN in $AMG_STRN_LIST
do
  for AMG_COARSE in $AMG_COARSE_LIST
  do

if [ "$AMG_STRN" = "0.5" ]; then
  AMG_STRN_STR="Strn05"
else
  AMG_STRN_STR="Strn07"
fi

if [ "$AMG_COARSE" = "0" ]; then
  AMG_COARSE_STR="CLJP"
elif [ "$AMG_COARSE" = "3" ]; then
  AMG_COARSE_STR="mRS"
elif [ "$AMG_COARSE" = "6" ]; then
  AMG_COARSE_STR="Falgout"
elif [ "$AMG_COARSE" = "8" ]; then
  AMG_COARSE_STR="PMIS"
else
  AMG_COARSE_STR="HMIS"
fi

OUTFILE="three_d_poisson_${AMG_STRN_STR}_${AMG_COARSE_STR}_${NPROC_STR}"

RUN_COMMAND="mpirun -np $NPROC"
PARAM="--amg_iter 1 --amg_smiter 2 --amg_sim_smoo 0 --amg_damp 0.8 --amg_strn ${AMG_STRN} --amg_coarse ${AMG_COARSE} --noel 101"

echo "$RUN_COMMAND ./$PROGRAM $PARAM > $OUTFILE 2>&1" >> $TESTLIST

 done
done

}

TESTLIST="testlist_np1"
gen_tests 1

TESTLIST="testlist_np2"
gen_tests 2


TESTLIST="testlist_np4"
gen_tests 4

TESTLIST="testlist_np8"
gen_tests 8


