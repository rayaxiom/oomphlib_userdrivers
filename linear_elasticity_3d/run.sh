#!/bin/bash

PROGRAM="periodic_load_3d"

NOEL="--noel 8"
AMG_ITER="--amg_iter 1"
AMG_SMITER="--amg_smiter 2"

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
AMG_SIM_SMOO="--amg_sim_smoo 1"

#   /// \short Complex smoothing methods used in BoomerAMG. Relaxation types
#   /// are:
#   ///  6 = Schwarz
#   ///  7 = Pilut
#   ///  8 = ParaSails
#   ///  9 = Euclid
#   /// To use these methods set AMG_using_simple_smoothing to false
AMG_COM_SMOO=""

AMG_DAMP="--amg_damp -3"
AMG_STRN="--amg_strn 0.7"

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
AMG_COARSE="--amg_coarse 1"


PARAM="$NOEL $AMG_ITER $AMG_SMITER $AMG_SIM_SMOO $AMG_COM_SMOO $AMG_DAMP "
PARAM+="$AMG_STRN $AMG_COARSE"

echo $PARAM
mpirun -np 1 ./$PROGRAM $PARAM


