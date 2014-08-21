#!/bin/bash


## These are handled by the NavierStokesProblemParameters namespace
DIST_PROB="--dist_prob"
PROB_ID="--prob_id 11" # RAYRAY SET
DOC_SOLN="--doc_soln tmp_soln" # RAYRAY SET
VISC="--visc 0" # RAYRAY SET
REY="--rey 100" # RAYRAY SET
MAX_SOLVER_ITER="--max_solver_iter 100"
ITSTIMEDIR="--itstimedir tmp_iter_times"
SOLVER_TYPE="--solver_type 2"
DT="--dt -1.0"
TIME_START="--time_start 0"
TIME_END="--time_end 1"
MESH_TYPE="--mesh_type 0"

## Concatenate the above.
NSPP="$DIST_PROB $PROB_ID $DOC_SOLN $VISC $REY $MAX_SOLVER_ITER $ITSTIMEDIR "
NSPP+="$SOLVER_TYPE $DT $TIME_START $TIME_END $MESH_TYPE"

#############################################################################
#############################################################################

## These are handled by LagrangianPreconditionerHelpers namespace
DOC_PREC="--doc_prec" # RAYRAY SET
LSC_ONLY="--lsc_only"
SIGMA="--sigma 42.0"
PRINT_HYPRE="--print_hypre"

## Settiing the W solver options ############################################
W_SOLVER_EXACT="--w_solver 0"
W_SOLVER_CG="--w_solver 1" # Using 4 iteration of CG as an inner conditioner

BDW="--bdw" # We do not use this.

# Settings for the Navier Stokes block
NS_SOLVER_EXACT="--ns_solver 0"
NS_SOLVER_LSC="--ns_solver 1"

P_SOLVER_EXACT="--p_solver 0"
P_SOLVER_AMG2D="--p_solver 1"
P_SOLVER_AMG3D="--p_solver 13"

F_SOLVER_EXACT="--f_solver 0"

F_SOLVER_AMG="--f_solver 96"

F_AMG_STRN_25="--f_amg_str 0.25"
F_AMG_STRN_668="--f_amg_str 0.668"
F_AMG_STRN_75="--f_amg_str 0.75"


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
F_AMG_COARSE_CLJP="--f_amg_coarse 0"
F_AMG_COARSE_RS="--f_amg_coarse 1"
F_AMG_COARSE_MRS="--f_amg_coarse 3"
F_AMG_COARSE_FALGOUT="--f_amg_coarse 6"
F_AMG_COARSE_PMIS="--f_amg_coarse 8"
F_AMG_COARSE_HMIS="--f_amg_coarse 10"
F_AMG_COARSE_1passRS="--f_amg_coarse 11"


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
F_AMG_SIM_SMOO_J="--f_amg_sim_smoo 0"
F_AMG_DAMP="--f_amg_damp -1.0"

F_AMG_SIM_SMOO_GS="--f_amg_sim_smoo 1"
F_AMG_SIM_SMOO_GSiterP="--f_amg_sim_smoo 2"
F_AMG_SIM_SMOO_SOR_forward="--f_amg_sim_smoo 3"
F_AMG_SIM_SMOO_SOR_backward="--f_amg_sim_smoo 4"
F_AMG_SIM_SMOO_SSOR="--f_amg_sim_smoo 6"

#   /// \short Complex smoothing methods used in BoomerAMG. Relaxation types
#   /// are:
#   ///  6 = Schwarz
#   ///  7 = Pilut
#   ///  8 = ParaSails
#   ///  9 = Euclid
#   /// To use these methods set AMG_using_simple_smoothing to false
F_AMG_COM_SMOO_SCHWARZ="--f_amg_com_smoo 6"
F_AMG_COM_SMOO_PILUT="--f_amg_com_smoo 7"
F_AMG_COM_SMOO_PARASAILS="--f_amg_com_smoo 8"
F_AMG_COM_SMOO_EUCLID="--f_amg_com_smoo 9"

#F_AMG_ITER="--f_amg_iter 1"
#F_AMG_SMITER="--f_amg_smiter 2"
F_AMG_CYCLE_1V22="--f_amg_iter 1 --f_amg_smiter 2"
F_AMG_CYCLE_2V22="--f_amg_iter 2 --f_amg_smiter 2"


## Concatenate the AMG parameters...

# Here are some preset ones:

# Full Exact
PREC_WLu_NSLu="$W_SOLVER_EXACT $NS_SOLVER_EXACT"

# LSC Exact
PREC_WLu_NSLsc_PLu_FLu="$W_SOLVER_EXACT $NS_SOLVER_LSC "
PREC_WLu_NSLsc_PLu_FLu+="$P_SOLVER_EXACT $F_SOLVER_EXACT"

# AMG for P solve
PREC_WLu_NSLsc_P2Da_FLu="$W_SOLVER_EXACT $NS_SOLVER_LSC "
PREC_WLu_NSLsc_P2Da_FLu+="$P_SOLVER_AMG2D $F_SOLVER_EXACT"

PREC_WLu_NSLsc_P3Da_FLu="$W_SOLVER_EXACT $NS_SOLVER_LSC "
PREC_WLu_NSLsc_P3Da_FLu+="$P_SOLVER_AMG3D $F_SOLVER_EXACT"

# AMG for F solve in 2D (simple form)
PREC_WLu_NSLsc_PLu_FAmg_2D_RS_GS_SIM="$W_SOLVER_EXACT $NS_SOLVER_LSC "
PREC_WLu_NSLsc_PLu_FAmg_2D_RS_GS_SIM+="$P_SOLVER_EXACT $F_SOLVER_AMG "
PREC_WLu_NSLsc_PLu_FAmg_2D_RS_GS_SIM+="$F_AMG_STRN_25 $F_AMG_COARSE_RS "
PREC_WLu_NSLsc_PLu_FAmg_2D_RS_GS_SIM+="$F_AMG_SIM_SMOO_GS --f_amg_damp -1.0 "
PREC_WLu_NSLsc_PLu_FAmg_2D_RS_GS_SIM+="$F_AMG_CYCLE_1V22 "

# and again for stress divergence form.
PREC_WLu_NSLsc_PLu_FAmg_2D_RS_GS_STR="$W_SOLVER_EXACT $NS_SOLVER_LSC "
PREC_WLu_NSLsc_PLu_FAmg_2D_RS_GS_STR+="$P_SOLVER_EXACT $F_SOLVER_AMG "
PREC_WLu_NSLsc_PLu_FAmg_2D_RS_GS_STR+="$F_AMG_STRN_668 $F_AMG_COARSE_RS "
PREC_WLu_NSLsc_PLu_FAmg_2D_RS_GS_STR+="$F_AMG_SIM_SMOO_GS --f_amg_damp -1.0 "
PREC_WLu_NSLsc_PLu_FAmg_2D_RS_GS_STR+="$F_AMG_CYCLE_1V22 "


# You get what to do...


#############################################################################
#############################################################################








