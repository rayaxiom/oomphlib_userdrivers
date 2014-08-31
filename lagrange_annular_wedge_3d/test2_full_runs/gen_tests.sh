#!/bin/bash
 
set -u

#############################################################################
#############################################################################
#############################################################################

## Get a handle on the directories I may require
ABS_TEST_DIR=`pwd` # This is the current directory

# Go one up to get the program directory
cd ..
ABS_PROGRAM_DIR=`pwd`

# Go one up to get the user drivers directory
cd ..
ABS_USERDRIVERS_DIR=`pwd`

# Go one up to get the oomph root directory
cd ..
ABS_OOMPHROOT_DIR=`pwd`


# Now take note of the user driver revisions and oomph-lib revisions, so
# we can reproduce these results
cd $ABS_OOMPHROOT_DIR
git log -1 > $ABS_TEST_DIR/git_rev_oomphlib

cd $ABS_USERDRIVERS_DIR
git log -1 > $ABS_TEST_DIR/git_rev_user_drivers


cd $ABS_TEST_DIR
### Now continue with the rest of the script

#############################################################################
#############################################################################
#############################################################################






###
# Prec: Exact, LSC Exact, AMG
# Visc: Sim, Str
# Ang 0, 67
# Rey 100, 200, 500
#
# Noel 4, 6, 8, 10, 12, 14, 16
# Noel 4, 6, 8, 10, 12, 14, 16, 18 20, 22, 24, 26, 28, 30, 32, 34, 35
#


# Initially an empty array. Fill this up later
Files_to_copy=()

RES_DIR="res_iter_times"
mkdir -p $RES_DIR
PROGRAM="annular_wedge_threed"
TESTLIST="testlist.list"
RUN_COMMAND="mpirun -np 1"

#################
TIMETYPE="--time_type 1"
SOLVERTYPE="--solver_type 1"
DISTPROB="--dist_prob"
MAXSOLVERITER="--max_solver_iter 1000"
DT="--dt 0.01"
TIMESTART="--time_start 0.0"
TIMEEND="--time_end 1.0"
DOCSOLN=""
ITSTIMEDIR="--itstimedir $RES_DIR"

GENHELPER="$TIMETYPE $SOLVERTYPE $DISTPROB $MAXSOLVERITER "
GENHELPER+="$DT $TIMESTART $TIMEEND "
GENHELPER+="$DOCSOLN $ITSTIMEDIR"
#################


WSOLVER_EXACT="--w_solver 0"
WSOLVER_CG="--w_solver 1"

NSSOLVER_EXACT="--ns_solver 0"
NSSOLVER_LSC="--ns_solver 1"

PSOLVER_EXACT="--p_solver 0"
PSOLVER_AMG="--p_solver 1"

###### 2D poisson problem
#   // Set iterations to 1
#   hypre_preconditioner_pt->set_amg_iterations(1);
#
#   // Use simple smoother
#   hypre_preconditioner_pt->amg_using_simple_smoothing();
#   
#   // Smoother types:
#   //           0=Jacobi
#   //           1=Gauss-Seidel
#   hypre_preconditioner_pt->amg_simple_smoother() = 1;
#   
#   // AMG preconditioner
#   hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;
#   
#   // Choose strength parameter for amg
#   hypre_preconditioner_pt->amg_strength() = 0.25;
#
#   // Coarsening type
#   hypre_preconditioner_pt->amg_coarsening() = 0; 

#####  3D poisson problem
#   // Set default settings as for 2D Poisson
#   set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
#   
#   // Change strength parameter for amg
#   hypre_preconditioner_pt->amg_strength() = 0.7;
P_ITER="--p_amg_iter 1"
P_SMITER="--p_amg_smiter 2"
P_SIM_SMOO="--p_amg_sim_smoo 1" # Gauss Seidel
P_COM_SMOO=""
P_DAMP=""
P_STR="--p_amg_str 0.7"
P_COARSE="--p_amg_coarse 0" # CLJP
P_PRINT="--print_p_hypre"

P_PREC="$P_ITER $P_SMITER "
P_PREC+="$P_SIM_SMOO $P_COM_SMOO $P_DAMP "
P_PREC+="$P_STR $P_COARSE $P_PRINT"

FSOLVER_EXACT="--f_solver 0"
FSOLVER_AMG="--f_solver 1"


#################


F_ITER="--f_amg_iter 2"
F_SMITER="--f_amg_smiter 2"
F_PRINT="--print_f_hypre"

F_STRN_25="--f_amg_str 0.25"
F_STRN_75="--f_amg_str 0.75"
F_STRN_90="--f_amg_str 0.90"

#//========================================================================
#\short Default AMG coarsening strategy. Coarsening types include:
# 0 = CLJP (parallel coarsening using independent sets)
# 1 = classical RS with no boundary treatment (not recommended
#     in parallel)
# 3 = modified RS with 3rd pass to add C points on the boundaries
# 6 = Falgout (uses 1 then CLJP using interior coarse points as
#     first independent set)
# 8 = PMIS (parallel coarsening using independent sets - lower
#     complexities than 0, maybe also slower convergence)
# 10= HMIS (one pass RS on each processor then PMIS on interior
#     coarse points as first independent set)
# 11= One pass RS on each processor (not recommended)
#//========================================================================
F_COARSE_CLJP="--f_amg_coarse 0"
F_COARSE_RS="--f_amg_coarse 1"
F_COARSE_MRS="--f_amg_coarse 3"
F_COARSE_FALGOUT="--f_amg_coarse 6"
F_COARSE_PMIS="--f_amg_coarse 8"
F_COARSE_HMIS="--f_amg_coarse 10"
F_COARSE_1passRS="--f_amg_coarse 11"

#\short Simple smoothing methods used in BoomerAMG. Relaxation types
#include:
# 0 = Jacobi 
# 1 = Gauss-Seidel, sequential
#     (very slow in parallel!)
# 2 = Gauss-Seidel, interior points in parallel, boundary sequential
#     (slow in parallel!)
# 3 = hybrid Gauss-Seidel or SOR, forward solve
# 4 = hybrid Gauss-Seidel or SOR, backward solve
# 6 = hybrid symmetric Gauss-Seidel or SSOR
#To use these methods set AMG_using_simple_smoothing to true
F_SIM_SMOO_J="--f_amg_sim_smoo 0"
F_DAMP="--f_amg_damp -1.0" # REMEMBER TO SET THIS IF USING Jacobi

F_SIM_SMOO_GS="--f_amg_sim_smoo 1"
F_SIM_SMOO_GSiterP="--f_amg_sim_smoo 2"
F_SIM_SMOO_SOR_forward="--f_amg_sim_smoo 3"
F_SIM_SMOO_SOR_backward="--f_amg_sim_smoo 4"
F_SIM_SMOO_SSOR="--f_amg_sim_smoo 6"

# \short Complex smoothing methods used in BoomerAMG. Relaxation types
# are:
#  6 = Schwarz
#  7 = Pilut
#  8 = ParaSails
#  9 = Euclid
# To use these methods set AMG_using_simple_smoothing to false
F_COM_SMOO_SCHWARZ="--f_amg_com_smoo 6"
F_COM_SMOO_PILUT="--f_amg_com_smoo 7"
F_COM_SMOO_PARASAILS="--f_amg_com_smoo 8"
F_COM_SMOO_EUCLID="--f_amg_com_smoo 9"


#  void set_defaults_for_navier_stokes_momentum_block(
#   HyprePreconditioner* hypre_preconditioner_pt)
#  {
#   // Set default settings as for 2D Poisson
#   set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
#   
#   // Change smoother type:
#   //           0=Jacobi
#   //           1=Gauss-Seidel
#   hypre_preconditioner_pt->amg_simple_smoother() = 0;
#    
#   // Set smoother damping
#   hypre_preconditioner_pt->amg_damping() = 0.5;
#   
#   // Change strength parameter for amg
#   hypre_preconditioner_pt->amg_strength() = 0.75;
#  }



#F_2v2cRSGs_Strn025="$F_ITER $F_SMITER $F_PRINT "
#F_2v2cRSGs_Strn025+="$F_SIM_SMOO_GS $F_DAMP $F_COARSE_RS $F_STRN_25"

F_2v2cRSGs_Strn090="$F_ITER $F_SMITER $F_PRINT "
F_2v2cRSGs_Strn090+="$F_SIM_SMOO_GS $F_DAMP $F_COARSE_RS $F_STRN_90"

##################################################
PREC_EXACT="$WSOLVER_EXACT $NSSOLVER_EXACT"

##
PREC_LSC_EXACT="$WSOLVER_EXACT $NSSOLVER_LSC "
PREC_LSC_EXACT+="$FSOLVER_EXACT $PSOLVER_EXACT"

##
PREC_LSC_AMG_SIM="$WSOLVER_EXACT $NSSOLVER_LSC "
PREC_LSC_AMG_SIM+="$FSOLVER_AMG $PSOLVER_AMG "
PREC_LSC_AMG_SIM+="$P_PREC $F_2v2cRSGs_Strn090"

PREC_LSC_AMG_STR="$WSOLVER_EXACT $NSSOLVER_LSC "
PREC_LSC_AMG_STR+="$FSOLVER_AMG $PSOLVER_AMG "
PREC_LSC_AMG_STR+="$P_PREC $F_2v2cRSGs_Strn090"

function generate_exact_4_14()
{

PRECLIST="0 1"
VISCLIST="0 1"
REYLIST="100 200 500"
NOELLIST="4 6 8 10 12 14"

for PREC in $PRECLIST
do
  for VISC in $VISCLIST
  do
      for REY in $REYLIST
      do
        for NOEL in $NOELLIST
        do

PRECPARAM=""
if [ "$PREC" = "0" ]; then
  PRECPARAM="$PREC_EXACT"
elif [ "$PREC" = "1" ]; then
  PRECPARAM="$PREC_LSC_EXACT"
elif [ "$PREC" = "2" ]; then
  if [ "$VISC" = "0" ]; then
    PRECPARAM="$PREC_LSC_AMG_SIM"
  elif [ "$VISC" = "1" ]; then
    PRECPARAM="$PREC_LSC_AMG_STR"
  else
    PRECPARAM="NULL"
  fi
else
  PRECPARAM="NULL"
fi



NSHELPER="--visc $VISC --rey $REY"
PROBHELPER="--prob_id 0 --noel $NOEL"
PARAM="$GENHELPER $NSHELPER $PROBHELPER $PRECPARAM"
echo "$RUN_COMMAND ./$PROGRAM $PARAM" >> $TESTLIST
        done
      done
  done
done
}

function generate_amg()
{

PRECLIST="2"
VISCLIST="0 1"
REYLIST="100 200 500"
NOELLIST="4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36"

for PREC in $PRECLIST
do
  for VISC in $VISCLIST
  do
      for REY in $REYLIST
      do
        for NOEL in $NOELLIST
        do

PRECPARAM=""
if [ "$PREC" = "0" ]; then
  PRECPARAM="$PREC_EXACT"
elif [ "$PREC" = "1" ]; then
  PRECPARAM="$PREC_LSC_EXACT"
elif [ "$PREC" = "2" ]; then
  if [ "$VISC" = "0" ]; then
    PRECPARAM="$PREC_LSC_AMG_SIM"
  elif [ "$VISC" = "1" ]; then
    PRECPARAM="$PREC_LSC_AMG_STR"
  else
    PRECPARAM="NULL"
  fi
else
  PRECPARAM="NULL"
fi



NSHELPER="--visc $VISC --rey $REY"
PROBHELPER="--prob_id 0 --noel $NOEL"
PARAM="$GENHELPER $NSHELPER $PROBHELPER $PRECPARAM"
echo "$RUN_COMMAND ./$PROGRAM $PARAM" >> $TESTLIST
        done
      done
  done
done
}


function generate_exact_16()
{

PRECLIST="0 1"
VISCLIST="0 1"
REYLIST="100 200 500"
NOELLIST="16"

for PREC in $PRECLIST
do
  for VISC in $VISCLIST
  do
      for REY in $REYLIST
      do
        for NOEL in $NOELLIST
        do

PRECPARAM=""
if [ "$PREC" = "0" ]; then
  PRECPARAM="$PREC_EXACT"
elif [ "$PREC" = "1" ]; then
  PRECPARAM="$PREC_LSC_EXACT"
elif [ "$PREC" = "2" ]; then
  if [ "$VISC" = "0" ]; then
    PRECPARAM="$PREC_LSC_AMG_SIM"
  elif [ "$VISC" = "1" ]; then
    PRECPARAM="$PREC_LSC_AMG_STR"
  else
    PRECPARAM="NULL"
  fi
else
  PRECPARAM="NULL"
fi



NSHELPER="--visc $VISC --rey $REY"
PROBHELPER="--prob_id 0 --noel $NOEL"
PARAM="$GENHELPER $NSHELPER $PROBHELPER $PRECPARAM"
echo "$RUN_COMMAND ./$PROGRAM $PARAM" >> $TESTLIST
        done
      done
  done
done
}


#############################################################################
# First compile the program and move it into this folder.
cd ..
make $PROGRAM
cd $ABS_TEST_DIR
cp ./../$PROGRAM .
Files_to_copy+=($PROGRAM)
Files_to_copy+=("git_rev_oomphlib")
Files_to_copy+=("git_rev_user_drivers")
#############################################################################

rm -rf $TESTLIST

generate_amg
generate_exact_4_14
#generate_exact_16


Files_to_copy+=($TESTLIST)
Files_to_copy+=("Aw3D_full.qsub")

#############################################################################
ABS_SCRATCH_DIR="/mnt/iusers01/mh01/mbax5ml3/scratch/mpi_optimized/user_drivers/lagrange_annular_wedge_3d/test2_full_runs/"

rm -rf $ABS_SCRATCH_DIR
mkdir -p $ABS_SCRATCH_DIR

QSUB_OUT_DIR="qsub_output"

mkdir -p $ABS_SCRATCH_DIR$QSUB_OUT_DIR




for i in "${Files_to_copy[@]}"
do
  rsync -av $i $ABS_SCRATCH_DIR
done





