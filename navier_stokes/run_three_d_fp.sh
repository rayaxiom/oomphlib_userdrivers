#!/bin/bash

PROGRAM="three_d_fp"

TIME_TYPE="--time_type 0" # Steady state
SOLVER_TYPE="--solver_type 2" # direct solve
DIST_PROB="--dist_prob"
MAX_SOLVER_ITER="--max_solver_iter 100"

DT=""
TIME_START=""
TIME_END=""

DOC_SOLN="--doc_soln tmp_soln_dir"
ITSTIMEDIR="--itstimedir tmp_time_dir"

##
GPH="$TIME_TYPE $SOLVER_TYPE $DIST_PROB $MAX_SOLVER_ITER "
GPH+="$DT $TIME_START $TIME_END "
GPH+="$DOC_SOLN $ITSTIMEDIR"

#################################

VISC="--visc 1"
REY="--rey 100"
REY_START=""
REY_INCRE=""
REY_END=""

##
NSH="$VISC $REY $REY_START $REY_INCRE $REY_END"
###############################

PROB_ID="--prob_id 0"
NOEL="--noel 4"

##
PSP="$NOEL $PROB_ID"
###############################

SIGMA=""
W_SOLVER=""
NS_SOLVER=""
F_SOLVER="--f_solver 1"
P_SOLVER="--p_solver 1"

LGR_PREC="$SIGMA $W_SOLVER $NS_SOLVER $F_SOLVER $P_SOLVER"



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


F_AMG_ITER="--f_amg_iter 1"
F_AMG_SMITER="--f_amg_smiter 2"
F_AMG_SIM_SMOO="--f_amg_sim_smoo 1"
F_AMG_COM_SMOO=""
F_AMG_DAMP=""
F_AMG_STR="--f_amg_str 0.668"
F_AMG_COARSE="--f_amg_coarse 1"
F_AMG_PRINT="--print_f_hypre"

##
F_AMG_PREC="$F_AMG_ITER $F_AMG_SMITER "
F_AMG_PREC+="$F_AMG_SIM_SMOO $F_AMG_COM_SMOO $F_AMG_DAMP "
F_AMG_PREC+="$F_AMG_STR $F_AMG_COARSE $F_AMG_PRINT"


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


P_AMG_ITER="--p_amg_iter 1"
P_AMG_SMITER="--p_amg_smiter 2"
P_AMG_SIM_SMOO="--p_amg_sim_smoo 1"
P_AMG_COM_SMOO=""
P_AMG_DAMP=""
P_AMG_STR="--p_amg_str 0.25"
P_AMG_COARSE="--p_amg_coarse 0"
P_AMG_PRINT="--print_f_hypre"


##
P_AMG_PREC="$P_AMG_ITER $P_AMG_SMITER "
P_AMG_PREC+="$P_AMG_SIM_SMOO $P_AMG_COM_SMOO $P_AMG_DAMP "
P_AMG_PREC+="$P_AMG_STR $P_AMG_COARSE $P_AMG_PRINT"

PH="$LGR_PREC $F_AMG_PREC $P_AMG_PREC"



PARAM="$GPH $NSH $PSP $PH"

echo $PARAM
make $PROGRAM && mpirun -np 2 ./$PROGRAM $PARAM


