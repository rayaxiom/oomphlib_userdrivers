#!/bin/bash
CURRENTDIR=`pwd`
BASE="/home/ray/oomphlib/mpi_debug_paranoid"
UDRI="/home/ray/oomphlib/mpi_debug_paranoid/user_drivers/lagrange_annular_wedge_3d"
PROGRAM="annular_wedge_threed"
RUNCOMMAND="mpirun -np 1"


function compileoomph {
CDIR=`pwd`
cd $BASE
./non_interactive_autogen.sh -j 2
cd $CDIR
}
function compileudri {
CDIR=`pwd`
cd $UDRI
make $PROGRAM
cd $CDIR
}

function runtest {
PARAM1="--time_type 1 --solver_type 2 --dist_prob --max_solver_iter 100"
PARAM2="--dt 0.01 --time_start 0.0 --time_end 0.5  --itstimedir res_iterations"
PARAM3="--visc 0 --rey 200 --prob_id 0 --noel 18"

LPREC="--w_solver 0 --ns_solver 1"
PPREC="--p_solver 1 --p_amg_iter 2 --p_amg_smiter 1 --p_amg_sim_smoo 4 --p_amg_damp 0.668 --p_amg_coarse 6 --p_amg_str 0.7 --print_p_hypre"
FPREC="--f_solver 1 --f_amg_iter 1 --f_amg_smiter 2 --f_amg_sim_smoo 4 --f_amg_damp 1 --f_amg_str 0.668 --print_f_hypre --f_amg_coarse 6"
PREC="$LPREC $PPREC $FPREC"

$RUNCOMMAND ./$PROGRAM $PARAM1 $PARAM2 $PARAM3 $PREC
}

#compileoomph && compileudri && runtest
runtest

