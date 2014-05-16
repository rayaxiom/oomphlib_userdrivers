#!/bin/bash
CURRENTDIR=`pwd`

# --f_amg_iter is generally 1, --f_amg_smiter is generally 2
# --f_amg_str is 0.668 for stress divergence or 0.25 for simple
#
# --f_coarse:
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
#
# --f_amg_sim_smoo:
# 0 - Jacobi (Need to set damping as well) (damping needed)
# 1 - Gauss-Seidel, sequential, very slow in parallel
# 2 - GS - interior parallel, serial on boundary.
# 3 - hybrid GS or SOR, forward solve (damping needed ?)
# 4 - hybrid GS or SOR, backwards solve (damping needed ?)
# 6 - hybrid symmetric GS or SSOR. (damping needed ?)
#
# --f_amg_com_smoo:
# 6 - Schwarz
# 7 - Pilut
# 8 - ParaSails
# 9 - Euclid
#
# --f_amg_damp should only be positive for Jacobi smoother, (f_amg_sim_smoo 0)
# if we are not using Jacobi for smoothing, then the damping value is ignored.
#
#
# --f_solver 69 is equivalent to these settings:
#FAMG_ITER="--f_amg_iter 1"
#FAMG_SMITER="--f_amg_smiter 2"
#FAMG_SSMOOTHER="--f_amg_sim_smoo 1"
#FAMG_CSMOOTHER=""
#FAMG_DAMP="--f_amg_damp -1"
#FAMG_STRN="--f_amg_str 0.668"
#FAMG_COARSE="--f_amg_coarse 1"


ANGLIST="0 30 67"
NOELLIST="4 8 16 32 64 128 256"


FAMG_ITER="--f_amg_iter 1"
FAMG_SMITER="--f_amg_smiter 2"
FAMG_SSMOOTHER="--f_amg_sim_smoo 1" # GS
FAMG_CSMOOTHER=""
FAMG_DAMP="--f_amg_damp -1"
FAMG_STRN="--f_amg_str 0.25" # REMEMBER TO CHANGE THIS FOR SIMPLE/STRESS VISCOUS FORMS.
FAMG_COARSE="--f_amg_coarse 1" #RS

FAMGPARAM="--f_solver 96 $FAMG_ITER $FAMG_SMITER $FAMG_SSMOOTHER $FAMG_CSMOOTHER $FAMG_DAMP $FAMG_STRN $FAMG_COARSE"
ITSTIMEDIR="It1Sit2_GSs_RSc_Strn0.25"
touch $ITSTIMEDIR
rm -rf $ITSTIMEDIR
mkdir $ITSTIMEDIR
for ANG in $ANGLIST
do
  for NOEL in $NOELLIST
  do
  mpirun -np 1 ./sq_lgr --dist_prob --trilinos_solver --prob_id 11 --w_solver 0 --ns_solver 1 --p_solver 1 $FAMGPARAM --visc 0 --ang $ANG --rey 200 --print_hypre --noel $NOEL --itstimedir $ITSTIMEDIR
  done
done

	

FAMG_ITER="--f_amg_iter 2"
FAMG_SMITER="--f_amg_smiter 2"
FAMG_SSMOOTHER="--f_amg_sim_smoo 1" # GS
FAMG_CSMOOTHER=""
FAMG_DAMP="--f_amg_damp -1"
FAMG_STRN="--f_amg_str 0.25" # REMEMBER TO CHANGE THIS FOR SIMPLE/STRESS VISCOUS FORMS.
FAMG_COARSE="--f_amg_coarse 1" #RS

FAMGPARAM="--f_solver 96 $FAMG_ITER $FAMG_SMITER $FAMG_SSMOOTHER $FAMG_CSMOOTHER $FAMG_DAMP $FAMG_STRN $FAMG_COARSE"
ITSTIMEDIR="It2Sit2_GSs_RSc_Strn0.25"
touch $ITSTIMEDIR
rm -rf $ITSTIMEDIR
mkdir $ITSTIMEDIR
for ANG in $ANGLIST
do
  for NOEL in $NOELLIST
  do
  mpirun -np 1 ./sq_lgr --dist_prob --trilinos_solver --prob_id 11 --w_solver 0 --ns_solver 1 --p_solver 1 $FAMGPARAM --visc 0 --ang $ANG --rey 200 --print_hypre --noel $NOEL --itstimedir $ITSTIMEDIR
  done
done


FAMG_ITER="--f_amg_iter 1"
FAMG_SMITER="--f_amg_smiter 2"
FAMG_SSMOOTHER=""
FAMG_CSMOOTHER="--f_amg_com_smoo 9" # Euclid
FAMG_DAMP="--f_amg_damp -1"
FAMG_STRN="--f_amg_str 0.25" # REMEMBER TO CHANGE THIS FOR SIMPLE/STRESS VISCOUS FORMS.
FAMG_COARSE="--f_amg_coarse 1" #RS

FAMGPARAM="--f_solver 96 $FAMG_ITER $FAMG_SMITER $FAMG_SSMOOTHER $FAMG_CSMOOTHER $FAMG_DAMP $FAMG_STRN $FAMG_COARSE"
ITSTIMEDIR="It1Sit2_Euclid_RSc_Strn0.25"
touch $ITSTIMEDIR
rm -rf $ITSTIMEDIR
mkdir $ITSTIMEDIR
for ANG in $ANGLIST
do
  for NOEL in $NOELLIST
  do
  mpirun -np 1 ./sq_lgr --dist_prob --trilinos_solver --prob_id 11 --w_solver 0 --ns_solver 1 --p_solver 1 $FAMGPARAM --visc 0 --ang $ANG --rey 200 --print_hypre --noel $NOEL --itstimedir $ITSTIMEDIR
  done
done


FAMG_ITER="--f_amg_iter 1"
FAMG_SMITER="--f_amg_smiter 2"
FAMG_SSMOOTHER=""
FAMG_CSMOOTHER="--f_amg_com_smoo 9" # Euclid
FAMG_DAMP="--f_amg_damp -1"
FAMG_STRN="--f_amg_str 0.668" # REMEMBER TO CHANGE THIS FOR SIMPLE/STRESS VISCOUS FORMS.
FAMG_COARSE="--f_amg_coarse 1" #RS

FAMGPARAM="--f_solver 96 $FAMG_ITER $FAMG_SMITER $FAMG_SSMOOTHER $FAMG_CSMOOTHER $FAMG_DAMP $FAMG_STRN $FAMG_COARSE"
ITSTIMEDIR="It1Sit2_Euclid_RSc_Strn0.668"
touch $ITSTIMEDIR
rm -rf $ITSTIMEDIR
mkdir $ITSTIMEDIR
for ANG in $ANGLIST
do
  for NOEL in $NOELLIST
  do
  mpirun -np 1 ./sq_lgr --dist_prob --trilinos_solver --prob_id 11 --w_solver 0 --ns_solver 1 --p_solver 1 $FAMGPARAM --visc 1 --ang $ANG --rey 200 --print_hypre --noel $NOEL --itstimedir $ITSTIMEDIR
  done
done



