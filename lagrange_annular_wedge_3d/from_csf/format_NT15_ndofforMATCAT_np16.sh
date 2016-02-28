#!/bin/bash
# First untar strong50m.tar.gz with:
# Note: This has been tar'd with:
# 23:16:51[mbax5ml3@login1.prv.csf.compute.estate:/dev/pts/30 +1] ~/scratch/mpi_optimized/user_drivers/concatenate_matrix_without_communication 
# $ tar cvpzf weak50m.tar.gz weak50m
# So we untar it with:
# tar xpvzf strong50m.tar.gz

############################################################################
############################################################################
############################################################################

# In fact, I'll just do it for you.
TESTDIR="NT15_ndofforMATCAT_np16"

# NOTE: NT10_weak2_hashwell is also taken as the RESFILE, if not, edit it 
# below

# Check if the above directory exists.
if [ ! -d "$TESTDIR" ]; then
  # So the above directory doesn't exists, now check if the tar file exists.

  TESTTAR="${TESTDIR}.tar.gz"
  echo "${TESTDIR} does not exist, now trying to untar ${TESTTAR}"
  
  if [ ! -e $TESTTAR ]; then
    echo "${TESTTAR} does not exist, exiting..."
    exit 1
  else
    tar xpzf ${TESTTAR}
  fi
fi

############################################################################
############################################################################
############################################################################


#AWKRESFILE=""
#BCOL=""
#TAG=""
#function awkavg()
#{
### Now determine which column to average.
#GREPOUTPUT=$(grep "Processor 0" $AWKRESFILE | grep "${TAG}")
##echo "${GREPOUTPUT}"
##echo "${GREPOUTPUT}" | awk -v awkcol="${BCOL}" '{ print $awkcol }'
#
## From the above, I have determined that 7 is the correct column
#echo "${GREPOUTPUT}" | awk -v awkcol="$BCOL" '{ sum += $awkcol } END { if (NR > 0) print sum / NR }'
#}
#############################################################################
#############################################################################
#############################################################################
## First get the correct qsub directory.
#TESTDIRNUM=""
#RESFILENUM=""
#NPROC=""
#function get_average_times()
#{
#TEST1DIR="test${TESTDIRNUM}"
#QSUBDIR="${TESTDIR}/${TEST1DIR}/qsub_output/"
#
## Inside the qsub directory, the files are:
##"strong50m.qsub.o69441.1"
##"strong50m.qsub.o69441.2"
##"strong50m.qsub.o69441.3"
##"strong50m.qsub.o69441.4"
##"strong50m.qsub.o69441.5"
##"strong50m.qsub.o69441.6"
#
#RESFILE="${TESTDIR}_np${NPROC}.qsub.o*.${RESFILENUM}"
#AWKRESFILE="${QSUBDIR}/${RESFILE}"
#
#awkavg
#}
#
##function loop_through_testdir()
##{
##TESTDIRNUM="1"
##get_average_times
##TESTDIRNUM="2"
##get_average_times
##TESTDIRNUM="3"
##get_average_times
##}
#
#TESTDIRNUM="XX" # test1
#RESFILENUM="XX" # .oxxx.1
#NPROC="XX"
#
#function get_timest1t2()
#{
#TESTDIRNUM="1"
#RES1=$(get_average_times)
#echo "${RES1}"
#TESTDIRNUM="2"
#RES2=$(get_average_times)
#echo "${RES2}"
#}
#
#function get_min_avg()
#{
#RES=$(get_timest1t2)
##echo "${RES}"
#echo "$RES" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'
#}
#
#function get_all_tag_times()
#{
#TAG="RAYRAY1 Time to generate Jacobian"
#BCOL="10"
#get_min_avg
#
#TAG="LGR: clean_up_memory"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP: t_gen_coarsen_mapping"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP: t_block_to_dof_map"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP: t_clear_block_preconditioner_base"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP: t_initialise_index_in_dof_and_dof_num_sparse"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP_MO: t_get_mesh_ndof_types"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP_MO_SENDSPARSE: t_create_dense_row_required"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP_MO_SENDSPARSE: t_sparse_global_rows_for_block_lookup"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP_MO_SENDSPARSE: t_copy_to_global_index_sparse"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP_MO_SENDSPARSE: t_sending_recv_requests_sparse"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP_MO: t_setting_and_sending_sparse_vecs"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP_MO: t_problem_distributed"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP_MO: t_loop_through_meshes"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP_MO: t_sparse_rows_for_proc"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP_MO: t_index_in_dof_and_dof_dimension"
#BCOL="5"
#get_min_avg
#
#TAG="BLKSETUP: t_master_prec_only"
#BCOL="5"
#get_min_avg
#
#TAG="LGR: block_setup"
#BCOL="5"
#get_min_avg
#
#TAG="LGR: t_get_v_aug"
#BCOL="5"
#get_min_avg
#
#TAG="LGR: t_norm_time"
#BCOL="5"
#get_min_avg
#
##TAG="RAYSIGMA"
##BCOL="4"
##get_min_avg
#
#TAG="LGR: big_loop"
#BCOL="5"
#get_min_avg
#
#TAG="LGR: turn_into_subsidairy"
#BCOL="5"
#get_min_avg
#
#TAG="LSC: clean_up_memory_time"
#BCOL="5"
#get_min_avg
#
#TAG="LSC: block_setup"
#BCOL="5"
#get_min_avg
#
#TAG="LSC: get block B get_B_time"
#BCOL="8"
#get_min_avg
#
#TAG="LSC: ivmm_assembly_time"
#BCOL="5"
#get_min_avg
#
#TAG="LSC: get block Bt"
#BCOL="7"
#get_min_avg
#
#TAG="LSC: t_QBt_time (matrix multiplicaton)"
#BCOL="7"
#get_min_avg
#
#TAG="LSC: t_p_time (matrix multiplication)"
#BCOL="7"
#get_min_avg
#
#TAG="LSC: QBt (setup MV product)"
#BCOL="8"
#get_min_avg
#
#TAG="LSC: get_block t_get_F_time"
#BCOL="6"
#get_min_avg
#
#TAG="LSC: MV product setup t_F_MV_time"
#BCOL="8"
#get_min_avg
#
#TAG="LSC: get_block t_get_Bt_time2"
#BCOL="6"
#get_min_avg
#
#TAG="LSC: MV product setup t_Bt_MV_time"
#BCOL="8"
#get_min_avg
#
#TAG="LSC: p_prec setup time"
#BCOL="7"
#get_min_avg
#
#TAG="LSC: f_prec setup time"
#BCOL="7"
#get_min_avg
#
#TAG="LGR: ns_setup"
#BCOL="5"
#get_min_avg
#
#TAG="LGR: delete_v_aug"
#BCOL="5"
#get_min_avg
#
#TAG="LGR: t_w_prec_time"
#BCOL="5"
#get_min_avg
#
#TAG="LGR: delete_w time"
#BCOL="6"
#get_min_avg
#
#TAG="RAYRAY2OOMPH Time for preconditioner setup"
#BCOL="9"
#get_min_avg
#
#TAG="RAYRAY3 Time to generate Trilinos matrix"
#BCOL="10"
#get_min_avg
#
#TAG="Linear solver iterations"
#BCOL="7"
#get_min_avg
#
#TAG="RAYRAY4 Time for trilinos solve itself"
#BCOL="10"
#get_min_avg
#
#TAG="RAYRAY5 Time for complete trilinos solve"
#BCOL="10"
#get_min_avg
#
#TAG="Time for linear solver"
#BCOL="13"
#get_min_avg
#
######################################
#
#TAG="LGRSOLVE: get_block_vector w0"
#BCOL="6"
#get_min_avg
#
#TAG="LGRSOLVE: get_block_vector w1"
#BCOL="6"
#get_min_avg
#
#TAG="LGRSOLVE: solve w0"
#BCOL="6"
#get_min_avg
#
#TAG="LGRSOLVE: solve w1"
#BCOL="6"
#get_min_avg
#
#TAG="LGRSOLVE: return_block_vector w0"
#BCOL="6"
#get_min_avg
#
#TAG="LGRSOLVE: return_block_vector w1"
#BCOL="6"
#get_min_avg
#
###############################################
#TAG="LSCSOLVE: get_block_vec_p"
#BCOL="5"
#get_min_avg
#
#TAG="LSCSOLVE: p_prec_solve"
#BCOL="5"
#get_min_avg
#
#TAG="LSCSOLVE: p_prec_solve2"
#BCOL="5"
#get_min_avg
#
#TAG="LSCSOLVE: MATVEC QBt_mat_vec"
#BCOL="6"
#get_min_avg
#
#TAG="LSCSOLVE: MATVEC F_mat_vec"
#BCOL="6"
#get_min_avg
#
#TAG="LSCSOLVE: MATVEC QBt2_mat_vec"
#BCOL="6"
#get_min_avg
#
#TAG="LSCSOLVE: return_p_vec"
#BCOL="5"
#get_min_avg
#
#TAG="LSCSOLVE: MATVEC Bt_mat_vec"
#BCOL="6"
#get_min_avg
#
#TAG="LSCSOLVE: get_block_vec_f"
#BCOL="5"
#get_min_avg
#
#TAG="LSCSOLVE: F_prec_solve"
#BCOL="5"
#get_min_avg
#
#TAG="LSCSOLVE: F_prec_return_block_vec"
#BCOL="5"
#get_min_avg
#
#TAG="LSCSOLVE: total"
#BCOL="5"
#get_min_avg
#
#TAG="LGRSOLVE: total"
#BCOL="5"
#get_min_avg
#
##Processor 0:   LGRSOLVE: get_block_vector w0: 3.38554e-05
##Processor 0:   LGRSOLVE: get_block_vector w1: 7.86781e-06
##Processor 0:   LGRSOLVE: solve w0: 0.00628686
##Processor 0:   LGRSOLVE: solve w1: 4.69685e-05
##Processor 0:   LGRSOLVE: return_block_vector w0: 6.91414e-06
##Processor 0:   LGRSOLVE: return_block_vector w1: 2.14577e-06
#
##Processor 0:   LSCSOLVE: get_block_vec_p: 1.71661e-05
##Processor 0:   LSCSOLVE: p_prec_solve: 0.0126171
##Processor 0:   LSCSOLVE: p_prec_solve2: 0.00732207
##Processor 0:   LSCSOLVE: MATVEC QBt_mat_vec: 0.00179505
##Processor 0:   LSCSOLVE: MATVEC F_mat_vec: 0.00602794
##Processor 0:   LSCSOLVE: MATVEC QBt2_mat_vec: 0.00165009
##Processor 0:   LSCSOLVE: return_p_vec: 3.69549e-05
##Processor 0:   LSCSOLVE: MATVEC Bt_mat_vec: 0.00167799
##Processor 0:   LSCSOLVE: get_block_vec_f: 0.00560594
##Processor 0:   LSCSOLVE: F_prec_solve: 0.148047
##Processor 0:   LSCSOLVE: F_prec_return_block_vec: 0.000440121
##Processor 0:   LSCSOLVE: total: 0.185418
##Processor 0:   LGRSOLVE: total: 0.194974
#
#}
#
#
#RESFILENUM="5" # .oxxx.1
#NPROC="1"
#RESNPROC1=$(get_all_tag_times)
#NPROC="2"
#RESNPROC2=$(get_all_tag_times)
#NPROC="4"
#RESNPROC4=$(get_all_tag_times)
#NPROC="8"
#RESNPROC8=$(get_all_tag_times)
#NPROC="16"
#RESNPROC16=$(get_all_tag_times)
#paste <(echo "$RESNPROC1") <(echo "$RESNPROC2") <(echo "$RESNPROC4") <(echo "$RESNPROC8") <(echo "$RESNPROC16") > ${TESTDIR}C${RESFILENUM}









#Processor 0:   RAYRAY6 Time to build epetra matrix [sec] : 0.0328221

#Processor 0:   RAYRAY6 Time to build epetra matrix [sec] : 0.155656

#Processor 0:   RAYRAY6 Time to build epetra matrix [sec] : 0.0263181
#Processor 0:   Setting up BoomerAMG, Processor 0:   time for setup [s] : 0.041656
#Processor 0:   Setting up BoomerAMG, Processor 0:   time for setup [s] : 0.622228
#Processor 0:   Setting up SuperLU (exact) preconditioner
#Processor 0:   Setting up SuperLU (exact) preconditioner
#Processor 0:   RAYRAY4 Time for trilinos solve itself                 : 6.45478s
#Processor 0:   RAYRAY5 Time for complete trilinos solve                  : 9.29377s
#Processor 0:   Time for linear solver ( ndof = 70717 ) [sec]: 16.6274















#RESFILENUM="1"
#OUTPUT2=$(loop_through_testdir)
#echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'

#RESFILENUM="2"
#OUTPUT2=$(loop_through_testdir)
#echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'

#RESFILENUM="3"
#OUTPUT2=$(loop_through_testdir)
#echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'

#RESFILENUM="4"
#OUTPUT2=$(loop_through_testdir)
#echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'

#RESFILENUM="5"
#OUTPUT2=$(loop_through_testdir)
#echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'

#RESFILENUM="6"
#OUTPUT2=$(loop_through_testdir)
#echo "$OUTPUT2" | awk 'NR == 1 || $1 < min {min = $1}END{print min}'




############################################################################
############################################################################
############################################################################



############################################################################
############################################################################
############################################################################

