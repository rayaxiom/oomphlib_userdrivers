#!/bin/bash

NPROC="16"
COARSE="1"

RSLTFILE="awlgr_ndof200000_np${NPROC}.qsub.o*.${COARSE}"
############################################################################
############################################################################

#TEMPLATE
#XX
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "XX")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

############################################################################
############################################################################

awkit()
{
grepstring=$1
bashcol=$2
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "$grepstring")
echo "${OUTPUT}"
echo "${OUTPUT}" | awk -v awkcol="$bashcol" '{ print $awkcol }'
echo "${OUTPUT}" | awk -v awkcol="$bashcol" '{ sum += $awkcol } END { if (NR > 0) print sum / NR }'
}

awkitremoves()
{
echo "Got in awkitremoves"
grepstring=$1
bashcol=$2
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "$grepstring")
echo "${OUTPUT}"
echo "${OUTPUT}" | awk -v awkcol="$bashcol" '{ gsub(/s/, "", $awkcol); print $awkcol }'
echo "${OUTPUT}" | awk -v awkcol="$bashcol" '{ gsub(/s/, "", $awkcol); sum += $awkcol } END { if (NR > 0) print sum / NR }'
}


# 1
awkit "Time to generate Jacobian" 9

# 2
awkit "LGR: clean_up_memory" 5

# 3
awkit "BLKSETUP: t_gen_coarsen_mapping" 5

# 4
awkit "BLKSETUP: t_block_to_dof_map" 5

# 5
awkit "BLKSETUP: t_clear_block_preconditioner_base" 5

# 6
awkit "BLKSETUP: t_initialise_index_in_dof_and_dof_num_sparse" 5

# 7
awkit "BLKSETUP_MO: t_get_mesh_ndof_types" 5

# 8
awkit "BLKSETUP_MO_SENDSPARSE: t_create_dense_row_required" 5

# 9
awkit "BLKSETUP_MO_SENDSPARSE: t_sparse_global_rows_for_block_lookup" 5

# 10
awkit "BLKSETUP_MO_SENDSPARSE: t_copy_to_global_index_sparse" 5

# 11
awkit "BLKSETUP_MO_SENDSPARSE: t_sending_recv_requests_sparse" 5

# 12
awkit "BLKSETUP_MO: t_setting_and_sending_sparse_vecs" 5

# 13
awkit "BLKSETUP_MO: t_problem_distributed" 5

# 14
awkit "BLKSETUP_MO: t_loop_through_meshes" 5

# 15
awkit "BLKSETUP_MO: t_sparse_rows_for_proc" 5

# 16
awkit "BLKSETUP_MO: t_index_in_dof_and_dof_dimension" 5

# 17
awkit "BLKSETUP: t_master_prec_only" 5

# 18
awkit "LGR: block_setup" 5

# 19
awkit "LGR: t_get_v_aug" 5

# 20
awkit "LGR: t_norm_time" 5

# 21
awkit "LGR: big_loop" 5

# 22
awkit "LGR: turn_into_subsidairy" 5

# 23
awkit "LSC: clean_up_memory_time" 5

# 24
awkit "LSC: block_setup" 5

# 25
awkit "LSC: get block B get_B_time" 8

# 26
awkit "LSC: ivmm_assembly_time" 5

# 27
awkit "LSC: get block Bt" 7

# 28
awkit "LSC: t_QBt_time (matrix multiplicaton)" 7

# 29
awkit "LSC: t_p_time (matrix multiplication)" 7

# 30
awkit "LSC: QBt (setup MV product)" 8

# 31
awkit "LSC: get_block t_get_F_time" 6

# 32
awkit "LSC: MV product setup t_F_MV_time" 8

# 33
awkit "LSC: get_block t_get_Bt_time2" 6

# 34
awkit "LSC: MV product setup t_Bt_MV_time" 8

# 35
awkit "LSC: p_prec setup time" 7

# 36
awkit "LSC: f_prec setup time" 7

# 37
awkit "LGR: ns_setup" 5

# 38
awkit "LGR: delete_v_aug" 5

# 39
awkit "LGR: t_w_prec_time" 5

# 40
awkit "LGR: delete_w time" 6

# 41
awkit "Time for preconditioner setup" 8

# 42
awkitremoves "Time to generate Trilinos matrix" 9

# 43
awkit "Linear solver iterations" 7

# 44
awkitremoves "Time for trilinos solve itself" 9

# 45
awkitremoves "Time for complete trilinos solve" 9

# 46
awkit "Time for linear solver" 13












#############################################################################
#OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP: t_gen_coarsen_mapping")
##echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
#############################################################################
##OUTPUT=$(grep "Processor 0" awlgr_ndof200000_np16.qsub.o*.1 | grep "LSC: f_prec setup time")
##echo "${OUTPUT}"
##echo "${OUTPUT}" | awk '{ sum += $7 } END { if (NR > 0) print sum / NR }'
#############################################################################
##Time to generate Trilinos matrix
#OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time to generate Trilinos matrix")
##echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ gsub(/s/, "", $9); sum+=$9 } END { if (NR > 0) print sum / NR }'
#############################################################################






