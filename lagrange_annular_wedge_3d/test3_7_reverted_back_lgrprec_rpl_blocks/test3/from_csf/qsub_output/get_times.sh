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
#awkit "" 5

# 14
#awkit "" 5

# 15
#awkit "" 5

# 16
#awkit "" 5

# 17
#awkit "" 5

# 18
#awkit "" 5

# 19
#awkit "" 5

# 20
#awkit "" 5

# 21
#awkit "" 5

# 22
#awkit "" 5

# 23
#awkit "" 5

# 24
#awkit "" 5

# 25
#awkit "" 5

# 26
#awkit "" 5

# 27
#awkit "" 5

# 28
#awkit "" 5

# 29
#awkit "" 5

# 30
#awkit "" 5

# 31
#awkit "" 5

# 32
#awkit "" 5

# 33
#awkit "" 5

# 34
#awkit "" 5

# 35
#awkit "" 5

# 36
#awkit "" 5

# 37
#awkit "" 5

# 38
#awkit "" 5

# 39
#awkit "" 5

# 40
#awkit "" 5

# 41
#awkit "" 5

# 42
#awkit "" 5

# 43
#awkit "" 5

# 44
#awkit "" 5

# 45
#awkit "" 5

# 46
#awkit "" 5










#13 - BLKSETUP_MO: t_problem_distributed
#echo "#13 - BLKSETUP_MO: t_problem_distributed"
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP_MO: t_problem_distributed")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#14 - BLKSETUP_MO: t_loop_through_meshes
#echo "#14 - BLKSETUP_MO: t_loop_through_meshes"
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP_MO: t_loop_through_meshes")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#15 - BLKSETUP_MO: t_sparse_rows_for_proc
#echo "15 - BLKSETUP_MO: t_sparse_rows_for_proc"
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP_MO: t_sparse_rows_for_proc")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#16 - BLKSETUP_MO: t_index_in_dof_and_dof_dimension
#echo "#16 - BLKSETUP_MO: t_index_in_dof_and_dof_dimension"
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP_MO: t_index_in_dof_and_dof_dimension")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#17 - BLKSETUP: t_master_prec_only
#echo "17 - BLKSETUP: t_master_prec_only"
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP: t_master_prec_only")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#18 - LGR: block_setup
#echo "18 - LGR: block_setup"
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: block_setup")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#19 - LGR: t_get_v_aug
#echo "19 - LGR: t_get_v_aug"
#OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: t_get_v_aug")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#20 - LGR: t_norm_time
#echo "#20 - LGR: t_norm_time"
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: t_norm_time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#21 - LGR: big_loop
#echo "#21 - LGR: big_loop"
#OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: big_loop")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#22 - LGR: turn_into_subsidairy
#echo "#22 - LGR: turn_into_subsidairy"
#OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: turn_into_subsidairy")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#23 - LSC: clean_up_memory_time
#echo "#23 - LSC: clean_up_memory_time"
#OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: clean_up_memory_time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'

#24 - LSC: block_setup
#echo "#24 - LSC: block_setup"
#OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: block_setup")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'

#25 - LSC: get block B get_B_time
#echo "#25 - LSC: get block B get_B_time"
#bashcol=8
#OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: get block B get_B_time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk -v awkcol="$bashcol" '{ print $awkcol }' ## This is to be deleted
#echo "${OUTPUT}" | awk -v awkcol="$bashcol" '{ sum += $awkcol } END { if (NR > 0) print sum / NR }'

#26 - LSC: ivmm_assembly_time
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: ivmm_assembly_time")
#echo "${OUTPUT}"

#27 - LSC: get block Bt
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: get block Bt")
#echo "${OUTPUT}"

#28 - LSC: t_QBt_time (matrix multiplicaton)
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: t_QBt_time (matrix multiplicaton)")
#echo "${OUTPUT}"

#29 - LSC: t_p_time (matrix multiplication)
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: t_p_time (matrix multiplication)")
#echo "${OUTPUT}"

#30 - LSC: QBt (setup MV product)
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: QBt (setup MV product)")
#echo "${OUTPUT}"

#31 - LSC: get_block t_get_F_time
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: get_block t_get_F_time")
#echo "${OUTPUT}"

#32 - LSC: MV product setup t_F_MV_time
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: MV product setup t_F_MV_time")
#echo "${OUTPUT}"

#33 - LSC: get_block t_get_Bt_time2
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: get_block t_get_Bt_time2")
#echo "${OUTPUT}"

#34 - LSC: MV product setup t_Bt_MV_time
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: MV product setup t_Bt_MV_time")
#echo "${OUTPUT}"

#35 - LSC: p_prec setup time
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: p_prec setup time")
#echo "${OUTPUT}"

#36 - LSC: f_prec setup time
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: f_prec setup time")
#echo "${OUTPUT}"

#37 - LGR: ns_setup
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: ns_setup")
#echo "${OUTPUT}"

#38 - LGR: delete_v_aug
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: delete_v_aug")
#echo "${OUTPUT}"

#39 - LGR: t_w_prec_time
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: t_w_prec_time")
#echo "${OUTPUT}"

#40 - LGR: delete_w time
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: delete_w time")
#echo "${OUTPUT}"

#41 - Time for preconditioner setup
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time for preconditioner setup")
#echo "${OUTPUT}"

#42 - Time to generate Trilinos matrix
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time to generate Trilinos matrix")
#echo "${OUTPUT}"

#43 - Linear solver iterations
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Linear solver iterations")
#echo "${OUTPUT}"

#44 - Time for trilinos solve itself
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time for trilinos solve itself")
#echo "${OUTPUT}"

#45 - Time for complete trilinos solve
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time for complete trilinos solve")
#echo "${OUTPUT}"

#46 - Time for linear solver
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time for linear solver")
#echo "${OUTPUT}"













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






