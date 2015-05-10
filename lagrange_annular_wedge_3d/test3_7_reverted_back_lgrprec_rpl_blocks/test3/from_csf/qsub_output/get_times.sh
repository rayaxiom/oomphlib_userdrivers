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


#1 - Time to generate Jacobian #############################################
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time to generate Jacobian")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $9 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $9 } END { if (NR > 0) print sum / NR }'
# CHECKED

#2 - LGR: clean_up_memory: #################################################
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: clean_up_memory")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#3 - BLKSETUP: t_gen_coarsen_mapping #######################################
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP: t_gen_coarsen_mapping")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#4 - BLKSETUP: t_block_to_dof_map ##########################################
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP: t_block_to_dof_map")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#5 - BLKSETUP: t_clear_block_preconditioner_base ###########################
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP: t_clear_block_preconditioner_base")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#6 - BLKSETUP: t_initialise_index_in_dof_and_dof_num_sparse ################
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP: t_initialise_index_in_dof_and_dof_num_sparse")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#7 - BLKSETUP_MO: t_get_mesh_ndof_types ####################################
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP_MO: t_get_mesh_ndof_types")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#8 - BLKSETUP_MO_SENDSPARSE: t_create_dense_row_required ###################
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP_MO_SENDSPARSE: t_create_dense_row_required")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#9 - BLKSETUP_MO_SENDSPARSE: t_sparse_global_rows_for_block_lookup #########
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP_MO_SENDSPARSE: t_sparse_global_rows_for_block_lookup")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#10 - BLKSETUP_MO_SENDSPARSE: t_copy_to_global_index_sparse 
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP_MO_SENDSPARSE: t_copy_to_global_index_sparse")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' ##### DOING HERE?
# CHECKED

#11 - BLKSETUP_MO_SENDSPARSE: t_sending_recv_requests_sparse
#echo "#11 - BLKSETUP_MO_SENDSPARSE: t_sending_recv_requests_sparse"
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP_MO_SENDSPARSE: t_sending_recv_requests_sparse")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

#12 - BLKSETUP_MO: t_setting_and_sending_sparse_vecs
#echo "#12 - BLKSETUP_MO: t_setting_and_sending_sparse_vecs"
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "BLKSETUP_MO: t_setting_and_sending_sparse_vecs")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $5 }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $5 } END { if (NR > 0) print sum / NR }'
# CHECKED

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
echo "#25 - LSC: get block B get_B_time"
bashcol=8
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: get block B get_B_time")
echo "${OUTPUT}"
echo "${OUTPUT}" | awk -v awkcol="$bashcol" '{ print $awkcol }' ## This is to be deleted
echo "${OUTPUT}" | awk -v awkcol="$bashcol" '{ sum += $awkcol } END { if (NR > 0) print sum / NR }'

#26 - LSC: ivmm_assembly_time
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: ivmm_assembly_time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#27 - LSC: get block Bt
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: get block Bt")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#28 - LSC: t_QBt_time (matrix multiplicaton)
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: t_QBt_time (matrix multiplicaton)")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#29 - LSC: t_p_time (matrix multiplication)
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: t_p_time (matrix multiplication)")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#30 - LSC: QBt (setup MV product)
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: QBt (setup MV product)")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#31 - LSC: get_block t_get_F_time
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: get_block t_get_F_time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#32 - LSC: MV product setup t_F_MV_time
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: MV product setup t_F_MV_time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#33 - LSC: get_block t_get_Bt_time2
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: get_block t_get_Bt_time2")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#34 - LSC: MV product setup t_Bt_MV_time
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: MV product setup t_Bt_MV_time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#35 - LSC: p_prec setup time
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: p_prec setup time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#36 - LSC: f_prec setup time
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LSC: f_prec setup time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#37 - LGR: ns_setup
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: ns_setup")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#38 - LGR: delete_v_aug
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: delete_v_aug")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#39 - LGR: t_w_prec_time
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: t_w_prec_time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#40 - LGR: delete_w time
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "LGR: delete_w time")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#41 - Time for preconditioner setup
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time for preconditioner setup")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#42 - Time to generate Trilinos matrix
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time to generate Trilinos matrix")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#43 - Linear solver iterations
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Linear solver iterations")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#44 - Time for trilinos solve itself
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time for trilinos solve itself")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#45 - Time for complete trilinos solve
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time for complete trilinos solve")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'

#46 - Time for linear solver
OUTPUT=$(grep "Processor 0" ${RSLTFILE} | grep "Time for linear solver")
#echo "${OUTPUT}"
#echo "${OUTPUT}" | awk '{ print $YY }' ## This is to be deleted
#echo "${OUTPUT}" | awk '{ sum += $YY } END { if (NR > 0) print sum / NR }'













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






