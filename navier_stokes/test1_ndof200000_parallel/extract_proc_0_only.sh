#!/bin/bash

set -u

NPROC="1"
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.1 > ns_ndof200000_np${NPROC}.qsub.1_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.2 > ns_ndof200000_np${NPROC}.qsub.2_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.3 > ns_ndof200000_np${NPROC}.qsub.3_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.4 > ns_ndof200000_np${NPROC}.qsub.4_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.5 > ns_ndof200000_np${NPROC}.qsub.5_proc0

NPROC="2"
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.1 > ns_ndof200000_np${NPROC}.qsub.1_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.2 > ns_ndof200000_np${NPROC}.qsub.2_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.3 > ns_ndof200000_np${NPROC}.qsub.3_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.4 > ns_ndof200000_np${NPROC}.qsub.4_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.5 > ns_ndof200000_np${NPROC}.qsub.5_proc0

NPROC="4"
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.1 > ns_ndof200000_np${NPROC}.qsub.1_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.2 > ns_ndof200000_np${NPROC}.qsub.2_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.3 > ns_ndof200000_np${NPROC}.qsub.3_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.4 > ns_ndof200000_np${NPROC}.qsub.4_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.5 > ns_ndof200000_np${NPROC}.qsub.5_proc0

NPROC="8"
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.1 > ns_ndof200000_np${NPROC}.qsub.1_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.2 > ns_ndof200000_np${NPROC}.qsub.2_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.3 > ns_ndof200000_np${NPROC}.qsub.3_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.4 > ns_ndof200000_np${NPROC}.qsub.4_proc0
grep --no-filename "Processor 0:" ns_ndof200000_np${NPROC}.qsub.*.5 > ns_ndof200000_np${NPROC}.qsub.5_proc0




