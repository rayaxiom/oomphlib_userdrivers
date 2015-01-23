#!/bin/bash

WULFPATH="/home/mly/mpi_optimized/user_drivers/milan_permutation/test0_dump_matrices_wulfling/raw_lin_system"




function dostuff()
{
  FILE="$1"
  WPATH="${WULFPATH}/${FILE}"
  scpfw ${WPATH} .
}

### 12 on np1
#dostuff "jac12_np1r0"
#
#dostuff "residual12_np1r0"
#
### 13 on np2
##dostuff "jac13_np2r0"
##dostuff "jac13_np2r1"
#
#dostuff "residual13_np2r0"
#dostuff "residual13_np2r1"
#
#

## 14 on np4
dostuff "jac14_np4r0"
dostuff "jac14_np4r1"
dostuff "jac14_np4r2"
dostuff "jac14_np4r3"

dostuff "residual14_np4r0"
dostuff "residual14_np4r1"
dostuff "residual14_np4r2"
dostuff "residual14_np4r3"







