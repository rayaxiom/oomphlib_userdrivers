#!/bin/bash

# Executable
PROG="unstructured_three_d_fluid"

# Run command
RUNCOMMAND="mpirun -np 1"

# Command line arguments:
# NB: For convenience, you can just comment out the undesired variables
# instead of creating an empty string. However, this is bad practice!
#------------------------

# Use an iterative linear solver? If empty, SuperLU is used.
# Comment out if not required.
LINSOLVER="--use_iterative_lin_solver"

# Use trilinos? If Empty, oomph-lib's GMRES is used.
TRILINOSSOLVER="--use_trilinos"

# Use LSC for the (augmented) Navier-Stokes block?
# If empty, SuperLUPreconditioner is used.
# Comment out if not required.
NSSOLVER="--use_lsc"

# Use AMG for the velocity block?
# If empty, SuperLUPreconditioner is used.
# Comment out if not required.
FSOLVER="--use_amg_for_f"

# Use AMG for the pressure block?
# If empty, SuperLUPreconditioner is used.
# Comment out if not required.
PSOLVER="--use_amg_for_p"

# Use stress divergence form of the viscous term?
# If empty, the simple form of the viscous term is used.
# Comment out if not required.
VISCTERM="--use_stress_div"

# Reynolds number. If empty, the default is 100.0
# Comment out if not required.
RE="--re 50.0"

# Directory for DocInfo. If empty, the default is RESLT
DIR="RESLT"
DOCDIR="--doc_dir ${DIR}"

# Number for DocInfo. If empty, the default is 0
# Comment out if not required.
DOCNUM="--doc_num 0"

# Label for DocInfo. If empty, the default is fluid_soln
# Comment out if not required.
DOCLABEL="--doc_label fluid_soln"

# Concatenate the above parameters
PARAM="${TRILINOSSOLVER} ${LINSOLVER}"
PARAM="${PARAM} ${NSSOLVER} ${FSOLVER} ${PSOLVER} "
PARAM="${PARAM} ${VISCTERM} ${RE} "
PARAM="${PARAM} ${DOCDIR} ${DOCNUM} ${DOCLABEL}"

# mkdir $DIR
touch ${DIR} && rm -rf ${DIR} && mkdir ${DIR}

# Make and run
make ${PROG} && \
$RUNCOMMAND ./$PROG $PARAM



