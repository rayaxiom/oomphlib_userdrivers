#!/bin/bash

# Executable
PROG="unstructured_original"

# Run command
RUNCOMMAND="mpirun -np 1"

# Command line arguments:
# NB: For convenience, you can just comment out the undesired variables
# instead of creating an empty string. However, this is bad practice!
#------------------------

# Use brick mesh? If empty, tetgen mesh is used.
# Comment out if not required.
USEBRICK="--use_brick"

# Use an iterative linear solver? If empty, SuperLU is used.
# Comment out if not required.
LINSOLVER="--use_iterative_lin_solver"

# Use trilinos? If Empty, oomph-lib's GMRES is used.
TRILINOSSOLVER="--use_trilinos"

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
RE="--re 100.0"

# Directory for DocInfo. If empty, the default is RESLT
DIR="RESLT"
DOCDIR="--doc_dir ${DIR}"

# Number for DocInfo. If empty, the default is 0
# Comment out if not required.
DOCNUM="--doc_num 0"

# Label for DocInfo. If empty, the default is fluid_soln
# Comment out if not required.
DOCLABEL="--doc_label fluid_soln"

TETGENLABEL="--tetgen_label tetgen_original/fsi_bifurcation_fluid"

#### Time stepping stuff
DO_UNSTEADY="--do_unsteady"
TSTART="--tstart 0.0"
TEND="--tend 2.0"
DT="--dt 0.04"
DO_ADAPTTIME="--do_adapt_time"
TIMETOL="--time_tol 0.0001"
TIMEPARAM="${DO_UNSTEADY} ${TSTART} ${TEND} ${DT} ${DO_ADAPTTIME} ${TIMETOL}"




# (Try go close to 500,000 for exact, and 1,000,000 for amg)
#     tet/brick
# 1 - 476/17219
# 2 - 466/17597
# 3 - 679/24901
# 4 - 1689/58793*
# 5 - 3179/104925
# 6 - 7394/235539
# 7 - 16795/502301
# 8 - 35438/1015277
# 9 - 75391*/
# 10- 154788
# 11- 317697/
# 12- 644667/
# 13- 1314335/ 
TETGENNUM="--tetgen_num 1"

# Concatenate the above parameters
PARAM="${USEBRICK}"
PARAM="${PARAM} ${TRILINOSSOLVER} ${LINSOLVER}"
PARAM="${PARAM} ${FSOLVER} ${PSOLVER} "
PARAM="${PARAM} ${VISCTERM} ${RE} "
PARAM="${PARAM} ${DOCDIR} ${DOCNUM} ${DOCLABEL}"
PARAM="${PARAM} ${TETGENLABEL} ${TETGENNUM}"
PARAM="${PARAM} ${TIMEPARAM}"

# mkdir $DIR
touch ${DIR} && rm -rf ${DIR} && mkdir ${DIR}

# Make and run
make ${PROG} && \
$RUNCOMMAND ./$PROG $PARAM



