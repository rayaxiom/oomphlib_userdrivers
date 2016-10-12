#!/bin/bash


PROG="unstructured_three_d_fluid"
RUNCOMMAND="mpirun -np 1"

LINSOLVER="--use_trilinos"
NSSOLVER="" #"--use_lsc"
FSOLVER="" #"--use_amg_for_f"
PSOLVER="" #"--use_amg_for_p"
VISCTERM="--use_stress_div"
RE="--re 50.0"
DOCDIR="--doc_dir RESLT"
DOCNUM="--doc_num 0"
DOCLABEL="--doc_label soln"

PARAM="${LINSOLVER}"
PARAM="${PARAM} ${NSSOLVER} ${FSOLVER} ${PSOLVER} "
PARAM="${PARAM} ${VISCTERM} ${RE} "
PARAM="${PARAM} ${DOCDIR} ${DOCNUM} ${DOCLABEL}"

make ${PROG} && \
$RUNCOMMAND ./$PROG $PARAM



