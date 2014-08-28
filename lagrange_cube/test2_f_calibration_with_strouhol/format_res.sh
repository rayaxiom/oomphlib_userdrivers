#!/bin/bash
set -u

TAG="RAYAVGAVGITS"
PREC="WeNlFa2v2Strn0-75cRSGsPe"
############################################

LINE=""
VISC="Sim"

RESFILE="CuPo${VISC}R200${PREC}A30N4*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N6*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N8*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N10*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N12*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N14*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N16*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N18*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


echo $LINE
LINE=""

VISC="Str"
RESFILE="CuPo${VISC}R200${PREC}A30N4*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N6*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N8*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N10*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N12*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N14*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N16*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="CuPo${VISC}R200${PREC}A30N18*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


echo $LINE



