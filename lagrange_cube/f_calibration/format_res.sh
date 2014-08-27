#!/bin/bash
set -u

TAG="RAYAVGAVGITS"
LINE=""

RESFILE="*N4NP*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"

RESFILE="*N6NP*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"

RESFILE="*N8NP*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"

RESFILE="*N10NP*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"

RESFILE="*N12NP*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"

RESFILE="*N14NP*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"

#RESFILE="*N16NP*"
#TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
#LINE="$LINE $TOKEN"

#RESFILE="*N18NP*"
#TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
#LINE="$LINE $TOKEN"

echo $LINE


