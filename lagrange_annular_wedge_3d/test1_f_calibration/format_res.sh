#!/bin/bash
set -u

CURRENTDIR=`pwd`

## Change this for each folder.
RESDIR=""
PREC=""
TAG=""

function get_results()
{

############################################
cd $RESDIR


LINE=""
VISC="Sim"

RESFILE="Aw3D${VISC}R200${PREC}N4*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N6*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N8*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N10*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N12*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N14*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N16*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N18*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


echo $LINE
LINE=""

VISC="Str"
RESFILE="Aw3D${VISC}R200${PREC}N4*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N6*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N8*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N10*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N12*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N14*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N16*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


RESFILE="Aw3D${VISC}R200${PREC}N18*"
TOKEN=$(grep "$TAG" $RESFILE | awk '{print $NF}')
LINE="$LINE $TOKEN"


echo $LINE

cd $CURRENTDIR


}


#RESDIR:
#test1_2v2Strn0-75cRSGs_oomph
#test2_2v2Strn0-75cRSEuclid_oomph
#test3_2v2Strn0-25cRSGs_oomph
#test4_2v2Strn0-50cRSGs_oomph
#test5_2v2Strn0-80cRSGs_oomph
#test6_2v2Strn0-90cRSGs_oomph
#test7_2v2Strn0-668cRSGs_oomph

#PREC:
#WeNlFa2v2Strn0-75cRSGsPe
#WeNlFa2v2Strn0-75cRSEuclidPe
#WeNlFa2v2Strn0-25cRSGsPe
#WeNlFa2v2Strn0-5cRSGsPe
#WeNlFa2v2Strn0-8cRSGsPe
#WeNlFa2v2Strn0-9cRSGsPe
#WeNlFa2v2Strn0-668cRSGsPe

#TAG:
#RAYAVGAVGITS
#RAYAVGAVGPRECSETUP
#RAYAVGAVGLINSOLVER


TAG="RAYAVGAVGITS"


#RESDIR="test3_2v2Strn0-25cRSGs_oomph"
#PREC="WeNlFa2v2Strn0-25cRSGsPe"
#get_results

#RESDIR="test4_2v2Strn0-50cRSGs_oomph"
#PREC="WeNlFa2v2Strn0-5cRSGsPe"
#get_results

RESDIR="test1_2v2Strn0-75cRSGs_oomph"
PREC="WeNlFa2v2Strn0-75cRSGsPe"
get_results

RESDIR="test5_2v2Strn0-80cRSGs_oomph"
PREC="WeNlFa2v2Strn0-8cRSGsPe"
get_results

RESDIR="test6_2v2Strn0-90cRSGs_oomph"
PREC="WeNlFa2v2Strn0-9cRSGsPe"
get_results
echo ""
##########################
TAG="RAYAVGAVGPRECSETUP"


#RESDIR="test3_2v2Strn0-25cRSGs_oomph"
#PREC="WeNlFa2v2Strn0-25cRSGsPe"
#get_results

#RESDIR="test4_2v2Strn0-50cRSGs_oomph"
#PREC="WeNlFa2v2Strn0-5cRSGsPe"
#get_results

RESDIR="test1_2v2Strn0-75cRSGs_oomph"
PREC="WeNlFa2v2Strn0-75cRSGsPe"
get_results

RESDIR="test5_2v2Strn0-80cRSGs_oomph"
PREC="WeNlFa2v2Strn0-8cRSGsPe"
get_results

RESDIR="test6_2v2Strn0-90cRSGs_oomph"
PREC="WeNlFa2v2Strn0-9cRSGsPe"
get_results

echo ""
##########################
TAG="RAYAVGAVGLINSOLVER"


#RESDIR="test3_2v2Strn0-25cRSGs_oomph"
#PREC="WeNlFa2v2Strn0-25cRSGsPe"
#get_results

#RESDIR="test4_2v2Strn0-50cRSGs_oomph"
#PREC="WeNlFa2v2Strn0-5cRSGsPe"
#get_results

RESDIR="test1_2v2Strn0-75cRSGs_oomph"
PREC="WeNlFa2v2Strn0-75cRSGsPe"
get_results

RESDIR="test5_2v2Strn0-80cRSGs_oomph"
PREC="WeNlFa2v2Strn0-8cRSGsPe"
get_results

RESDIR="test6_2v2Strn0-90cRSGs_oomph"
PREC="WeNlFa2v2Strn0-9cRSGsPe"
get_results




