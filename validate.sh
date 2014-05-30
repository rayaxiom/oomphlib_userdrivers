#!/bin/bash

CURRENTDIR=`pwd`
folders=("lagrange_square" \
"lagrange_step")

touch validation.log
rm -rf validation.log
touch validation.log

for i in "${folders[@]}"
do
  echo "Doing $i" >> validation.log
  cd $i
  ./validate.sh 2>&1 | tee validate.output
  cat ./Validate/validation.log >> $CURRENTDIR/validation.log
  cd $CURRENTDIR
done

cd $CURRENTDIR
cat validation.log


