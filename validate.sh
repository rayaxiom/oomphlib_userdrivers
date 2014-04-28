#!/bin/bash

CURRENTDIR=`pwd`
folders=("lagrange_square")

touch validation.log
rm -rf validation.log
touch validation.log

for i in "${folders[@]}"
do
  echo "Doing $i" >> validation.log
  cd $i
  ./validate.sh
  cat ./Validate/validation.log >> ./../validation.log
  cd ..
done

cd $CURRENTDIR
cat validation.log


