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

  ## Now, if the ./Validate/validation.log file is empty, we know that
  ## all the tests pass, thus we can safely delete all the associated files
  ## of this self test.
  if [ ! -s ./Validate/validation.log ]
  then
    rm -rf Validate validate.output
  fi

  cd $CURRENTDIR
done

cd $CURRENTDIR
cat validation.log


