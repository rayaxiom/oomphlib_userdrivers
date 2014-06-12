#!/bin/bash

USER_DRIVERS_DIR=`pwd`
#folders=("lagrange_square" \
#"lagrange_step")

touch validation.log
rm -rf validation.log
touch validation.log

for d in */ ;
do
  if [[ -d "$d" && -f "$d/validate.sh" ]];
  then
    echo "$d/Validate/validate.log:" >> validation.log

    ## Go into the directory
    cd $d

    ## Run the validate script
    ./validate.sh 2>&1 | tee validate.output

    ## Now, if the ./Validate/validation.log file is empty, we know that
    ## all the tests pass, thus we can safely delete all the associated files
    ## of this self test. -s is true if file is not zero size.
    if [ -s ./Validate/validation.log ]
    then
      cat ./Validate/validation.log >> $USER_DRIVERS_DIR/validation.log
    else
      echo "OK" >> $USER_DRIVERS_DIR/validation.log
      rm -rf Validate validate.output
    fi
  
    cd $USER_DRIVERS_DIR

  fi # If the folder contains a validate.sh script.
done # for loop through folders.

cd $USER_DRIVERS_DIR
cat validation.log


