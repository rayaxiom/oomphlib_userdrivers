#!/bin/bash

USER_DRIVERS_DIR=`pwd`

touch validation.log
rm -rf validation.log
touch validation.log

# Get version of oomph-lib
# The OOMPHROOT_DIR is relative to the PROGRAM_DIR
# So we have to go into PROGRAM_DIR first (to be safe).
echo "OOMPH-LIB version:" >> validation.log
cd ..
git log -1 >> $USER_DRIVERS_DIR/validation.log
echo -e "\n" >> $USER_DRIVERS_DIR/validation.log
cd $USER_DRIVERS_DIR
echo "user_drivers version:" >> validation.log
git log -1 >> $USER_DRIVERS_DIR/validation.log
echo -e "\n" >> $USER_DRIVERS_DIR/validation.log


for d in */ ;
do
  if [[ -d "$d" && -f "$d/validate.sh" ]];
  then
    echo "$d/Validate/validation.log:" >> validation.log

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
      rm -rf ./Validate
      echo "OK" >> $USER_DRIVERS_DIR/validation.log
      rm -rf Validate validate.output
    fi
  
    cd $USER_DRIVERS_DIR

  fi # If the folder contains a validate.sh script.
done # for loop through folders.

cd $USER_DRIVERS_DIR
echo -e "\n\n"
echo "========================================================================"
echo "========================================================================"
echo "cat validation.log: "
echo "========================================================================"
cat validation.log


