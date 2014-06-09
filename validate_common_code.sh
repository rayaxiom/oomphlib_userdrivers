## Extract the required data from validata ####################
cd $PROGRAM_DIR
cd $VALIDATA_DIR
EXTRACT_COMMAND="tar -zxf $VALIDATA_TAR "
for i in "${files[@]}"
do
  EXTRACT_COMMAND="$EXTRACT_COMMAND $i"
done

$EXTRACT_COMMAND

cd $PROGRAM_DIR && cd $VALIDATE_DIR
#######################################

touch validation.log

for i in "${files[@]}"
do
	grep "RAYITS" $TEMPRES_DIR/$i > RAYITS_new
  grep "RAYITS" ./../$VALIDATA_DIR/$i > RAYITS_old
  
  DIFF=$(diff RAYITS_new RAYITS_old)
  if [ "$DIFF" != "" ]
  then
    echo "File not the same: $i" >> validation.log
  else
    rm -rf ./../$VALIDATA_DIR/$i
  fi

  rm -rf RAYITS_new RAYITS_old
done

