#! /bin/ksh

SCR_DIR=$1
ODIR1=$2
EXP1=$3
ODIR2=$4
EXP2=$5

# Expects environment variable QUICKDIFF_LINES to be set

# Tempfile
BUFFER1=`mktemp`
BUFFER2=`mktemp`
trap "rm -f $BUFFER1 $BUFFER2" 0

# comparison of outputfiles with cdo diff. Results are in ${ODIR}[12]/${EXP}[12]
OUTFILETYPE=`cat ${SCR_DIR}/outfiletype.dat`
EXPDIR1=${ODIR1}/${EXP1}
cd ${EXPDIR1}
if [ ${OUTFILETYPE} == 1 ]; then
FILES1T=`ls ${EXP1}*`
FILES1=''
for FF in ${FILES1T}; do
if [ ${FF%.codes} == ${FF} -a ${FF%.err} == ${FF} -a ${FF%.log} == ${FF} ]; then
FILES1=${FILES1}' '${FF}
fi
done
fi
if [ ${OUTFILETYPE} != 1 ]; then
FILES1=`ls ${EXP1}*.nc`
fi
NFILES1=`echo ${FILES1} |wc -w `
if [ ${NFILES1} == 0 ]; then
 echo "${EXP1} does not contain any result files in ${EXPDIR1}"
 exit
fi
EXPDIR2=${ODIR2}/${EXP2}
cd ${EXPDIR2}
if [ ${OUTFILETYPE} == 1 ]; then
FILES2T=`ls ${EXP2}*`
FILES2=''
for FF in ${FILES2T}; do
if [ ${FF%.codes} == ${FF} -a ${FF%.err} == ${FF} -a ${FF%.log} == ${FF} ]; then
FILES2=${FILES2}' '${FF}
fi
done
fi
if [ ${OUTFILETYPE} != 1 ]; then
FILES2=`ls ${EXP2}*.nc`
fi
NFILES2=`echo ${FILES2} |wc -w `
if [ ${NFILES2} == 0 ]; then
 echo "${EXP2} does not contain any result files in ${EXPDIR2}"
 exit
fi
FILES1X=''
FILES2X=''
FILESC=''
NFILESC=0
for FILE1 in ${FILES1}; do
 n=0
 for FILE2 in ${FILES2}; do
   if [ ${FILE1#$EXP1} == ${FILE2#$EXP2} ]; then
    n=$(( n + 1 ))
    FILESC=${FILESC}' '${FILE1#$EXP1}
    (( NFILESC += 1 ))
   fi
 done
 if [ $n -eq 0 ]; then
  FILES1X=${FILES1X}' '${FILE1}
 fi
done
for FILE2 in ${FILES2}; do
 n=0
 for FILE1 in ${FILES1}; do
  if [ ${FILE1#$EXP1} == ${FILE2#$EXP2} ]; then
   n=$(( n + 1 ))
  fi
 done
 if [ $n -eq 0 ]; then
  FILES2X=${FILES2X}' '${FILE2}
 fi
done
if [ `echo ${FILES1X}| wc -w ` == 0 ]; then
 echo "${EXPDIR1} contains ${NFILES1} files, all present in ${EXPDIR2} also"
else
 echo "${EXPDIR1} contains ${NFILES1} files, the following are not present in ${EXPDIR2}:"
 echo ${FILES1X}
fi
if [ `echo ${FILES2X}| wc -w ` == 0 ]; then
 echo "${EXPDIR2} contains ${NFILES2} files, all present in ${EXPDIR1} also"
else
 echo "${EXPDIR2} contains ${NFILES2} files, the following are not present in ${EXPDIR1}:"
 echo ${FILES2X}
fi
# comparison of files
DATAFILE=${EXPDIR1}/diff${FILE}_${EXP1}_${EXP2}.dat
rm -f $DATAFILE
NERROR=0
for FILE in ${FILESC}; do
 if cmp ${EXPDIR1}/${EXP1}${FILE} ${EXPDIR2}/${EXP2}${FILE} > /dev/null 2>&1
 then :
 else
     if cdo diffn ${EXPDIR1}/${EXP1}${FILE} ${EXPDIR2}/${EXP2}${FILE} 1> $BUFFER1 2> $BUFFER2
     then
         if [[ ! -s $BUFFER1 ]] || grep '^ *0 of [0-9][0-9]* records differ$' $BUFFER1 > /dev/null 2>&1
         then :
         else
	     (( NERROR += 1 ))
	     echo "> cdo diffn ${EXPDIR1}/${EXP1}${FILE} ${EXPDIR2}/${EXP2}${FILE}" | tee -a $DATAFILE
	     if [ $QUICKDIFF_LINES -gt 0 ]
             then
                 grep -v -e ' 0.0000      0.0000$' -e 'records differ' $BUFFER1 | head -n $QUICKDIFF_LINES
             fi
	     grep -v -e ' 0.0000      0.0000$' -e 'records differ' $BUFFER1 >> $DATAFILE
             grep 'records differ' $BUFFER1
             cat $BUFFER2 >> $DATAFILE
         fi
     else
	 (( NERROR += 1 ))
	 echo "> cdo diffn ${EXPDIR1}/${EXP1}${FILE} ${EXPDIR2}/${EXP2}${FILE}" | tee -a $DATAFILE
         sed 's/^/[CDO ERROR] /' $BUFFER2
         cat $BUFFER2 >> $DATAFILE
     fi
 fi
done
echo '--'
echo "$NERROR of $NFILESC files differ or cannot be compared"
