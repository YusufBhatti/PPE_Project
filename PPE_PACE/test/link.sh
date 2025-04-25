#!/bin/ksh
TARGET=$1
LNAME=$2
if [ -L ${LNAME} ]; then
 rm -f ${LNAME}
fi
if [ -f ${TARGET} ]; then
 ln -s ${TARGET} ${LNAME}
else
 echo "file ${TARGET} does not exist"
fi
