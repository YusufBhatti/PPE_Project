#!/bin/ksh
MODEL_DIR=$1
MODEL_BRANCH=$2
MODEL_REVISION=$3
MODEL_SVN=$4
MODEL_ODIR=$5
MODEL=$6
while [ ! -d $MODEL_DIR ]; do
echo '############################################################'
echo 'The '$MODEL' model will be stored in '$MODEL'_dir/'$MODEL'_branch_rev'$MODEL_REVISION
echo $MODEL'_dir is set to '$MODEL_DIR
echo 'This directory does not exist, give another directory or interrupt (0)'
echo '############################################################'
read MODEL_DIR
if [ $MODEL_DIR == 0 ]; then
exit 1
fi
done
#
REV_DIR=${MODEL_DIR}/${MODEL_BRANCH}_rev${MODEL_REVISION}
if [ ! -d ${REV_DIR} ]; then
echo '############################################################'
echo 'No checked out revision '$REV_DIR' can be found on your computer'
echo 'Try to find revision '$MODEL_REVISION ' at '$MODEL'_model URL '$MODEL_SVN
TEST_URL=`svn info $MODEL_SVN | grep URL`
TEST_URL=${TEST_URL#URL:}
while [ -z $TEST_URL ]; do
echo 'The '$MODEL'_model URL '$MODEL_SVN' does not exist'
echo 'Give a correct URL or interrupt (0)'
echo '############################################################'
read MODEL_SVN
if [ $MODEL_SVN == 0 ]; then
exit 1
else
TEST_URL=`svn info $MODEL_SVN | grep URL`
TEST_URL=${TEST_URL#URL:}
fi
done
TEST_REV=`svn info -r $MODEL_REVISION $MODEL_SVN | grep Revision:`
TEST_REV=${TEST_REV#Revision:}
while [ -z $TEST_REV ]; do
echo '############################################################'
echo 'The '$MODEL'-model revision '$MODEL_REVISION' does not exist'
LAST_REV=`svn info $MODEL_SVN | grep Rev:`
echo $LAST_REV
echo 'Give model revision here or interrupt (0)'
read MODEL_REVISION
if [ $MODEL_REVISION == 0 ]; then
exit 1
else
TEST_REV=`svn info -r $MODEL_REVISION $MODEL_SVN | grep Revision:`
TEST_REV=${TEST_REV#Revision:}
fi
done
fi
#
while [ ! -d $MODEL_ODIR ]; do
echo '############################################################'
echo 'The '$MODEL' model output will be stored in '$MODEL_ODIR
echo 'This directory does not exist, give another directory or interrupt (0)'
echo '############################################################'
read MODEL_ODIR
if [ $MODEL_ODIR == 0 ]; then
exit 1
fi
done
echo $MODEL_DIR $MODEL_BRANCH $MODEL_REVISION $MODEL_SVN $MODEL_ODIR > test_checkdirectories.dat
