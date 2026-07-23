#!/bin/ksh
SCR_DIR=$1
TEST_ODIR=$2
TEST_REVISION=$3
TEST_MODEL=$4
# nproma test on test model
echo ''
echo '-------------------------------------------------------------------'
echo "Now running test model rev. ${TEST_REVISION} - test case 1: nproma=17, nproca=$NPROCA, rerun "
echo '-------------------------------------------------------------------'
#1.1) T31L39, nproma=17, n CPUs, lrerun=.false., lco2=.true., 12 time steps
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00001rev${TEST_REVISION} ${TEST_MODEL} 17 $NPROCA .false. 11 1 .true.
if ! grep "Experiment checkpointed" ${TEST_ODIR}/00001rev${TEST_REVISION}/*err
then
  echo 'single run did not succeed'
  exit 1
fi
rm -f ${TEST_ODIR}/00001rev${TEST_REVISION}/*.err   # it is necessary to remove *.err otherwise model will not start on existing directory
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00001rev${TEST_REVISION} ${TEST_MODEL} 17 $NPROCA .true. 11 1 .true.
echo '-------------------------------------------------------------------'
grep "Experiment finished" ${TEST_ODIR}/00001rev${TEST_REVISION}/*err
