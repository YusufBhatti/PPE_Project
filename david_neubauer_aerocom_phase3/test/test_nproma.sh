#!/bin/ksh
SCR_DIR=$1
TEST_ODIR=$2
TEST_REVISION=$3
TEST_MODEL=$4
# nproma test on test model
echo ''
echo '-------------------------------------------------------------------'
echo "Now running test model rev. ${TEST_REVISION} - test case 1: nproma=17, no rerun"
echo '-------------------------------------------------------------------'
#1.1) T31L39, nproma=17, n CPUs, lrerun=.false., lco2=.true., 12 time steps
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00001rev${TEST_REVISION} ${TEST_MODEL} 17 $NPROCA .false. 11 32 .true.
#1>${SCR_DIR}/test_nproma.log 2>&1
echo ''
echo '-------------------------------------------------------------------'
echo "Now running test model rev. ${TEST_REVISION} - test case 2: nproma=23, no rerun"
echo '-------------------------------------------------------------------'
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00002rev${TEST_REVISION} ${TEST_MODEL} 23 $NPROCA .false. 11 32 .true.
echo ''
echo '-------------------------------------------------------------------'
echo "Comparison of results for nproma test on revision ${TEST_REVISION}"
echo '-------------------------------------------------------------------'
${SCR_DIR}/test_diff.sh ${SCR_DIR} ${TEST_ODIR} 00001rev${TEST_REVISION} ${TEST_ODIR} 00002rev${TEST_REVISION}
