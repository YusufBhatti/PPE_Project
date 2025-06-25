#!/bin/ksh
SCR_DIR=$1
TEST_ODIR=$2
TEST_REVISION=$3
TEST_MODEL=$4
# nproma test on test model
echo ''
echo '-------------------------------------------------------------------'
echo "Now running test model rev. ${TEST_REVISION} - test case 1: nproma=17, nproca=$NPROCA, no rerun"
echo '-------------------------------------------------------------------'
#1.1) T31L39, nproma=17, n CPUs, lrerun=.false., lco2=.true., 12 time steps
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00001rev${TEST_REVISION} ${TEST_MODEL} 17 $NPROCA .false. 11 32 .true.
#1>${SCR_DIR}/test_nproma.log 2>&1
echo ''
echo '-------------------------------------------------------------------'
echo "Now running test model rev. ${TEST_REVISION} - test case 2: nproma=17, nproca=$NPROCA, rerun"
echo '-------------------------------------------------------------------'
EXP_REMOVE=true
if [ -d  ${TEST_ODIR}/00002rev${TEST_REVISION} ]; then
EXP_REMOVE=false
fi 
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00002rev${TEST_REVISION} ${TEST_MODEL} 17 $NPROCA .false. 11 1 .true.
if [ ${EXP_REMOVE} == "true" ]; then
rm -f ${TEST_ODIR}/00002rev${TEST_REVISION}/*.err   # it is necessary to remove *.err otherwise model will not start on existing directory
fi
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00002rev${TEST_REVISION} ${TEST_MODEL} 17 $NPROCA .true. 11 1 .true.
echo ''
echo '-------------------------------------------------------------------'
echo "Comparison of results for rerun test on revision ${TEST_REVISION}"
echo '-------------------------------------------------------------------'
${SCR_DIR}/test_diff.sh ${SCR_DIR} ${TEST_ODIR} 00001rev${TEST_REVISION} ${TEST_ODIR} 00002rev${TEST_REVISION}
echo ''
echo '-------------------------------------------------------------------'
echo "Now running test model rev. ${TEST_REVISION} - test case 3: nproma=23, nproca=1, no rerun, but write restart file, lforcererun=.false."
echo '-------------------------------------------------------------------'
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00003rev${TEST_REVISION} ${TEST_MODEL} 23 1 .false. 11 1 .false.
echo ''
echo '-------------------------------------------------------------------'
echo "Comparison of results for parallelnproma test on revision ${TEST_REVISION}"
echo '-------------------------------------------------------------------'
${SCR_DIR}/test_diff.sh ${SCR_DIR} ${TEST_ODIR} 00001rev${TEST_REVISION} ${TEST_ODIR} 00003rev${TEST_REVISION}
