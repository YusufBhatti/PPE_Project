#!/bin/ksh
SCR_DIR=$1
REF_ODIR=$2
REF_REVISION=$3
REF_MODEL=$4
TEST_ODIR=$5
TEST_REVISION=$6
TEST_MODEL=$7
# reference simulation
echo ''
echo '-------------------------------------------------------------------'
echo "Now running reference model rev. ${REF_REVISION}, nproma=17, nproca=$NPROCA, no rerun"
echo '-------------------------------------------------------------------'
#1.1) T31L39, nproma=17, n CPUs, lrerun=.false., lco2=.true., 12 time steps
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} reference ${REF_ODIR} 00001rev${REF_REVISION} ${REF_MODEL} 17 $NPROCA .false. 11 32 .true.
#1>${SCR_DIR}/test_nproma.log 2>&1
# test model simulations
echo ''
echo '-------------------------------------------------------------------'
echo "Now running test model rev. ${TEST_REVISION} - test case 1: nproma=17, nproca=$NPROCA, no rerun"
echo '-------------------------------------------------------------------'
#1.1) T31L39, nproma=17, n CPUs, lrerun=.false., lco2=.true., 12 time steps
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00001rev${TEST_REVISION} ${TEST_MODEL} 17 $NPROCA .false. 11 32 .true.
#1>${SCR_DIR}/test_nproma.log 2>&1
#update test
echo ''
echo '-------------------------------------------------------------------'
echo "Comparison of results for update test on revision ${TEST_REVISION} compared to reference revision ${REF_REVISION}"
echo '-------------------------------------------------------------------'
${SCR_DIR}/test_diff.sh ${SCR_DIR} ${TEST_ODIR} 00001rev${TEST_REVISION} ${REF_ODIR} 00001rev${REF_REVISION}
