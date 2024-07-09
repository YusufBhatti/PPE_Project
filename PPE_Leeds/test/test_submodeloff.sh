#!/bin/ksh
SCR_DIR=$1
TEST_ODIR=$2
TEST_REVISION=$3
TEST_MODEL=$4
# submodeloff test on test model
echo ''
echo '-------------------------------------------------------------------'
echo "Now running test model rev. ${TEST_REVISION} - configuration 1: submodel=.false., no rerun"
echo '-------------------------------------------------------------------'
#1.1) T31L39, nproma=17, n CPUs, lrerun=.false., lco2=.true., 12 time steps
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00001rev${TEST_REVISION} ${TEST_MODEL} 17 $NPROCA .false. 00 32 .true.
#1>${SCR_DIR}/test_nproma.log 2>&1
echo '-------------------------------------------------------------------'
echo "Submodel off test: test whether program runs and finishes successfully"
echo "This test is sucessful only if \"Experiment finished\" appears underneath"
echo '-------------------------------------------------------------------'
grep "Experiment finished" ${TEST_ODIR}/00001rev${TEST_REVISION}/*err
