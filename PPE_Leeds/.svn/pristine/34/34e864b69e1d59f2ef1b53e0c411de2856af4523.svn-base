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
echo "Now running reference model rev. ${REF_REVISION} - configuration 1, nproma=17, nproca=$NPROCA, no rerun"
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
rm -f ${TEST_ODIR}/00002rev${TEST_REVISION}/*.err  # it is necessary to remove *.err otherwise model will not start on existing directory
fi
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00002rev${TEST_REVISION} ${TEST_MODEL} 17 $NPROCA .true. 11 1 .true.
#rerun
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
#parallelnproma
echo ''
echo '-------------------------------------------------------------------'
echo "Comparison of results for parallelnproma test on revision ${TEST_REVISION}"
echo '-------------------------------------------------------------------'
${SCR_DIR}/test_diff.sh ${SCR_DIR} ${TEST_ODIR} 00001rev${TEST_REVISION} ${TEST_ODIR} 00003rev${TEST_REVISION}
#submodeloff
echo ''
echo '-------------------------------------------------------------------'
echo "Now running reference model rev. ${REF_REVISION} - configuration 2: submodel off (nproma=17, $NPROCA processors), no rerun"
echo '-------------------------------------------------------------------'
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} reference ${REF_ODIR} 00002rev${REF_REVISION} ${REF_MODEL} 17 $NPROCA .false. 00 32 .true.
echo ''
echo '-------------------------------------------------------------------'
echo "Now running test model rev. ${TEST_REVISION} - test case 4: submodel off, no rerun"
echo '-------------------------------------------------------------------'
${SCR_DIR}/test_echam6_run.sh ${SCR_DIR} test ${TEST_ODIR} 00004rev${TEST_REVISION} ${TEST_MODEL} 17 1 .false. 00 32 .true.
echo '-------------------------------------------------------------------'
echo "Submodel off test: test whether program runs and finishes successfully"
echo "This test is sucessful only if \"Experiment finished\" appears underneath"
echo '-------------------------------------------------------------------'
grep "Experiment finished" ${TEST_ODIR}/00004rev${TEST_REVISION}/*err
echo ''
echo '-------------------------------------------------------------------'
echo "Comparison of results for submodel off test on revision ${TEST_REVISION}"
echo '-------------------------------------------------------------------'
${SCR_DIR}/test_diff.sh ${SCR_DIR} ${TEST_ODIR} 00002rev${REF_REVISION} ${TEST_ODIR} 00004rev${TEST_REVISION}
