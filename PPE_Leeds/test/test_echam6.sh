#!/bin/ksh
BINDIR=`dirname $0`
SCR_DIR=`cd $BINDIR && pwd` # directory with scripts
OUTFILETYPE=2               # 1=GRIB, 2=netcdf
FORTRANCOMPILER=nag/5.3.928
CCOMPILER=gcc/4.6.3
MPI_MODULE=mpich2/1.4.1p1-static-nag53
# Command line template for running MPI processes.
# When running, %n is replaced by number of processes, %x by name of executable
MPIRUN='mpiexec -n %n %x'
# Use this setting on p249
#MPIRUN="yes `hostname` | head -n %n > host.list && poe %x -procs %n"
#
##########    Test model ####################################
#
# set base directory, branch, and revision for test model
# source code of test model will be put into 
# ${TEST_DIR}/${TEST_BRANCH}_rev${TEST_REVISION}
# if the model already exists, this is used (and not overwritten).
#
TEST_DIR=/scratch/local2/$USER/echam6/
TEST_BRANCH=echam-dev
TEST_REVISION=3467
# source code location in svn system:
TEST_SVN=https://svn.zmaw.de/svn/echam6/trunk/${TEST_BRANCH}
#
# set base directory for test model output
#
TEST_ODIR=/scratch/local2/$USER/data/
#
# set whether compilation is forced (LCOMP=.true.) 
# or only done if no executable exists (LCOMP=.false.)
#
###########   Reference model ################################
#
# set base directory, branch, and revision for reference model
# source code of reference model will be put into 
# ${REF_DIR}/${REF_BRANCH}_rev${REF_REVISION}
# if the model already exists, this is used (and not overwritten).
#
REF_DIR=/scratch/local2/$USER/echam6/
REF_BRANCH=echam-dev
REF_REVISION=3289
# source code location in svn system:
REF_SVN=https://svn.zmaw.de/svn/echam6/trunk/${REF_BRANCH}
#
# set base directory for reference model output
#
REF_ODIR=/scratch/local2/$USER/data/
#
##########   Test mode       #################################
#
# set a test mode (only test mode 'all' compares reference and test version)
#MODE='compile'    # only compilation
MODE='single'     # a single model run on two processors 
#MODE='debug'      # single run on single processor
#MODE='parallel'   # parallel test: nproca=nprocb=1 versus nproca=2, nprocb=1
#MODE='nproma'     # nproma test: nproma=23 versus nproma=17
#MODE='rerun'      # rerun test: uninterrupted run against rerun
#MODE='update'     # update test: test two revisions against each other
#MODE='submodeloff' # submodel off test: does it run without submodels?
#MODE='parallelnproma' # parallel and nproma test
#MODE='parallelnpromarerun' # parallel, nproma, and rerun test
#MODE='parallelnpromarerunsubmodeloff' # parallel, nproma, rerun, and submodel
#MODE='all' # parallel, nproma, rerun, submodel off and update test
#
##########   Compilation     #################################
# Read mode from command line, if given
MODE=${1:-$MODE}
# Read revision of test model from command line when given
TEST_REVISION=${2:-$TEST_REVISION}
# Read revision of reference  model from command line when given
REF_REVISION=${3:-$REF_REVISION}
#
echo '############################################################'
echo "             test mode is: ${MODE}                          "
echo '############################################################'
##########   Compilation     #################################
#
#LCOMP=.true.
LCOMP=.false.
#
##########   Run-time        #################################
#
NPROCA=${NPROCA:-2} # Default may be overriden by environment variable
NPROCB=${NPROCB:-1} # Default may be overriden by environment variable
#
##########   Comparison      #################################
#
#-----------------------------------------------------------------------
QUICKDIFF_LINES=${QUICKDIFF_LINES:-0} # Default overriden by environment
#
# Check directories for test model
${SCR_DIR}/test_checkdirectories.sh $TEST_DIR $TEST_BRANCH $TEST_REVISION $TEST_SVN $TEST_ODIR test
read TEST_DIR TEST_BRANCH TEST_REVISION TEST_SVN TEST_ODIR < test_checkdirectories.dat
rm -f test_checkdirectories.dat
# Check directories for reference model
if [ $MODE == 'update' -o $MODE == 'all' ]; then
${SCR_DIR}/test_checkdirectories.sh $REF_DIR $REF_BRANCH $REF_REVISION $REF_SVN $REF_ODIR reference
read REF_DIR REF_BRANCH REF_REVISION REF_SVN REF_ODIR < test_checkdirectories.dat
rm -f test_checkdirectories.dat
fi
# write outfiletype into file
cat > ${SCR_DIR}/outfiletype.dat <<EOF
$OUTFILETYPE
EOF
# write fortran compiler version into file
cat > ${SCR_DIR}/fortran.dat <<EOF
$FORTRANCOMPILER
EOF
# write C compiler version into file
cat > ${SCR_DIR}/c.dat <<EOF
$CCOMPILER
EOF
# store MPI module in environment
export MPI_MODULE
# write mpirun command into file
cat > ${SCR_DIR}/mpirun.dat <<EOF
$MPIRUN
EOF
# get current directory for protocol files
CUR_DIR=`pwd`
# get source of reference model and compile it 
if [ $MODE == 'all' -o $MODE == 'update' ]; then
${SCR_DIR}/compile_echam6.sh reference $LCOMP $REF_DIR $REF_BRANCH $REF_REVISION $REF_SVN $SCR_DIR
STATUS=$?
if [ ${STATUS} != 0 ]; then
 echo 'test_echam6.sh: compilation of '$REF_BRANCH' failed'
 exit
fi
REF_MODEL=${REF_DIR}/${REF_BRANCH}_rev${REF_REVISION}/${REF_BRANCH}/bin/echam6
fi
# get source of test model and compile it
 ${SCR_DIR}/compile_echam6.sh test $LCOMP $TEST_DIR $TEST_BRANCH $TEST_REVISION $TEST_SVN $SCR_DIR
STATUS=$?
if [ ${STATUS} != 0 ]; then
 echo 'test_echam6.sh: compilation of '$TEST_BRANCH' failed'
 exit
fi
 TEST_MODEL=${TEST_DIR}/${TEST_BRANCH}_rev${TEST_REVISION}/${TEST_BRANCH}/bin/echam6

# set run-time environment

export NPROCA NPROCB
export QUICKDIFF_LINES

# run tests

case "$MODE" in
debug)
 # Run explicitly with a single process
 NPROCA=1 NPROCB=1 ${SCR_DIR}/test_single.sh ${SCR_DIR} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL}
;;
single)
 ${SCR_DIR}/test_single.sh ${SCR_DIR} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL}
;;
parallel)
 ${SCR_DIR}/test_parallel.sh ${SCR_DIR} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL}
;;
nproma)
 ${SCR_DIR}/test_nproma.sh ${SCR_DIR} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL} 
;;
rerun)
 ${SCR_DIR}/test_rerun.sh ${SCR_DIR} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL}
;;
submodeloff)
 ${SCR_DIR}/test_submodeloff.sh ${SCR_DIR} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL} 
;;
update)
 ${SCR_DIR}/test_update.sh ${SCR_DIR} ${REF_ODIR} ${REF_REVISION} ${REF_MODEL} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL}
;;
parallelnproma)
 ${SCR_DIR}/test_parallelnproma.sh ${SCR_DIR} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL}
;;
parallelnpromarerun)
 ${SCR_DIR}/test_parallelnpromarerun.sh ${SCR_DIR} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL}
;;
parallelnpromarerunsubmodeloff)
 ${SCR_DIR}/test_parallelnpromarerunsubmodeloff.sh ${SCR_DIR} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL}
;;
all)
 if [ ${REF_REVISION} == ${TEST_REVISION} ]; then
  echo "test_echam6.sh: test mode ''all'' not possible if REF_REVISION==TEST_REVISION: ${REF_REVISION}"
  exit 1
 fi
 ${SCR_DIR}/test_all.sh ${SCR_DIR} ${REF_ODIR} ${REF_REVISION} ${REF_MODEL} ${TEST_ODIR} ${TEST_REVISION} ${TEST_MODEL}
;;
compile)
 echo "only compilation was required, no other test"
;;
*)
 echo "test_echam6.sh: test mode '$MODE' does not exist" >&2
 exit 1
;;
esac
