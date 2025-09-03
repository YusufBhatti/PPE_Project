#!/bin/bash
module purge
cwd=${PWD}/

if_clean="yes"
#if_clean="no"
compile="yes"

#########################################################################

# Important
# Note in order to change the dependencies of the modules, you need to modify 
# /projects/0/prjs1474/aarifi/ECHAM-HAMMOZ/RUNDIR/COMPILE/mh-linux
# also very important is that you need to set the correct compiler verision for the packages to align with your echam version, you need to check
# /projects/0/prjs1474/aarifi/ECHAM-HAMMOZ/presetups_aarifi/install.sh

# Load modules
ml 2024
ml GCC/13.3.0
ml GCCcore/13.3.0
ml OpenMPI/5.0.3-GCC-13.3.0
ml netCDF/4.9.2-gompi-2024a
ml netCDF-Fortran/4.6.1-gompi-2024a 
#ml 2024
#ml gompi/2024a
#ml HDF5/1.14.5-gompi-2024a
#ml netCDF/4.9.2-gompi-2024a
#ml netCDF-Fortran/4.6.1-gompi-2024a

export ECHAM_MAIN_DIR="/projects/0/prjs1474/ybhatti/PPE_Project"
 
DIR_ECHAM_MAIN="$ECHAM_MAIN_DIR/PPE_PACE_Project"

#########################################################################

cd ${DIR_ECHAM_MAIN}
#cp -r ${cwd}COMPILE/mh-linux ${DIR_ECHAM_MAIN}config/mh-linux
printf "Hello, %s!\n" "world"
pwd
#sed -i "\
#s|\#IF|!\#IF|g; \
#s|\#ENDIF|!\#ENDIF|g; \
#" src/mo_jsbach.f90
#printf "Hello, %s!\n" "world"

#########################################################################
#module load 2021
#module load gompi/2021a
#module load netCDF/4.8.0-gompi-2021a
#module load netCDF-Fortran/4.5.3-gompi-2021a
#module load YAXT/0.9.1-gompi-2021a
#module load HDF5/1.10.7-gompi-2021a
#module load M4/1.4.18-GCCcore-10.3.0
#module load Autotools/20210128-GCCcore-10.3.0

#module load 2023
#module load GCC/12.3.0
#module load GCCcore/12.3.0
#module load MPICH/4.2.1-GCC-12.3.0
#module load Autotools/20220317-GCCcore-12.3.0
#module load M4/1.4.19-GCCcore-12.3.0
#module load netCDF/4.9.2-gompi-2023a
#module load netCDF-Fortran/4.6.1-gompi-2023a
#module load HDF5/1.14.0-gompi-2023a


printf "Hello, %s!\n" "Module loading"
	


##cd ${DIR_ECHAM_MAIN}
##
if [[ ${if_clean} == "yes" ]]; then
	make clean
fi
printf "Cleaning time HAS COMPLETED!!"

if [[ ${compile} == "yes" ]]; then

	./configure --with-fortran=gcc
	printf "./configure --with-fortran=gcc  HAS COMPLETED"
	./config/createMakefiles.pl
	printf "./config/createMakefiles.pl HAS COMPLETED"
	libtoolize --force --install

	autoreconf -i --no-recursive
	printf "autoreconf -i --no-recursive HAS COMPLETED"

fi
make -j 128
printf "make HAS COMPLETED!!"

make install

printf "make install HAS COMPLETED!!"
printf "Compiling HAS COMPLETED!!"

cd ${cwd}
