%--------------------------------------------------------------------------------------------
Preliminary note from S. Ferrachat (2016-03):

This directory contains experiment setups / scripts-related files, inherited from the echam6 base model that was used to build this echam-hammoz version. None of that is necessary for running echam-hammoz (not even for running echam), as long as you use the jobscript toolkit to build your own experiments. These files are kept as reference only.

The following text is the original README.txt that is contained in [echam base code]/run.

%--------------------------------------------------------------------------------------------

=====================================================================================
 Using ECHAM at MPIMet or DKRZ on mistral, thunder, or CIS desktops
=====================================================================================

+---------------------------------------------------------------------------+
|Some part of the input below is variable and marked by                     |
|``<angle_brackets>``. Remember to replace these markers by actual          |
|experiment name, project name, version tag, etc. before trying any of the  |
|examples.                                                                  |
|e.g. version_tag=6.3.02p1; project_name=mh0081                           |
+---------------------------------------------------------------------------+

More detailed information about generated scripts and used tools can be found
on ECHAM's wiki pages in redmine.


*mistral*
==========

1. Build the ECHAM model
------------------------

1.1 Download the current ECHAM source to your home directory:

     cd
     svn checkout https://svn.zmaw.de/svn/echam6/tags/echam-<version_tag>
     
1.2 load the necessary/recommended modules for building and scripting

     module load intel/16.0
     module load mxm/3.3.3002
     module load fca/2.5.2393 
     module load bullxmpi_mlx/bullxmpi_mlx-1.2.8.3
     module load python
 
     * Remark: the numbers of mxm/fca modules might change.

1.3 Change to your ECHAM source directory:

     cd ~/echam-<version_tag>

1.4 Configure ECHAM compilation, including parallel I/O:

     ./configure --with-fortran=intel --enable-cdi-pio

1.5 Compile using a number of parallel processes, e.g.:

     make -j 12
     make install

The model executable will be generated and stored as 'echam6' in the 'bin' directory.


2. Build the scripts
--------------------

2.1 Go to the 'run'-directory and create a copy of the 'amiptest.config' example:

     cd run
     cp examples/amiptest.config <experiment_name>.config

2.2 Edit this file and change the model settings to:

     MODEL_SUBDIR = echam-<version_tag>
 
     ACCOUNT = <project_name>

2.3 Create scripts and experiment directories:

     ../util/mkexp/mkexp <experiment_name>.config

    This prints the 'script directory' and 'data directory' needed in later steps.

2.4 Change to the script directory written by the previous step:

     cd <script_directory>

2.5 Submit the first experiment job 

     sbatch <experiment_name>.run_start (restart from control experiment) or 
     sbatch <experiment_name>.run_init  (initialize model)

2.6 To check the status of your jobs, use one of :

     squeue
     squeue -u $USER

=====================================================================================

*thunder*
=========

1. Build the ECHAM model
------------------------

1.1 Download the current ECHAM source to your project directory:

     cd /scratch/mpi/<department_name>/<project_name>/$USER
     svn checkout https://svn.zmaw.de/svn/echam6/tags/echam-<version_tag>
     
1.2 load the necessary/recommended modules for building and scripting

     module load intel/15.0.2
     module load gcc/4.8.2
     module load python
     module load openmpi/1.6.5-static-intel15

1.3 Change to your ECHAM source directory:

     cd /scratch/mpi/<department_name>/<project_name>/$USER/echam-<version_tag>

1.4 Configure ECHAM compilation:

     ./configure --with-fortran=intel

1.5 Compile using a number of parallel processes, e.g.:

     make -j4
     make install

The model executable will be generated and stored as 'echam6' in the 'bin' directory.

2. Build the scripts
--------------------

2.1 Go to the 'run'-directory and create a copy of the 'amiptest.config' example:

     cd run
     cp examples/amiptest.config <experiment_name>.config

2.2 Edit this file; change the model and environment settings to

     ENVIRONMENT = thunder

     MODEL_SUBDIR = echam-<version_tag>

     ACCOUNT = <project_name>

    Remove the lines 'DATA_ROOT' and 'WORK_ROOT' and add the project setting:

     PROJECT_SUBDIR = <department_name>/<project_name>

2.3 Create scripts and experiment directories:

     ../util/mkexp/mkexp <experiment_name>.config

    This prints the 'script directory' and 'data directory' needed in later steps.

2.4 Change to the script directory written by the previous step:

     cd <script_directory>

2.5 Submit the first experiment job

     sbatch <experiment_name>.run_start (restart from control experiment) or              
     sbatch <experiment_name>.run_init  (initialize model)
   
2.6 To check the status of your jobs, use:

     squeue -u $USER ### Add -l for more details

=====================================================================================

On *CIS desktops*
=================


1. Build the ECHAM model
------------------------

1.1 Download the current ECHAM source to your local scratch directory:

     cd /scratch/local1/$USER
     svn checkout https://svn.zmaw.de/svn/echam6/tags/echam-<version_tag>

1.2 Load the necessary/recommended modules for building and scripting

     module load nag/6.0.1038
     module load gcc/4.8.2
     module load python 
     module load openmpi/1.6.5-static-nag60

1.3 Change to your ECHAM source directory:

     cd /scratch/local1/$USER/echam-<version_tag>

1.4 Configure ECHAM to use the default compiler (NAG Fortran):

     ./configure

1.5 Compile:

     make
     make install

The model executable will be generated and stored as 'echam6' in the 'bin' directory.

2. Build the scripts
--------------------

2.1 Go to the 'run'-directory and create a copy of the 'amiptest.config' example:

     cp examples/amiptest.config <experiment_name>.config

2.2 Edit this file; change the model and environment settings to

     ENVIRONMENT = DEFAULT ### you may also remove the entire line

     MODEL_SUBDIR = echam-<version_tag>
 
    Remove the lines 'DATA_ROOT' and 'WORK_ROOT' and add the model root:

     MODEL_ROOT = /scratch/local1/$USER

    On a desktop you should use the coarse model resolution (CR), so change

     EXP_TYPE = amip-CR

     Note that this version is thought for technical tests only, not for scientific use.

2.3 Create scripts and experiment directories:

     ../util/mkexp/mkexp <experiment_name>.config

    This prints the 'script directory' and 'data directory' needed in later steps.

2.4 Change to the script directory written by the previous step:

     cd <script_directory>

2.5 Run the first experiment job in background

     ./<experiment_name>.run_init &
   
2.6 To check the status of your jobs, use:

     ps -fu $USER

=====================================================================================
