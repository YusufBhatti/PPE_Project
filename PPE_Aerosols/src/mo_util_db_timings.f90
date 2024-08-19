!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_util_db_timings

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_DOUBLE, C_NULL_CHAR

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    FUNCTION open_db_timings(connection_info) RESULT(iret) &
         BIND(C,NAME='open_db_timings')
      IMPORT :: C_CHAR,C_INT
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: connection_info
    END FUNCTION open_db_timings
  END INTERFACE
  
  INTERFACE
    FUNCTION begin_transaction_block() RESULT(iret) &
         BIND(C,NAME='begin_transaction_block')
      IMPORT :: C_INT
      INTEGER(C_INT) :: iret
    END FUNCTION begin_transaction_block
  END INTERFACE

  INTERFACE
    FUNCTION commit_transaction_block() RESULT(iret) &
         BIND(C,NAME='commit_transaction_block')
      IMPORT :: C_INT
      INTEGER(C_INT) :: iret
    END FUNCTION commit_transaction_block
  END INTERFACE

  INTERFACE
    FUNCTION insert_experiment(experiment,outpath,host,os,user) RESULT(iret) &
         BIND(C,NAME='insert_experiment')
      IMPORT :: C_CHAR,C_INT
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: experiment
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: outpath
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: host
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: os
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: user
    END FUNCTION insert_experiment
  END INTERFACE

  INTERFACE
    FUNCTION insert_job(label,experiment,nproca,nprocb, &
         nprocio,nthreads,job_run_time_unit,job_run_time,timestep,     &
         truncation,levels,ncycle) RESULT(iret) & 
         BIND(C,NAME='insert_job')
      IMPORT :: C_CHAR,C_INT
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: label
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: experiment
      INTEGER(C_INT), VALUE,           INTENT(in) :: nproca
      INTEGER(C_INT), VALUE,           INTENT(in) :: nprocb
      INTEGER(C_INT), VALUE,           INTENT(in) :: nprocio
      INTEGER(C_INT), VALUE,           INTENT(in) :: nthreads
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: job_run_time_unit
      INTEGER(C_INT), VALUE,           INTENT(in) :: job_run_time
      INTEGER(C_INT), VALUE,           INTENT(in) :: timestep
      INTEGER(C_INT), VALUE,           INTENT(in) :: truncation
      INTEGER(C_INT), VALUE,           INTENT(in) :: levels
      INTEGER(C_INT), VALUE,           INTENT(in) :: ncycle
    END FUNCTION insert_job
  END INTERFACE

  INTERFACE
    FUNCTION insert_timing(label,job_label, &
         minval,avgval,maxval,sumval,efficiency) RESULT(iret) &
         BIND(C,NAME='insert_timing')
      IMPORT :: C_CHAR,C_DOUBLE,C_INT
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: label
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: job_label
      REAL(C_DOUBLE), VALUE,           INTENT(in) :: minval
      REAL(C_DOUBLE), VALUE,           INTENT(in) :: avgval
      REAL(C_DOUBLE), VALUE,           INTENT(in) :: maxval
      REAL(C_DOUBLE), VALUE,           INTENT(in) :: sumval
      REAL(C_DOUBLE), VALUE,           INTENT(in) :: efficiency
    END FUNCTION insert_timing
  END INTERFACE

  INTERFACE
    SUBROUTINE close_db_timings() BIND(C,NAME='close_db_timings')
    END SUBROUTINE close_db_timings
  END INTERFACE

  PUBLIC :: set_db_connection
  PUBLIC :: set_experiment_name
  PUBLIC :: set_job_name
  PUBLIC :: set_job_data

  PUBLIC :: open_db_timings
  PUBLIC :: insert_db_experiment
  PUBLIC :: insert_db_job
  PUBLIC :: insert_db_timing
  PUBLIC :: close_db_timings

  PUBLIC :: begin_transaction_block
  PUBLIC :: commit_transaction_block

  PUBLIC :: connection_info

  ! DB connection parameter
  CHARACTER(len=128,kind=C_CHAR) :: connection_info = ''

  ! Experiment table entries
  CHARACTER(len=32,kind=C_CHAR)   :: experiment_name   = ''
  CHARACTER(len=1024,kind=C_CHAR) :: output_path       = ''
  CHARACTER(len=256,kind=C_CHAR)   :: compute_host      = ''
  CHARACTER(len=256,kind=C_CHAR)   :: operating_system  = ''
  CHARACTER(len=256,kind=C_CHAR)   :: user_name         = ''

  ! Job table entries
  CHARACTER(len=32,kind=C_CHAR) :: job_name          = ''
  INTEGER(C_INT)                :: nproca            = -1_C_INT
  INTEGER(C_INT)                :: nprocb            = -1_C_INT
  INTEGER(C_INT)                :: nprocio           = -1_C_INT
  INTEGER(C_INT)                :: nthreads          = -1_C_INT
  CHARACTER(len=32,kind=C_CHAR) :: job_run_time_unit = ''
  INTEGER(C_INT)                :: job_run_time      = -1_C_INT
  INTEGER(C_INT)                :: timestep          = -1_C_INT
  INTEGER(C_INT)                :: truncation        = -1_C_INT
  INTEGER(C_INT)                :: levels            = -1_C_INT
  INTEGER(C_INT)                :: ncycle            = -1_C_INT

CONTAINS

  SUBROUTINE set_db_connection(dbhost, dbname, dbuser)
    CHARACTER(len=*), INTENT(in) :: dbhost
    CHARACTER(len=*), INTENT(in) :: dbname
    CHARACTER(len=*), INTENT(in) :: dbuser
    
    connection_info = " host = "               // TRIM(dbhost) // &
                      " port = 5432 dbname = " // TRIM(dbname) // &
                      " user = "               // TRIM(dbuser) // &
                      " password = '"          // TRIM(dbuser) // &
                      "!' "                    // C_NULL_CHAR 

  END SUBROUTINE set_db_connection

  SUBROUTINE set_experiment_name(expname, outpath)
    CHARACTER(len=*), INTENT(in) :: expname
    CHARACTER(len=*), INTENT(in) :: outpath

    CHARACTER (len=256)   :: ytmp
    INTEGER :: nlena, nlenb, nlenc

    EXTERNAL :: util_os_system, util_user_name, util_node_name

    experiment_name   = TRIM(expname)
    output_path       = TRIM(outpath)

    operating_system = ''
    user_name        = ''
    compute_host     = ''
    
    ytmp = ''
    CALL util_os_system (ytmp, nlena)
    compute_host = ytmp(1:nlena)

    ytmp = ''
    CALL util_user_name (ytmp, nlenb)
    user_name = ytmp(1:nlenb)

    ytmp = ''
    CALL util_node_name (ytmp, nlenc)
    operating_system = ytmp(1:nlenc)    

  END SUBROUTINE set_experiment_name

  SUBROUTINE set_job_name(nlcycle)
    INTEGER, INTENT(in), OPTIONAL :: nlcycle
    CHARACTER(len= 8) :: date    = ''
    CHARACTER(len=10) :: time    = ''
    CHARACTER(len=32) :: jobName = ''
    CHARACTER(len=32) :: jobId   = ''
    CHARACTER(len= 4) :: clcycle = ''
    INTEGER :: istat
    job_name = ''
    IF (PRESENT(nlcycle)) THEN
      WRITE(clcycle,'(a1,i0)') '.', nlcycle
    ENDIF
    ! check for LoadLeveler
    CALL get_environment_variable(name='LOADL_STEP_ID', value=job_name, status=istat)            
    job_name = TRIM(job_name)//TRIM(clcycle)
    IF (istat == 0) RETURN
    ! check for Sun GridEngine
    CALL get_environment_variable(name='JOB_NAME', value=jobName, status=istat)            

    CALL get_environment_variable(name='JOB_ID', value=jobId, status=istat)            
    job_name = TRIM(jobName)//'.'//TRIM(jobId)//TRIM(clcycle)
    IF (istat == 0) RETURN
    ! have to build our own job_label
    CALL date_and_time(date=date, time=time)
    job_name = TRIM(experiment_name)//'.'//TRIM(date)//TRIM(time)//TRIM(clcycle)
  END SUBROUTINE set_job_name
  
  SUBROUTINE set_job_data(in_nproca, in_nprocb, in_nprocio, in_nthreads, &
       in_job_run_time_unit, in_job_run_time,                            &
       in_timestep, in_truncation, in_levels, in_ncycle)     
    INTEGER,          INTENT(in) :: in_nproca
    INTEGER,          INTENT(in) :: in_nprocb
    INTEGER,          INTENT(in) :: in_nprocio
    INTEGER,          INTENT(in) :: in_nthreads
    CHARACTER(len=*), INTENT(in) :: in_job_run_time_unit
    INTEGER,          INTENT(in) :: in_job_run_time
    INTEGER,          INTENT(in) :: in_timestep
    INTEGER,          INTENT(in) :: in_truncation
    INTEGER,          INTENT(in) :: in_levels
    INTEGER,          INTENT(in) :: in_ncycle
    nproca            = in_nproca
    nprocb            = in_nprocb
    nprocio           = in_nprocio
    nthreads          = in_nthreads
    job_run_time_unit = in_job_run_time_unit
    job_run_time      = in_job_run_time     
    timestep          = in_timestep         
    truncation        = in_truncation       
    levels            = in_levels           
    ncycle            = in_ncycle          
  END SUBROUTINE set_job_data

  FUNCTION insert_db_experiment() RESULT(iret)
    INTEGER(C_INT) :: iret
    iret = insert_experiment(TRIM(experiment_name)//C_NULL_CHAR,  &
                             TRIM(output_path)//C_NULL_CHAR,      &
                             TRIM(compute_host)//C_NULL_CHAR,     &
                             TRIM(operating_system)//C_NULL_CHAR, &
                             TRIM(user_name)//C_NULL_CHAR)
  END FUNCTION insert_db_experiment
  
  FUNCTION insert_db_job() RESULT(iret)
    INTEGER(C_INT) :: iret
    iret = insert_job(TRIM(job_name)//C_NULL_CHAR,          &
                      TRIM(experiment_name)//C_NULL_CHAR,   &
                      nproca,nprocb, nprocio,nthreads,      &
                      TRIM(job_run_time_unit)//C_NULL_CHAR, &
                      job_run_time,timestep,                &
                      truncation,levels,ncycle)
  END FUNCTION insert_db_job

  FUNCTION insert_db_timing(timer_name, &
       minval,avgval,maxval,sumval,efficiency) RESULT(iret) 
    INTEGER(C_INT) :: iret
    CHARACTER(len=*,kind=C_CHAR), INTENT(in) :: timer_name
    REAL(C_DOUBLE),               INTENT(in) :: minval
    REAL(C_DOUBLE),               INTENT(in) :: avgval
    REAL(C_DOUBLE),               INTENT(in) :: maxval
    REAL(C_DOUBLE),               INTENT(in) :: sumval
    REAL(C_DOUBLE),               INTENT(in) :: efficiency
    iret = insert_timing(TRIM(timer_name)//C_NULL_CHAR,     &
                         TRIM(job_name)//C_NULL_CHAR, &
                         minval,avgval,maxval,sumval,efficiency)
  END FUNCTION insert_db_timing

END MODULE mo_util_db_timings

