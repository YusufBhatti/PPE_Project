#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE subjob(kjob)

  ! Description:
  !
  ! Transfer information from model into a shell-script.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, April 1995, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, March 2001, date and time control revision
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_filename,        ONLY: yexp, ypath, syr, smo, sdy, shr, smn, sse, &
                                find_next_free_unit
  USE mo_time_control,    ONLY: previous_date, current_date, next_date,    &
                                get_time_step, get_date_components
  USE mo_exception,       ONLY: finish, message_text
  USE mo_memory_base,     ONLY: nofiles, cfilenames
  USE mo_mpi,             ONLY: p_parallel_io

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: kjob

  !  Local scalars: 
  INTEGER :: idx, i
  LOGICAL :: loex
  INTEGER :: njin, njout
  
  ! define keywords for the script filter

  CHARACTER(1), PARAMETER  :: str_pre = '#'
  CHARACTER(1), PARAMETER  :: str_suf = '='

  CHARACTER(LEN=*), PARAMETER :: str_files = 'FILENAMES'

  ! elements for date/time from present output file
  CHARACTER(LEN=*), PARAMETER :: str_yr = 'YEAR'
  CHARACTER(LEN=*), PARAMETER :: str_mo = 'MONTH'
  CHARACTER(LEN=*), PARAMETER :: str_dy = 'DAY'
  CHARACTER(LEN=*), PARAMETER :: str_hr = 'HOUR'
  CHARACTER(LEN=*), PARAMETER :: str_mn = 'MINUTE'
  CHARACTER(LEN=*), PARAMETER :: str_se = 'SECOND'
  
  ! the next three blocks are related to the time window

  ! elements for date/time from previus date
  CHARACTER(LEN=*), PARAMETER :: str_o_yr = 'PREV_YEAR'
  CHARACTER(LEN=*), PARAMETER :: str_o_mo = 'PREV_MONTH'
  CHARACTER(LEN=*), PARAMETER :: str_o_dy = 'PREV_DAY'
  CHARACTER(LEN=*), PARAMETER :: str_o_hr = 'PREV_HOUR'
  CHARACTER(LEN=*), PARAMETER :: str_o_mn = 'PREV_MINUTE'
  CHARACTER(LEN=*), PARAMETER :: str_o_se = 'PREV_SECOND'
  CHARACTER(LEN=*), PARAMETER :: str_o_no = 'PREV_STEP'
  
  ! elements for date/time from present date
  CHARACTER(LEN=*), PARAMETER :: str_p_yr = 'PRES_YEAR'
  CHARACTER(LEN=*), PARAMETER :: str_p_mo = 'PRES_MONTH'
  CHARACTER(LEN=*), PARAMETER :: str_p_dy = 'PRES_DAY'
  CHARACTER(LEN=*), PARAMETER :: str_p_hr = 'PRES_HOUR'
  CHARACTER(LEN=*), PARAMETER :: str_p_mn = 'PRES_MINUTE'
  CHARACTER(LEN=*), PARAMETER :: str_p_se = 'PRES_SECOND'
  CHARACTER(LEN=*), PARAMETER :: str_p_no = 'PRES_STEP'
  
  ! elements for date/time from next date
  CHARACTER(LEN=*), PARAMETER :: str_n_yr = 'NEXT_YEAR'
  CHARACTER(LEN=*), PARAMETER :: str_n_mo = 'NEXT_MONTH'
  CHARACTER(LEN=*), PARAMETER :: str_n_dy = 'NEXT_DAY'
  CHARACTER(LEN=*), PARAMETER :: str_n_hr = 'NEXT_HOUR'
  CHARACTER(LEN=*), PARAMETER :: str_n_mn = 'NEXT_MINUTE'
  CHARACTER(LEN=*), PARAMETER :: str_n_se = 'NEXT_SECOND'
  CHARACTER(LEN=*), PARAMETER :: str_n_no = 'NEXT_STEP'
  
  ! keywords
  CHARACTER(LEN=*), PARAMETER :: str_ex = 'EXP'
  CHARACTER(LEN=*), PARAMETER :: str_dp = 'DPATH'

  CHARACTER(LEN=30) :: word

  INTEGER :: oyr, omo, ody, ohr, omn, ose       ! o* previous (old) date/time
  INTEGER :: pyr, pmo, pdy, phr, pmn, pse, pno  ! p* present date/time
  INTEGER :: nyr, nmo, ndy, nhr, nmn, nse       ! n* next date/time

  CHARACTER (4)   :: yoifil   ! input script file name
  CHARACTER (28)  :: yoofil   ! output script file name

  CHARACTER (512) :: yoline   ! script line

  !  Executable statements 

  IF (.NOT. p_parallel_io) RETURN

  CALL get_date_components(previous_date,oyr,omo,ody,ohr,omn,ose)
  CALL get_date_components( current_date,pyr,pmo,pdy,phr,pmn,pse)
  CALL get_date_components(    next_date,nyr,nmo,ndy,nhr,nmn,nse)

  pno = get_time_step()

  WRITE(yoifil,'(a3,i0)') 'job',    kjob
  WRITE(yoofil,'(a6,i0)') 'subjob', kjob

  INQUIRE (file=yoifil,exist=loex)
  IF (loex) THEN
    njin = find_next_free_unit(90,99)
    OPEN (unit=njin,file=yoifil)
  ELSE
    WRITE (message_text,*) 'Input file ' // yoifil // ' not found!'
    CALL finish('subjob',message_text)
  END IF

  INQUIRE (file=yoofil,exist=loex)
  njout = find_next_free_unit(90,99)
  IF (loex) THEN
    OPEN (unit=njout,file=yoofil,status='OLD')
  ELSE
    OPEN (unit=njout,file=yoofil,status='NEW')
  END IF

  ! filter the job script

  DO
    READ (njin,'(a)', END=999) yoline

    WRITE (njout,'(a)') TRIM(yoline)
    ! search keywords, ignore leading blanks
    idx = INDEX(yoline,str_pre)
    WRITE(word,'(a)') TRIM(yoline(idx+1:idx+30))

    SELECT CASE(word)

      ! filter experiment and path
    CASE(str_ex);    WRITE(njout,'(a,a,a)') TRIM(str_ex),str_suf,TRIM(yexp)
    CASE(str_dp);    WRITE(njout,'(a,a,a)') TRIM(str_dp),str_suf,TRIM(ypath)
      
      ! filename generation
    CASE(str_files)
      DO i=1,nofiles
        WRITE(njout,'(a)') TRIM(cfilenames(i))
      END DO

      ! filter date/time -- file label information
    CASE(str_yr);    WRITE(njout,'(a,a,i4.4)') TRIM(str_yr),str_suf,syr
    CASE(str_mo);    WRITE(njout,'(a,a,i2.2)') TRIM(str_mo),str_suf,smo
    CASE(str_dy);    WRITE(njout,'(a,a,i2.2)') TRIM(str_dy),str_suf,sdy
    CASE(str_hr);    WRITE(njout,'(a,a,i2.2)') TRIM(str_hr),str_suf,shr
    CASE(str_mn);    WRITE(njout,'(a,a,i2.2)') TRIM(str_mn),str_suf,smn
    CASE(str_se);    WRITE(njout,'(a,a,i2.2)') TRIM(str_se),str_suf,sse

      ! filter present date/time elements
    CASE(str_o_yr);  WRITE(njout,'(a,a,i4.4)') TRIM(str_o_yr),str_suf,oyr
    CASE(str_o_mo);  WRITE(njout,'(a,a,i2.2)') TRIM(str_o_mo),str_suf,omo
    CASE(str_o_dy);  WRITE(njout,'(a,a,i2.2)') TRIM(str_o_dy),str_suf,ody
    CASE(str_o_hr);  WRITE(njout,'(a,a,i2.2)') TRIM(str_o_hr),str_suf,ohr
    CASE(str_o_mn);  WRITE(njout,'(a,a,i2.2)') TRIM(str_o_mn),str_suf,omn
    CASE(str_o_se);  WRITE(njout,'(a,a,i2.2)') TRIM(str_o_se),str_suf,ose
    CASE(str_o_no);WRITE(njout,'(a,a,i20.20)') TRIM(str_o_no),str_suf,pno-1

      ! filter present date/time elements
    CASE(str_p_yr);  WRITE(njout,'(a,a,i4.4)') TRIM(str_p_yr),str_suf,pyr
    CASE(str_p_mo);  WRITE(njout,'(a,a,i2.2)') TRIM(str_p_mo),str_suf,pmo
    CASE(str_p_dy);  WRITE(njout,'(a,a,i2.2)') TRIM(str_p_dy),str_suf,pdy
    CASE(str_p_hr);  WRITE(njout,'(a,a,i2.2)') TRIM(str_p_hr),str_suf,phr
    CASE(str_p_mn);  WRITE(njout,'(a,a,i2.2)') TRIM(str_p_mn),str_suf,pmn
    CASE(str_p_se);  WRITE(njout,'(a,a,i2.2)') TRIM(str_p_se),str_suf,pse
    CASE(str_p_no);WRITE(njout,'(a,a,i20.20)') TRIM(str_p_no),str_suf,pno

      ! filter next date/time elements
    CASE(str_n_yr);  WRITE(njout,'(a,a,i4.4)') TRIM(str_n_yr),str_suf,nyr
    CASE(str_n_mo);  WRITE(njout,'(a,a,i2.2)') TRIM(str_n_mo),str_suf,nmo
    CASE(str_n_dy);  WRITE(njout,'(a,a,i2.2)') TRIM(str_n_dy),str_suf,ndy
    CASE(str_n_hr);  WRITE(njout,'(a,a,i2.2)') TRIM(str_n_hr),str_suf,nhr
    CASE(str_n_mn);  WRITE(njout,'(a,a,i2.2)') TRIM(str_n_mn),str_suf,nmn
    CASE(str_n_se);  WRITE(njout,'(a,a,i2.2)') TRIM(str_n_se),str_suf,nse
    CASE(str_n_no);WRITE(njout,'(a,a,i20.20)') TRIM(str_n_no),str_suf,pno+1

    END SELECT

  END DO

999 CONTINUE

  CLOSE (njin)
  CLOSE (njout)

END SUBROUTINE subjob
