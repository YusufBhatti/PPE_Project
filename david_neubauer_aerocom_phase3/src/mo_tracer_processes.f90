!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mo_tracer_processes: collection of generic routines affecting the xt arrays:
!!      xt_initialize      : intialisation of tracer concentrations
!!      xt_burden_init_mem : initialisation of burden diagnostics for tracers
!!      xt_burden          : burden diagnostics
!!      xt_conv_massfix    : mass fixer for convective transport (borrowing scheme)
!!      xt_borrow          : borrowing scheme to prevent negative tracer mass
!!
!!
!! @author 
!! <ol>
!! <li>ECHAM5 developers 
!! <li>M. Schultz (FZ-Juelich)
!! </ol>
!!
!! $Id: 1423$
!!
!! @par Responsible coder
!! m.schultzfz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE mo_tracer_processes

  USE mo_kind,             ONLY: wp 
  USE mo_filename,         ONLY: trac_filetype, find_next_free_unit
  IMPLICIT NONE

  
  PRIVATE

  !! Interface routines         ! purpose                   ! called by
                                                      
  PUBLIC :: xt_initialize       ! set initial values        ! ioinitial, iorestart   (### former xtini)
  PUBLIC :: xt_burden_init_mem  ! set up burden diagnostics
  PUBLIC :: xt_burden           ! calculate tracer burdens
  PUBLIC :: xt_conv_massfix     ! adjust tracer mass in convection
  PUBLIC :: xt_borrow           ! another mass fixer

  !--- module variables
  REAL(wp), ALLOCATABLE :: zxtte_old(:,:,:)          ! for xt_convmassfix
  !$OMP THREADPRIVATE(zxtte_old)

  CONTAINS


  SUBROUTINE xt_initialize (xt, xtm1)

  ! Description:
  !
  ! Initialize the tracer fields
  !
  ! Method:
  !
  ! This routine is called from *ioinitial* and from *iorestart*.
  ! Tracers not yet initialised are set to constant values.
  ! Tracer input from file will be possible in the future.
  !
  ! Authors:
  !
  ! J. Feichter, MI, August 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, rewrite
  ! M. Schultz, FZ Juelich - adaptation to ECHAM6 (2009-09-25)
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_mpi,                 ONLY: p_parallel_io
  USE mo_semi_impl,           ONLY: eps
  USE mo_exception,           ONLY: finish,         &       
                                    message,        &       
                                    message_text
  USE mo_util_string,         ONLY: separator          ! format string (----)
  USE mo_tracdef,             ONLY: trlist,         &  ! tracer info variable
!!mgs!!                                    ln, ll, nf,     &  ! len of char components of trlist 
                                    ntrac,          &  ! number of tracers defined
                                    INITIAL,        &  ! 
                                    RESTART,        & 
                                    CONSTANT

  !
  !  Arguments 
  !
  REAL(wp),INTENT(inout)           :: xt  (:,:,:,:) ! tracer array
  REAL(wp),INTENT(inout) ,OPTIONAL :: xtm1(:,:,:,:) ! tracer array at t-dt
  !
  !  Local scalars:
  !
  INTEGER          :: jt    ! tracer index variable
  CHARACTER(len=9) :: cini  ! initialisation method to print
  LOGICAL          :: fault ! indicates that some tracer was not initialised
  !
  ! Loop over all tracers , flag tracers read from restart file.
  !
  DO jt = 1, trlist% ntrac
    trlist% ti(jt)% init = 0
    IF (trlist% mi(jt)% xt%   restart_read .AND. &
        trlist% mi(jt)% xtm1% restart_read )     &
        trlist% ti(jt)% init = RESTART
  END DO
  !
  ! If restart flag is not set, ignore values read so far.
  !
  trlist% ti(1:ntrac)% init = IAND (trlist% ti(1:ntrac)% init, &
                                    trlist% ti(1:ntrac)% ninit)
  !
  ! Set 'lrerun' flag to .false. for XT and XTM1 on output.
  ! Old restart file format is not written any more.
  !
  trlist% mixt  % lrerun = .FALSE.    
  trlist% mixtm1% lrerun = .FALSE.
  !
  ! read initial concentrations from file
  !
  CALL xt_init_file(xt,xtm1)
  !
  ! loop over tracers not initialised so far
  !
  DO jt = 1, trlist% ntrac
    IF (trlist% ti(jt)% init > 0) CYCLE ! skip if already initialised
    !
    ! Set to constant value
    !
    IF (IAND (trlist% ti(jt)% ninit, CONSTANT) /= 0) THEN
      xt(:,:,jt,:)         = trlist% ti(jt)% vini
      trlist% ti(jt)% init = CONSTANT
      IF (PRESENT(xtm1)) xtm1(:,:,jt,:) = (1._wp - eps) * xt(:,:,jt,:)
    ENDIF
  END DO
  !
  ! loop over all tracers, print initialisation method and min,max value
  !
  IF (ntrac > 0 .AND. p_parallel_io) THEN
    CALL message('',separator)
    CALL message('','')
    CALL message('','')
    CALL message('xt_initialize','initial values of tracers:')
    CALL message('','')
    CALL message('','  tracer          source      minval    maxval   minval,maxval(xtm1)')
    fault = .FALSE.
    DO jt = 1, trlist% ntrac
      SELECT CASE (trlist% ti(jt)% init)
      CASE (CONSTANT)
        cini = 'constant'
      CASE (INITIAL)
        cini = 'initial '
      CASE (RESTART)
        cini = 'restart '
      CASE (0)
        cini = 'NOT INIT.'
        fault=.TRUE.
      CASE DEFAULT
        cini = 'unknown '
      END SELECT  
      IF(PRESENT(xtm1)) THEN
        WRITE(message_text,'(2x,a16,a8,2x,4g10.3)')       &
             trlist%ti(jt)%fullname, cini,                &
             MINVAL(xt(:,:,jt,:)),MAXVAL(xt(:,:,jt,:)),   &
             MINVAL(xtm1(:,:,jt,:)),MAXVAL(xtm1(:,:,jt,:))
        CALL message('',message_text)
      ELSE
        WRITE(message_text,'(2x,a16,a8,2x,4g10.3)')       &
             trlist%ti(jt)%fullname, cini,                &
             MINVAL(xt(:,:,jt,:)),MAXVAL(xt(:,:,jt,:))
        CALL message('',message_text)
      ENDIF
    END DO
    CALL message('','')
    CALL message('',separator)
    !
    ! Abort if any tracer is not initialised
    !
    IF(fault) CALL finish ('xt_initialize','tracer not initialised')
  ENDIF

  END SUBROUTINE xt_initialize



  SUBROUTINE xt_init_file(xt, xtm1)

  !
  ! This routine is called from xt_initialize to allow initialization of tracer
  ! variables besides initialization with constant values or from the
  ! rerun file.
  !
  ! Base the decision whether to initialize on the following conditions:
  !
  ! trlist%ti(jt)%init   == 0       (not initialized so far)
  ! trlist%ti(jt)%ninit: >= INITIAL (initialisation requested)
  !
  ! Set trlist%ti(jt)%init to a value =/0 afterwards.
  !

  USE mo_control,                   ONLY: ngl, nlon, nlev, ldebugio
  USE mo_mpi,                       ONLY: p_parallel_io, p_io, p_bcast
  USE mo_exception,                 ONLY: message, message_text, em_error, em_warn, em_info, em_debug
  USE mo_netcdf
  USE mo_transpose,                 ONLY: scatter_gp
  USE mo_decomposition,             ONLY: dcg => global_decomposition
  USE mo_read_netcdf77,             ONLY: read_var_nf77_3d, read_var_hs_nf77_3d
  USE mo_tracdef,                   ONLY: trlist, INITIAL
  USE mo_physical_constants,        ONLY: amd
  USE mo_filename,                  ONLY: find_next_free_unit

  REAL(wp),INTENT(inout)           :: xt  (:,:,:,:) ! tracer array
  REAL(wp),INTENT(inout) ,OPTIONAL :: xtm1(:,:,:,:) ! tracer array at t-dt

  CHARACTER(len=16), PARAMETER  :: trinifile = 'tracer_ic.nc'    ! file name for tracer file
  CHARACTER(len=16)   :: clon, clat, clev, ctime                 ! variable names for coordinates
  CHARACTER(len=16)   :: varname                                 ! variable name for tracer variable
  CHARACTER(len=80)   :: cunits                                  ! units string of variable
  CHARACTER(len=32)   :: cfactor                                 ! character string of conversion factor
  INTEGER             :: jt, ierr, ierr0, iunit
  INTEGER             :: fileID, varID, timeID, nrec             ! netcdf indices
  INTEGER             :: ndims                                   ! number of variable dimensions
  INTEGER             :: nnodef                                  ! number of tracers not found in file
  CHARACTER(len=194)  :: cnodef                                  ! names of tracers not found in file
  LOGICAL             :: lfirst = .TRUE.,  &                     ! first attempt to read
                         lfound = .FALSE.                        ! netcdf file trinifile valid
  REAL(wp)            :: zfactor                                 ! factor for unit conversion
  REAL(wp), POINTER   :: zptr(:,:,:) ! for reading variables from netCDF-file
       
  !-- Initialisation
  IF (p_parallel_io) THEN
    ALLOCATE(zptr(nlon, nlev, ngl))
  ELSE
    NULLIFY(zptr)
  ENDIF
  clon = 'lon'
  clat = 'lat'
  clev = 'mlev'
  ctime = 'time'    ! optional 4th dimension: last index will be read
  nnodef = 0
  cnodef = ''
  cunits = ''
  cfactor = ''

  !-- Find tracer file
  !-- Find tracer file
  ! check, if there is a file existing (probably with 0 bytes -- see issue #225)
  iunit = find_next_free_unit (30, 100)
  OPEN(iunit,file=TRIM(trinifile),status='OLD',iostat=ierr0)
  ierr = NF_OPEN(TRIM(trinifile), NF_NOWRITE, fileID)
  IF (ierr == NF_NOERR) THEN
     lfound = .TRUE.
     CALL nf_check(NF_CLOSE(fileID), fname=TRIM(trinifile))
  ELSE
     IF (ierr0 == 0) CALL message('xt_init_file', 'damaged or empty tracer_ic.nc!', level=em_error)
  END IF  
  IF(ldebugio) THEN
    ! Note: if ldebugio=false we give a warning only when tracers are actually expected to be read
    IF (lfound) THEN
      CALL message('xt_init_file', 'Successfully opened tracer IC file '//TRIM(trinifile), &
                   level=em_debug)
    ELSE
      CALL message('xt_init_file', 'Failed to open tracer IC file '//TRIM(trinifile), &
                   level=em_debug)
    END IF
  END IF

  !-- Loop over all tracers and read all of those that have not yet been initialized and
  !   which allow initialisation (from file)
  DO jt = 1, trlist%ntrac
    IF (trlist%ti(jt)%init == 0 .AND. IAND(trlist%ti(jt)%ninit,INITIAL) /= 0 ) THEN
      ierr=1       ! default: return error
      !-- report warning if no initial tracer file was found
      IF (lfirst .AND. .NOT. lfound) THEN
        CALL message('xt_init_file', 'Cannot find file '//TRIM(trinifile)//'!', level=em_warn)
        CALL message('', 'Problem occured while trying to initialize '//TRIM(trlist%ti(jt)%basename), &
                     level=em_info)
        RETURN      ! no sense to continue
      END IF
      lfirst = .FALSE.
      !-- try to read variable from file
      IF (p_parallel_io) THEN
        ! look for variable containing tracer mixing ratio
        CALL nf_check(NF_OPEN(TRIM(trinifile), NF_NOWRITE, fileID), fname=TRIM(trinifile))
        CALL get_nf_varname(fileID, trlist%ti(jt)%fullname, varname, varID)
        ! ... if not found try again with tracer basename
        IF (varID < 0) CALL get_nf_varname(fileID, trlist%ti(jt)%basename, varname, varID)
        IF (ldebugio) THEN
          WRITE(message_text,*) 'Variable ID for '//TRIM(varname)//' : ', varID,   &
                                '   (tracer name='//TRIM(trlist%ti(jt)%basename)//')'
          CALL message('', message_text, level=em_debug)
        END IF
        IF (varID >= 0) THEN
          ! find out number of dimensions
          CALL nf_check(nf_inq_varndims(fileID, varID, ndims), fname=TRIM(trinifile))
          IF (ndims == 4) THEN 
            ! find maximum index of time dimension
            CALL nf_check(nf_inq_dimid(fileID, ctime, timeID), fname=TRIM(trinifile))
            CALL nf_check(nf_inq_dimlen(fileID, timeID, nrec), fname=TRIM(trinifile))
            CALL nf_check(NF_CLOSE(fileID), fname=TRIM(trinifile))
            CALL read_var_hs_nf77_3d (trinifile, clon, clev, clat, ctime,    &
                                      nrec, varname, zptr, ierr)
            IF (ldebugio) THEN
              WRITE(message_text,*) 'Read in hyperslab for variable '//TRIM(varname)//   &
                                    ': nrec=', nrec, ', ierr = ', ierr
              CALL message('', message_text, level=em_debug)
            END IF
          ELSE IF (ndims == 3) THEN
            CALL nf_check(NF_CLOSE(fileID), fname=TRIM(trinifile))
            CALL read_var_nf77_3d (trinifile, clon, clev, clat,        &
                                   trlist%ti(jt)%basename, zptr, ierr)
            IF (ldebugio) THEN
              WRITE(message_text,*) 'Read in variable '//TRIM(varname)//   &
                                    ': ierr = ', ierr
              CALL message('', message_text, level=em_debug)
            END IF
          ELSE
            WRITE(message_text,*) 'Invalid dimensions of variable '//TRIM(trlist%ti(jt)%basename)//   &
                                  ' in file '//TRIM(trinifile)//'! ndims = ',ndims
            CALL message('xt_init_file', message_text, level=em_error)
            CYCLE
          END IF
          ! report error if reading failed
          IF (ierr /= 0) THEN
            WRITE(message_text,*) 'Error while reading variable '//TRIM(trlist%ti(jt)%basename)//   &
                                  ' from file '//TRIM(trinifile)//'! ierr = ',ierr
            CALL message('xt_init_file', message_text, level=em_error)
            CYCLE
          END IF
          ! obtain text of units attribute and appropriate scaling factor
          CALL get_units_factor (trinifile, varname, &
                                 trlist%ti(jt)%moleweight, amd, zfactor, cunits, ierr)
          IF (ierr == 0) THEN
            WRITE(cfactor,*) zfactor
            WRITE(message_text,*) 'Read '//TRIM(trlist%ti(jt)%basename)//         &
                                  ' from '//TRIM(trinifile)//                     &
                                  '. Units = '//TRIM(cunits)//                    &
                                  ', conversion factor = '//TRIM(cfactor)
            CALL message('xt_init_file', message_text, level=em_info)
            zptr(:,:,:) = zptr(:,:,:) * zfactor
          ELSE
            WRITE(message_text,*) 'Error while trying to convert units for variable '//    &
                                  TRIM(trlist%ti(jt)%basename)//' in file '//TRIM(trinifile)
            CALL message('xt_init_file', message_text, level=em_error)
          END IF

        END IF    ! variable found in file
      END IF
      CALL p_bcast (ierr, p_io)
      IF (ierr == 0) THEN
        CALL scatter_gp(zptr,xt(:,:,jt,:),dcg)
        IF (PRESENT(xtm1)) xtm1(:,:,jt,:) = xt(:,:,jt,:)
        ! set flag that shows that tracer is already initialized
        trlist%ti(jt)%init = INITIAL
      ELSE
        IF (nnodef == 0) THEN
          cnodef = TRIM(trlist%ti(jt)%fullname)
        ELSE
          IF (len(TRIM(cnodef)//', '//TRIM(trlist%ti(jt)%fullname)) <= 194) THEN
            cnodef = TRIM(cnodef)//', '//TRIM(trlist%ti(jt)%fullname)
          ELSE
! already write undefined tracers to logfile
! make sure message_text is not longer than defined in subroutine message (len=256)
! subroutine "message" changes public variable "message_text" in extending it with
! (in this case!) "WARNING: xt_init_file: "!!!
! ==> make sure that extended message_text will not exceed 256 characters!
! string "WARNING: xt_init_file: Tracers not found in TRIM(trinifile): " already accounts for 62 characters
! The following if statement is just for security reasons if max length of trinifile is modified.
            IF (len(TRIM('WARNING: xt_init_file: Tracers not found in '//  &
                    TRIM(trinifile)//': '//TRIM(cnodef))) <= 256) THEN
              message_text='Tracers not found in '//TRIM(trinifile)//': '//TRIM(cnodef)
              CALL message('xt_init_file', message_text, level=em_warn)
            ELSE
              message_text='Tracers not found in trinifile: '//TRIM(cnodef)
              CALL message('xt_init_file', message_text, level=em_warn)
            END IF
            cnodef = "...(continued) "//TRIM(trlist%ti(jt)%fullname)
         END IF
        END IF
        nnodef = nnodef+1
      END IF
    END IF    ! tracer needs initialisation from file
  END DO

  IF (nnodef > 0) THEN
! make sure message_text is not longer than defined in subroutine message (len=256)
! (see explanation above)
    IF (len(TRIM('WARNING: xt_init_file: Tracers not found in '//  &
            TRIM(trinifile)//': '//TRIM(cnodef))) <= 256) THEN
      message_text='Tracers not found in '//TRIM(trinifile)//': '//TRIM(cnodef)
      CALL message('xt_init_file', message_text, level=em_warn)
    ELSE
      message_text='Tracers not found in trinifile: '//TRIM(cnodef)
      CALL message('xt_init_file', message_text, level=em_warn)
    END IF
  END IF

  IF (ASSOCIATED(zptr)) DEALLOCATE(zptr)

  END SUBROUTINE xt_init_file


  SUBROUTINE get_nf_varname(fileID, trname, varname, varID)

  ! find a tracer variable in the initial condition file. The search is greedy and will
  ! allow for upper and lower case characters or prepended or appended 'vmr_' or 'mmr_' strings
  ! The netcdf file must be open and fileID must be valid

  USE mo_netcdf
  USE mo_util_string,         ONLY: tolower

  INTEGER, INTENT(in)           :: fileID
  CHARACTER(len=*), INTENT(in)  :: trname
  CHARACTER(len=*), INTENT(out) :: varname
  INTEGER, INTENT(out)          :: varID

  CHARACTER(len=256)   :: ctmp
  INTEGER              :: jv, nvars, ierr

  ! test original tracer name first
  varname = TRIM(trname)
  ierr=nf_inq_varid(fileID, TRIM(varname), varID)
  ! if not successful, try pattern match for each variable
  IF (ierr /= 0) THEN
    varID = -1
    CALL nf_check(nf_inq_nvars(fileID, nvars))
    DO jv=1,nvars
      CALL nf_check(nf_inq_varname(fileID, jv, ctmp))
      IF (TRIM(tolower(varname)) == TRIM(tolower(ctmp))) THEN
        varID = jv
        varname = TRIM(ctmp)      ! copy to preserve upper/lower case spelling
      ELSE IF ('vmr_'//TRIM(tolower(varname)) == TRIM(tolower(ctmp))) THEN
        varID = jv
        varname = TRIM(ctmp)
      ELSE IF (TRIM(tolower(varname))//'_vmr' == TRIM(tolower(ctmp))) THEN
        varID = jv
        varname = TRIM(ctmp)
      ELSE IF ('mmr_'//TRIM(tolower(varname)) == TRIM(tolower(ctmp))) THEN
        varID = jv
        varname = TRIM(ctmp)
      ELSE IF (TRIM(tolower(varname))//'_mmr' == TRIM(tolower(ctmp))) THEN
        varID = jv
        varname = TRIM(ctmp)
      END IF
    END DO
  END IF
    
  END SUBROUTINE get_nf_varname


  SUBROUTINE get_units_factor(file_name, var_name, var_mwght, air_mwght, &
                              pfactor, punits, ierr)
  ! Description:
  !
  ! Obtain conversion factor for initial tracer concentrations.
  ! Currently handles only mass mixing ratio (factor = 1) and volume mixing
  ! ratio (factor = var_mwght/air_mwght)
  !
  ! Authors:
  !
  ! J.S. Rast, MPI, December 2003, original source
  ! M.G. Schultz, FZ Juelich, May 2010 -- adapted for ECHAM6-HAMMOZ

  USE mo_netcdf
  USE mo_kind,             ONLY: wp
  USE mo_util_string,      ONLY: tolower
  USE mo_exception,        ONLY: message, em_warn
  !
  ! Arguments:
  !
  CHARACTER (LEN = *), INTENT (in) :: file_name
  CHARACTER (LEN = *), INTENT (in) :: var_name
  REAL(wp), INTENT(in)             :: var_mwght
  REAL(wp), INTENT(in)             :: air_mwght
  REAL(wp), INTENT(out)            :: pfactor
  CHARACTER (len=80), INTENT(out)  :: punits
  INTEGER, INTENT(out)             :: ierr
  !
  ! Local variables:
  !
  INTEGER                          :: zncid, zvarid, zierr, zlen
  LOGICAL                          :: lfound
  CHARACTER (len=80)               :: ctmp, zunits

  !-- Initialize
  pfactor = 1.0_wp
  ierr = 0
  lfound = .FALSE.

  !-- Open netCDF file
  CALL nf_check(nf_open(TRIM(file_name), nf_nowrite, zncid), &
          fname=TRIM(file_name))
  !-- inquire variable name
  zierr=nf_inq_varid(zncid, TRIM(var_name), zvarid)
  IF (zierr /= NF_NOERR) THEN
    ierr = zierr
    write(0,*) '**** mo_tracer_processes: get_units_factor :: Should never be here! ***'
    CALL nf_check(nf_close(zncid), fname=TRIM(file_name))
    RETURN
  END IF
  !-- inquire units attribute
  zierr = nf_inq_attlen(zncid, zvarid, 'units', zlen)
  IF (zierr /= NF_NOERR) THEN
    ierr = zierr
    CALL nf_check(nf_close(zncid), fname=TRIM(file_name))
    RETURN         ! failed to read units attribute
  END IF
  !-- Get text of units attribute
  ctmp = '   '
  zierr = nf_get_att_text(zncid, zvarid, 'units', ctmp)
  IF (zierr /= NF_NOERR) THEN
    ierr = zierr
    zunits = 'unknown'
    CALL nf_check(nf_close(zncid), fname=TRIM(file_name))
    RETURN         ! failed to read units attribute
  END IF
  CALL nf_check(nf_close(zncid), fname=TRIM(file_name))
  !-- Process units string
  zunits = tolower(TRIM(ctmp))
  !-- Check for mass mixing ratio units (safety check)
  IF (zunits(1:7) == 'kg kg-1' .OR. zunits(1:5) == 'kg/kg' .OR.    &
      zunits(1:3) == 'mmr') THEN
    lfound = .TRUE.
    pfactor = 1.0_wp      ! no conversion necessary
  !-- Check for volume mixing ratio units
  ELSE IF (zunits(1:11) == 'mole mole-1' .OR. zunits(1:9) == 'mole/mole' .OR.   &
           zunits(1:3) == 'vmr') THEN
    lfound = .TRUE.
    pfactor = var_mwght / air_mwght
  END IF

  IF (.NOT. lfound) THEN
    CALL message('get_units_factor',       &
                 'Cannot find mass or volume mixing ratio unit. Units = '//TRIM(zunits),    &
                 level=em_warn)
  END IF

  punits = zunits     ! for reporting

  END SUBROUTINE get_units_factor

!-----  xt_burden  ---------------------------------------------------------------
 
  ! Description:
  !
  ! Calculation of tracer burdens (column integrals)
  ! The mo_tracer module comprises the formerly independent routine xt_burden and the 
  ! burden memory management which was formerly handled in mo_aero_mem
  !
  ! Authors:
  !
  ! P. Stier, MPI, 2001, original source
  ! M. Schultz, FZ Juelich - adaptation to ECHAM6 (2009-09-25)
  ! 
  ! for more details see file AUTHORS
  !

!### tstream: option to add burden diag to existing stream. May be obsolete, but leave for now.
  SUBROUTINE xt_burden_init_mem ( tstream )

  USE mo_submodel,       ONLY: lburden     ! switch burden diagnostics on/off
  USE mo_tracer,         ONLY: nburden, d_burden
  USE mo_linked_list,    ONLY: t_stream, SURFACE
  USE mo_memory_base,    ONLY: new_stream, add_stream_element,     &
                               default_stream_setting, add_stream_reference, AUTO

  ! Arguments
  TYPE (t_stream), TARGET, INTENT(inout), OPTIONAL :: tstream    ! diag stream reference
   
  ! Local variables
  TYPE (t_stream), POINTER     :: diagstream
  CHARACTER(len=32)            :: cunit
  INTEGER                      :: jn


  IF (.NOT. lburden) RETURN     ! don't do anything if burden diagnostics is inactive 

  ! Set stream pointer to existing stream if given
  IF (PRESENT(tstream)) THEN
     diagstream => tstream
  ELSE
  ! Otherwise open a new diagnostic stream for burden diagnostics
  ! file type taken from tracer files.
    CALL new_stream (diagstream ,'burden',filetype=trac_filetype)       
    !--- Add standard fields for post-processing:
    CALL add_stream_reference (diagstream, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (diagstream, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (diagstream, 'aps'     ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (diagstream, 'gboxarea','geoloc',lpost=.TRUE.)
    !--- Default stream element settings
    CALL default_stream_setting (diagstream, lrerun    = .TRUE. , &
                                 contnorest= .TRUE. ,      &
                                 laccu     = .TRUE. ,      &
                                 lpost     = .TRUE. ,      &
                                 leveltype = SURFACE,      &
                                 table     = 199,          &
                                 code      = AUTO     )
  END IF

  ! Add burden fields for all d_burden entries
  ! total column mass
  DO jn=1, nburden
    IF (IAND(d_burden(jn)%itype, 1) /= 0) THEN
      cunit = 'kg m-2'
      CALL add_stream_element(diagstream, 'burden_'//TRIM(d_burden(jn)%name),  &
                              d_burden(jn)%ptr1, units=cunit,                  &
                              longname='atmospheric burden of '//TRIM(d_burden(jn)%name))
    END IF
  END DO
  ! tropospheric column mass
  DO jn=1, nburden
    IF (IAND(d_burden(jn)%itype, 2) /= 0) THEN
      cunit = 'kg m-2'
      CALL add_stream_element(diagstream, 'trburden_'//TRIM(d_burden(jn)%name),  &
                              d_burden(jn)%ptr2, units=cunit,                  &
                              longname='tropospheric burden of '//TRIM(d_burden(jn)%name))
    END IF
  END DO
  ! stratospheric column mass
  DO jn=1, nburden
    IF (IAND(d_burden(jn)%itype, 4) /= 0) THEN
      cunit = 'kg m-2'
      CALL add_stream_element(diagstream, 'stburden_'//TRIM(d_burden(jn)%name),  &
                              d_burden(jn)%ptr4, units=cunit,                  &
                              longname='stratospheric burden of '//TRIM(d_burden(jn)%name))
    END IF
  END DO
  ! total column density
  DO jn=1, nburden
    IF (IAND(d_burden(jn)%itype, 8) /= 0) THEN
      cunit = 'm-2'
      CALL add_stream_element(diagstream, 'coldensity_'//TRIM(d_burden(jn)%name),  &
                              d_burden(jn)%ptr8, units=cunit,                  &
                              longname='atmospheric column density of '//TRIM(d_burden(jn)%name))
    END IF
  END DO
  ! tropospheric column density
  DO jn=1, nburden
    IF (IAND(d_burden(jn)%itype, 16) /= 0) THEN
      cunit = 'm-2'
      CALL add_stream_element(diagstream, 'trcoldensity_'//TRIM(d_burden(jn)%name),  &
                              d_burden(jn)%ptr16, units=cunit,                  &
                              longname='tropospheric column density of '//TRIM(d_burden(jn)%name))
    END IF
  END DO
  ! stratospheric column density
  DO jn=1, nburden
    IF (IAND(d_burden(jn)%itype, 32) /= 0) THEN
      cunit = 'm-2'
      CALL add_stream_element(diagstream, 'stcoldensity_'//TRIM(d_burden(jn)%name),  &
                              d_burden(jn)%ptr32, units=cunit,                  &
                              longname='stratospheric column density of '//TRIM(d_burden(jn)%name))
    END IF
  END DO

  END SUBROUTINE xt_burden_init_mem


  SUBROUTINE xt_burden(kproma, kbdim, klev, klevp1, krow, & ! Local dimensions
                       papp1,  paphp1,                    & ! Thermodynamic quantities
                       pxtm1,  pxtte                      ) ! Tracer fields

  !  *xt_burden* calculates the atmospheric burden and column densities for tracers
  !  for which a burden diagnostics was defined. Note that several tracers can share 
  !  the same diagnostics to allow output for example by aerosol species or mode.
  !
  !  Author:
  !  -------
  !  Martin Schultz, FZ-Juelich     (2009-09-25)
  !  based on code from P. Stier (2002)
  !
  !  Method:
  !  -------
  !  Calculate burden [kg m-2] from mixing ratio [kg kg-1] by multiplying 
  !  with the term dp/g and summing up the vertical integral.
  !  The burden is multiplied by dt and divided by the output interval 
  !  in the stream management to give a mean over the output interval. 
  !
  !  ToDo:
  !  -----
  !  ### implement column density calculation (burden%itype >= 8)
  !  ### implement global sum diagnostics (itype=64) - use trastat. ??

  USE mo_kind,         ONLY: wp
  USE mo_exception,    ONLY: finish, message, em_info, em_warn
  USE mo_time_control, ONLY: delta_time, time_step_len
  USE mo_physical_constants,    ONLY: grav 
  USE mo_submodel,     ONLY: lburden
  USE mo_vphysc,       ONLY: vphysc
  USE mo_tracdef,      ONLY: ntrac, trlist
  USE mo_tracer,       ONLY: nburden, d_burden


  !--- Dummy Variables:

  INTEGER,  INTENT(in) :: kproma, kbdim, klev, klevp1, krow

  REAL(wp), INTENT(in) :: paphp1(kbdim,klevp1),    papp1(kbdim,klev)
  REAL(wp), INTENT(in) :: pxtm1(kbdim,klev,ntrac), pxtte(kbdim,klev,ntrac)


  !--- Local Variables:

  INTEGER              :: jl, jk, jt, iburden, itype

  REAL(wp)             :: zxtp1(kbdim)
  REAL(wp)             :: zdpg(kbdim,klev)

  LOGICAL, SAVE        :: lfirst=.TRUE.



  !--- 0) Initializations:

  IF (ntrac == 0 .OR. .NOT. lburden) RETURN    ! don't do anything if burden diagnostics is inactive

  !--- Calculate burden only for tracers that are given the unit of mixing ratio,
  !    i.e. assume that the unit is mixing ratio occurs if the string kg-1 occurs
  !    in the tracer units:
  
  IF (lfirst) THEN
    CALL message('xt_burden', '----------------------------------------------------', level=em_info)
    ! test for validity of burdenid
    DO jt = 1, ntrac
      IF (trlist%ti(jt)%burdenid > nburden) CALL finish('xt_burden', 'Invalid burden id for tracer '// &
                                                        trlist%ti(jt)%fullname)
    END DO
  END IF

  !--- 1) Calculate auxiliary variable dp/g :

  !--- Uppermost level:
  zdpg(1:kproma,1)=2._wp*(paphp1(1:kproma,2)-papp1(1:kproma,1))/grav
  !--- Other levels:
  DO jk=2, klev
     zdpg(1:kproma,jk)=(paphp1(1:kproma,jk+1)-paphp1(1:kproma,jk))/grav
  END DO


  !--- 2) Calculate and store d_burden:

  DO jt=1, ntrac
     iburden = trlist%ti(jt)%burdenid
     IF (iburden <= 0) CYCLE    ! nothing to be done for this tracer

     itype = d_burden(iburden)%itype

     ! total column mass
     IF (IAND(itype, 1) /= 0) THEN
        DO jk=1, klev
          zxtp1(1:kproma)=pxtm1(1:kproma,jk,jt)+pxtte(1:kproma,jk,jt)*time_step_len
          d_burden(iburden)%ptr1(1:kproma,krow) = d_burden(iburden)%ptr1(1:kproma,krow)   &
                                            + zxtp1(1:kproma)*zdpg(1:kproma,jk)*delta_time
        END DO
     END IF
     ! tropospheric column mass
     IF (IAND(itype, 2) /= 0) THEN
        DO jk=1, klev
          DO jl=1, kproma
            IF (vphysc%trpwmo(jl,krow) < jk) THEN
              zxtp1(jl)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*time_step_len
              d_burden(iburden)%ptr2(jl,krow) = d_burden(iburden)%ptr2(jl,krow)   &
                                                + zxtp1(jl)*zdpg(jl,jk)*delta_time
            END IF
          END DO
        END DO
     END IF
     ! stratospheric column mass
     IF (IAND(itype, 4) /= 0) THEN
        DO jk=1, klev
          DO jl=1, kproma
            IF (vphysc%trpwmo(jl,krow) > jk) THEN
              zxtp1(jl)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*time_step_len
              d_burden(iburden)%ptr4(jl,krow) = d_burden(iburden)%ptr4(jl,krow)   &
                                                + zxtp1(jl)*zdpg(jl,jk)*delta_time
            END IF
          END DO
        END DO
     END IF

!### MUST ADD UNIT CONVERSION / CHECK HERE !!!
!### the following code is just providing the structure but will not work properly!!
     ! total column density
     IF (IAND(itype, 8) /= 0) THEN
IF (lfirst) CALL message('xt_burden','burdentype 8 not properly implemented yet!', level=em_warn)
        DO jk=1, klev
          zxtp1(1:kproma)=pxtm1(1:kproma,jk,jt)+pxtte(1:kproma,jk,jt)*time_step_len
          d_burden(iburden)%ptr8(1:kproma,krow) = d_burden(iburden)%ptr8(1:kproma,krow)   &
                                            + zxtp1(1:kproma)*zdpg(1:kproma,jk)*delta_time
        END DO
     END IF
     ! tropospheric column density
     IF (IAND(itype, 16) /= 0) THEN
IF (lfirst) CALL message('xt_burden','burdentype 16 not properly implemented yet!', level=em_warn)
        DO jk=1, klev
          DO jl=1, kproma
            IF (vphysc%trpwmo(jl,krow) < jk) THEN
              zxtp1(jl)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*time_step_len
              d_burden(iburden)%ptr16(jl,krow) = d_burden(iburden)%ptr16(jl,krow)   &
                                                + zxtp1(jl)*zdpg(jl,jk)*delta_time
            END IF
          END DO
        END DO
     END IF
     ! stratospheric column density
     IF (IAND(itype, 32) /= 0) THEN
IF (lfirst) CALL message('xt_burden','burdentype 32 not properly implemented yet!', level=em_warn)
        DO jk=1, klev
          DO jl=1, kproma
            IF (vphysc%trpwmo(jl,krow) > jk) THEN
              zxtp1(jl)=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*time_step_len
              d_burden(iburden)%ptr32(jl,krow) = d_burden(iburden)%ptr32(jl,krow)   &
                                                + zxtp1(jl)*zdpg(jl,jk)*delta_time
            END IF
          END DO
        END DO
     END IF
  END DO


  lfirst=.FALSE.

  END SUBROUTINE xt_burden


  SUBROUTINE xt_conv_massfix (kproma,        kbdim,             klev,         &
                              klevp1,        ktrac,             krow,         &
                              papp1,         paphp1,            pxtte,        &
#ifdef HAMMOZ
                              loini, pxtbound)
#else
                              loini)
#endif

    ! *xt_massfix* corrects the tendencies of each column to
    !              conserve mass
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-Met,         11/2003
    !
    ! Method:
    ! -------
    ! To keep the model as closely as possible to the intended physical
    ! tendency changes, the correction of the tendencies is
    ! imposed proportionally to the absolute values of the
    ! applied tendency in each layer.
    !
    ! Total mass error tendency [kg m-2 s-1]:
    !
    !    dxtdt(jl)=SUM(xtte(jl,jk)*dpg(jl,jk))+xtbound(jl,jk)
    !
    !    where xtbound is a boundary condition, i.e. deposition flux
    !
    ! Correction proportional to tendency [kg m-2 s-1]:
    !
    !                         - ABS(xtte(jl,jk)*dpg(jl,jk))
    !    dxtfix(jl,jk)   = --------------------------------- * dxdt(jl)
    !                       SUM(ABS(xtte(jl,jk)*dpg(jl,jk)))
    !
    ! Resulting in a corrective tendency [kg kg-1 s-1]:
    !
    !    zxttefix(jl,jk) = dxtfix(jl,jk) / dpg(jl,jk)
    !
    !
    ! Arguments:
    ! ----------
    !
    ! pxtte     = cumulative tendency for all tracers [kg kg-1 s-1]
    !
    ! Usage:
    ! ------
    ! Call twice: first with loini=.TRUE. to store the old tendency
    !             then to fix the mass of processes that have modified
    !             pxtte meanwhile with the current pxtte and loini=.FALSE.

  USE mo_tracdef,             ONLY: trlist             ! tracer info variable
  USE mo_physical_constants,  ONLY: grav

    IMPLICIT NONE


    !--- Arguments:
    INTEGER,      INTENT(IN)    :: kproma, kbdim, klev, klevp1, ktrac, krow
    LOGICAL,      INTENT(IN)    :: loini
    REAL(wp),     INTENT(IN)    :: papp1(kbdim,klev),       &
                                   paphp1(kbdim,klevp1)
    REAL(wp),     INTENT(INOUT) :: pxtte(kbdim,klev,ktrac)
!>>SF
#ifdef HAMMOZ
    REAL(wp),     INTENT(IN)    :: pxtbound(kbdim,ktrac) ! boundary condition (wet deposition) [kg m-2 s-1]
#endif
!<<SF


    !--- Local Variables:
    INTEGER               :: jl, jk, jt
    REAL(wp)              :: zeps, zxttefix
    REAL(wp)              :: zdxtdt(kbdim), zdxtdtsum(kbdim)
    REAL(wp)              :: zxtte(kbdim,klev), zdpg(kbdim,klev)

    !--- 1) Initialization mode:

    IF (.NOT. ALLOCATED(zxtte_old)) ALLOCATE(zxtte_old(kbdim,klev,ktrac))
    IF (loini) THEN
       zxtte_old(1:kproma,:,:) = pxtte(1:kproma,:,:)
       RETURN             ! initialisation done: return from subroutine
    END IF

    !--- 2) Mass fix mode:

    zeps=EPSILON(1.0_wp)

    !--- 2.1) Calculate auxiliary variable dp/g :
    !--- Uppermost level:
    zdpg(1:kproma,1)=2._wp*(paphp1(1:kproma,2)-papp1(1:kproma,1))/grav
    !--- Other levels:
    DO jk=2, klev
       zdpg(1:kproma,jk)=(paphp1(1:kproma,jk+1)-paphp1(1:kproma,jk))/grav
    END DO

    !--- 2.2) Apply mass fixer
    DO jt=1, ktrac
!!     IF(trlist%ti(jt)%nwetdep > 0) THEN
       IF(trlist%ti(jt)%nconvmassfix > 0) THEN

          zdxtdt(1:kproma)   =0.0_wp
          zdxtdtsum(1:kproma)=0.0_wp

          !--- Accumulated tendency since initialization:
          zxtte(1:kproma,:)=pxtte(1:kproma,:,jt)-zxtte_old(1:kproma,:,jt)

          !--- Calculate vertically integrated mass error tendency [kg m-2 s-1]:
          DO jk=1, klev
             DO jl=1, kproma

                zdxtdt(jl)=zdxtdt(jl)      +zxtte(jl,jk)*zdpg(jl,jk)

                zdxtdtsum(jl)=zdxtdtsum(jl)+ABS(zxtte(jl,jk)*zdpg(jl,jk))

             END DO
          END DO

!>>SF
#ifdef HAMMOZ
          zdxtdt(1:kproma)=zdxtdt(1:kproma)+pxtbound(1:kproma,jt) !SF restore boundary cond. calc.
#endif
!<<SF

          DO jk=1, klev
             DO jl=1, kproma
                IF (ABS(zdxtdtsum(jl)) > zeps) THEN

                   zxttefix=-((ABS(zxtte(jl,jk)*zdpg(jl,jk))/zdxtdtsum(jl))*zdxtdt(jl))/zdpg(jl,jk)

                   pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zxttefix

                END IF
             END DO ! kproma
          END DO ! klev

       END IF ! nwetdep
    END DO ! ktrac


  END SUBROUTINE xt_conv_massfix

  SUBROUTINE xt_borrow(kproma, kbdim,  klev,  klevp1, ktrac,       &
                       papp1,  paphp1,                             &
                       pxtm1,  pxtte                               )

    ! *xt_borrow* borrowing scheme to correct for negative
    !             tracer masses conserveing mass within the
    !             column
    ! Note by SF: actually, the scheme is also good for non-mass tracers, therefore it is
    ! now applied to all transported tracers.
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-Met,         06/2004
    !
    ! Method:
    ! -------
    ! If negative tracer values occur the scheme is iteratively
    ! borrowing tracer mixing ratios from the grid-box below.
    ! If the lowest layer can not compensate for accumulated remaining
    ! corrections, the procedure is repeated from bottom to top.
    !
    ! Restrictions:
    ! -------------
    ! Columns with a negative total integrated tracer
    ! content are set to zero despite the associated mass
    ! error.

    USE mo_physical_constants,    ONLY: grav
    USE mo_time_control, ONLY: time_step_len
    USE mo_tracdef,      ONLY: trlist, AEROSOLMASS, GAS, AEROSOLNUMBER
    USE mo_advection,    ONLY: no_advection !SF #246
    USE mo_submodel,     ONLY: lmoz  !csld #330

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: kproma, kbdim, klev, klevp1, ktrac
    REAL(wp),INTENT(IN)    :: papp1(kbdim,klev), paphp1(kbdim,klevp1)
    REAL(wp),INTENT(IN)    :: pxtm1(kbdim,klev,ktrac)
    REAL(wp),INTENT(INOUT) :: pxtte(kbdim,klev,ktrac)

    !--- Local variables:

    INTEGER :: jl, jk, jt

    REAL(wp):: zxtp1, ztmst, zeps

    REAL(wp):: zxtbor(kbdim)

    REAL(wp):: zdpg(kbdim,klev)

    LOGICAL :: lborrtrac !csld #330 : determines if the tracer can be borrowed or not

    !--- 0) Initializations:

    ztmst=time_step_len
    zeps=10._wp*EPSILON(1.0_wp)

    !--- 1) Calculate auxiliary variable dp/g :

    !--- Uppermost level:

    zdpg(1:kproma,1)=2._wp*(paphp1(1:kproma,2)-papp1(1:kproma,1))/grav

    !--- Other levels:

    DO jk=2, klev

       zdpg(1:kproma,jk)=(paphp1(1:kproma,jk+1)-paphp1(1:kproma,jk))/grav

    END DO

    !--- 2) Borrowing scheme:
    
    DO jt=1, ktrac

       !>>csld #330
       ! use xt_borrow on all aerosols species
       lborrtrac = trlist%ti(jt)%nphase == AEROSOLMASS .OR. trlist%ti(jt)%nphase == AEROSOLNUMBER

       IF (.NOT. lmoz) THEN    ! in case of "ham only" runs, use of xt_borrow on gaz species as well
          lborrtrac = lborrtrac .OR. (trlist%ti(jt)%nphase == GAS)
       END IF

       IF ( (trlist%ti(jt)%ntran /= no_advection)       .AND. &  !SF #246 exclude non-advected tracers
            lborrtrac) THEN
       !<<csld #330

          !--- Start borrowing scheme:
    
          !--- Integrate from top to bottom:

          zxtbor(1:kproma)=0.0_wp

          DO jk=1, klev
             DO jl=1, kproma

                !--- Check if corrected mass, including the fix of the layer above,
                !    yields a negative mixing ratio:
                !    (Convert mass correction converted to [kg kg-1] in current layer)

                zxtp1=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst + zxtbor(jl)/zdpg(jl,jk)

                IF ( zxtp1 > 0.0_wp ) THEN

                   !--- Subtract corrected tracer mass from current layer:
                   !    (zxtbor is negative)
       
                   pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zxtbor(jl)/(zdpg(jl,jk)*ztmst)
    
                   !--- Reset mass correction:
    
                   zxtbor(jl)=0.0_wp

                ELSE

                   !--- Adjust tendency to yield zero:

                   pxtte(jl,jk,jt)=-pxtm1(jl,jk,jt)/ztmst

                   !--- Add correcting mass mixing ratio and convert to [kg m-2]:
                   !    (implicit summation due to the inclusion
                   !     of zxtbor(jk-1) in zxtp1)

                   zxtbor(jl)=zxtp1*zdpg(jl,jk)

                END IF

             END DO
          END DO

          !--- If surface layer cannot compensate accumulated correction:
          !    Iterate from bottom to top

          DO jk=klev, 1, -1
             DO jl=1, kproma
                IF (zxtbor(jl) < -zeps) THEN

                   !--- Check if corrected mass, including the fix of the layer below,
                   !    yields a negative mixing ratio:
                   !    (Convert mass correction converted to [kg kg-1] in current layer)

                   zxtp1=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst + zxtbor(jl)/zdpg(jl,jk)

                   IF ( zxtp1 > 0.0_wp ) THEN

                      !--- Subtract corrected tracer mass from current layer:
                      !    (zxtbor is negative)

                      pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zxtbor(jl)/(zdpg(jl,jk)*ztmst)

                      !--- Reset mass correction:

                      zxtbor(jl)=0.0_wp

                   ELSE

                      !--- Adjust tendency to yield zero:

                      pxtte(jl,jk,jt)=-pxtm1(jl,jk,jt)/ztmst
                   
                      !--- Add correcting mass mixing ratio and convert to [kg m-2]:
                      !    (implicit summation due to the inclusion
                      !     of zxtbor(jk-1) in zxtp1)
                   
                      zxtbor(jl)=zxtp1*zdpg(jl,jk) 
                   
                   END IF
                   
                END IF
             END DO 
          END DO
       
       END IF 
    END DO
          
  END SUBROUTINE xt_borrow
          




  
!!mgs!!   -------   obsolete code  (perhaps re-use some parts of this to enhance burden diag? ------  
!!mgs!!   SUBROUTINE trastat
!!mgs!! 
!!mgs!!   ! Description:
!!mgs!!   !
!!mgs!!   ! Prints out accumulated mass budgets for tracers at
!!mgs!!   ! the end of a run
!!mgs!! 
!!mgs!!   USE mo_control,        ONLY: ngl
!!mgs!!   USE mo_mpi,            ONLY: p_sum, p_communicator_d, p_pe, p_io
!!mgs!! 
!!mgs!!   !  Local scalars: 
!!mgs!!   REAL(wp) :: zmglob, zmnhk, zmshk, zmstrat, zmtrop, zqcount
!!mgs!!   INTEGER ::  jt
!!mgs!! 
!!mgs!!   !  Local arrays: 
!!mgs!!   REAL(wp) :: zmstratn(ntrac+1), zmstrats(ntrac+1), zmtropn(ntrac+1),    &
!!mgs!!               zmtrops(ntrac+1)
!!mgs!! 
!!mgs!!   !  Intrinsic functions 
!!mgs!!   INTRINSIC SUM
!!mgs!! 
!!mgs!! 
!!mgs!!   !  Executable statements 
!!mgs!! 
!!mgs!!   zqcount = 1.0_wp/(REAL(icount,wp))
!!mgs!!   IF (p_pe == p_io) THEN
!!mgs!!     CALL message('',' Tracer mass budget:')
!!mgs!!     CALL message('',separator)
!!mgs!!     CALL message('',' Averaged mass budgets in [kg] ')
!!mgs!!     CALL message('',' global   n-hem  s-hem  tropo  strat  n-tro s-tro  n-str s-str ')
!!mgs!!   ENDIF
!!mgs!!   DO jt = 1, ntrac + 1
!!mgs!!      zmtropn(jt)  = SUM(tropm(1:ngl/2,jt))  * zqcount
!!mgs!!      zmstratn(jt) = SUM(stratm(1:ngl/2,jt)) * zqcount
!!mgs!!      zmtrops(jt)  = SUM(tropm(ngl/2+1:ngl,jt))  * zqcount
!!mgs!!      zmstrats(jt) = SUM(stratm(ngl/2+1:ngl,jt)) * zqcount
!!mgs!!   END DO
!!mgs!!   zmtropn  = p_sum (zmtropn,  p_communicator_d)
!!mgs!!   zmstratn = p_sum (zmstratn, p_communicator_d)
!!mgs!!   zmtrops  = p_sum (zmtrops,  p_communicator_d)
!!mgs!!   zmstrats = p_sum (zmstrats, p_communicator_d)
!!mgs!! 
!!mgs!!   IF (p_pe == p_io) THEN    
!!mgs!!     DO jt = 1, ntrac + 1
!!mgs!!        zmnhk   = zmtropn(jt)  + zmstratn(jt)
!!mgs!!        zmshk   = zmtrops(jt)  + zmstrats(jt)
!!mgs!!        zmtrop  = zmtropn(jt)  + zmtrops(jt)
!!mgs!!        zmstrat = zmstratn(jt) + zmstrats(jt)
!!mgs!!        zmglob  = zmnhk + zmshk
!!mgs!! 
!!mgs!!        IF (jt <= ntrac) THEN
!!mgs!!          WRITE (message_text,'(a,i2,9e9.2)')                    &
!!mgs!!               ' Tracer: ', jt, zmglob, zmnhk,                   &
!!mgs!!               zmshk, zmtrop, zmstrat, zmtropn(jt), zmtrops(jt), &
!!mgs!!               zmstratn(jt), zmstrats(jt)
!!mgs!!          CALL message('',message_text)
!!mgs!!        ELSE
!!mgs!!          WRITE (message_text,'(a   ,9e9.2)')                    &
!!mgs!!               ' Air mass: ', zmglob, zmnhk,                     &
!!mgs!!               zmshk, zmtrop, zmstrat, zmtropn(jt), zmtrops(jt), &
!!mgs!!               zmstratn(jt), zmstrats(jt)
!!mgs!!          CALL message('',message_text)
!!mgs!!        END IF
!!mgs!!     END DO
!!mgs!!   END IF
!!mgs!!   
!!mgs!!   END SUBROUTINE trastat
!!mgs!! 
  
END MODULE mo_tracer_processes

