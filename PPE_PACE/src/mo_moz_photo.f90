!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_photo
!!
!! \brief
!!  Module containing general set-up and output definition for photolysis schemes
!!  Can be used by waccm photo and fastj
!!
!! \author Martin G. Schultz (FZ-Juelich)
!! \author Luca Pozzoli (JRC Ispra)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - MGS and LP: original version (2008)
!!  - MGS: adaptation for ECHAM6-HAMMOZ (2009-08-31)
!!
!! \limitations
!!  fast-J not yet implemented
!!
!! \details
!!
!! \bibliographic_references
!!  none
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_moz_photo

  USE mo_kind,           ONLY: dp
  USE mo_memory_base,    ONLY: t_stream
  USE mo_time_event,     ONLY: io_time_event
  USE mo_submodel_diag,  ONLY: vmem3d
  USE mo_moz_mods,       ONLY: phtcnt, rxt_tag_lst    ! number and names of photo reactions defined in MOZART

  IMPLICIT NONE

  PRIVATE

  PUBLIC     :: photo            ! output stream for photolysis frequencies and aux. information
  PUBLIC     :: photo_initialize, init_photo_stream
  PUBLIC     :: lpost_photo, tinterval_photo,   &
                nphotovars, nphotoreac, cvarn_photo, cphotovars, clonglabels, jdefined
  PUBLIC     :: o3_rad           ! save ozone field from ECHAM radiation calculation
  PUBLIC     :: od_w, od_i, lwp, iwp, erw, eri,  &
                od_aer, od_aerm1, od_aerm2, od_aerm3, od_aerm4, od_aerm5, &
                od_aerm6, od_aerm7, aclcac, photo_mem



  ! photo stream, choice of output variables via mozctl namelist photovarn
  ! copied locally into cvarn_photo
  TYPE (t_stream), POINTER :: photo   ! output stream
  LOGICAL                  :: lpost_photo
  TYPE(io_time_event), SAVE:: tinterval_photo
  INTEGER, PARAMETER       :: nphotoextra=15
  INTEGER, PARAMETER       :: nphotoreac=phtcnt, maxphotovars=nphotoextra+nphotoreac
  INTEGER, SAVE            :: nphotovars                 ! number of photovars defined in cphotovars
  CHARACTER(LEN=32)        :: cphotovars(maxphotovars)   ! list of all possible output names
  CHARACTER(LEN=32)        :: cvarn_photo(maxphotovars)  ! list of selected output names
  CHARACTER(LEN=60)        :: clonglabels(maxphotovars)  ! list of long names for output
  LOGICAL                  :: jdefined(nphotoreac)       ! manage list of photorates actually computed


! to save optinal fields in the stream 'photo' 
  REAL(dp), POINTER :: od_w(:,:,:)      ! optical depth of water clouds
  REAL(dp), POINTER :: od_i(:,:,:)      ! optical depth of ice clouds
  REAL(dp), POINTER :: lwp(:,:,:)      ! liquid water path
  REAL(dp), POINTER :: iwp(:,:,:)      ! ice    water path
  REAL(dp), POINTER :: erw(:,:,:)      ! effective radius water clouds
  REAL(dp), POINTER :: eri(:,:,:)      ! effective radius ice   clouds
  REAL(dp), POINTER :: od_aer(:,:,:)     ! total aerosol optical depth from HAM
  REAL(dp), POINTER :: od_aerm1(:,:,:)     ! optical depth for aerosol mode 1 (HAM)
  REAL(dp), POINTER :: od_aerm2(:,:,:)     ! optical depth for aerosol mode 2 (HAM)
  REAL(dp), POINTER :: od_aerm3(:,:,:)     ! optical depth for aerosol mode 3 (HAM)
  REAL(dp), POINTER :: od_aerm4(:,:,:)     ! optical depth for aerosol mode 4 (HAM)
  REAL(dp), POINTER :: od_aerm5(:,:,:)     ! optical depth for aerosol mode 5 (HAM)
  REAL(dp), POINTER :: od_aerm6(:,:,:)     ! optical depth for aerosol mode 6 (HAM)
  REAL(dp), POINTER :: od_aerm7(:,:,:)     ! optical depth for aerosol mode 7 (HAM)
  REAL(dp), POINTER :: aclcac(:,:,:)    ! cloud cover fraction profiles
  REAL(dp), POINTER :: o3_rad(:,:,:)    ! ozone field from radiation
  TYPE (vmem3d)     :: photo_mem(nphotoreac)


  CONTAINS 

!--------------------------------------------------------------------------------------------

  SUBROUTINE photo_initialize (photovarn, nphotovarn)

  USE mo_exception,             ONLY : message, message_text, em_error, em_debug
  USE mo_util_string,           ONLY : tolower
  USE mo_time_control,          ONLY : putdata

  !###debug###
  USE mo_mpi,               ONLY : p_parallel_io

  !-- Arguments
  CHARACTER (LEN=32), INTENT(in)   :: photovarn(:)
  INTEGER, INTENT(in)              :: nphotovarn

  !-- Local variables
  INTEGER          :: jr, jmax

  !-- Set list of possible variable names in photo stream
  cphotovars(:) = ''
  clonglabels(:) = ''
  jdefined(:) = .false.     ! This list must be changed by mo_fastj or equivalent...
!sschr:
  !-- Add general properties from ECHAM
  cphotovars( 1) = 'od_w'       ! optical depth of water clouds
  cphotovars( 2) = 'od_i'       ! optical depth of ice clouds
  cphotovars( 3) = 'lwp'       ! liquid water path
  cphotovars( 4) = 'iwp'       ! ice    water path
  cphotovars( 5) = 'erw'       ! effective radius water clouds
  cphotovars( 6) = 'eri'       ! effective radius ice   clouds
  cphotovars( 7) = 'aclcac'       ! cloud cover fraction profiles
  cphotovars( 8) = 'od_aer'       ! total aerosol optical depth
  nphotovars = 8

!sschr: fastj not fully implemeted (lham is a questionable switch)
  !-- Add HAM specific variables for fastj
! IF ( lham ) THEN
!   cphotovars(nphotovars+1) = 'od_aerm1'       ! optical depth for aerosol mode 1 (HAM)
!   cphotovars(nphotovars+2) = 'od_aerm2'       ! optical depth for aerosol mode 2 (HAM)
!   cphotovars(nphotovars+3) = 'od_aerm3'       ! optical depth for aerosol mode 3 (HAM)
!   cphotovars(nphotovars+4) = 'od_aerm4'       ! optical depth for aerosol mode 4 (HAM)
!   cphotovars(nphotovars+5) = 'od_aerm5'       ! optical depth for aerosol mode 5 (HAM)
!   cphotovars(nphotovars+6) = 'od_aerm6'       ! optical depth for aerosol mode 6 (HAM)
!   cphotovars(nphotovars+7) = 'od_aerm7'       ! optical depth for aerosol mode 7 (HAM)
!   nphotovars = nphotovars + 7
! END IF

  !-- Add MOZ specific variables (photolysis frequencies)
  DO jr = 1,nphotoreac
     cphotovars(nphotovars+jr) = TRIM(rxt_tag_lst(jr))
  END DO
  nphotovars = nphotovars + nphotoreac
 
  ! initialize cvarn_photo list of actual output variables in photo stream
  cvarn_photo(:) = ''
  ! copy requested variables from photovarn argument
  IF (nphotovarn > maxphotovars) THEN
    WRITE (message_text,'(a,i0,a,i0)') 'Number of requested photovariables (', nphotovarn, &
                            ') exceeds maximum allowed value of ',maxphotovars
    CALL message('photo_initialize', message_text, level=em_error)
  END IF
  jmax = MAX(maxphotovars, nphotovarn)
  cvarn_photo(1:jmax) = photovarn(1:jmax)
  
! --- interprete cvarn_photo string: if first entry is 'all', then expand list
  IF (tolower(TRIM(cvarn_photo(1))) == 'all') cvarn_photo(:) = cphotovars(:)
! --- if first entry is 'default' then set default entries
  IF (tolower(TRIM(cvarn_photo(1))) == 'default') THEN
!sschr: remove od_ from default output
!   cvarn_photo(1) = 'od_w'
!   cvarn_photo(2) = 'od_i'
!   cvarn_photo(3) = 'od_aer'
    cvarn_photo(1) = 'jo3_a'
    cvarn_photo(2) = 'jo3_b'
    cvarn_photo(3) = 'jno2'
  END IF

! DO jr=1,jmax
!    WRITE (message_text,'(2a)') 'photo_initialize: cvarn_photo = ',TRIM(cvarn_photo(jr))
!    CALL message('photo_initialize', message_text, level=em_debug)
! ENDDO

  !-- set time interval for output stream photo
  tinterval_photo = putdata

  DO jr=1,nphotovars 
     WRITE (message_text,'(2a)') 'photo_initialize : cphotovars = ', TRIM(cphotovars(jr))
     CALL message('photo_initialize', message_text, level=em_debug)
  ENDDO

  END SUBROUTINE photo_initialize


  SUBROUTINE init_photo_stream 
  !--------------------------------------------------------------------------
  ! This subroutine defines a new output streams containing different
  ! variables in order to analyse the output of the coupled echam/mozart/fastj 
  ! model.
  !--------------------------------------------------------------------------
  ! stream _photo:
  ! name        description                               dimension module
  !
  ! od_w, od_i  cloud optical depth                                  rad_int
  ! fjo1d       photorate O3+hv->O1D+O2                     1/s      photo_fast
  ! fjo3p       photorate O3+hv->O3P+O2                     1/s      photo_fast
  ! fjno2       photorate NO2+hv->NO+O                      1/s      photo_fast
  !
  !--------------------------------------------------------------------------
  ! Authors:
  ! --------
  ! Luca Pozzoli, EPFL                       01/2004
  ! M. Schultz,  FZ Juelich,  Jan 2009
  !
  !


  USE mo_memory_base,    ONLY: new_stream, add_stream_element, &
                               default_stream_setting, AUTO
  USE mo_exception,      ONLY: message, message_text, em_warn, em_error, em_param, em_debug
  USE mo_string_utls,    ONLY: st1_in_st2_proof

  IMPLICIT NONE
  !---------------
  !local variables
  !---------------

  INTEGER    :: ii, ierr, ivar, iphotoextra
  LOGICAL    :: lpost


  !create photo stream and add variables
  lpost_photo = LEN_TRIM(cvarn_photo(1)) > 1    ! only output stream if at least 
                                                ! one variable requested
  CALL new_stream (photo,'photo', lpost=lpost_photo, lrerun=.FALSE., &
                   interval=tinterval_photo)

  CALL default_stream_setting (photo ,lrerun    = .FALSE.,   &
                               contnorest = .TRUE.,   &
                               laccu = .FALSE.,   &
                               lpost     = .TRUE.,   &
                               table      = 199,     &
                               code = AUTO )

   !add variables to stream
   IF (st1_in_st2_proof( cvarn_photo, cphotovars(1:nphotovars), ierr=ierr )) THEN
     IF (ierr > 0) THEN
       WRITE (message_text,'(a,i0,a,a)') 'Invalid variable in photo stream. ierr = ', ierr, &
                               ', variable = ',cvarn_photo(ierr)
       CALL message('init_photo_stream', message_text, level=em_error)
     END IF
     DO ivar=1,nphotovars
       WRITE (message_text,'(a,i0,a,a)') 'cvarn_photo(',ivar,') = ',cvarn_photo(ivar)
       CALL message('', message_text, level=em_debug)
     END DO
   END IF

   ! ozone field from radiation (no file output)
   lpost = .false.
   CALL add_stream_element (photo, 'o3_rad', o3_rad, &
                            longname = 'ozone mass mixing ratio from last radiation step', &
                            units = '1', lpost = lpost)

   lpost = st1_in_st2_proof( 'od_w', cvarn_photo)
   CALL add_stream_element (photo, 'od_w', od_w, &                                 
                            longname = 'optical depth of water clouds', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'od_i', cvarn_photo)
   CALL add_stream_element (photo, 'od_i', od_i, &                                 
                            longname = 'optical depth of ice clouds', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'lwp', cvarn_photo)
   CALL add_stream_element (photo, 'lwp', lwp, &                                 
                            longname = 'liquid water path', &
 !!                         standardname = '', &
                            units = 'kg m-2', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'iwp', cvarn_photo)
   CALL add_stream_element (photo, 'iwp', iwp, &                                 
                            longname = 'ice water path', &
 !!                         standardname = '', &
                            units = 'kg m-2', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'erw', cvarn_photo)
   CALL add_stream_element (photo, 'erw', erw, &                                 
                            longname = 'effective droplet radius water clouds', &
 !!                         standardname = '', &
                            units = 'm', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'eri', cvarn_photo)
   CALL add_stream_element (photo, 'eri', eri, &                                 
                            longname = 'effective droplet radius ice clouds', &
 !!                         standardname = '', &
                            units = 'm', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'od_aer', cvarn_photo)
   CALL add_stream_element (photo, 'od_aer', od_aer, &                                 
                            longname = 'optical depth of aerosol', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'od_aerm1', cvarn_photo)
   CALL add_stream_element (photo, 'od_aerm1', od_aerm1, &                                 
                            longname = 'optical depth of nucleation mode aerosol', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'od_aerm2', cvarn_photo)
   CALL add_stream_element (photo, 'od_aerm2', od_aerm2, &                                 
                            longname = 'optical depth of aitken soluble aerosol', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'od_aerm3', cvarn_photo)
   CALL add_stream_element (photo, 'od_aerm3', od_aerm3, &                                 
                            longname = 'optical depth of accumulation soluble aerosol', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'od_aerm4', cvarn_photo)
   CALL add_stream_element (photo, 'od_aerm4', od_aerm4, &                                 
                            longname = 'optical depth of coarse soluble aerosol', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'od_aerm5', cvarn_photo)
   CALL add_stream_element (photo, 'od_aerm5', od_aerm5, &                                 
                            longname = 'optical depth of aitken insoluble aerosol', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'od_aerm6', cvarn_photo)
   CALL add_stream_element (photo, 'od_aerm6', od_aerm6, &                                 
                            longname = 'optical depth of accumulation insoluble aerosol', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'od_aerm7', cvarn_photo)
   CALL add_stream_element (photo, 'od_aerm7', od_aerm7, &                                 
                            longname = 'optical depth of coarse insoluble aerosol', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   lpost = st1_in_st2_proof( 'aclcac', cvarn_photo)
   CALL add_stream_element (photo, 'aclcac', aclcac, &                                 
                            longname = 'cloud fraction', &
 !!                         standardname = '', &
                            units = '1', lpost = lpost)
 
   iphotoextra = nphotovars-nphotoreac    !!! bug fix. Used to be ...+1
   DO ii=1,nphotoreac
      ivar=iphotoextra+ii
      lpost = st1_in_st2_proof( cphotovars(ivar), cvarn_photo)
      IF ( LEN_TRIM(clonglabels(ivar)) > 1 ) THEN
         CALL add_stream_element (photo, cphotovars(ivar), photo_mem(ii)%ptr, & 
                                  longname = TRIM(clonglabels(ivar)), &
                                  units = 's-1', lpost = lpost)
      ELSE
         CALL add_stream_element (photo, cphotovars(ivar), photo_mem(ii)%ptr, & 
                                  units = 's-1', lpost = lpost)
      END IF
   END DO
    
END SUBROUTINE init_photo_stream

END MODULE mo_moz_photo

