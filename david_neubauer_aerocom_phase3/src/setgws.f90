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
SUBROUTINE setgws

  !- Description:
  !
  !  Modify preset variables of module mo_gwspectrum which control
  !  the configuration of the Hines gravity waves parameterization 
  !
  !-  Method:
  !
  !  Read the gwsclt namelist and modify constants.
  !
  !  *setgws* is called from *initialize*.
  !
  !- Authors:
  !
  !   E. Manzini, MPI,  January 2002
  !     (re-write from M. Charron original)
  !   H. Schmidt/S. Misios, MPI, April 2010
  !     (latitude dependent GW source enabled)
  !

  USE mo_kind           ,ONLY: dp
  USE mo_mpi            ,ONLY: p_parallel, p_parallel_io, p_bcast, &
                               p_io
  USE mo_exception      ,ONLY: finish, message, message_text
  USE mo_gwspectrum     ,ONLY: lextro,lfront,lozpr,iheatcal,       &
                               rmscon,emiss_lev,kstar,m_min,       &
                               rms_front,front_thres,pcrit,pcons,  &
                               naz,slope,smco, lrmscon_lat,        &
                               rmscon_lat, rmscon_lo, rmscon_hi,   &
                               lat_rmscon_lo, lat_rmscon_hi
  USE mo_control        ,ONLY: nn, ngl, nlev
  USE mo_namelist       ,ONLY: open_nml, position_nml, POSITIONED
  USE mo_gaussgrid      ,ONLY: philat

  IMPLICIT NONE
 
  ! Local variable
  INTEGER :: ierr, jl, jstart, jlast, inml, iunit

  INCLUDE 'gwsctl.inc'

  ! Executable statements 

  ! Set latitude dependent GW source
  !
  ! aim: In particular in coupled (atmosphere/ocean) model runs wave forcing
  ! is sometimes insufficient to produce a QBO. A general increase of the GW
  ! source (rmscon) may have adverse effects on the stratospheric
  ! mid-latitude circulation. This option allows to increase GW forcing only
  ! in low latitudes. A value of rmscon_lo is used for latitudes below
  ! +/-lat_rmscon_lo, a value of rmscon_hi above +/- lat_rmscon_hi. Between
  ! these latitude boundaries an interpolation is performed.
  ! Here, default values are given that may be overwritten by setting parameters
  ! in the namelist gwsctl.
  ! This option will only be used if lrmscon_lat=.true. (default is .false.)
  !  attention: may be overwritten if lfront or lozpr true

!>>SF Adaptation for running @ T63L31. <-- Validated by Hauke Schmidt on 2016-05-25
  !SForig IF (nn==31 .OR. (nn==63 .AND. nlev==47))  THEN
  IF (nn==31 .OR. (nn==63 .AND. (nlev==31 .OR. nlev==47)))  THEN
!<<SF
    lrmscon_lat   = .false.
    lat_rmscon_lo =  5.0_dp
    lat_rmscon_hi = 10.0_dp
    rmscon_lo     =  1.0_dp
    rmscon_hi     =  1.0_dp
  ELSE IF (nn==63 .AND. nlev==95)  THEN
    ! note: values had been tuned for echam-6.1.01 (pre mo_midatm-bugfix)
    lrmscon_lat   = .true.
    lat_rmscon_lo =  5.0_dp
    lat_rmscon_hi = 10.0_dp
    rmscon_lo     =  1.2_dp
    rmscon_hi     =  1.0_dp
  ELSE IF (nn==127) THEN
    lrmscon_lat   = .true.
    lat_rmscon_lo =  5.0_dp
    lat_rmscon_hi = 10.0_dp
    rmscon_lo     =  1.07_dp
    rmscon_hi     =  1.0_dp
  ELSE IF (nn==255) THEN
    lrmscon_lat   = .true.
    lat_rmscon_lo =  5.0_dp
    lat_rmscon_hi = 10.0_dp
    rmscon_lo     =  1.07_dp
    rmscon_hi     =  1.0_dp
  ELSE
    CALL finish('setgws','Truncation not supported')
  END IF

  ! preset emiss_lev for L199
  ! ================================
  IF (nlev==199) THEN
     emiss_lev = 26
  ELSE                ! default
     emiss_lev = 10
  ENDIF
  !!IF ( (nlev /= 199) .AND. (nlev /= 47) .AND. (nlev /= 95)) THEN
  !!  WRITE(message_text,'(a,i0)') &
  !!       'emiss_lev not implemented for nlev=', nlev 
  !!  CALL message('',message_text)
  !!  CALL finish('setgws','Run terminated')
  !!ENDIF

  ! Read gwsctl namelist to modify mo_gwspectrum
  ! ============================================

  IF (p_parallel_io) THEN
    inml = open_nml('namelist.echam')
     iunit = position_nml ('GWSCTL', inml, status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
       READ (iunit, gwsctl)
     END SELECT
  ENDIF
  IF (p_parallel) THEN
     CALL p_bcast (lextro, p_io)
     CALL p_bcast (lfront, p_io)
     CALL p_bcast (lozpr, p_io)
     CALL p_bcast (iheatcal, p_io)
     CALL p_bcast (rmscon, p_io)
     CALL p_bcast (emiss_lev, p_io)
     CALL p_bcast (kstar, p_io)
     CALL p_bcast (m_min, p_io)
     CALL p_bcast (rms_front, p_io)
     CALL p_bcast (front_thres, p_io)
     CALL p_bcast (pcrit, p_io)
     CALL p_bcast (pcons, p_io)
     CALL p_bcast (lrmscon_lat,p_io )
     CALL p_bcast (lat_rmscon_lo,p_io )
     CALL p_bcast (lat_rmscon_hi,p_io )
     CALL p_bcast (rmscon_lo,p_io )
     CALL p_bcast (rmscon_hi,p_io )
  ENDIF

  ! Check lextro
  ! ============
  IF (.NOT.lextro) THEN
       CALL message('','lextro = .FALSE. --> no hines gw spectrum')
  ENDIF

  IF (lextro) THEN

     ! Initialization for Hines gravity wave parameterization
     ! =======================================================

     !  check that things are setup correctly and log error if not


     IF (lfront) THEN
           CALL message('','lfront = .TRUE.  --> gravity waves from fronts')
        IF (nn==31)  THEN
           front_thres  = 0.10_dp
        ENDIF
        IF (nn /= 31) THEN
          WRITE(message_text,'(a,i0)') &
               'lfront = .TRUE. not implemented for nn=', nn 
          CALL message('',message_text)
          CALL finish('setgws','Run terminated')
        ENDIF
     ENDIF

     IF (lozpr) THEN
        CALL message('','lozpr = .TRUE. not implemented for echam6')
        CALL finish('setgws','Run terminated') 
     ENDIF

     SELECT CASE(iheatcal)
     CASE(0)
       CALL message('','standard echam6: momentum flux deposition ONLY')
     CASE(1)
       CALL message('','iheatcal = 1 : momentum deposion, heating and diff coeff')
        IF ( ABS(slope-1._dp) > EPSILON(1._dp) ) THEN
          CALL message('','Heat and coeff for slope not 1 not implemented')
          CALL finish('setgws','Run terminated')
        ENDIF
     END SELECT

     IF ( naz /= 8 )  THEN
       CALL message('','naz /= 8 not implemented')
       CALL finish('setgws','Run terminated')
     ENDIF

     IF (ABS(slope-1._dp) > EPSILON(1._dp)  ) THEN
       CALL message('','slope /= 1 not implemented')
       CALL finish('setgws','Run terminated')
     ENDIF

     IF ( (m_min > EPSILON(1._dp)) .AND. (ABS(slope-1._dp) > EPSILON(1._dp)) ) THEN
       CALL message('','m_min /= 0 for slope /= 1 not implemented')
       CALL finish('setgws','Run terminated')
     ENDIF

     IF ( smco < 1._dp ) THEN
       CALL message('','smco less than 1 not implemented') 
       CALL finish('setgws','Run terminated')
     ENDIF
     !

     ! Set latitude dependent GW source
     !
     ! aim: see above

     IF ( lrmscon_lat ) THEN

       IF (p_parallel_io) THEN
         WRITE(message_text,*)"Latidudinal dependence of rmscon is enabled."
         CALL message('setgws',message_text)
         CALL message('','lrmscon_lat =.TRUE.  --> increase GW forcing in low lat.')
         WRITE (message_text, '(a,f7.2)') 'lat_rmscon_lo = ', lat_rmscon_lo
         CALL message('setgws', message_text)
         WRITE (message_text, '(a,f7.2)') 'lat_rmscon_hi = ', lat_rmscon_hi
         CALL message('setgws', message_text)
         WRITE (message_text, '(a,f7.2)') 'rmscon_lo = ', rmscon_lo
         CALL message('setgws', message_text)
         WRITE (message_text, '(a,f7.2)') 'rmscon_hi = ', rmscon_hi
         CALL message('setgws', message_text)
       END IF

       IF (lfront .OR. lozpr) THEN
         CALL finish ('setgws', 'Latitude dependent GW source not supported'// &
                                ' in connection with lfront or lozpr set true')
       END IF

       IF (.NOT. ALLOCATED(rmscon_lat)) ALLOCATE(rmscon_lat(ngl))
      
       DO jl=1,ngl/2

         IF ( philat(jl) >= lat_rmscon_hi ) THEN
           rmscon_lat(jl) = rmscon_hi
         ELSE IF ( philat(jl) < lat_rmscon_hi .AND. philat(jl) > lat_rmscon_lo ) THEN
           rmscon_lat(jl) = rmscon_hi + (philat(jl)-lat_rmscon_hi)* &
                           (rmscon_lo-rmscon_hi)/(lat_rmscon_lo-lat_rmscon_hi)
         ELSE
           rmscon_lat(jl) = rmscon_lo
         END IF
  
         !copy the values to the other hemisphere
         rmscon_lat(ngl+1-jl)=rmscon_lat(jl)

       END DO
     
       IF (p_parallel_io) THEN
         WRITE(message_text,*)"The lat dependent rmscon array is:"
         CALL message('setgws',message_text)
         DO jl = 1,ngl,10
            jstart = jl
            jlast  = jl+9
            IF (jlast > ngl) jlast = ngl     
            WRITE (message_text, '(10f6.3)') rmscon_lat(jstart:jlast)
            CALL message ('',message_text)
         ENDDO
       END IF
      
     END IF

     !------------------------------------------------------------
     ! set up for Hammonia:
     !  rmscon       = 1.0
     !  m_minimum    = 3.141593e-4     ! = 2*pi/(20 km)
     !  iheatcal    = 1
     !------------------------------------------------------------

  ENDIF
  !
END SUBROUTINE setgws

