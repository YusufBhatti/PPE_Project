MODULE mo_activ

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  USE mo_submodel_diag, ONLY: vmem3d

  IMPLICIT NONE

  PUBLIC activ_initialize
  PUBLIC activ_updraft
  PUBLIC activ_lin_leaitch
  PUBLIC construct_activ_stream

  PRIVATE

  INTEGER,         PUBLIC :: idt_cdnc, idt_icnc, nfrzmod

  TYPE (t_stream), PUBLIC, POINTER :: activ

  INTEGER,         PUBLIC          :: nw ! actual number of updraft velocity (w) bins 
                                         ! (can be 1 if characteristic updraft is used)

  REAL(dp),        PUBLIC, POINTER :: swat(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: w_large(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: w_turb(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: w_cape(:,:)
  REAL(dp),        PUBLIC, POINTER :: w_sigma(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: reffl(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: reffi(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: na(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qnuc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qaut(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qacc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qfre(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qeva(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qmel(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: icnc_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: icnc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: icnc_instantan(:,:,:) ! Ice crystal number concentration (ICNC), actual instantaneous value [1/m3]
  REAL(dp),        PUBLIC, POINTER :: lwc_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: iwc_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cloud_time(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cliwc_time(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc_burden_acc(:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc_burden(:,:)
  REAL(dp),        PUBLIC, POINTER :: icnc_burden_acc(:,:)
  REAL(dp),        PUBLIC, POINTER :: icnc_burden(:,:)
  REAL(dp),        PUBLIC, POINTER :: burden_time(:,:)
  REAL(dp),        PUBLIC, POINTER :: burdic_time(:,:)
  REAL(dp),        PUBLIC, POINTER :: reffl_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: reffi_acc(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: cloud_cover_duplic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: sice(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: reffl_ct(:,:)
  REAL(dp),        PUBLIC, POINTER :: reffl_time(:,:)
  REAL(dp),        PUBLIC, POINTER :: cdnc_ct(:,:)
  REAL(dp),        PUBLIC, POINTER :: reffi_tovs(:,:)
  REAL(dp),        PUBLIC, POINTER :: reffi_time(:,:)
  REAL(dp),        PUBLIC, POINTER :: iwp_tovs(:,:)

  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: w(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: w_pdf(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: swat_max_strat(:)
  TYPE(vmem3d), PUBLIC, ALLOCATABLE :: swat_max_conv(:)

  REAL(dp)            :: w_min = 0.0_dp       ! minimum characteristic w for activation [m s-1]
  REAL(dp), PARAMETER :: w_sigma_min = 0.1_dp ! minimum value of w standard deviation [m s-1]

  !--- Subroutines:

CONTAINS

  SUBROUTINE activ_updraft(kproma,   kbdim,  klev,    krow, &
                           ptkem1,   pwcape, pvervel, prho, &
                           pw,       pwpdf                  )

    ! *activ_updraft* calculates the updraft vertical velocity
    !                 as sum of large scale and turbulent velocities
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford                 2008
    !
    ! References:
    ! -----------
    ! Lohmann et al., ACP, (2008)
    !

    USE mo_physical_constants, ONLY: grav
    !>>SF #345
    USE mo_cloud_utils,        ONLY: fact_tke
    USE mo_param_switches,     ONLY: ncd_activ
    !<<SF #345
    USE mo_param_switches, ONLY: nactivpdf !ZK

    IMPLICIT NONE

    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow

    REAL(dp), INTENT(out) :: pw(kbdim,klev,nw)        ! stratiform updraft velocity bins, large-scale+TKE (>0.0) [m s-1]
    REAL(dp), INTENT(out) :: pwpdf(kbdim,klev,nw)     ! stratiform updraft velocity PDF

    REAL(dp), INTENT(in)  :: prho(kbdim,klev),      & ! air density
                             ptkem1(kbdim,klev),    & ! turbulent kinetic energy
                             pvervel(kbdim,klev),   & ! large scale vertical velocity [Pa s-1]
                             pwcape(kbdim)            ! CAPE contribution to convective vertical velocity [m s-1]

    REAL(dp) :: zwlarge(kbdim, klev), & ! large-scale vertical velocity [m s-1]
                zwturb(kbdim, klev)     ! TKE-derived vertical velocity or st. dev. thereof [m s-1]

    !--- Large scale vertical velocity in SI units:

    zwlarge(1:kproma,:)      = -1._dp* pvervel(1:kproma,:)/(grav*prho(1:kproma,:))
    w_large(1:kproma,:,krow) = zwlarge(1:kproma,:)

    !--- Turbulent vertical velocity:

    w_turb(1:kproma,:,krow)  = fact_tke*SQRT(ptkem1(1:kproma,:))

    !>>SF #345: correction for the TKE prefactor, in case of Lin & Leaitch scheme only
    IF (ncd_activ == 1) THEN ! Lin & Leaitch scheme
       w_turb(1:kproma,:,krow)  = 1.33_dp*SQRT(ptkem1(1:kproma,:))
    ENDIF
    !<<SF #345

    !--- Convective updraft velocity from CAPE:

    w_cape(1:kproma,krow)  = pwcape(1:kproma) !SF although this is no longer used as a contribution to the
                                              ! convective updraft velocity, this is just kept here
                                              ! for recording it into the activ stream

    !--- Total stratiform updraft velocity:

    IF (nactivpdf == 0) THEN
       !--- Turbulent vertical velocity:
       pw(1:kproma,:,1) = MAX(w_min,w_large(1:kproma,:,krow)+w_turb(1:kproma,:,krow))
       w(1)%ptr(1:kproma,:,krow) = pw(1:kproma,:,1)
       ! Only one "bin", with probability of 1. The actual value doesn't
       ! matter so long as it's finite, since it cancels out of the CDNC
       ! calculation.
       pwpdf(1:kproma,:,1) = 1.0_dp
    ELSE
       CALL aero_activ_updraft_sigma(kproma,   kbdim,   klev,    krow, &
                                     ptkem1,  zwturb                   )

       CALL aero_activ_updraft_pdf(kproma,  kbdim,  klev,  krow, &
                                   zwlarge, zwturb, pw,    pwpdf )
    END IF

  END SUBROUTINE activ_updraft

  SUBROUTINE aero_activ_updraft_sigma(kproma,   kbdim,   klev,    krow, &
                                      ptkem1,   pwsigma                 )

    ! *aero_activ_updraft_sigma* calculates the standard deviation of the pdf of uf 
    !                            updraft vertical velocity
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford                 2013
    !
    ! References:
    ! -----------
    ! West et al., ACP, 2013. 
    !

    IMPLICIT NONE

    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow

    REAL(dp), INTENT(in)  :: ptkem1(kbdim,klev)  ! turbulent kinetic energy

    REAL(dp), INTENT(out) :: pwsigma(kbdim,klev) ! st. dev. of vertical velocity

    !--- Large scale vertical velocity in SI units:

    pwsigma(1:kproma,:)      = MAX(w_sigma_min, ((2.0_dp/3.0_dp)*ptkem1(1:kproma,:))**0.5_dp) ! m/s
    w_sigma(1:kproma,:,krow) = pwsigma(1:kproma,:)

  END SUBROUTINE aero_activ_updraft_sigma

  SUBROUTINE aero_activ_updraft_pdf(kproma,  kbdim,   klev, krow, &
                                    pwlarge, pwsigma, pw,   pwpdf )

    ! *aero_activ_updraft_* calculates Gaussian pdf of  
    !                       updraft vertical velocity
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford                 2013
    !
    ! References:
    ! -----------
    ! West et al., ACP, 2013. 
    !

    USE mo_math_constants, ONLY: pi
    USE mo_param_switches, ONLY: nactivpdf

    IMPLICIT NONE

    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow

    REAL(dp), INTENT(in)  :: pwlarge(kbdim,klev), & ! large-scale vertical velocity [m s-1]
                             pwsigma(kbdim,klev)    ! st. dev. of vertical velocity [m s-1]

    REAL(dp), INTENT(out) :: pw(kbdim,klev,nw),   & ! vertical velocity bins [m s-1]
                             pwpdf(kbdim,klev,nw)   ! vettical velocity PDF [s m-1]

    INTEGER               :: jl, jk, jw

    REAL(dp)              :: zw_width(kbdim,klev), &
                             zw_min(kbdim,klev), &
                             zw_max(kbdim,klev)

    zw_min(1:kproma,:)   = 0.0_dp
    zw_max(1:kproma,:)   = 4.0_dp*pwsigma(1:kproma,:)
    zw_width(1:kproma,:) = (zw_max(1:kproma,:) - zw_min(1:kproma,:)) / DBLE(nw)

    DO jw=1, nw
      pw(1:kproma,:,jw) = zw_min(1:kproma,:) + (DBLE(jw) - 0.5_dp) * zw_width(1:kproma,:)

      pwpdf(1:kproma,:,jw) = (1.0_dp / ((2.0_dp*pi)**0.5_dp))                         &
                             * (1.0_dp / pwsigma(1:kproma,:))                         &
                             * EXP( -((pw(1:kproma,:,jw) - pwlarge(1:kproma,:))**2_dp &
                                    / (2.0_dp*pwsigma(1:kproma,:)**2.0_dp)) )
    END DO

    IF (nactivpdf < 0) THEN
       DO jw=1, nw
          w(jw)%ptr(1:kproma,:,krow) = pw(1:kproma,:,jw)
          w_pdf(jw)%ptr(1:kproma,:,krow) = pwpdf(1:kproma,:,jw)
       END DO
    END IF

  END SUBROUTINE aero_activ_updraft_pdf

!-------------------------------------------

  SUBROUTINE activ_lin_leaitch(kproma,  kbdim,    klev,     krow, &
                               pw, pcdncact                       )

    ! *activ_lin_leaitch* calculates the number of activated aerosol 
    !                     particles from the aerosol number concentration
    !SF now independent of HAM, since HAM-specific calculation are computed in mo_ham_activ
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-MET                       2004
    !
    ! Method:
    ! -------
    ! The parameterisation follows the simple empirical relations of 
    ! Lin and Leaitch (1997).
    ! Updraft velocity is parameterized following Lohmann et al. (1999).
    !

    USE mo_kind,       ONLY: dp
    USE mo_conv,       ONLY: na_cv, cdncact_cv

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kproma, kbdim, klev, krow

    REAL(dp), INTENT(IN)  :: pw(kbdim,klev)  ! stratiform updraft velocity, large-scale+TKE (>0.0) [m s-1]
    REAL(dp), INTENT(out) :: pcdncact(kbdim,klev)  ! number of activated particles

    REAL(dp), PARAMETER :: c2=2.3E-10_dp, & ! [m4 s-1]
                           c3=1.27_dp       ! [1]

    INTEGER  :: jl, jk
    REAL(dp) :: zNmax, zeps

    zeps=EPSILON(1.0_dp)

    pcdncact(:,:) = 0._dp
    cdncact_cv(:,:,krow) = 0._dp

    !--- Aerosol activation:

    DO jk=1, klev
       DO jl=1, kproma

          !--- Stratiform clouds:

          ! Activation occurs only in occurrence of supersaturation

          !>>SF note: 
          !     The previous temperature restriction (temp > homogeneous freezing temp)
          !     has been removed because it was preventing to diagnose the number of
          !     dust and BC particules in soluble modes where temp < hom. freezing.
          !     The rationale behind is that diagnosing this allows further
          !     devel to implement concurrent homogeneous vs heterogenous freezing processes
          !     (which is not yet part of this version, though).
          !
          !IMPORTANT: 
          !     This temperature condition removal is completely transparent for the sanity 
          !     of the current code, since relevant temperature ranges are now safely checked
          !     directly in cloud_cdnc_icnc
          !<<SF

          IF(pw(jl,jk)>zeps .AND. na(jl,jk,krow)>zeps) THEN

             !--- Maximum number of activated particles [m-3]:

             zNmax=(na(jl,jk,krow)*pw(jl,jk))/(pw(jl,jk)+c2*na(jl,jk,krow))

             ! Average number of activated particles [m-3]:
             ! zNmax need to be converted to [cm-3] and the
             ! result to be converted back to [m-3].

             pcdncact(jl,jk)=0.1E6_dp*(1.0E-6_dp*zNmax)**c3

          END IF

          !--- Convective clouds:

          IF(pw(jl,jk)>zeps .AND. na_cv(jl,jk,krow)>zeps) THEN

             zNmax=(na_cv(jl,jk,krow)*pw(jl,jk))/(pw(jl,jk)+c2*na_cv(jl,jk,krow))
             cdncact_cv(jl,jk,krow)=0.1E6_dp*(1.0E-6_dp*zNmax)**c3

          ENDIF

       END DO
    END DO

  END SUBROUTINE activ_lin_leaitch

 SUBROUTINE activ_initialize

    USE mo_control,            ONLY: nlev, nn, &
                                     lcouple !SF
    USE mo_exception,          ONLY: message, em_param
    USE mo_submodel,           ONLY: print_value, lham, lhammoz, lccnclim
    USE mo_echam_cloud_params, ONLY: ccsaut, ccraut
    USE mo_param_switches,     ONLY : icover, nauto, &          !++mgs
                                      ncd_activ, nactivpdf, nic_cirrus, lcdnc_progn, &
                                      cdnc_min_fixed
    USE mo_tracer,             ONLY: get_tracer

    CHARACTER(len=24)      :: csubmname

    !--- Set number of updraft bins: 

    SELECT CASE(ABS(nactivpdf))
      CASE(0)
        nw = 1
      CASE(1)
        nw = 20
      CASE DEFAULT
        nw = ABS(nactivpdf)
    END SELECT

    IF (nactivpdf <= 0) THEN
      ! These are used either if not using a PDF, or if per-bin
      ! diagnostics are requested.
      ALLOCATE(w(nw))
      ALLOCATE(w_pdf(nw))
      ALLOCATE(swat_max_strat(nw))
      ALLOCATE(swat_max_conv(nw))
    END IF

    !
    !-- overwrite values for coupled CDNC/ICNC cloud scheme
    !
    IF (lcdnc_progn)  THEN
      IF (nlev == 31) THEN
         IF (nn == 63) THEN
            SELECT CASE (ncd_activ)
               CASE(1) ! LL activtion
                  !SF: updated on 2015.02.25 (David Neubauer / Katty Huang, pure atm run, HAM-M7, LL activation)
                  ccsaut = 1200._dp
                  ccraut = 3.5_dp
               CASE(2) !AR&G activation
                  SELECT CASE(cdnc_min_fixed)
                     CASE(10)
                        !SF: updated on 2017.02.14 (David Neubauer, pure atm run, HAM-M7)
                        ccsaut = 900._dp
                        ccraut = 2.8_dp
                     CASE(40)
                        !SF: updated on 2017.02.14 (David Neubauer, pure atm run, HAM-M7)
                        ccsaut = 900._dp
                        ccraut = 10.6_dp
                  END SELECT
            END SELECT
         ENDIF
      ENDIF

      IF (nlev == 47) THEN
         IF (nn == 63) THEN
            SELECT CASE (ncd_activ)
               CASE(1) ! LL activtion
                  !SF: updated on 2015.02.19 (David Neubauer, pure atm run, HAM-M7, LL activation)
                  ccsaut = 800._dp
                  ccraut = 5._dp
               CASE(2) ! AR&G activtion
                  SELECT CASE(cdnc_min_fixed)
                     CASE(10)
                        !SF: updated on 2017.02.14 (David Neubauer, pure atm run, HAM-M7)
                        ccsaut = 900._dp
                        ccraut = 2.8_dp
                     CASE(40)
                        !SF: updated on 2017.02.14 (David Neubauer, pure atm run, HAM-M7)
                        ccsaut = 900._dp
                        ccraut = 10.6_dp
                  END SELECT
            END SELECT
         ENDIF
      ENDIF
    ENDIF

!>>SF
    !-- Define the cdnc and icnc tracer index to point to the correct tracer:
    CALL get_tracer('CDNC',idx=idt_cdnc)
    CALL get_tracer('ICNC',idx=idt_icnc)
!<<SF

    !
    !-- Write out new parameters
    !
    IF (ncd_activ>0 .OR. nic_cirrus>0) THEN

      csubmname = 'UNKNOWN'
      IF (lham) csubmname = 'HAM'
      IF (lhammoz) csubmname = 'HAMMOZ'
      IF (lccnclim) csubmname = 'CCNCLIM'

      CALL message('','')
      CALL message('','----------------------------------------------------------')
      CALL message('activ_initialize','Parameter settings for the ECHAM-'//TRIM(csubmname)  &
                   //' cloud microphysics scheme')
      CALL message('','---')
      CALL print_value('              ncd_activ                       = ', ncd_activ)
      CALL print_value('              nic_cirrus                       = ', nic_cirrus)
      CALL message('', ' => Parameter adjustments in mo_activ:', level=em_param)
      CALL print_value('              ccsaut =', ccsaut)
      CALL print_value('              ccraut =', ccraut)
      CALL message('','---')
      CALL message('','----------------------------------------------------------')

    ENDIF
  END SUBROUTINE activ_initialize

  SUBROUTINE construct_activ_stream

    ! *construct_stream_activ* allocates output streams
    !                          for the activation schemes
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-MET                       2004
    !

  USE mo_memory_base,    ONLY: new_stream, add_stream_element, AUTO,  &
                               default_stream_setting, add_stream_reference
  USE mo_filename,       ONLY: trac_filetype
  USE mo_linked_list,    ONLY: HYBRID
  USE mo_param_switches, ONLY: ncd_activ, nactivpdf, nic_cirrus !SF

  IMPLICIT NONE

  INTEGER           :: jw
  CHARACTER(len=10) :: cbin


  !--- Create new stream:

  CALL new_stream (activ ,'activ',filetype=trac_filetype)


  !--- Add standard fields for post-processing:

  CALL add_stream_reference (activ, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
  CALL add_stream_reference (activ, 'lsp'     ,'sp'    ,lpost=.TRUE.)
  CALL add_stream_reference (activ, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
  CALL add_stream_reference (activ, 'gboxarea','geoloc',lpost=.TRUE.)

  CALL default_stream_setting (activ, lpost     = .TRUE. , &
                                      lrerun    = .TRUE. , &
                                      leveltype = HYBRID , &
                                      table     = 199,     &
                                      code      = AUTO     )
  !--- 1) Cloud Properties:

  CALL add_stream_element (activ,   'SWAT',       swat,                                   &
                           longname='ECHAM supersaturation over water',   units='% [0-1]' )

  IF (ncd_activ==2) THEN

     IF (nactivpdf == 0) THEN
        CALL add_stream_element (activ,   'SWAT_MAX_STRAT', swat_max_strat(1)%ptr, &
                                 longname='maximum supersaturation stratiform', units='% [0-1]' )

        CALL add_stream_element (activ,   'SWAT_MAX_CONV',  swat_max_conv(1)%ptr, &
                                 longname='maximum supersaturation convective', units='% [0-1]' )
     ELSE IF (nactivpdf < 0) THEN
        DO jw=1,nw
           WRITE (cbin, "(I2.2)") jw
           CALL add_stream_element (activ,   'SWAT_MAX_STRAT_'//TRIM(cbin), swat_max_strat(jw)%ptr, &
                                    longname='maximum supersaturation stratiform, vertical velocity bin '//TRIM(cbin), &
                                    units='% [0-1]' )

           CALL add_stream_element (activ,   'SWAT_MAX_CONV_'//TRIM(cbin), swat_max_conv(jw)%ptr, &
                                    longname='maximum supersaturation convective, vertical velocity bin '//TRIM(cbin), &
                                    units='% [0-1]' )
        END DO
     END IF
  ENDIF

  IF (nactivpdf == 0) THEN
     CALL add_stream_element (activ,   'W',          w(1)%ptr, &
                              longname='total vertical velocity for activation',units='m s-1')
  ELSE IF (nactivpdf < 0) THEN
     DO jw=1, nw
       WRITE (cbin, "(I2.2)") jw
       CALL add_stream_element (activ,   'W_'//TRIM(cbin), w(jw)%ptr, &
                                longname='Vertical velocity bin '//TRIM(cbin)//' for activation', &
                                units='m s-1')

       CALL add_stream_element (activ,   'W_PDF_'//TRIM(cbin), w_pdf(jw)%ptr, &
                                longname='Vertical velocity PDF in bin '//TRIM(cbin)//' for activation', &
                                units='s m-1')
     END DO
  END IF

  CALL add_stream_element (activ,   'W_LARGE',    w_large,                                &
                           longname='large scale vertical velocity',      units='m s-1'   )

  IF (nactivpdf == 0) THEN
     CALL add_stream_element (activ, 'W_TURB',     w_turb,                                 &
                              longname='turbulent vertical velocity',      units='m s-1'   )
  ELSE
     CALL add_stream_element (activ, 'W_SIGMA',    w_sigma,                                    &
                              longname='sub-grid st. dev. of vertical velocity', units='m s-1' )
  END IF

  CALL add_stream_element (activ,   'W_CAPE',     w_cape,                                 &
                           longname='convective updraft velocity from CAPE', units='m s-1')

  CALL add_stream_element (activ,   'REFFL',      reffl,                                  &
                           longname='cloud drop effectiv radius',         units='um'      )

  IF (nic_cirrus>0) THEN

  CALL add_stream_element (activ,   'REFFI',      reffi,                                  &
                           longname='ice crystal effectiv radius',        units='um'      )
  END IF

  CALL add_stream_element (activ,   'NA',         na,                                     &
                           longname='aerosol number for activation',      units='m-3'     )

  CALL default_stream_setting (activ, laccu=.TRUE.)

  CALL add_stream_element (activ,   'QNUC',       qnuc,                                   &
                           longname='CD nucleation rate',                 units='m-3 s-1' )

  CALL add_stream_element (activ,   'QAUT',       qaut,                                   &
                           longname='CD autoconversion rate',             units='m-3 s-1' )

  CALL add_stream_element (activ,   'QACC',       qacc,                                   &
                           longname='CD accretion rate',                  units='m-3 s-1' )

  CALL add_stream_element (activ,   'QFRE',       qfre,                                   &
                           longname='CD freezing rate',                   units='m-3 s-1' )
  !>>dod deleted QEVA, not used anywhere
  !  CALL add_stream_element (activ,   'QEVA',       qeva,                                   &
  !                           longname='CD evaporation rate',                units='m-3 s-1' )

  CALL add_stream_element (activ,   'QMEL',       qmel,                                   &
                           longname='CD source rate from melting ice',    units='m-3 s-1' )

  CALL add_stream_element (activ,   'CDNC_ACC',   cdnc_acc,                               &
                           longname='CDNC occurence acc.+ cloud weighted',units='m-3'     )

  CALL add_stream_element (activ,   'CDNC',       cdnc,                                   &
                           longname='CDNC',units='m-3'                                    )

  CALL add_stream_element (activ,   'CDNC_BURDEN_ACC',cdnc_burden_acc,                    &
                           longname='CDNC burden occurence accumulated',  units='m-2'     )

  CALL add_stream_element (activ,   'CDNC_BURDEN',cdnc_burden,                            &
                           longname='CDNC burden',                        units='m-2'     )

  CALL add_stream_element (activ,   'BURDEN_TIME',burden_time,                            &
                           longname='acc. cdnc burden occ.time fraction', units='1'       )

  CALL add_stream_element (activ,   'LWC_ACC',    lwc_acc,                                &
                           longname='liq wat cont acc.+ cloud weighted',  units='kg m-3'  )

  CALL add_stream_element (activ,   'CLOUD_TIME', cloud_time,                             &
                           longname='acc. cloud occurence time fraction', units='1'       )

  CALL add_stream_element (activ,   'REFFL_ACC',  reffl_acc,                              &
                           longname='cloud drop effectiv radius weighted',units='um'      )

  CALL add_stream_element (activ,   'REFFL_CT',  reffl_ct,                                &
                           longname='cloud top effectiv radius weighted',units='um'       )

  CALL add_stream_element (activ,   'REFFL_TIME',  reffl_time,                            &
                           longname='cloud top effectiv radius occ.time',units='1'        )

  CALL add_stream_element (activ,   'CDNC_CT',  cdnc_ct,                                  &
                           longname='cloud top cloud droplet number conc.',units='cm-3'   )

  CALL add_stream_element (activ,   'IWC_ACC',    iwc_acc,                                &
                           longname='ice wat cont acc.+ cloud weighted',  units='kg m-3'  )

  CALL add_stream_element (activ,   'CLIWC_TIME', cliwc_time,                             &
                           longname='acc. cloud occurence time fraction', units='1'       )

  CALL default_stream_setting (activ, laccu=.FALSE., lpost=.FALSE.)

  CALL add_stream_element (activ,   'CLOUD_COVER_DUPLIC', cloud_cover_duplic,             &
                           longname='cloud cover duplicate for record at t+1', units='1'  )


  IF (nic_cirrus>0) THEN

  CALL add_stream_element (activ, 'ICNC_instantaneous', icnc_instantan, &
                           longname='ICNC instantaneous', units='m-3',  &
                           laccu=.FALSE., lpost=.TRUE., lrerun=.TRUE.)

  CALL default_stream_setting (activ, laccu=.TRUE., lpost=.TRUE.)

  CALL add_stream_element (activ,   'ICNC_ACC',   icnc_acc,                               &
                           longname='ICNC occurence acc.+ cloud weighted',units='m-3'     )

  CALL add_stream_element (activ,   'ICNC',       icnc,                                   &
                           longname='ICNC',units='m-3'                                    )

  CALL add_stream_element (activ,   'ICNC_BURDEN_ACC',icnc_burden_acc,                    &
                           longname='ICNC burden occurence accumulated',  units='m-2'     )

  CALL add_stream_element (activ,   'ICNC_BURDEN',icnc_burden,                            &
                           longname='ICNC burden',                        units='m-2'     )

  CALL add_stream_element (activ,   'BURDIC_TIME',burdic_time,                            &
                           longname='acc. icnc burden occ.time fraction', units='1'       )

  CALL add_stream_element (activ,   'REFFI_ACC',  reffi_acc,                              &
                           longname='ice crystal effectiv radius weighted',units='um'     )

  CALL add_stream_element (activ,   'REFFI_TOVS',  reffi_tovs,                            &
                           longname='semi-transparent cirrus effectiv radius',units='um'  )

  CALL add_stream_element (activ,   'REFFI_TIME',  reffi_time,                            &
                           longname='accumulted semi-transp. cirrus time',units='1'       )

  CALL add_stream_element (activ,   'IWP_TOVS',  iwp_tovs,                                &
                           longname='IWP sampled a la TOVS',units='kg m-2'                )

  CALL default_stream_setting (activ, laccu=.FALSE.)

  CALL add_stream_element (activ,   'SICE',       sice,                                   &
                           longname='ECHAM supersaturation over ice',     units='% [0-1]' )

  END IF

  CALL default_stream_setting (activ, laccu=.FALSE.)

END SUBROUTINE construct_activ_stream

END MODULE mo_activ
