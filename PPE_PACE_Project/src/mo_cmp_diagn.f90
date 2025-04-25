MODULE mo_cmp_diagn
! Based on mo_active
! Author: UP
! Date: 07 2020

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream

  IMPLICIT NONE

  PUBLIC construct_cmp_diagn_stream
  PUBLIC diags_cc_bytemp

  PRIVATE

  TYPE (t_stream), PUBLIC, POINTER :: cmpdiagn

  REAL(dp),        PUBLIC, POINTER :: tend_sacl(:,:,:) ! accretion of snow flakes with cloud droplets [kg/kg], tendency
  REAL(dp),        PUBLIC, POINTER :: tend_inucl(:,:,:) ! ice nucleation within section 1.3 of cloud_micro_interface
  REAL(dp),        PUBLIC, POINTER :: icnc_bfrs238(:,:,:) ! ICNC before freezing at T < 238 K routine
  REAL(dp),        PUBLIC, POINTER :: icnc_pfrs238(:,:,:) ! ICNC after/past freezing at T < 238 K routine
  REAL(dp),        PUBLIC, POINTER :: icnc_dfrs238(:,:,:) ! ICNC difference from freezing at T < 238 K routine
  REAL(dp),        PUBLIC, POINTER :: icnc_ctrfrs238(:,:,:) ! counter for going inside if statemenet of the freezing at T < 238 K routine
  REAL(dp),        PUBLIC, POINTER :: icnc_bfrl238(:,:,:) ! ICNC before freezing at T > 238 K routine
  REAL(dp),        PUBLIC, POINTER :: icnc_pfrl238(:,:,:) ! ICNC after/past freezing at T > 238 K routine
  REAL(dp),        PUBLIC, POINTER :: icnc_dfrl238(:,:,:) ! ICNC difference from freezing at T > 238 K routine
  REAL(dp),        PUBLIC, POINTER :: icnc_ctrfrl238(:,:,:) ! counter for going inside if statemenet of the freezing at T > 238 K routine
  REAL(dp),        PUBLIC, POINTER :: tend_testdummy_1(:,:,:) ! testing stuff
  REAL(dp),        PUBLIC, POINTER :: tend_testdummy_2(:,:,:) ! testing stuff

  !>>UP #821
  REAL(dp),        PUBLIC, POINTER :: diag_delta_sci(:,:,:) ! diagnose the Delta
                                      ! x inflicted by self-collection of ice,
                                      ! which can be phased with eta_sci
  REAL(dp),        PUBLIC, POINTER :: diag_cntr_sci(:,:,:) ! count how often
                                      ! something gets added as delta sci
  REAL(dp),        PUBLIC, POINTER :: diag_delta_icnucl(:,:,:) ! diagnose the Delta
                                      ! x inflicted by nucleation of ice (cirrus),
                                      ! which can be phased with eta_icnucl
  REAL(dp),        PUBLIC, POINTER :: diag_cntr_icnucl(:,:,:) ! count how often
                                      ! something gets added as delta icnucl
  REAL(dp),        PUBLIC, POINTER :: diag_delta_rime(:,:,:) ! diagnose the Delta
                                      ! x inflicted by riming,
                                      ! which can be phased with eta_riming
  REAL(dp),        PUBLIC, POINTER :: diag_cntr_rime(:,:,:) ! count how often
                                      ! something gets added as delta rime
  REAL(dp),        PUBLIC, POINTER :: diag_delta_subfs(:,:,:) ! diagnose the Delta
                                      ! x inflicted by sublimation of falling
                                      ! snow, which can be phased with eta_subfis
  REAL(dp),        PUBLIC, POINTER :: diag_cntr_subfs(:,:,:) ! count how often
                                      ! something gets added as delta subfs
  REAL(dp),        PUBLIC, POINTER :: diag_delta_subfi(:,:,:) ! diagnose the Delta
                                      ! x inflicted by sublimation of falling
                                      ! ice, which can be phased with eta_subfis
  REAL(dp),        PUBLIC, POINTER :: diag_cntr_subfi(:,:,:) ! count how often
                                      ! something gets added as delta subfi
  REAL(dp),        PUBLIC, POINTER :: diag_delta_icaccr(:,:,:) ! diagnose the Delta
                                      ! x inflicted by IC acccretion,
                                      ! which can be phased with eta_icaccr
  REAL(dp),        PUBLIC, POINTER :: diag_cntr_icaccr(:,:,:) ! count how often
                                      ! something gets added as delta icaccr
  REAL(dp),        PUBLIC, POINTER :: diag_delta_mlt(:,:,:) ! diagnose the Delta
                                      ! x inflicted by melting,
                                      ! which can be phased with eta_mlt
  REAL(dp),        PUBLIC, POINTER :: diag_cntr_mlt(:,:,:) ! count how often
                                      ! something gets added as delta mlt
  !<<UP #821
  !<<UP #821

  REAL(dp), PARAMETER :: thres_cld = 0.001 !UP #783

  !--- Subroutines:

CONTAINS

  SUBROUTINE construct_cmp_diagn_stream

    ! *construct_stream_cmp_diagn* allocates output streams
    !                          for the cmp scheme for tendency/pathway analysis
    !
    ! Author:
    ! -------
    ! Ulrike Proske, ETHZ                       2020
    !

  USE mo_memory_base,    ONLY: new_stream, add_stream_element, AUTO,  &
                               default_stream_setting, add_stream_reference
  USE mo_filename,       ONLY: out_filetype
  USE mo_linked_list,    ONLY: HYBRID

  IMPLICIT NONE

  !--- Create new stream:

  CALL new_stream (cmpdiagn ,'cmpdiagn',filetype=out_filetype) !out_filetype is different than in mo_activ.f90
                                         ! because I guess that this could be the correct type


  !--- Add standard fields for post-processing:

  CALL add_stream_reference (cmpdiagn, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
  CALL add_stream_reference (cmpdiagn, 'lsp'     ,'sp'    ,lpost=.TRUE.)
  CALL add_stream_reference (cmpdiagn, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
  CALL add_stream_reference (cmpdiagn, 'gboxarea','geoloc',lpost=.TRUE.)

  CALL default_stream_setting (cmpdiagn, lpost     = .TRUE. , &
                                      lrerun    = .TRUE. , &
                                      leveltype = HYBRID , &
                                      table     = 199,     &
                                      code      = AUTO     )
  !--- 1) Nucleation:

  CALL default_stream_setting (cmpdiagn, laccu=.TRUE.)

  CALL add_stream_element (cmpdiagn,   'TEND_SACL',    tend_sacl, &
                           longname='accretion of snow flakes with cloud droplets',      units='kg kg-1'   )
  CALL add_stream_element (cmpdiagn,   'TEND_INUCL',    tend_inucl, &
                           longname='number concentration of newly nucleated IC',      units='m-3'   )
  CALL add_stream_element (cmpdiagn,   'ICNC_BFRS238',  icnc_bfrs238, &
                           longname='icnc before the freezing colder 238K routine', units='m-3')
  CALL add_stream_element (cmpdiagn,   'ICNC_PFRS238',  icnc_pfrs238, &
                           longname='icnc past the freezing colder 238K routine', units='m-3')
  CALL add_stream_element (cmpdiagn,   'ICNC_DFRS238',  icnc_dfrs238, &
                           longname='icnc difference through the freezing colder 238K routine', units='m-3')
  CALL add_stream_element (cmpdiagn,   'ICNC_CTRFRS238',  icnc_ctrfrs238, &
                           longname='counter for the freezing colder 238K routine', units='m-3')
  CALL add_stream_element (cmpdiagn,   'ICNC_BFRL238',  icnc_bfrl238, &
                           longname='icnc before the freezing warmer 238K routine', units='m-3')
  CALL add_stream_element (cmpdiagn,   'ICNC_PFRL238',  icnc_pfrl238, &
                           longname='icnc past the freezing warmer 238K routine', units='m-3')
  CALL add_stream_element (cmpdiagn,   'ICNC_DFRL238',  icnc_dfrl238, &
                           longname='icnc difference through the freezing warmer 238K routine', units='m-3')
  CALL add_stream_element (cmpdiagn,   'ICNC_CTRFRL238',  icnc_ctrfrl238, &
                           longname='counter for the freezing warmer 238K routine', units='m-3')
  CALL add_stream_element (cmpdiagn, 'TEST-DUMMY-1', tend_testdummy_1, &
                           longname='test dummy 1', units = '')
  CALL add_stream_element (cmpdiagn, 'TEST-DUMMY-2', tend_testdummy_2, &
                           longname='test dummy 2', units = '')

  !>>UP #821
  CALL add_stream_element (cmpdiagn, 'DIAG_DELTA_SCI', diag_delta_sci, &
                           longname='Delta inflicted by self-collection of ice', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_CNTR_SCI', diag_cntr_sci, &
                           longname='How often the delta inflicted by sci gets added', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_DELTA_ICNUCL', diag_delta_icnucl, &
                           longname='Delta inflicted by nucleation of ice (cirrus)', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_CNTR_ICNUCL', diag_cntr_icnucl, &
                           longname='How often the delta inflicted by icnucl gets added', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_DELTA_RIME', diag_delta_rime, &
                           longname='Delta inflicted by riming', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_CNTR_RIME', diag_cntr_rime, &
                           longname='How often the delta inflicted by riming gets added', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_DELTA_SUBFS', diag_delta_subfs, &
                           longname='Delta inflicted by sublimation of falling snow', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_CNTR_SUBFS', diag_cntr_subfs, &
                           longname='How often the delta inflicted by subfs gets added', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_DELTA_SUBFI', diag_delta_subfi, &
                           longname='Delta inflicted by sublimation of falling ice', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_CNTR_SUBFI', diag_cntr_subfi, &
                           longname='How often the delta inflicted by subfi gets added', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_DELTA_ICACCR', diag_delta_icaccr, &
                           longname='Delta inflicted by icaccr', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_CNTR_ICACCR', diag_cntr_icaccr, &
                           longname='How often the delta inflicted by icaccr gets added', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_DELTA_MLT', diag_delta_mlt, &
                           longname='Delta inflicted by mlt', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  CALL add_stream_element (cmpdiagn, 'DIAG_CNTR_MLT', diag_cntr_mlt, &
                           longname='How often the delta inflicted by melting gets added', &
                           units = '', lrerun=.TRUE., laccu=.FALSE.)
  !<<UP #821
  !CALL default_stream_setting (activ, laccu=.FALSE., lpost=.FALSE.)

  !CALL default_stream_setting (activ, laccu=.TRUE., lpost=.TRUE.)

  !CALL default_stream_setting (activ, laccu=.FALSE.)

  !CALL default_stream_setting (activ, laccu=.FALSE.)

END SUBROUTINE construct_cmp_diagn_stream

SUBROUTINE diags_cc_bytemp(kproma, kbdim, klev, krow, paclc, ptm1, ptte)

    !UP: #783
    !Note: for the cloud top properties, the max random overlap assumption is
    !chosen

    USE mo_time_control,   ONLY: time_step_len
    USE mo_physical_constants,      ONLY: tmelt
    USE mo_echam_cloud_params,          ONLY: cthomi

    USE mo_activ,           ONLY: lcc_mod, mcc_mod, icc_mod
    !-- arguments
    INTEGER, INTENT(IN) :: kproma, kbdim, klev, krow

    REAL(dp), INTENT(IN) :: paclc(kbdim,klev) !cloud fraction
    REAL(dp), INTENT(IN) :: ptm1(kbdim,klev)  !temperature (t-1)
    REAL(dp), INTENT(IN) :: ptte(kbdim,klev)  !temperature tendency

    !-- local vars
    INTEGER :: jk

    REAL(dp) :: ztemp(kbdim,klev)  ! temperature (t)
    REAL(dp) :: zacltot_liq(kbdim) ! liq cloud cover seen from above down to the current level
    REAL(dp) :: zacltot_mix(kbdim) ! mixed-p cloud cover seen from above down to the current level
    REAL(dp) :: zacltot_ice(kbdim) ! ice cloud cover seen from above down to the current level
    REAL(dp) :: zdeltacc(kbdim)    ! utility var 

    REAL(dp) :: ztmp1_1d(kbdim), ztmp2_1d(kbdim), ztmp3_1d(kbdim)

    LOGICAL :: ll_liq(kbdim,klev)
    LOGICAL :: ll_mix(kbdim,klev)
    LOGICAL :: ll_ice(kbdim,klev)

    !-- Cloud top properties (plus some related 3D diags):

    !Note: the following 2 cod variables can't be directly copied by means of
    !add_stream_reference 
    !      from the COSP stream in the construct_stream_aerocom_ic subroutine
    !      because the COSP stream
    !      is built AFTER the aerocom_ic stream

    ztemp(1:kproma,:) = ptm1(1:kproma,:) + time_step_len*ptte(1:kproma,:)

    ll_liq(1:kproma,:) = (ztemp(1:kproma,:) >= tmelt)
    ll_mix(1:kproma,:) = (ztemp(1:kproma,:) < tmelt) .AND. (ztemp(1:kproma,:) >= cthomi)
    ll_ice(1:kproma,:) = .NOT. ll_liq(1:kproma,:) .AND. .NOT. ll_mix(1:kproma,:)

    lcc_mod(1:kproma,krow) = 1._dp !1 instead of zero because max random overlap calc. needs it this way
    mcc_mod(1:kproma,krow) = 1._dp !1 instead of zero because max random overlap calc. needs it this way
    icc_mod(1:kproma,krow) = 1._dp !1 instead of zero because max random overlap calc. needs it this way

    DO jk = 2,klev
       ztmp1_1d(1:kproma) = MAX(paclc(1:kproma,jk), paclc(1:kproma,jk-1))
       ztmp2_1d(1:kproma) = MIN(paclc(1:kproma,jk-1), 1._dp - thres_cld)
       !UP comment: I think above is because you don't want to divide by 0 down
       !below

       ztmp3_1d(1:kproma) = lcc_mod(1:kproma,krow) * (1._dp - ztmp1_1d(1:kproma)) &
                                                   / (1._dp - ztmp2_1d(1:kproma))

       zacltot_liq(1:kproma) = MERGE(ztmp3_1d(1:kproma), zacltot_liq(1:kproma), ll_liq(1:kproma,jk))

       ztmp3_1d(1:kproma) = mcc_mod(1:kproma,krow) * (1._dp - ztmp1_1d(1:kproma)) &
                                                   / (1._dp - ztmp2_1d(1:kproma))

       zacltot_mix(1:kproma) = MERGE(ztmp3_1d(1:kproma), zacltot_mix(1:kproma), ll_mix(1:kproma,jk))

       ztmp3_1d(1:kproma) = icc_mod(1:kproma,krow) * (1._dp - ztmp1_1d(1:kproma)) &
                                                   / (1._dp - ztmp2_1d(1:kproma))

       zacltot_ice(1:kproma) = MERGE(ztmp3_1d(1:kproma), zacltot_ice(1:kproma), ll_ice(1:kproma,jk))

       !-- final:
       lcc_mod(1:kproma,krow) = MERGE(zacltot_liq(1:kproma),lcc_mod(1:kproma,krow), ll_liq(1:kproma,jk))
       mcc_mod(1:kproma,krow) = MERGE(zacltot_mix(1:kproma),mcc_mod(1:kproma,krow), ll_mix(1:kproma,jk))
       icc_mod(1:kproma,krow) = MERGE(zacltot_ice(1:kproma),icc_mod(1:kproma,krow), ll_ice(1:kproma,jk))
       !UP comment: I think it's the random overlap that requires you to
       !calculate the cloud free fraction before

    ENDDO

    lcc_mod(1:kproma,krow) = 1._dp - lcc_mod(1:kproma,krow) !now this is really a cloud cover
    mcc_mod(1:kproma,krow) = 1._dp - mcc_mod(1:kproma,krow) !now this is really a cloud cover
    icc_mod(1:kproma,krow) = 1._dp - icc_mod(1:kproma,krow) !now this is really a cloud cover

  END SUBROUTINE diags_cc_bytemp

END MODULE mo_cmp_diagn
