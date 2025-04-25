MODULE mo_remos

  USE mo_kind,                    ONLY: dp
  USE mo_memory_base,             ONLY: new_stream, add_stream_element, &
                                        default_stream_setting, AUTO, t_stream, &
                                        add_stream_reference
  USE mo_time_event,              ONLY: io_time_event
  USE mo_time_control,            ONLY: p_bcast_event
  USE mo_linked_list,             ONLY: HYBRID
  USE mo_netcdf,                  ONLY: HYBRID_H

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: init_remos

! Remos essential diagnostics
  REAL(dp), PUBLIC, POINTER :: dslf(:,:,:), dfr(:,:,:), dfr_het(:,:,:), dfr_liq(:,:,:), dfr_oliq(:,:,:), &
                               dfr_sosi(:,:,:), dfr_socm(:,:,:), dfr_nihet(:,:,:), dfr_nihom(:,:,:), &
                               dfr_ninuc(:,:,:), dfr_nidet(:,:,:)

! Remos 2 moment properties
  REAL(dp), PUBLIC, POINTER :: dqsnow(:,:,:), dqsflx(:,:,:), dsacl(:,:,:), dsacln(:,:,:), dsaci(:,:,:), &
                               dsacin(:,:,:), dsaut(:,:,:), dsautn(:,:,:), driv_2m(:,:,:), &
                               drieff_2m(:,:,:), dqssub(:,:,:), dqsmlt(:,:,:), dvtim_2m(:,:,:), dvtin_2m(:,:,:), &
                               dqismlt(:,:,:), dnismlt(:,:,:), dqissub(:,:,:), dnissub(:,:,:), dqisflx(:,:,:), &
                               dnisflx(:,:,:), dqifal(:,:,:), dnifal(:,:,:), dvtimfal(:,:,:), dvtinfal(:,:,:)
  
! Remos process rates
  REAL(dp), PUBLIC, POINTER :: dqiflx(:,:,:), dniflx(:,:,:), dqrimflx(:,:,:), dbrimflx(:,:,:), &
                       dfrl_het(:,:,:), dfrl_hom(:,:,:), dfrln_hom(:,:,:), dfrln_het(:,:,:), &
                       ddep_wbf(:,:,:), dqimlt(:,:,:), dnimlt(:,:,:), drpr(:,:,:), dqrflx(:,:,:), &
                       dqrevp(:,:,:), dnislf(:,:,:), dnifrz_ci(:,:,:), ddep_ci1(:,:,:), ddep_ci2(:,:,:), &
                       dxite_ls(:,:,:), dxite_cv(:,:,:), ddep_co(:,:,:), ddep_aj(:,:,:), dni_aj(:,:,:), &
                       dxlte_ls(:,:,:), dxlte_cv(:,:,:), dcnd_co(:,:,:), dcnd_aj(:,:,:), dnc_aj(:,:,:), &
                       dxisub_cv(:,:,:), dxlevp_cv(:,:,:), dxisub_tp(:,:,:), dxlevp_tp(:,:,:), &
                       dqccol(:,:,:), dnccol(:,:,:), dncnuc(:,:,:), dicnc_cv(:,:,:), &
                       dcdnc_cv(:,:,:), dxisub_cc(:,:,:), dxlevp_cc(:,:,:), dximlt_cv(:,:,:), &
                       dcdnc_wbf(:,:,:), dqisub(:,:,:), dni_lkp(:,:,:), dqimlt_rain(:,:,:), drprn(:,:,:), &
                       dnisub(:,:,:), dnisub_cc(:,:,:), dnisub_tp(:,:,:), dncevp_cc(:,:,:), dncevp_tp(:,:,:), &
                       dnevp_co(:,:,:), dnsub_co(:,:,:), ddeltaqf(:,:,:), ddeltaqs(:,:,:), dqcdif(:,:,:), &
                       dqisten(:,:,:), dnisten(:,:,:), dbgsten(:,:,:), dqristen(:,:,:), dnimlt_rain(:,:,:), &
                       daut(:,:,:), dacc(:,:,:), ddep_a(:,:,:), dqhsten(:,:,:), dcdncte_ls(:,:,:), dicncte_ls(:,:,:)

! Remos ice properties
  REAL(dp), PUBLIC, POINTER :: dvtim(:,:,:), dvtin(:,:,:), drim(:,:,:), drcm(:,:,:), &
                       driv(:,:,:), drieff(:,:,:), drieff_optics(:,:,:), &
                       dqirim(:,:,:), dbirim(:,:,:), dqihet(:,:,:), dqiliq(:,:,:), drhop(:,:,:), &
                       dnihet(:,:,:), dnihom(:,:,:), dninuc(:,:,:), dnidet(:,:,:), &
                       drhoice(:,:,:), dqioliq(:,:,:), &
                       dsedfracm(:,:,:), dsedfracn(:,:,:), daclci(:,:,:), dicnc(:,:,:), &
                       dcdnc(:,:,:), dxlb(:,:,:), dxib(:,:,:), dicncb(:,:,:), dcdncb(:,:,:), &
                       dxi_cloud(:,:,:), dxi_snow(:,:,:), dxivi_cloud(:,:), dxivi_snow(:,:), &
                       dice_mu(:,:,:), dice_lam(:,:,:)

! Remos liquid properties
  REAL(dp), PUBLIC, POINTER :: dqr(:,:,:), dcdncact(:,:,:)

! Remos budgets
  REAL(dp), PUBLIC, POINTER :: dmete(:,:), dncdiff(:,:), dnidiff(:,:), dnisurf(:,:), &
                               dncdiff2(:,:), dnidiff2(:,:), dnictrl(:,:), dncctrl(:,:)

! Remos scavenging values
  REAL(dp), PUBLIC, POINTER :: dfrain(:,:,:), dfevapr(:,:,:), dmratepr(:,:,:), &
                               dfsnow(:,:,:), dfsubls(:,:,:), dmrateps(:,:,:), &
                               dmlwc(:,:,:), dmiwc(:,:,:), dmsnowacl(:,:,:)

! Remos dynamics
  REAL(dp), PUBLIC, POINTER :: dupdraft(:,:,:), dupdraftmax(:,:,:), dtke(:,:,:), drhoair(:,:,:), &
                               dcfl(:,:,:), ddpg(:,:,:), dsedtm(:,:,:), dsedtn(:,:,:), &
                               dnsedi(:,:), dnmicro(:,:), ddz(:,:,:)

! Remos cloud types
  REAL(dp), PUBLIC, POINTER :: ddp_cld(:,:,:), dt_top(:,:,:), diclnb(:,:), dclcol_time(:,:)

  CONTAINS

  SUBROUTINE init_remos

    USE mo_mpi,                   ONLY: p_parallel_io, p_barrier, p_bcast, p_io
    USE mo_namelist,              ONLY: open_nml, position_nml, POSITIONED
    USE mo_exception,             ONLY: finish
    USE mo_p3_fields,             ONLY: lfull_diag

    CHARACTER(LEN=32) :: ichar
    INTEGER           :: i, j, ierr, inml, iunit, ilength, nlength
    TYPE (t_stream), POINTER      :: remos

    CALL new_stream (remos,'remos',lrerun=.false.)

    !--- Add standard fields for post-processing:

    CALL add_stream_reference (remos, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (remos, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (remos, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (remos, 'gboxarea','geoloc',lpost=.TRUE.)

    CALL default_stream_setting (remos, lpost     = .TRUE. , &
                                        lrerun    = .FALSE. , &
                                        leveltype = HYBRID , &
                                        table     = 199,     &
                                        code      = AUTO     )

    CALL default_stream_setting (remos, laccu=.TRUE.)

    ! Remos essential diagnostics
    CALL add_stream_element (remos, 'dslf', dslf,                   &
                             longname='Supercooled liquid fraction', units='-')
    CALL add_stream_element (remos, 'dfr', dfr,                   &
                             longname='riming fraction', units='(0-1)')
    CALL add_stream_element (remos, 'dfr_het', dfr_het,                   &
                             longname='heterogeneosly formed ice fraction', units='(0-1)')
    CALL add_stream_element (remos, 'dfr_liq', dfr_liq,                   &
                             longname='fraction of ice formed through liquid', units='(0-1)')
    CALL add_stream_element (remos, 'dfr_oliq', dfr_oliq,                   &
                             longname='fraction of ice initiated through liquid', units='(0-1)')
    CALL add_stream_element (remos, 'dfr_sosi', dfr_sosi,                   &
                             longname='fraction sink over source terms', units='(0-1)')
    CALL add_stream_element (remos, 'dfr_socm', dfr_socm,                   &
                             longname='fraction source over cloud mass', units='(0-10)')

    CALL add_stream_element (remos, 'dfr_nihet', dfr_nihet,                   &
                             longname='fraction of het. formed ice number', units='(0-1)')
    CALL add_stream_element (remos, 'dfr_nihom', dfr_nihom,                   &
                             longname='fraction of hom. formed ice number', units='(0-1)')
    CALL add_stream_element (remos, 'dfr_ninuc', dfr_ninuc,                   &
                             longname='fraction of nucleated ice number', units='(0-1)')
    CALL add_stream_element (remos, 'dfr_nidet', dfr_nidet,                   &
                             longname='fraction detrained ice number', units='(0-1)')

    ! anything below this line is ONLY OUTPUT IF lfull_diag == .TRUE.
    CALL default_stream_setting (remos, lpost=lfull_diag)

    ! Remos 2 moment properties 
    CALL add_stream_element (remos, 'dqsnow', dqsnow,                   &
                             longname='snow mass', units='kg m-2 s-1')
    CALL add_stream_element (remos, 'dqsflx', dqsflx,                   &
                             longname='snow flux', units='kg m-2 s-1')

    CALL add_stream_element (remos, 'dqisflx', dqisflx,                   &
                             longname='falling ice mass flux', units='kg m-2 s-1')
    CALL add_stream_element (remos, 'dnisflx', dnisflx,                   &
                             longname='falling ice number flux', units='m-2 s-1')
    CALL add_stream_element (remos, 'dnifal', dnifal,                   &
                             longname='ice number sedimentation tendency', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dqifal', dqifal,                   &
                             longname='ice mass sedimentation tendency', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dvtimfal', dvtimfal,                   &
                             longname='fall-speed of falling ice mass', units='m s-1')
    CALL add_stream_element (remos, 'dvtinfal', dvtinfal,                   &
                             longname='fall-speed of falling ice number', units='m s-1')

    CALL add_stream_element (remos, 'dsacl', dsacl,                   &
                             longname='snow accretion liquid', units='-')
    CALL add_stream_element (remos, 'dsacln', dsacln,                   &
                             longname='snow accretion liquid number', units='-')
    CALL add_stream_element (remos, 'dsaci', dsaci,                   &
                             longname='snow accretion ice', units='-')
    CALL add_stream_element (remos, 'dsacin', dsacin,                   &
                             longname='snow accretion ice number', units='-')
    CALL add_stream_element (remos, 'dsaut', dsaut,                   &
                             longname='snow autoconversion', units='-')
    CALL add_stream_element (remos, 'dsautn', dsautn,                   &
                             longname='snow autoconversion number', units='-')
    CALL add_stream_element (remos, 'dqssub', dqssub,                   &
                             longname='snow sublimation', units='-')
    CALL add_stream_element (remos, 'dqsmlt', dqsmlt,                   &
                             longname='snow melting', units='-')

    CALL add_stream_element (remos, 'dqismlt', dqismlt,                   &
                             longname='falling ice melting', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dnismlt', dnismlt,                   &
                             longname='falling ice melting number', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dqissub', dqissub,                   &
                             longname='falling ice sublimation', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dnissub', dnissub,                   &
                             longname='falling ice sublimation number', units='kg-1 s-1')

    CALL add_stream_element (remos, 'driv_2m', driv_2m,                   &
                             longname='volume mean ice crystal radius', units='-')
    CALL add_stream_element (remos, 'drieff_2m', drieff_2m,                   &
                             longname='effective ice crystal radius', units='m')
    CALL add_stream_element (remos, 'dvtim_2m', dvtim_2m,                   &
                             longname='mass-weighted crystal fallspeed', units='m s-1')
    CALL add_stream_element (remos, 'dvtin_2m', dvtin_2m,                   &
                             longname='number-weighted crystal fallspeed', units='m s-1')

  ! Remos liquid properties
    CALL add_stream_element( remos, 'dqr', dqr,                           &
                             longname='rain water mixing ratio', units='kg kg-1')
    CALL add_stream_element( remos, 'dcdncact', dcdncact,                           &
                             longname='number of activated aerosols', units='kg-1')

  ! Remos debug streams
    CALL add_stream_element (remos, 'dqiflx', dqiflx,                   &
                             longname='ice mass flux', units='units')
    CALL add_stream_element (remos, 'dniflx', dniflx,                   &
                             longname='ice number flux', units='units')
    CALL add_stream_element (remos, 'dqrimflx', dqrimflx,                   &
                             longname='rimed mass flux', units='units')
    CALL add_stream_element (remos, 'dbrimflx', dbrimflx,                   &
                             longname='rimed volume flux', units='units')
    CALL add_stream_element (remos, 'dnifrz_ci', dnifrz_ci,                   &
                             longname='number from cirrus scheme', units='kg-1 s-1')
    CALL add_stream_element (remos, 'ddep_ci1', ddep_ci1,                   &
                             longname='mass dep from cirrus scheme', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'ddep_ci2', ddep_ci2,                   &
                             longname='explicit mass dep in cirrus regime', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dfrl_het', dfrl_het,                   &
                             longname='heterogeneous freezing mass', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dfrl_hom', dfrl_hom,                   &
                             longname='homogeneous freezing mass', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dfrln_hom', dfrln_hom,                   &
                             longname='homogeneous freezing number', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dfrln_het', dfrln_het,                   &
                             longname='heterogeneous freezing number', units='kg-1 s-1')
    CALL add_stream_element (remos, 'ddep_wbf', ddep_wbf,                   &
                             longname='WBF process deposition', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dcdnc_wbf', dcdnc_wbf,                   &
                             longname='WBF process number', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dqimlt', dqimlt,                   &
                             longname='mass melting', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dqimlt_rain', dqimlt_rain,                   &
                             longname='mass melting to rain', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dnimlt', dnimlt,                   &
                             longname='number melting', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dnimlt_rain', dnimlt_rain,                   &
                             longname='number melting to rain', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dnislf', dnislf,                   &
                             longname='self-collection ice', units='kg-1 s-1')
    CALL add_stream_element (remos, 'drpr', drpr,                   &
                             longname='rain production', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'drprn', drprn,                   &
                             longname='rain production number', units='kg-1 s-1')
    CALL add_stream_element (remos, 'daut', daut,                   &
                             longname='rain production (autoconversion)', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dacc', dacc,                   &
                             longname='rain production (accretion)', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dqrflx', dqrflx,                   &
                             longname='rain flux', units='units')
    CALL add_stream_element (remos, 'dqrevp', dqrevp,                   &
                             longname='evaporation of rain', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dxite_ls', dxite_ls,                   &
                             longname='ice from transport', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dicncte_ls', dicncte_ls,                   &
                             longname='ice number from transport', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dxite_cv', dxite_cv,                   &
                             longname='ice from convection', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'ddep_co', ddep_co,                   &
                             longname='sub gridscale deposition', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'ddep_a', ddep_a,                   &
                             longname='sub gridscale required deposition', units='kg kg-1')
    CALL add_stream_element (remos, 'ddep_aj', ddep_aj,                   &
                             longname='deposition mass from adjustment', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dni_aj', dni_aj,                   &
                             longname='deposition number from adjustment', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dxlte_ls', dxlte_ls,                   &
                             longname='liquid water from transport', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dcdncte_ls', dcdncte_ls,                   &
                             longname='liquid number water from transport', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dxlte_cv', dxlte_cv,                   &
                             longname='liquid water from convection', units='kg kg-1 s-1')    
    CALL add_stream_element (remos, 'dcnd_co', dcnd_co,                   &
                             longname='sub gridscale condensation', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dcnd_aj', dcnd_aj,                   &
                             longname='condensation mass from adjustment', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dnc_aj', dnc_aj,                   &
                             longname='condensation number from adjustment', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dxisub_cv', dxisub_cv,                   &
                             longname='clear-sky sub from cv', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dximlt_cv', dximlt_cv,                   &
                             longname='melting of ice from cv', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dxlevp_cv', dxlevp_cv,                   &
                             longname='clear-sky evp from cv', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dxisub_tp', dxisub_tp,                   &
                             longname='clear-sky sub from transport', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dnisub_tp', dnisub_tp,                   &
                             longname='clear-sky sub from transport (number)', units='kg-1')
    CALL add_stream_element (remos, 'dxlevp_tp', dxlevp_tp,                   &
                             longname='clear-sky evp from transport', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dncevp_tp', dncevp_tp,                   &
                             longname='clear-sky evp from transport (number)', units='kg-1')
    CALL add_stream_element (remos, 'dqccol', dqccol,                   &
                             longname='riming cloud droplets mass', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dnccol', dnccol,                   &
                             longname='riming cloud droplets number', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dncnuc', dncnuc,                   &
                             longname='cloud droplet nucleation', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dicnc_cv', dicnc_cv,                   &
                             longname='ICNC from convection', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dcdnc_cv', dcdnc_cv,                   &
                             longname='CDNC from convection', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dxisub_cc', dxisub_cc,                   &
                             longname='Sublimation adjustment to cover', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dnisub_cc', dnisub_cc,                   &
                             longname='Sublimation adjustment to cover (number)', units='kg-1')
    CALL add_stream_element (remos, 'dxlevp_cc', dxlevp_cc,                   &
                             longname='Evaporation adjustment to cover', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dncevp_cc', dncevp_cc,                   &
                             longname='Evaporation adjustment to cover (number)', units='kg-1')
    CALL add_stream_element (remos, 'dqisub', dqisub,                   &
                             longname='Sublimation', units='kg kg-1 s-1')
    CALL add_stream_element (remos, 'dnisub', dnisub,                   &
                             longname='Sublimation number', units='1/kg')
    CALL add_stream_element (remos, 'dni_lkp', dni_lkp,                   &
                             longname='Lookuptable overflow', units='-')

  ! Remos ice properties
    CALL add_stream_element (remos, 'dvtim', dvtim,                   &
                             longname='ice fall vel. mass-weighted', units='m s-1')
    CALL add_stream_element (remos, 'dvtin', dvtin,                   &
                             longname='ice fall vel. number-weighted', units='m s-1')
    CALL add_stream_element (remos, 'daclci', daclci,                   &
                             longname='ice cloud cover', units='(0-1)')
    CALL add_stream_element (remos, 'dcdnc', dcdnc,                   &
                             longname='cloud droplet number concentration', units='kg-1')
    CALL add_stream_element (remos, 'dicnc', dicnc,                   &
                             longname='ice crystal number concentration', units='kg-1')
    CALL add_stream_element (remos, 'dcdncb', dcdncb,                   &
                             longname='cloud droplet number concentration (in-cloud)', units='kg-1')
    CALL add_stream_element (remos, 'dicncb', dicncb,                   &
                             longname='ice crystal number concentration (in-cloud)', units='kg-1')
    CALL add_stream_element (remos, 'dxib', dxib,                   &
                             longname='ice mass mixing ratio (in-cloud)', units='kg kg-1')
    CALL add_stream_element (remos, 'dxlb', dxlb,                   &
                             longname='liquid mass mixing ratio (in-cloud)', units='kg kg-1')
    CALL add_stream_element (remos, 'drcm', drcm,                   &
                             longname='liquid radius mass-weighted', units='units')
    CALL add_stream_element (remos, 'drim', drim,                   &
                             longname='ice radius mass-weighted', units='units')
    CALL add_stream_element (remos, 'driv', driv,                   &
                             longname='ice volume radius', units='m')
    CALL add_stream_element (remos, 'drieff', drieff,                   &
                             longname='ice effective radius (in micro)', units='m')
    CALL add_stream_element (remos, 'drieff_optics', drieff_optics,                   &
                             longname='ice effective radius (in optics)', units='m', &
                             laccu=.FALSE., lrerun=.TRUE.)
    CALL add_stream_element (remos, 'dqirim', dqirim,                   &
                             longname='riming mass', units='kg kg-1')
    CALL add_stream_element (remos, 'dbirim', dbirim,                   &
                             longname='riming volume', units='m3')
    CALL add_stream_element (remos, 'dqihet', dqihet,                   &
                             longname='heterogeneously formed ice', units='kg kg-1')
    CALL add_stream_element (remos, 'dnihet', dnihet,                   &
                             longname='heterogeneously formed ice number', units='kg-1')
    CALL add_stream_element (remos, 'dnihom', dnihom,                   &
                             longname='homogeneously formed ice number', units='kg-1')
    CALL add_stream_element (remos, 'dninuc', dninuc,                   &
                             longname='nucleated ice number (cirrus-scheme)', units='kg-1')
    CALL add_stream_element (remos, 'dnidet', dnidet,                   &
                             longname='detrained ice number', units='kg-1')
    CALL add_stream_element (remos, 'dqiliq', dqiliq,                   &
                             longname='ice formed through liquid', units='kg kg-1')
    CALL add_stream_element (remos, 'dqioliq', dqioliq,                   &
                             longname='liquid origin ice', units='kg kg-1')
    CALL add_stream_element (remos, 'drhop', drhop,                   &
                             longname='riming density', units='kg m-3')
    CALL add_stream_element (remos, 'drhoice', drhoice,                   &
                             longname='mean ice density mass-weighted', units='kg m-3')
    CALL add_stream_element (remos, 'dxi_cloud', dxi_cloud,                   &
                             longname='Cloud portion of ice', units='kg kg-1')
    CALL add_stream_element (remos, 'dxi_snow', dxi_snow,                   &
                             longname='Snow portion of ice', units='kg kg-1')
    CALL add_stream_element (remos, 'dxivi_cloud', dxivi_cloud,                   &
                             longname='Cloud portion of v. integrated ice', units='kg m-2')
    CALL add_stream_element (remos, 'dxivi_snow', dxivi_snow,                   &
                             longname='Snow portion of v. integrated ice', units='kg m-2')
    CALL add_stream_element (remos, 'dnsub_co', dnsub_co,                   &
                             longname='Sublimation due to divergence', units='s-1')
    CALL add_stream_element (remos, 'dnevp_co', dnevp_co,                   &
                             longname='Evaporation due to divergence', units='s-1')
    CALL add_stream_element (remos, 'ddeltaqf', ddeltaqf,                   &
                             longname='Humidity forcing term: transport', units='kg-1s-1')
    CALL add_stream_element (remos, 'ddeltaqs', ddeltaqs,                   &
                             longname='Humidity forcing term: heating', units='kg-1s-1')
    CALL add_stream_element (remos, 'dqcdif', dqcdif,                   &
                             longname='Humidity forcing term: Sum', units='kg-1s-1')
    CALL add_stream_element (remos, 'dqisten', dqisten,                   &
                             longname='Ice sedimentation tendency: mass', units='kg-1s-1')
    CALL add_stream_element (remos, 'dnisten', dnisten,                   &
                             longname='Ice sedimentation tendency: number', units='s-1')
    CALL add_stream_element (remos, 'dqristen', dqristen,                   &
                             longname='Ice sedimentation tendency: rimed mass', units='kg-1s-1')
    CALL add_stream_element (remos, 'dbgsten', dbgsten,                   &
                             longname='Ice sedimentation tendency: rimed volume', units='m3s-1')
    CALL add_stream_element (remos, 'dqhsten', dqhsten,                   &
                             longname='Ice sedimentation tendency: heterogeneous ice', units='kg-1 s-1')
    CALL add_stream_element (remos, 'dice_mu', dice_mu,                   &
                             longname='Parameter mu of ice PSD', units='-')
    CALL add_stream_element (remos, 'dice_lam', dice_lam,                   &
                             longname='Parameter lambda of ice PSD', units='m-1')

  ! Remos budgets
    CALL add_stream_element (remos, 'dmete', dmete, &
                             longname='Microphysical energy change', units='W/m2')
    CALL add_stream_element (remos, 'dncdiff', dncdiff, &
                             longname='Microphysical number change', units='#')
    CALL add_stream_element (remos, 'dnidiff', dnidiff, &
                             longname='Microphysical number change', units='#')
    CALL add_stream_element (remos, 'dncdiff2', dncdiff2, &
                             longname='Microphysical number change 2', units='#')
    CALL add_stream_element (remos, 'dnidiff2', dnidiff2, &
                             longname='Microphysical number change 2', units='#')
    CALL add_stream_element (remos, 'dnisurf', dnisurf, &
                             longname='Ice number loss to ground', units='#')
    CALL add_stream_element (remos, 'dnictrl', dnictrl, &
                             longname='difference between ni 1 and 2', units='#')
    CALL add_stream_element (remos, 'dncctrl', dncctrl, &
                             longname='difference between nc 1 and 2', units='#')

  ! Remos scavenging
    CALL add_stream_element (remos, 'ain', dfrain, &
                             longname='rain flux before evaporation', units='kg/m2/s')
    CALL add_stream_element (remos, 'dfsnow', dfsnow, &
                             longname='snow flux before evaporation', units='kg/m2/s')
    CALL add_stream_element (remos, 'dfevapr', dfevapr, &
                             longname='evaporation of rain', units='kg/m2/s')
    CALL add_stream_element (remos, 'dfsubls', dfsubls, &
                             longname='sublimation of snow', units='kg/m2/s')
    CALL add_stream_element (remos, 'dmsnowacl', dmsnowacl, &
                             longname='accretion rate of snow with cloud droplets', units='kg/m2/s')
    CALL add_stream_element (remos, 'dmlwc', dmlwc, &
                             longname='cloud liquid content before rain', units='kg/kg') 
    CALL add_stream_element (remos, 'dmiwc', dmiwc, &
                             longname='cloud ice content before rain', units='kg/kg')
    CALL add_stream_element (remos, 'dmratepr', dmratepr, &
                             longname='rain formation rate in cloudy part', units='kg/m2/s')
    CALL add_stream_element (remos, 'dmrateps', dmrateps, &
                             longname='snow formation rate in cloudy part', units='kg/m2/s')

  ! Remos dynamics
    CALL add_stream_element (remos, 'dtke', dtke,                   &
                             longname='Turbulent Kinetic Energy', units='units')
    CALL add_stream_element (remos, 'dupdraft', dupdraft,                   &
                             longname='Updraft velocity', units='cm/s')
    CALL add_stream_element (remos, 'dupdraftmax', dupdraftmax,                   &
                             longname='Maximal updraft velocity', units='cm/s')
    CALL add_stream_element (remos, 'drhoair', drhoair,                   &
                             longname='Density of air', units='kg m-3')
    CALL add_stream_element (remos, 'dcfl', dcfl,                   &
                             longname='CFL-number for ice', units='')
    CALL add_stream_element (remos, 'ddpg', ddpg,                   &
                             longname='Delta p / g', units='kg m-2')
    CALL add_stream_element (remos, 'ddz', ddz,                   &
                             longname='Height difference across level', units='m')
    CALL add_stream_element (remos, 'dsedtm', dsedtm,                   &
                             longname='Sedimentation time (mass)', units='s')
    CALL add_stream_element (remos, 'dsedtn', dsedtn,                   &
                             longname='Sedimentation time (number)', units='s')
    CALL add_stream_element (remos, 'dnsedi', dnsedi,                   &
                             longname='Number of substeps (sedimentation)', units='#')  
    CALL add_stream_element (remos, 'dnmicro', dnmicro,                   &
                             longname='Number of substeps (micro)', units='#')

  ! Remos cloud types
    CALL add_stream_element (remos, 'ddp_cld', ddp_cld,                   &
                             longname='Pressure difference across cloud', units='Pa')
    CALL add_stream_element (remos, 'dt_top', dt_top,                   &
                             longname='Temperature at cloud top', units='K')
    CALL add_stream_element (remos, 'diclnb', diclnb,                   &
                             longname='Number of clouds in column', units='')
    CALL add_stream_element (remos, 'dclcol_time', dclcol_time,                   &
                             longname='Time with at least 1 cloud in column', units='')

  END SUBROUTINE init_remos
END MODULE mo_remos
