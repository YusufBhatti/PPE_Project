MODULE mo_ham_m7_processing

  ! *mo_aero_processing* contains subroutines for aerosol processing 
  ! in cloud droplets and ice crystals.
  ! Corinna Hoose, ETH Zurich, October 2005 - April 2008

  ! Basis: Aerosol processing in mixed-phase clouds in ECHAM5-HAM: 
  !        Model description and comparison to observations  !        Hoose et al, JGR 113, D07210, 2008

  ! Declan O'Donnell ETH Zuerich, may 2010: rewritten completely
  ! David Neubauer, ETH Zurich, October 2012 - Februar 2013: porting to ECHAM6-HAM

  USE mo_kind,         ONLY: dp
  USE mo_linked_list,  ONLY: t_stream
  USE mo_submodel_diag,     ONLY: vmem2d, vmem3d
  USE mo_math_constants,    ONLY: pi

  IMPLICIT NONE
!>>DN
  PRIVATE
!<<DN
  !---public member functions
  PUBLIC :: aeroproc_initialize, construct_stream_aeroproc, aeroproc_interface

  !---private member functions
  PRIVATE :: aeroproc_nuc, aeroproc_nuci, aeroproc_freeze, aeroproc_evap, aeroproc_melt, aeroproc_sub, aeroproc_coll

  REAL(dp), PARAMETER, PUBLIC :: sigma_cd = 2._dp             ! geometric standard deviation of in-cloud droplet size distribution
  REAL(dp), PARAMETER, PUBLIC :: sigma_ic = 2._dp             ! geometric standard deviation of in-ice crystal size distribution

  TYPE (t_stream), PUBLIC, POINTER :: aeroproc

  REAL(dp),        PUBLIC, POINTER :: quptk_incd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: quptk_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: quptk_inic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: quptk_ms4ic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qevap_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qfrz_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qevap_incd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcoll_incd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcoll_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcoll_inic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcoll_ms4ic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qfrz_incd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qmelt_incd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qmelt_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsub_inic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qsub_ms4ic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qwdep_incd(:,:,:)    !to be diagnosed in xt_wetdep
  REAL(dp),        PUBLIC, POINTER :: qwdep_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qwdep_inic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qwdep_ms4ic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qprod_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcor_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qcor_ms4ic(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qloss1_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qloss2_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: qloss3_ms4cd(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: pnprecip(:,:,:)
  REAL(dp),        PUBLIC, POINTER :: revap_vol(:,:,:)
!>>DN [additional diagnostics]
  REAL(dp),        PUBLIC, POINTER :: devap_incd(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsub_inic(:,:)
  REAL(dp),        PUBLIC, POINTER :: duptk_incd(:,:)
  REAL(dp),        PUBLIC, POINTER :: duptk_inic(:,:)
  REAL(dp),        PUBLIC, POINTER :: devap_ncd(:,:)
  REAL(dp),        PUBLIC, POINTER :: dsub_nic(:,:)
  REAL(dp),        PUBLIC, POINTER :: dnuc_ncd(:,:)
  REAL(dp),        PUBLIC, POINTER :: dnuc_nic(:,:)
!<<DN
!  REAL(dp),        PUBLIC, POINTER :: rdry_precip(:,:,:)        !radius of particles released from evap. precip              
!!$  !variables which are set in cloud_cdnc_icnc
!!$  REAL(dp),        PUBLIC, POINTER :: pnucnic(:,:,:)
!!$  REAL(dp),        PUBLIC, POINTER :: pnucncd(:,:,:)
!!$  REAL(dp),        PUBLIC, POINTER :: pcdncactevap(:,:,:)    
!!$  REAL(dp),        PUBLIC, POINTER :: pevapncd(:,:,:)
!!$  REAL(dp),        PUBLIC, POINTER :: pevapncdbf(:,:,:)
!!$  REAL(dp),        PUBLIC, POINTER :: pcdncactfrz(:,:,:)     
!!$  REAL(dp),        PUBLIC, POINTER :: pfrzncd(:,:,:)
!!$  REAL(dp),        PUBLIC, POINTER :: picncactmelt(:,:,:)    
!!$  REAL(dp),        PUBLIC, POINTER :: pmeltncd(:,:,:)
!!$  REAL(dp),        PUBLIC, POINTER :: picncactsub(:,:,:)     
!!$  REAL(dp),        PUBLIC, POINTER :: psubnic(:,:,:)

  !vertically integrated rates for all in-cloud species
  TYPE (vmem2d), PUBLIC, POINTER :: devap(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dsubl(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: duptkcd(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: duptkic(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dcollcd(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dcollic(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dfrz(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dmelt(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dloss1(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dloss2(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dloss3(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dcorrcd(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dcorric(:)              
!>>DN [additional diagnostics]
  TYPE (vmem2d), PUBLIC, POINTER :: dreevap(:)              
  TYPE (vmem2d), PUBLIC, POINTER :: dresubl(:)              
!<<DN

  TYPE(vmem3d), PUBLIC, POINTER :: qbceva(:)

  !---dry radius of aerosol in cloud droplets / ice crystals
  REAL(dp), PUBLIC, POINTER :: rdry_incd(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rdry_inic(:,:,:)

  REAL(dp), PARAMETER, PRIVATE :: zeps=EPSILON(1._dp)
  REAL(dp), PARAMETER, PRIVATE :: zeps_massratio=1.206e-11_dp   !m/N for d=5mum dust particle
  REAL(dp), PARAMETER, PRIVATE :: zeps_mass=1.e-30_dp           !small number for aerosol mass
  REAL(dp), PARAMETER, PRIVATE :: zeps_num=1._dp                !threshold for evaporation (#/kg)

  REAL(dp), PARAMETER, PRIVATE :: astoso4 = 98.0734_dp/32.07_dp
  REAL(dp), PARAMETER, PRIVATE :: a4piover3 = 4._dp*pi/3._dp

  CHARACTER(len=20), PRIVATE :: cunit




!--------------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------
  SUBROUTINE aeroproc_initialize
    
    ! aeroproc_initialize allocates arrays that hold the tracer identities
    ! for in-cloud droplet and in-ice crystal aerosol

    USE mo_ham_species,   ONLY: idt_cd, idt_ic
    USE mo_species,        ONLY: naerospec,aero_idx

    !---allocate the tracer id arrays
    ALLOCATE(idt_cd(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(idt_ic(aero_idx(1):aero_idx(naerospec)))
    
  END SUBROUTINE aeroproc_initialize

!-----------------------------------------------------------------------
  SUBROUTINE construct_stream_aeroproc

    ! *construct_stream_aeroproc* allocates output streams
    ! for the aerosol processing schemes

    USE mo_memory_base,   ONLY: new_stream, add_stream_element, AUTO,  &
                                default_stream_setting, add_stream_reference
    USE mo_linked_list,   ONLY: HYBRID, SURFACE
    USE mo_filename,      ONLY: trac_filetype
!>>DN
!    USE mo_aero_species,  ONLY: naerospec, aerospec
    USE mo_species,       ONLY: naerospec, speclist, aero_idx
    USE mo_ham_species,   ONLY: id_wat
!<<DN
!    USE mo_outctl,        ONLY: lproc_detail
    USE mo_ham_m7ctl,    ONLY: lcoll
    
    INTEGER :: jn

    !--- Create new stream:

    CALL new_stream (aeroproc ,'aerproc',filetype=trac_filetype)


    !--- Add standard fields for post-processing:

    CALL add_stream_reference (aeroproc, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (aeroproc, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (aeroproc, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (aeroproc, 'gboxarea','geoloc',lpost=.TRUE.)

    CALL default_stream_setting (aeroproc, lpost     = .TRUE. , &
                                           lrerun    = .TRUE. , & 
                                           laccu     = .TRUE. , &
                                           leveltype = HYBRID , &
                                           table     = 199,     &
                                           code      = AUTO     )
    CALL add_stream_element (aeroproc,   'QUPTK_INCD',     quptk_incd,                   &
                             longname='Uptake rate of aerosol number into CD', &
                             units='m-3 s-1' )
    CALL add_stream_element (aeroproc,   'QUPTK_MS4CD',     quptk_ms4cd,                   &
                             longname='Uptake rate of sulfate aerosol mass into CD', &
                             units='kg m-3 s-1' )
    CALL add_stream_element (aeroproc,   'QUPTK_INIC',     quptk_inic,                   &
                             longname='Uptake rate of aerosol number into IC', &
                             units='m-3 s-1' )
    CALL add_stream_element (aeroproc,   'QUPTK_MS4IC',     quptk_ms4ic,                   &
                             longname='Uptake rate of sulfate aerosol mass into IC', &
                             units='kg m-3 s-1' )
!pr    CALL add_stream_element (aeroproc,   'RDRY_PRECIP',      rdry_precip,      &
!pr                             longname='Radius of particle released from evaporating precipitation',         &
!pr                             units='m' )


!    IF (lproc_detail) THEN
!!>>DN [SO4 in kg(S) kg-1]
       CALL add_stream_element (aeroproc,   'QEVAP_MS4CD',     qevap_ms4cd,                   &
                                longname='Loss rate of sulfate aerosol from CD evaporation', &
                                units='kg m-3 s-1' )
!!$       CALL add_stream_element (aeroproc,   'QEVAP_MS4CD',     qevap_ms4cd,                   &
!!$                                longname='Loss rate of sulfate aerosol from CD evaporation', &
!!$                                units='kg(S) m-3 s-1' )
!!<<DN
       CALL add_stream_element (aeroproc,   'QFRZ_MS4CD',     qfrz_ms4cd,                   &
                                longname='Loss rate of sulfate aerosol from CD freezing', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QEVAP_INCD',      qevap_incd,                   &
                                longname='Loss rate of aerosol inside CD by evaporation', &
                                units='m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QCOLL_INCD',      qcoll_incd,                   &
                                longname='Uptake rate of aerosol inside CD by collisions', &
                                units='m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QCOLL_MS4CD',      qcoll_ms4cd,                   &
                                longname='Uptake rate of aerosol inside CD by collisions', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QFRZ_INCD',      qfrz_incd,                   &
                                longname='Loss rate of aerosol inside CD by freezing', &
                                units='m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QMELT_INCD',      qmelt_incd,                   &
                                longname='Gain rate of aerosol inside CD by melting of IC', &
                                units='m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QMELT_MS4CD',      qmelt_ms4cd,                   &
                                longname='Gain rate of aerosol inside CD by melting of IC', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QCOLL_INIC',      qcoll_inic,                   &
                                longname='Uptake rate of aerosol into IC by collisions', &
                                units='m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QCOLL_MS4IC',      qcoll_ms4ic,                   &
                                longname='Uptake rate of sulfate into IC by collisions', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QSUB_INIC',      qsub_inic,                   &
                               longname='Loss rate of aerosol inside IC by sublimation', &
                               units='m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QSUB_MS4IC',      qsub_ms4ic,                   &
                                longname='Loss rate of aerosol inside IC by sublimation', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QWDEP_INCD',      qwdep_incd,                   &
                                longname='Loss rate of aerosol inside CD by wet deposition', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QWDEP_MS4CD',      qwdep_ms4cd,                   &
                                longname='Loss rate of aerosol inside CD by wet deposition', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QWDEP_INIC',      qwdep_inic,                   &
                                longname='Loss rate of aerosol inside IC by wet deposition', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QWDEP_MS4IC',      qwdep_ms4ic,                   &
                                longname='Loss rate of aerosol inside IC by wet deposition', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QPROD_MS4CD',      qprod_ms4cd,                   &
                                longname='Production rate of in-droplet sulfate', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QCOR_MS4CD',      qcor_ms4cd,                   &
                                longname='Removal of in-droplet sulfate for consistency', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QCOR_MS4IC',      qcor_ms4ic,                   &
                                longname='Removal of in-droplet sulfate for consistency', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QLOSS1_MS4CD',      qloss1_ms4cd,                   &
                                longname='Loss of in-droplet sulfate during evaporation (wrong diameter)', &
                               units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QLOSS2_MS4CD',      qloss2_ms4cd,                   &
                                longname='Loss of in-droplet sulfate during evaporation (wrong diameter)', &
                                units='kg m-3 s-1' )
       CALL add_stream_element (aeroproc,   'QLOSS3_MS4CD',      qloss3_ms4cd,                   &
                                longname='Loss of in-droplet sulfate during evaporation (wrong diameter)', &
                                units='kg m-3 s-1' )

       CALL default_stream_setting (aeroproc, laccu     = .FALSE. )

       CALL add_stream_element (aeroproc,   'pnprecip',      pnprecip,                   &
                                longname='estimated number of precipitating hydrometeors', &
                                units='[N] kg-1(air)' )
       CALL add_stream_element (aeroproc,   'revap_vol',      revap_vol,                   &
                                longname='volume released from evaporation of precipitating hydrometeors', &
                                units='m3 kg-1(air)' )
!>>DN [additional diagnostics]

       CALL default_stream_setting (aeroproc, leveltype = SURFACE, lpost=.TRUE., lrerun=.TRUE., laccu=.TRUE. )
       CALL add_stream_element (aeroproc,   'DEVAP_INCD',       devap_incd,       &
                                longname='Burden of number of evap. from CD',          &
                                units='[N] m-2 s-1' )
       CALL add_stream_element (aeroproc,   'DSUB_INIC',       dsub_inic,       &
                                longname='Burden of number of subl. from IC',          &
                                units='[N] m-2 s-1' )
       CALL add_stream_element (aeroproc,   'DUPTK_INCD',       duptk_incd,       &
                                longname='Burden of number uptake into CD',          &
                                units='[N] m-2 s-1' )
       CALL add_stream_element (aeroproc,   'DUPTK_INIC',       duptk_inic,       &
                                longname='Burden of number uptake into IC',          &
                                units='[N] m-2 s-1' )
       CALL add_stream_element (aeroproc,   'DEVAP_NCD',       devap_ncd,       &
                                longname='Burden of number of evap. CD',          &
                                units='[N] m-2 s-1' )
       CALL add_stream_element (aeroproc,   'DSUB_NIC',       dsub_nic,       &
                                longname='Burden of number of subl. IC',          &
                                units='[N] m-2 s-1' )
       CALL add_stream_element (aeroproc,   'DNUC_NCD',       dnuc_ncd,       &
                                longname='Burden of number of nucl. CD',          &
                                units='[N] m-2 s-1' )
       CALL add_stream_element (aeroproc,   'DNUC_NIC',       dnuc_nic,       &
                                longname='Burden of number of nucl. IC',          &
                                units='[N] m-2 s-1' )
!<<DN

!    END IF

    !---fields used for online calculations in cloud_cdnc_icnc...which should maybe be moved here...
!>>DN [laccu must be .FALSE. for these stream elements]
!!    CALL default_stream_setting (aeroproc, lpost=.FALSE., lrerun=.FALSE., laccu=.FALSE.)
    CALL default_stream_setting (aeroproc, laccu=.FALSE., leveltype = HYBRID)
!<<DN
!!$    CALL add_stream_element (aeroproc, 'NUC_NCD', pnucncd)
!!$    CALL add_stream_element (aeroproc, 'NUC_NIC', pnucnic)
!!$    CALL add_stream_element (aeroproc, 'EVAP_NCD', pevapncd)
!!$    CALL add_stream_element (aeroproc, 'EVAP_NCD_BF', pevapncdbf)
!!$    CALL add_stream_element (aeroproc, 'FRZ_NCD', pfrzncd)
!!$    CALL add_stream_element (aeroproc, 'MELT_NCD', pmeltncd)
!!$    CALL add_stream_element (aeroproc, 'SUB_NIC', psubnic)
!!$    CALL add_stream_element (aeroproc, 'CDNC_ACT_EVAP', pcdncactevap)
!!$    CALL add_stream_element (aeroproc, 'CDNC_ACT_FRZ', pcdncactfrz)
!!$    CALL add_stream_element (aeroproc, 'ICNC_ACT_MELT', picncactmelt)
!!$    CALL add_stream_element (aeroproc, 'ICNC_ACT_SUB', picncactsub)
!!    CALL add_stream_element (aeroproc, 'IN-CLOUD PARTICLE DRY RADIUS', rdry_incd)
!>>DN [added long name and units]
    CALL add_stream_element (aeroproc, 'INCLOUDPARTICLEDRYRADIUS', rdry_incd, longname='Dry particle radius in cloud droplets', &
                             units='m')
    CALL add_stream_element (aeroproc, 'INCRYSTALPARTICLEDRYRADIUS', rdry_inic, longname='Dry particle radius in ice crystals', &
                             units='m')
!<<DN

    !---more diagnostics: vertically integrated rates for all in-cloud species
    ALLOCATE(duptkcd(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(duptkic(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dcollcd(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dcollic(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dfrz(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dmelt(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(devap(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dsubl(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dloss1(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dloss2(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dloss3(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dcorrcd(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dcorric(aero_idx(1):aero_idx(naerospec)))
!>>DN [additional diagnostics]   
    ALLOCATE(dreevap(aero_idx(1):aero_idx(naerospec)))
    ALLOCATE(dresubl(aero_idx(1):aero_idx(naerospec)))
!<<DN

    ALLOCATE(qbceva(aero_idx(1):aero_idx(naerospec)))
    
    CALL default_stream_setting (aeroproc, leveltype = SURFACE, lpost=.TRUE., lrerun=.TRUE., laccu=.TRUE.) 
    cunit='kg m-2 s-1'

    DO jn=aero_idx(1),aero_idx(naerospec)
!>>DN [no water tracer in cloud particles]
       IF (jn==id_wat) CYCLE
!<<DN

       CALL add_stream_element(aeroproc, 'D_UPTKCD_'//TRIM(speclist(jn)%shortname), duptkcd(jn)%ptr, units=cunit,         &
                               longname='uptake into cloud droplets'      )
       CALL add_stream_element(aeroproc, 'D_UPTKIC_'//TRIM(speclist(jn)%shortname), duptkic(jn)%ptr, units=cunit,         &
                               longname='uptake into ice crystals'      )
       CALL add_stream_element(aeroproc, 'D_COLLCD_'//TRIM(speclist(jn)%shortname),   dcollcd(jn)%ptr, units=cunit,         &
                               longname='collisions with cloud droplets'    )
       CALL add_stream_element(aeroproc, 'D_COLLIC_'//TRIM(speclist(jn)%shortname),   dcollic(jn)%ptr, units=cunit,         &
                               longname='collisions with ice crystals'    )
       CALL add_stream_element(aeroproc, 'D_FRZ_'//TRIM(speclist(jn)%shortname),      dfrz(jn)%ptr, units=cunit,         &
                               longname='freezing of cloud droplets'  )
       CALL add_stream_element(aeroproc, 'D_MELT_'//TRIM(speclist(jn)%shortname),     dmelt(jn)%ptr, units=cunit,         &
                               longname='melting of crystals'         )
       CALL add_stream_element(aeroproc, 'D_EVAP_'//TRIM(speclist(jn)%shortname),     devap(jn)%ptr, units=cunit,         &
                               longname='evaporation from cloud droplets' )
       CALL add_stream_element(aeroproc, 'D_SUB_'//TRIM(speclist(jn)%shortname),     dsubl(jn)%ptr, units=cunit,          &
                               longname='evaporation from ice crystals' )
       CALL add_stream_element(aeroproc, 'D_LOSS1_'//TRIM(speclist(jn)%shortname),    dloss1(jn)%ptr, units=cunit,         &
                               longname='loss of particles smaller than nucleation mode (???)'    )
       CALL add_stream_element(aeroproc, 'D_LOSS2_'//TRIM(speclist(jn)%shortname),    dloss2(jn)%ptr, units=cunit,         &
                               longname='loss of particles smaller than Aitken mode'        )
       CALL add_stream_element(aeroproc, 'D_LOSS3_'//TRIM(speclist(jn)%shortname),    dloss3(jn)%ptr, units=cunit,         &
                               longname='loss of particles larger than coarse mode'         )
       CALL add_stream_element(aeroproc, 'D_CORRCD_'//TRIM(speclist(jn)%shortname),   dcorrcd(jn)%ptr, units=cunit,         &
                               longname='correction for negative tracer values CD'    )
       CALL add_stream_element(aeroproc, 'D_CORRIC_'//TRIM(speclist(jn)%shortname),   dcorric(jn)%ptr, units=cunit,         &
                               longname='correction for negative tracer values IC'    )
!>>DN [additional diagnostics]   
       CALL add_stream_element(aeroproc, 'D_REEVAP_'//TRIM(speclist(jn)%shortname),     dreevap(jn)%ptr, units=cunit,         &
                               longname='re-evaporation from cloud droplets' )
       CALL add_stream_element(aeroproc, 'D_RESUB_'//TRIM(speclist(jn)%shortname),     dresubl(jn)%ptr, units=cunit,          &
                               longname='re-evaporation from ice crystals' )
!<<DN

!       IF (lproc_detail) &
            CALL add_stream_element(aeroproc,   'QBCEVA_'//TRIM(speclist(jn)%shortname),    qbceva(jn)%ptr,  &
                                    longname='Below-cloud evaporation of precipitating in-droplet '//&
                                    &TRIM(speclist(jn)%shortname),  &
                                    units='kg m-3 s-1' )
    END DO

  END SUBROUTINE construct_stream_aeroproc

  !------------------------------------------------------------------------

  ! SUBROUTINES FOR MICROPHYSICAL PROCESSES THAT CHANGE THE MASS/NUMBER OF 
  ! CLOUDBORNE PARTICLES

  !------------------------------------------------------------------------
  SUBROUTINE aeroproc_interface(kbdim,  kproma,  klev,  krow,     & 
                                paclc,  pxtm1,  pxtte,            &
!davidn                                paphp1,     papp1)!davidn
                                paphp1,     papp1, ptm1)!davidn
  ! 
  ! Microphysical processes that change the mass/number of cloudborne particles
  !
  ! This routine is called from cloud_cdnc_icnc...or xt_wetdep...not yet decided
  !
    USE mo_kind,            ONLY: dp 
    USE mo_tracdef,         ONLY: ntrac, trlist, AEROSOLMASS, AEROSOLNUMBER, &
                                  IN_CLOUD_LIQUID, IN_CLOUD_ICE
    USE mo_ham_species,     ONLY: id_so4, idt_cd, idt_ic, id_bc, id_ss, id_du, id_wat, id_oc !davidn
    USE mo_time_control,    ONLY: time_step_len,  get_time_step, delta_time
    USE mo_ham_m7ctl,       ONLY: lcoll, icoas, iaits, iaccs
    USE mo_ham,             ONLY: sizeclass
    USE mo_control,         ONLY: ltimer
    USE mo_timer,           ONLY: timer_start, timer_stop,                    &
                                  timer_aeroproc_nuc, timer_aeroproc_evap,    &
                                  timer_aeroproc_melt, timer_aeroproc_freeze, &
                                  timer_aeroproc_subl, timer_aeroproc_coll,   &
!>>DN
                                  timer_aeroproc_nuci
!<<DN
    USE mo_activ,           ONLY: idt_cdnc, idt_icnc
    USE mo_tracer_processes,ONLY: xt_borrow !davidn
    USE mo_activ,           ONLY: pnucncd,pcdncactevap,pnucnic,                    &
                                  pevapncd,pevapncdbf,                             &
                                  pcdncactfrz,pfrzncd,                             &
                                  picncactmelt,pmeltncd,                           &
                                  picncactsub,psubnic, pcdncactap,                 &
                                  picncactap, dpg, cloud_cover_begin_micro_2m
  USE mo_physical_constants,     ONLY: tmelt !davidn
  USE mo_echam_cloud_params,         ONLY: cthomi !davidn
  USE mo_ham_streams, ONLY: nact_strat !davidn
  USE mo_vphysc,          ONLY: vphysc !davidn

    !---subroutine argument list
    INTEGER, INTENT(IN)  :: kbdim                        ! geographic block maximum number of locations
    INTEGER, INTENT(IN)  :: kproma                       ! geographic block number of locations on this row
    INTEGER, INTENT(IN)  :: krow                         ! geographic block number
    INTEGER, INTENT(IN)  :: klev                         ! number of vertical levels
    REAL(dp), INTENT(IN) :: paclc(kbdim,klev)            ! fractional cloud cover
    REAL(dp), INTENT(IN) :: pxtm1(kbdim,klev,ntrac)      ! tracer concentrations at t-dt 
    REAL(dp), INTENT(INOUT) :: pxtte(kbdim,klev,ntrac)   ! tracer concentrations tendency
    REAL(dp):: papp1(kbdim,klev),paphp1(kbdim,klev+1) !davidn
  REAL(dp), INTENT(in)    :: ptm1     (kbdim,klev)       !davidn

    !---local data
    REAL(dp) :: zxtp1(kbdim,klev,ntrac)                  ! tracer concentraion at t+dt
    REAL(dp) :: pxtp1c(kbdim,klev,ntrac)                 ! tracer concentrations in the cloudy fraction at t+dt
    REAL(dp) :: pxtp10(kbdim,klev,ntrac)                 ! tracer concentrations in the cloud-free fraction at t+dt
    REAL(dp) :: pcdnc    (kbdim,klev)        ! Cloud droplet number concentration (CDNC) [1/m3]
    REAL(dp) :: picnc    (kbdim,klev)        ! Ice crystal number concentration (ICNC) [1/m3]
    REAL(dp) :: pdpg(kbdim,klev)             ! air mass per unit area
    INTEGER  :: jt,lev,jtcoas, jtaccs, jtaits,idt_no, jn
    INTEGER :: ispec,jk,jl !davidn
    REAL(dp) :: zscavnuc(kbdim,klev,3,3,2),zscavnuci(kbdim,klev,3,3,2),zscavcoll(kbdim,klev,3,3,2),&
         zaerm(kbdim,klev,3,3,2),zsratio(kbdim,klev,3,3,2),zam(kbdim,klev,3,3,2)!davidn
    REAL(dp) :: zaclc(kbdim,klev)
    !--------------------------------------------------------------------------

    jtaits = sizeclass(iaits)%idt_no
    jtaccs = sizeclass(iaccs)%idt_no
    jtcoas = sizeclass(icoas)%idt_no

    pdpg(1:kproma,:)=dpg(1:kproma,:,krow)
    pxtp1c(1:kproma,:,idt_cdnc)=pcdncactap(1:kproma,:,krow)
    pxtp1c(1:kproma,:,idt_icnc)=picncactap(1:kproma,:,krow)
    pxtp10(1:kproma,:,idt_cdnc)=0._dp
    pxtp10(1:kproma,:,idt_icnc)=0._dp
    
    zscavnuc(1:kproma,:,:,:,:)=0._dp !davidn
    zscavnuci(1:kproma,:,:,:,:)=0._dp !davidn
    zscavcoll(1:kproma,:,:,:,:)=0._dp !davidn
    zaerm(1:kproma,:,:,:,:)=0._dp !davidn
    zsratio(1:kproma,:,:,:,:)=0._dp !davidn
    zam(1:kproma,:,:,:,:)=0._dp !davidn
    zaclc(1:kproma,:)=cloud_cover_begin_micro_2m(1:kproma,:,krow)

    WHERE(pxtp1c(1:kproma,:,idt_cdnc)<1E-9)
       pxtp1c(1:kproma,:,idt_cdnc)=0.0_dp
    END WHERE
    WHERE(pxtp1c(1:kproma,:,idt_icnc)<1E-9)
       pxtp1c(1:kproma,:,idt_icnc)=0.0_dp
    END WHERE

    !---set tracer in cloudy and cloud-free according to cloud cover
    DO jt=1,ntrac
       IF (trlist%ti(jt)%nphase == AEROSOLMASS .OR. &
           trlist%ti(jt)%nphase == AEROSOLNUMBER ) THEN

          zxtp1(1:kproma,:,jt) = pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt) * time_step_len
          !---weight the tracers by the cloud fraction 
          pxtp1c(1:kproma,:,jt) = zxtp1(1:kproma,:,jt) * zaclc(1:kproma,:)
          pxtp10(1:kproma,:,jt) = zxtp1(1:kproma,:,jt) * (1._dp - zaclc(1:kproma,:))

       !---in-cloud tracers:
       ELSE IF (trlist%ti(jt)%nphase == IN_CLOUD_LIQUID .OR. trlist%ti(jt)%nphase == IN_CLOUD_ICE) THEN
          zxtp1(1:kproma,:,jt) = pxtm1(1:kproma,:,jt) + pxtte(1:kproma,:,jt) * time_step_len

   !---2-D diagnostics for each species
          ispec=trlist%ti(jt)%spid
          DO jk=1,klev
             IF (trlist%ti(jt)%nphase == IN_CLOUD_LIQUID)THEN
                dcorrcd(ispec)%ptr(1:kproma,krow) = dcorrcd(ispec)%ptr(1:kproma,krow)    -           &
                     MAX(0._dp,-zxtp1(1:kproma,jk,jt))*pdpg(1:kproma,jk)*delta_time
             ELSE
                dcorric(ispec)%ptr(1:kproma,krow) = dcorric(ispec)%ptr(1:kproma,krow)    -           &
                     MAX(0._dp,-zxtp1(1:kproma,jk,jt))*pdpg(1:kproma,jk)*delta_time
             END IF
          END DO

         WHERE(zxtp1(1:kproma,:,jt)<0._dp)!davidn
             pxtte(1:kproma,:,jt)=-pxtm1(1:kproma,:,jt)/time_step_len!davidn
          END WHERE!davidn
          zxtp1(1:kproma,:,jt) = MAX(0._dp,zxtp1(1:kproma,:,jt))!davidn

          pxtp1c(1:kproma,:,jt) = zxtp1(1:kproma,:,jt)
          pxtp10(1:kproma,:,jt) = 0._dp    

       END IF
    END DO

    !---cloud droplet nucleation
    IF (ltimer) CALL timer_start(timer_aeroproc_nuc)
!>>DNdebug
!    CALL aeroproc_nuc(kbdim, kproma, klev, krow, pdpg, zxtp1, & 
    CALL aeroproc_nuc(kbdim, kproma, klev, krow, pdpg, pxtm1, & 
                      pxtp1c, pnucncd(:,:,krow))
!<<DNdebug
    IF (ltimer) CALL timer_stop(timer_aeroproc_nuc)
!>>DN [separate routine for ice crystal nucleation]
    !---ice crystal nucleation
    IF (ltimer) CALL timer_start(timer_aeroproc_nuci)
!    CALL aeroproc_nuci(kbdim, kproma, klev, krow, pdpg, zxtp1, & 
    CALL aeroproc_nuci(kbdim, kproma, klev, krow, pdpg, pxtm1, & 
         pxtp1c, pnucnic(:,:,krow)) 
    IF (ltimer) CALL timer_stop(timer_aeroproc_nuci)
!<<DN

    !---freezing of cloud droplets
    IF (ltimer) CALL timer_start(timer_aeroproc_freeze)
    CALL aeroproc_freeze(kbdim, kproma, klev, krow, pdpg,       & 
         pxtp1c, pcdncactfrz(:,:,krow), pfrzncd(:,:,krow))
    IF (ltimer) CALL timer_stop(timer_aeroproc_freeze)

  !---evaporation of cloud droplets
    IF (ltimer) CALL timer_start(timer_aeroproc_evap)
    CALL aeroproc_evap(kbdim, kproma, klev, krow, pdpg,       & 
         pxtp1c, pxtp10, pcdncactevap(:,:,krow), pevapncd(:,:,krow)) 
    IF (ltimer) CALL timer_stop(timer_aeroproc_evap)

    !---melting of ice crystals
    IF (ltimer) CALL timer_start(timer_aeroproc_melt)
    CALL aeroproc_melt(kbdim, kproma, klev, krow, pdpg,       & 
                       pxtp1c, picncactmelt(:,:,krow), pmeltncd(:,:,krow)) 
    IF (ltimer) CALL timer_stop(timer_aeroproc_melt)

    !---sublimation of ice crystals
    IF (ltimer) CALL timer_start(timer_aeroproc_subl)
    CALL aeroproc_sub(kbdim, kproma, klev, krow, pdpg,       & 
         pxtp1c, pxtp10, picncactsub(:,:,krow), psubnic(:,:,krow)) 
    IF (ltimer) CALL timer_stop(timer_aeroproc_subl)

    !---in cloud collisions
    IF (lcoll) THEN
       IF (ltimer) CALL timer_start(timer_aeroproc_coll)
       CALL aeroproc_coll  (kbdim, kproma, klev, krow, pdpg, pxtp1c) 
       IF (ltimer) CALL timer_stop(timer_aeroproc_coll)
    END IF

    !---update tracer tendencies
    DO jt=1,ntrac
       IF (trlist%ti(jt)%nphase == AEROSOLMASS .OR. &
           trlist%ti(jt)%nphase == AEROSOLNUMBER .OR. &
           trlist%ti(jt)%nphase == IN_CLOUD_LIQUID .OR. &
           trlist%ti(jt)%nphase == IN_CLOUD_ICE ) THEN
!>>DN [1/time_step_len was missing]
!          pxtte(1:kproma,:,jt) = pxtte(1:kproma,:,jt) + pxtp1c(1:kproma,:,jt) + &
!                                 pxtp10(1:kproma,:,jt) - zxtp1(1:kproma,:,jt)
          pxtte(1:kproma,:,jt) = pxtte(1:kproma,:,jt) + (pxtp1c(1:kproma,:,jt) + &
                                 pxtp10(1:kproma,:,jt) - zxtp1(1:kproma,:,jt))/time_step_len
!davidn          pxtte(1:kproma,:,jt) = pxtte(1:kproma,:,jt) + (pxtp1c(1:kproma,:,jt)*paclc(1:kproma,:) + &
!davidn                                 pxtp10(1:kproma,:,jt)*(1._dp - paclc(1:kproma,:))- zxtp1(1:kproma,:,jt))/time_step_len
!<<DN          
!       END IF
       END IF

    END DO

  !--- Mass conserving correction of negative tracer values:

    CALL xt_borrow(kproma, kbdim,  klev,  klev+1, ntrac, &!davidn
                   papp1,  paphp1,                       &
                   pxtm1,  pxtte                         )

  END SUBROUTINE aeroproc_interface
  !--------------------------------------------------------------------------

  SUBROUTINE aeroproc_nuc(kbdim, kproma, klev, krow, pdpg,             &
                          pxtp1, pxtp1c, pqnuc)

    USE mo_kind,            ONLY: dp
    USE mo_tracdef,         ONLY: ntrac
    USE mo_time_control,    ONLY: time_step_len, delta_time
    USE mo_ham_streams,     ONLY: nact, densaer, rwet
    USE mo_ham,             ONLY: aerocomp, sizeclass, naerocomp, aerowater, &
                                  nclass, nsol
    USE mo_ham_m7ctl,       ONLY: sigmaln, &
                                  iaits, iaccs, icoas, crdiv, cmr2ram
    USE mo_ham_species,     ONLY: idt_cd, id_so4
    USE mo_vphysc,          ONLY: vphysc
    USE mo_time_control,    ONLY: current_date, get_time_step
    USE mo_activ,           ONLY: idt_cdnc!PR testing
    USE mo_species,         ONLY: naerospec
    USE mo_physical_constants, ONLY: tmelt !davidn
    USE mo_echam_cloud_params, ONLY: cthomi !davidn
    USE mo_ham_species,     ONLY: id_bc, id_oc, id_ss, id_du !davidn

    !---subroutine argument list
    INTEGER, INTENT(IN)  :: kbdim                        ! geographic block maximum number of locations
    INTEGER, INTENT(IN)  :: kproma                       ! geographic block number of locations on this row
    INTEGER, INTENT(IN)  :: krow                         ! geographic block number
    INTEGER, INTENT(IN)  :: klev                         ! number of vertical levels
    REAL(dp), INTENT(IN) :: pdpg(kbdim,klev)             ! air mass per unit area
    REAL(dp), INTENT(IN) :: pqnuc(kbdim,klev)            ! number of nucleated cloud droplets this time step [# m-3]
    REAL(dp), INTENT(IN) :: pxtp1(kbdim,klev,ntrac)      ! tracer concentrations at t+dt [mmr]
    REAL(dp), INTENT(INOUT) :: pxtp1c(kbdim,klev,ntrac)  ! in-cloud tracer concentrations [kg m-3]

    !PR
    REAL(dp) :: pxtm1(kbdim,klev,ntrac)      ! tracer concentrations at t-dt
    REAL(dp) :: zwork(kbdim,klev)
    REAL(dp) :: zworkacc(kbdim,klev,ntrac)
    REAL(dp) :: pncd(kbdim,klev)
!    INTEGER  :: jspec,idt_aer,idt_num,idt_no,naerospec,imod,imod_crdiv
    INTEGER  :: jspec,idt_aer,idt_num,idt_no,imod,imod_crdiv
    !PRend
    !---local data    
    REAL(dp) :: zmodmass(kbdim,klev,nclass)                 ! total aerosol mass per mode
    REAL(dp) :: znactsum(kbdim,klev),      &              ! number of activated particles in all modes
                zmnratio(kbdim,klev),      &              ! mass to number conversion 
                znact(kbdim,klev,nclass)                    ! number of activated particles in each mode

    INTEGER :: jl,jk,jt, jmod, jmod_crdiv, lev

    REAL(dp), POINTER :: duptk_p(:,:) 

    REAL(dp) :: zmodefac(nclass)
    INTEGER :: jn, jtcl, jtaer, jtnum, ispec
    INTEGER :: jtcoas, jtaccs, jtaits
    LOGICAL :: lnumupdated_cd(nclass)

    REAL(dp) :: zcdworknum(kbdim,klev)
    REAL(dp) :: zcdworkmass(kbdim,klev)
    REAL(dp) :: zworkdiagnum(kbdim,klev)
    REAL(dp) :: zworkdiagmass(kbdim,klev)
    LOGICAL :: llwork(kbdim,klev)
 
    !-----------------------------------------------------------------------------------
    !---initialisations
    znactsum(1:kproma,:) = 0._dp
    zcdworkmass(1:kproma,:) = 0._dp
    zworkdiagnum(1:kproma,:) = 0._dp
    zworkdiagmass(1:kproma,:) = 0._dp
!>>DNdebug
    lnumupdated_cd(:) = .FALSE.
!<<DNdebug
    !pr 110209 zmodmass was never set....
    zmodmass(1:kproma,:,:)=0._dp
    zcdworknum(1:kproma,:)=0._dp
    znact(1:kproma,:,:)=0._dp
    zworkacc(1:kproma,:,:)=0._dp
    zwork(1:kproma,:)=0._dp
    !pr
    !---changes of mass of the in-cloud droplet numbers
    !   (terms with Qnuc in equations (3) and (4) in Hoose et al.) 
    
    !---calculate the total number of activated particles
!davidn
     !let big modes nucleate first, like for ice
    jtaits = sizeclass(iaits)%idt_no
    jtaccs = sizeclass(iaccs)%idt_no
    jtcoas = sizeclass(icoas)%idt_no
!!$    WHERE (pxtp1c(1:kproma,:,jtcoas) > zeps_num)
!!$       znact(1:kproma,:,icoas) = MIN(1._dp,pqnuc(1:kproma,:)/pxtp1c(1:kproma,:,jtcoas))*pxtp1c(1:kproma,:,jtcoas)
!!$    END WHERE
!!$    
!!$    WHERE (pxtp1c(1:kproma,:,jtaccs) > zeps_num)
!!$       znact(1:kproma,:,iaccs) = MIN(1._dp,MAX(0._dp,(pqnuc(1:kproma,:)-pxtp1c(1:kproma,:,jtcoas))        &
!!$                                                     /pxtp1c(1:kproma,:,jtaccs)))*pxtp1c(1:kproma,:,jtaccs)
!!$    END WHERE
!!$
!!$    WHERE (pxtp1c(1:kproma,:,jtaits) > zeps_num)
!!$       znact(1:kproma,:,iaits) = MIN(1._dp,MAX(0._dp,(pqnuc(1:kproma,:) - pxtp1c(1:kproma,:,jtcoas)   &
!!$                                                                        - pxtp1c(1:kproma,:,jtaccs))  &
!!$                                                       / pxtp1c(1:kproma,:,jtaits)))*pxtp1c(1:kproma,:,jtaits)
!!$    END WHERE
    DO jmod=1,nclass!nsol !davidn
       znact(1:kproma,:,jmod) = MAX(0._dp,nact(jmod)%ptr(1:kproma,:,krow))
!davidn
       znactsum(1:kproma,:) = znactsum(1:kproma,:) + znact(1:kproma,:,jmod)
       jtnum = sizeclass(jmod)%idt_no !davidn
    END DO
    ! if no explicit nucleation, but cloud cover still increased, assign everything
    ! to the accumulation mode
    llwork(1:kproma,:) = znactsum(1:kproma,:) > zeps
    znact(1:kproma,:,iaccs) = MERGE( znact(1:kproma,:,iaccs), pqnuc(1:kproma,:), llwork(1:kproma,:) )
    znactsum(1:kproma,:) = MERGE( znactsum(1:kproma,:), pqnuc(1:kproma,:), llwork(1:kproma,:) )

    DO jn=1,naerocomp
       jspec = aerocomp(jn)%spid
       jmod = aerocomp(jn)%iclass
       idt_aer = aerocomp(jn)%idt             ! tracer index, aerosol species jspec in mode jmod
       idt_num  = sizeclass(jmod)%idt_no         ! tracer index, aerosol number for mode jmod


!       IF (sizeclass(jmod)%lsoluble) THEN           ! treat only soluble modes !davidn
          !PR: calculation of total mass zmodmass
          zmodmass(1:kproma,:,jmod) = zmodmass(1:kproma,:,jmod) + pxtp1(1:kproma,:,idt_aer)
          
          !PR: Mail from DOD, 110208
          !Delete the local integer idt_no from this subroutine and replace it with idt_num. 
          !NB Do NOT replace idt_no where it is part of the sizeclass(jmod) structure, i.e. any reference to
          !sizeclass(jmod)%idt_no should be left alone.

          !---calculate (Mx,j/Nj)*Nact,j (x=species, j=mode) 
          WHERE (pxtp1c(1:kproma,:,idt_num) > zeps)
             zwork(1:kproma,:) = ( pxtp1c(1:kproma,:,idt_aer) / pxtp1c(1:kproma,:,idt_num) ) &
                  *  znact(1:kproma,:,jmod)
             zworkacc(1:kproma,:,jspec) = zworkacc(1:kproma,:,jspec) + zwork(1:kproma,:)                          
          END WHERE

!!$          !---loss of interstitial mass (term in Qnuc in eqn(18))
!!$          WHERE (llwork(1:kproma,:))
!!$             pxtp1c(1:kproma,:,idt_aer) = pxtp1c(1:kproma,:,idt_aer) - zwork(1:kproma,:) &
!!$                  * pqnuc(1:kproma,:) / znactsum(1:kproma,:)
!!$          END WHERE
          
!       END IF
    END DO

!>>DN [zmodmass must include aerosol water]
    DO jmod=1,nsol   
       zmodmass(1:kproma,:,jmod) = zmodmass(1:kproma,:,jmod) + pxtp1(1:kproma,:,aerowater(jmod)%idt)          
    END DO
!<<DN

!!$    !---loss of interstitial number 
!!$    DO jmod=1,nclass
!!$       idt_num  = sizeclass(jmod)%idt_no         ! tracer index, aerosol number for mode jmod
!!$       ! davidn: should be calculated for soluble modes only?
!!$       WHERE (llwork(1:kproma,:))
!!$          pxtp1c(1:kproma,:,idt_num) = pxtp1c(1:kproma,:,idt_num) - znact(1:kproma,:,jmod) &
!!$               * pqnuc(1:kproma,:) / znactsum(1:kproma,:)
!!$       END WHERE
!!$    END DO
!!$                  
!!$    !---update in-cloud droplet concentration for each species
!!$    DO jn=1,naerospec
!!$       IF (idt_cd(jn) > 0) THEN
!!$            WHERE (llwork(1:kproma,:))
!!$             !---complete the calculation: (Mx,j/Nj)*(Nact,j/Nact,total)*Qnuc
!!$             pxtp1c(1:kproma,:,idt_cd(jn)) = pxtp1c(1:kproma,:,idt_cd(jn)) + & 
!!$                  zworkacc(1:kproma,:,jn) * pqnuc(1:kproma,:) / znactsum(1:kproma,:) 
!!$          END WHERE
!!$       END IF
!!$    END DO

!>>DN [additional diagnostics]
    DO jk=1,klev
!davidn       dnuc_ncd(1:kproma,krow) = dnuc_ncd(1:kproma,krow) + pqnuc(1:kproma,jk)&
       dnuc_ncd(1:kproma,krow) = dnuc_ncd(1:kproma,krow) + MAX(0._dp,pqnuc(1:kproma,jk))&
                                 *pdpg(1:kproma,jk)*delta_time
    END DO
!<<DN

     !---cloud droplets:
    DO jn=1,naerocomp !davidn: computation is additional to previous computation/done twice
!>>DNdebug
       zmnratio(1:kproma,:)=0._dp
!<<DNdebug
       imod = aerocomp(jn)%iclass
       ispec = aerocomp(jn)%spid
       
       !---set mode size limits
       IF (imod .LE. nsol) THEN
          imod_crdiv = imod
       ELSE
          imod_crdiv = imod-(nclass-nsol)
       END IF

       !---tracer index of aerosol mass
       jtaer = aerocomp(jn)%idt
       
       !---tracer index of aerosol number
       jtnum = sizeclass(imod)%idt_no
       
       !---tracer index of in-cloud species
       jtcl = idt_cd(ispec)

       !   WHERE (zmodmass(1:kproma,:,imod) > zeps_mass)
       WHERE (zmodmass(1:kproma,:,imod) > 0._dp)
          zmnratio(1:kproma,:) = a4piover3 * densaer(imod)%ptr(1:kproma,:,krow)               &
               !*(MAX(0.01*crdiv(imod_crdiv),rwet(imod)%ptr(:,:,krow)))**3       &
!>>DNdebug
!               *(cmr2ram(imod)*rwet(imod)%ptr(:,:,krow))**3       &
               *(cmr2ram(imod)*MAX(0.01_dp*crdiv(imod_crdiv),rwet(imod)%ptr(1:kproma,:,krow)))**3       &
               *MAX(0._dp,pxtp1(1:kproma,:,jtaer)/zmodmass(1:kproma,:,imod))
!               * pxtp1(1:kproma,:,jtaer)/zmodmass(1:kproma,:,imod)
!<<DNdebug
       END WHERE
      
       llwork(1:kproma,:) = (zmnratio(1:kproma,:) <= zeps_massratio)

       WHERE (znactsum(1:kproma,:) > 0._dp) 
          !---for aerosol number, in-droplet
!>>DNdebug
!          zcdworknum(1:kproma,:) = znact(1:kproma,:,imod) * pqnuc(1:kproma,:) / znactsum(1:kproma,:)
          zcdworknum(1:kproma,:) = MAX(0._dp,znact(1:kproma,:,imod)) *&
                                   pqnuc(1:kproma,:) / znactsum(1:kproma,:)
!<<DNdebug
       ELSEWHERE
          zcdworknum(1:kproma,:) = 0._dp
       END WHERE

       !---for aerosol mass, in-droplet 
       zcdworkmass(1:kproma,:) = zcdworknum(1:kproma,:) * zmnratio(1:kproma,:)
!>>DN [correction if number/mass transfer would be too large]
       zcdworknum(1:kproma,:) = MAX(0._dp,MIN(zcdworknum(1:kproma,:),pxtp1c(1:kproma,:,jtnum)))
       zcdworkmass(1:kproma,:) = MAX(0._dp,MIN(zcdworkmass(1:kproma,:),pxtp1c(1:kproma,:,jtaer)))
!<<DN

       !---update the in-droplet tracer
       pxtp1c(1:kproma,:,jtcl) = pxtp1c(1:kproma,:,jtcl)                                              &
                                 +MERGE(zcdworkmass(1:kproma,:),0._dp,llwork(1:kproma,:))

       !---mass loss from interstitial modes due to uptake in droplets
       pxtp1c(1:kproma,:,jtaer) = pxtp1c(1:kproma,:,jtaer)                                             &
                                 -MERGE(zcdworkmass(1:kproma,:),0._dp,llwork(1:kproma,:))
                
       !---diagnostics: sulfate mass uptake in droplets
       IF (ispec == id_so4) THEN
          zworkdiagmass(1:kproma,:) = zcdworkmass(1:kproma,:)*vphysc%rhoam1(1:kproma,:,krow)*delta_time/time_step_len
          quptk_ms4cd(1:kproma,:,krow) = quptk_ms4cd(1:kproma,:,krow) + &
                                         MERGE(zworkdiagmass(1:kproma,:), 0._dp, llwork(1:kproma,:))
       END IF

       !---number loss from interstitial modes due to uptake in droplets
       IF (.NOT. lnumupdated_cd(imod)) THEN
          pxtp1c(1:kproma,:,jtnum)=pxtp1c(1:kproma,:,jtnum)-MERGE(zcdworknum(1:kproma,:),0._dp,llwork(1:kproma,:))

          !---diagnostics
          zworkdiagnum(1:kproma,:) = zcdworknum(1:kproma,:)*vphysc%rhoam1(1:kproma,:,krow)*delta_time/time_step_len
          quptk_incd(1:kproma,:,krow) =quptk_incd(1:kproma,:,krow) + &
                                       MERGE(zworkdiagnum(1:kproma,:),0._dp,llwork(1:kproma,:))
!>>DN [additional diagnostics]
          zworkdiagnum(1:kproma,:) = zcdworknum(1:kproma,:)*pdpg(1:kproma,:)*delta_time   
          DO jk=1,klev
             duptk_incd(1:kproma,krow) = duptk_incd(1:kproma,krow) + &
                                         MERGE(zworkdiagnum(1:kproma,jk),0._dp,llwork(1:kproma,jk))
          END DO
!<<DN
          lnumupdated_cd(imod) = .TRUE.

       END IF

       !--vertically integrated diagnostic per species
       DO jk=1,klev
          duptkcd(ispec)%ptr(1:kproma,krow) = duptkcd(ispec)%ptr(1:kproma,krow)        &
                                              + zcdworkmass(1:kproma,jk) * pdpg(1:kproma,jk) &
                                              * delta_time/time_step_len
       END DO

    END DO
   
  END SUBROUTINE aeroproc_nuc

  !------------------------------------------------------------------------

  SUBROUTINE aeroproc_nuci(kbdim, kproma, klev, krow, pdpg,             &
                          pxtp1, pxtp1c, pqnuci)

    USE mo_kind,            ONLY: dp
    USE mo_tracdef,         ONLY: ntrac
    USE mo_time_control,    ONLY: time_step_len, delta_time
    USE mo_ham_streams,     ONLY: nact, densaer, rwet
    USE mo_ham,             ONLY: aerocomp, aerowater, sizeclass, naerocomp, &
                                  nclass, nsol
    USE mo_ham_m7ctl,       ONLY: sigmaln, &
                                  iaits, iaccs, icoas, crdiv, cmr2ram
    USE mo_ham_species,     ONLY: idt_ic, id_so4
    USE mo_vphysc,          ONLY: vphysc
    USE mo_time_control,    ONLY: current_date, get_time_step
    USE mo_activ,           ONLY: idt_icnc!PR testing
    USE mo_species,         ONLY: speclist
    USE mo_physical_constants,       ONLY: tmelt !davidn
    USE mo_echam_cloud_params,           ONLY: cthomi !davidn
    USE mo_ham_species,     ONLY: id_bc,id_oc,id_ss,id_du !davidn

    !---subroutine argument list
    INTEGER, INTENT(IN)  :: kbdim                        ! geographic block maximum number of locations
    INTEGER, INTENT(IN)  :: kproma                       ! geographic block number of locations on this row
    INTEGER, INTENT(IN)  :: krow                         ! geographic block number
    INTEGER, INTENT(IN)  :: klev                         ! number of vertical levels
    REAL(dp), INTENT(IN) :: pdpg(kbdim,klev)             ! air mass per unit area
    REAL(dp), INTENT(IN) :: pqnuci(kbdim,klev)           ! number of nucleated ice crystals this time step [# m-3] 
    REAL(dp), INTENT(IN) :: pxtp1(kbdim,klev,ntrac)      ! tracer concentrations at t+dt [mmr]
    REAL(dp), INTENT(INOUT) :: pxtp1c(kbdim,klev,ntrac)  ! in-cloud tracer concentrations [kg m-3]

    !PR
    REAL(dp) :: pxtm1(kbdim,klev,ntrac)      ! tracer concentrations at t-dt
    REAL(dp) :: zwork(kbdim,klev)
    REAL(dp) :: zworkacc(kbdim,klev,ntrac)
!    INTEGER  :: jspec,idt_aer,idt_num,idt_no,naerospec,imod,imod_crdiv
    INTEGER  :: jspec,idt_aer,idt_num,idt_no,imod,imod_crdiv
    !PRend
    !---local data    
    REAL(dp) :: zmodmass(kbdim,klev,nclass)                 ! total aerosol mass per mode
    REAL(dp) :: znactsum(kbdim,klev),      &              ! number of activated particles in all modes
                zmnratio(kbdim,klev),      &              ! mass to number conversion 
                znact(kbdim,klev,nclass)                    ! number of activated particles in each mode

    INTEGER :: jl,jk,jt, jmod, jmod_crdiv, lev

    REAL(dp), POINTER :: duptk_p(:,:) 

    REAL(dp) :: zmodefac(nclass)
    INTEGER :: jn, jtcl, jtic, jtaer, jtnum, ispec
    INTEGER :: jtcoas, jtaccs, jtaits
    LOGICAL :: lnumupdated_ic(nclass)

    REAL(dp) :: zfracn(kbdim,klev,nsol)
    REAL(dp) :: zicworknum(kbdim,klev)
    REAL(dp) :: zicworkmass(kbdim,klev)
    REAL(dp) :: zworkdiagnum(kbdim,klev)
    REAL(dp) :: zworkdiagmass(kbdim,klev)
!>>DNdebug
    REAL(dp) :: zxtp1c(kbdim,klev,ntrac)
!<<DNdebug
    LOGICAL :: llwork(kbdim,klev)
    LOGICAL :: lnuci
  
    !-----------------------------------------------------------------------------------
    !---initialisations
    znactsum(1:kproma,:) = 0._dp
    zfracn(1:kproma,:,:) = 0._dp
    zicworknum(1:kproma,:) = 0._dp
    zicworkmass(1:kproma,:) = 0._dp
    zworkdiagnum(1:kproma,:) = 0._dp
    zworkdiagmass(1:kproma,:) = 0._dp
!>>DNdebug
    lnumupdated_ic(:) = .FALSE.
    zxtp1c(1:kproma,:,:)=pxtp1c(1:kproma,:,:)
!<<DNdebug
    
    !pr 110209 zmodmass was never set....
    zmodmass(1:kproma,:,:)=0._dp
    znact(1:kproma,:,:)=0._dp
    zworkacc(1:kproma,:,:)=0._dp
    zwork(1:kproma,:)=0._dp
    !pr
    !---changes of mass of the in-cloud droplet and in-ice crystal numbers
    !   (terms with Qnuci in equations (3) and (4) in Hoose et al.) 
    
    DO jn=1,naerocomp
       jspec = aerocomp(jn)%spid
       jmod = aerocomp(jn)%iclass
       idt_aer = aerocomp(jn)%idt             ! tracer index, aerosol species jspec in mode jmod
       idt_num  = sizeclass(jmod)%idt_no         ! tracer index, aerosol number for mode jmod

       IF (sizeclass(jmod)%lsoluble) THEN           ! treat only soluble modes
          !PR: calculation of total mass zmodmass
          zmodmass(1:kproma,:,jmod) = zmodmass(1:kproma,:,jmod) + pxtp1(1:kproma,:,idt_aer)          
       END IF
    END DO
    
!>>DN [zmodmass must include aerosol water]
    DO jmod=1,nsol       
       zmodmass(1:kproma,:,jmod) = zmodmass(1:kproma,:,jmod) + pxtp1(1:kproma,:,aerowater(jmod)%idt)          
    END DO
!<<DN

    !---ice crystal nucleation: take aerosols first from coarse mode, then accumulation mode and 
    !   finally accumulation mode
    !   calculate fraction of nucleated crystals for aitken, accumulation and coarse soluble modes
    !---tracers for aerosol number in aitken, accumulation and coarse modes
    jtaits = sizeclass(iaits)%idt_no
    jtaccs = sizeclass(iaccs)%idt_no
    jtcoas = sizeclass(icoas)%idt_no
    WHERE (pxtp1c(1:kproma,:,jtcoas) > zeps)
       zfracn(1:kproma,:,icoas) = MIN(1._dp,pqnuci(1:kproma,:)/pxtp1c(1:kproma,:,jtcoas))    
    END WHERE
    
    WHERE (pxtp1c(1:kproma,:,jtaccs) > zeps)
       zfracn(1:kproma,:,iaccs) = MIN(1._dp,MAX(0._dp,(pqnuci(1:kproma,:)-pxtp1c(1:kproma,:,jtcoas))        &
                                                     /pxtp1c(1:kproma,:,jtaccs)) )
    END WHERE

    WHERE (pxtp1c(1:kproma,:,jtaits) > zeps)
       zfracn(1:kproma,:,iaits) = MIN(1._dp,MAX(0._dp,(pqnuci(1:kproma,:) - pxtp1c(1:kproma,:,jtcoas)   &
                                                                        - pxtp1c(1:kproma,:,jtaccs))  &
                                                       / pxtp1c(1:kproma,:,jtaits)) )
    END WHERE
    
!>>DN [additional diagnostics]
    DO jk=1,klev
       dnuc_nic(1:kproma,krow) = dnuc_nic(1:kproma,krow) + pqnuci(1:kproma,jk)&
                                 *pdpg(1:kproma,jk)*delta_time
    END DO
!<<DN

    DO jn=1,naerocomp
!>>DNdebug
       zmnratio(1:kproma,:)=0._dp
!<<DNdebug
       imod = aerocomp(jn)%iclass
       ispec = aerocomp(jn)%spid

       !---set mode size limits
       IF (imod .LE. nsol) THEN
          imod_crdiv = imod
       ELSE
          imod_crdiv = imod-(nclass-nsol)
       END IF

       !---ice crystal processing only for aitken, accumulation and coarse soluble modes
       lnuci = (imod == iaits .OR. imod == iaccs .OR. imod == icoas)

       !---tracer index of aerosol mass
       jtaer = aerocomp(jn)%idt
       
       !---tracer index of aerosol number
       jtnum = sizeclass(imod)%idt_no
            
       !---tracer index of in-ice crystal species
       jtic = idt_ic(ispec)

       !   WHERE (zmodmass(1:kproma,:,imod) > zeps_mass)
       WHERE (zmodmass(1:kproma,:,imod) > 0._dp)
          zmnratio(1:kproma,:) = a4piover3 * densaer(imod)%ptr(1:kproma,:,krow)               &
               !*(MAX(0.01*crdiv(imod_crdiv),rwet(imod)%ptr(:,:,krow)))**3       &
!>>DNdebug
!               *(cmr2ram(imod)*rwet(imod)%ptr(:,:,krow))**3       &
               *(cmr2ram(imod)*MAX(0.01_dp*crdiv(imod_crdiv),rwet(imod)%ptr(1:kproma,:,krow)))**3       &
               *MAX(0._dp,pxtp1(1:kproma,:,jtaer)/zmodmass(1:kproma,:,imod))
!               * pxtp1(1:kproma,:,jtaer)/zmodmass(1:kproma,:,imod)
!<<DNdebug
       END WHERE
     
       llwork(1:kproma,:) = (zmnratio(1:kproma,:) <= zeps_massratio)
       IF (lnuci) THEN
       !---aerosol number, in-crystal
!>>DNdebug
!          zicworknum(1:kproma,:) = zfracn(1:kproma,:,imod) * pxtp1c(1:kproma,:,jtnum)
          zicworknum(1:kproma,:) = zfracn(1:kproma,:,imod) * zxtp1c(1:kproma,:,jtnum)
!<<DNdebug
          !---aerosol mass, in-crystal
          zicworkmass(1:kproma,:) = zicworknum(1:kproma,:) * zmnratio(1:kproma,:)
       
!>>DN [correction if number/mass transfer would be too large]
          zicworknum(1:kproma,:) = MAX(0._dp,MIN(zicworknum(1:kproma,:),zxtp1c(1:kproma,:,jtnum)))
          zicworkmass(1:kproma,:) = MAX(0._dp,MIN(zicworkmass(1:kproma,:),pxtp1c(1:kproma,:,jtaer)))
!<<DN

       !---update the in-crystal tracer
          pxtp1c(1:kproma,:,jtic) = pxtp1c(1:kproma,:,jtic)                                              &
                                           +MERGE(zicworkmass(1:kproma,:),0._dp,llwork(1:kproma,:))
       !---mass loss from interstitial modes due to uptake in crystals
          pxtp1c(1:kproma,:,jtaer) = pxtp1c(1:kproma,:,jtaer)                                             &
                                            -MERGE(zicworkmass(1:kproma,:),0._dp,llwork(1:kproma,:))
          
       !---diagnostics: sulfate mass uptake in crystals
          IF (ispec == id_so4) THEN
             zworkdiagmass(1:kproma,:) = zicworkmass(1:kproma,:)*vphysc%rhoam1(1:kproma,:,krow)*delta_time/time_step_len
             quptk_ms4ic(1:kproma,:,krow) = quptk_ms4ic(1:kproma,:,krow) + &
                                            MERGE(zworkdiagmass(1:kproma,:), 0._dp, llwork(1:kproma,:))
          END IF

       !---loss of interstitial aerosol number due to uptake in crystals
          IF (.NOT. lnumupdated_ic(imod)) THEN
             pxtp1c(1:kproma,:,jtnum) = pxtp1c(1:kproma,:,jtnum) - zicworknum(1:kproma,:)
!             pxtp1c(1:kproma,:,jtnum) = MAX(0._dp,pxtp1c(1:kproma,:,jtnum))!davidn
             !---diagnostics
             zworkdiagnum(1:kproma,:) = zicworknum(1:kproma,:)*vphysc%rhoam1(1:kproma,:,krow)*delta_time/time_step_len
             quptk_inic(1:kproma,:,krow) = quptk_inic(1:kproma,:,krow) + &
                                           MERGE(zworkdiagnum(1:kproma,:), 0._dp, llwork(1:kproma,:))
!>>DN [additional diagnostics]
             zworkdiagnum(1:kproma,:) = zicworknum(1:kproma,:)*pdpg(1:kproma,:)*delta_time  
             DO jk=1,klev
                duptk_inic(1:kproma,krow) = duptk_inic(1:kproma,krow) + &
                                            MERGE(zworkdiagnum(1:kproma,jk),0._dp,llwork(1:kproma,jk))
             END DO
!<<DN
             lnumupdated_ic(imod) = .TRUE.
          END IF

       !--vertically integrated diagnostic per species
          DO jk=1,klev
             duptkic(ispec)%ptr(1:kproma,krow) = duptkic(ispec)%ptr(1:kproma,krow)        &
                                                 + zicworkmass(1:kproma,jk) * pdpg(1:kproma,jk) &
                                                 * delta_time/time_step_len
          END DO

       END IF

    END DO
    
  END SUBROUTINE aeroproc_nuci

  !------------------------------------------------------------------------
  
  SUBROUTINE aeroproc_freeze(kbdim, kproma, klev, krow, pdpg, pxtp1c, pnactfrz, pnfrz)

    USE mo_kind,            ONLY: dp
!>>DN
!    USE mo_aero_species,    ONLY: naerospec, idt_cd, idt_ic, id_so4
    USE mo_ham_species,     ONLY: idt_cd, idt_ic, id_so4, id_wat
    USE mo_species,         ONLY: naerospec, aero_idx
!<<DN
    USE mo_time_control,    ONLY: time_step_len, delta_time
    USE mo_tracdef,         ONLY: ntrac
    USE mo_activ,           ONLY: idt_cdnc, idt_icnc
    USE mo_vphysc,          ONLY: vphysc
    
    INTEGER, INTENT(IN)  :: kbdim                            ! geographic block maximum number of locations
    INTEGER, INTENT(IN)  :: kproma                           ! geographic block number of locations on this row
    INTEGER, INTENT(IN)  :: krow                             ! geographic block number
    INTEGER, INTENT(IN)  :: klev                             ! number of vertical levels
    REAL(dp), INTENT(IN) :: pdpg(kbdim,klev)                 ! air mass per unit area
    REAL(dp), INTENT(IN) :: pnfrz(kbdim,klev)                ! cloud droplet freezing rate ?
    REAL(dp), INTENT(IN) :: pnactfrz(kbdim,klev)             ! 
    REAL(dp), INTENT(INOUT) :: pxtp1c(kbdim,klev,ntrac)      ! tracer concentration in the cloudy part

    !---local data
    REAL(dp) :: zaerodm(kbdim,klev)
    
    INTEGER :: jk, jn, jtcl, jtic 

    REAL(dp), POINTER :: dfrz_p(:,:)
    !---executable procedure

    zaerodm(1:kproma,:) = 0._dp



    DO jn=aero_idx(1),aero_idx(naerospec)
!>>DN [no water tracer in cloud particles]
       IF (jn==id_wat) CYCLE
!<<DN
       jtcl = idt_cd(jn)
       jtic = idt_ic(jn)
       dfrz_p     => dfrz(jn)%ptr
       
       WHERE (pnfrz(1:kproma,:) > zeps .AND. pnactfrz(1:kproma,:) > zeps)
!>>DNdebug
!          zaerodm(1:kproma,:) = pxtp1c(1:kproma,:,jtcl) / pnactfrz(1:kproma,:) * pnfrz(1:kproma,:)
          zaerodm(1:kproma,:) = MAX(0._dp,MIN(pxtp1c(1:kproma,:,jtcl),&
               pxtp1c(1:kproma,:,jtcl) / pnactfrz(1:kproma,:) * pnfrz(1:kproma,:)))
!<<DNdebug
          !---reduce mass of in-cloud droplet aerosol
          pxtp1c(1:kproma,:,jtcl) = pxtp1c(1:kproma,:,jtcl) - zaerodm(1:kproma,:)

          !---add to mass of in-ice crystal aerosol
          pxtp1c(1:kproma,:,jtic) = pxtp1c(1:kproma,:,jtic) + zaerodm(1:kproma,:)

       END WHERE

       !---vertically integrated diagnostics per species
       DO jk=1,klev
          dfrz_p(1:kproma,krow) = dfrz_p(1:kproma,krow)  + zaerodm(1:kproma,jk)       &
                                  * pdpg(1:kproma,jk)*delta_time/time_step_len
       END DO

       !---3-dimensional diagnostics for sulfate mass
       IF (jn == id_so4)  THEN
          qfrz_ms4cd(1:kproma,:,krow) = qfrz_ms4cd(1:kproma,:,krow)                          &
                                        + zaerodm(1:kproma,:)  *  vphysc%rhoam1(1:kproma,:,krow)         &
                                        *delta_time/time_step_len
       END IF

    END DO

    !---3-dimensional diagnostics (number)
    WHERE (pnfrz(1:kproma,:) > zeps .AND. pnactfrz(1:kproma,:) > zeps)
       qfrz_incd(1:kproma,:,krow) =  qfrz_incd(1:kproma,:,krow)                &
                                     + pnfrz(1:kproma,:) * vphysc%rhoam1(1:kproma,:,krow)      &
                                     *delta_time/time_step_len   
       
    END WHERE

  END SUBROUTINE aeroproc_freeze

  !------------------------------------------------------------------------

  SUBROUTINE aeroproc_melt(kbdim, kproma, klev, krow, pdpg, pxtp1c, pnactmelt, pnmelt)

    USE mo_kind,            ONLY: dp
    USE mo_tracdef,         ONLY: ntrac
    USE mo_activ,           ONLY: idt_icnc
!>>DN
!    USE mo_aero_species,    ONLY: naerospec, idt_cd, idt_ic, id_so4
    USE mo_ham_species,     ONLY: idt_cd, idt_ic, id_so4, id_wat
    USE mo_species,         ONLY: naerospec, aero_idx
!<<DN
    USE mo_time_control,    ONLY: delta_time, time_step_len
    USE mo_vphysc,          ONLY: vphysc

    INTEGER, INTENT(IN)  :: kbdim                            ! geographic block maximum number of locations
    INTEGER, INTENT(IN)  :: kproma                           ! geographic block number of locations on this row
    INTEGER, INTENT(IN)  :: krow                             ! geographic block number
    INTEGER, INTENT(IN)  :: klev                             ! number of vertical levels
    REAL(dp), INTENT(IN) :: pdpg(kbdim,klev)                 ! air mass per unit area
    REAL(dp), INTENT(IN) :: pnmelt(kbdim,klev)               ! ice crystal melting rate ?
    REAL(dp), INTENT(IN) :: pnactmelt(kbdim,klev)            ! 
    REAL(dp), INTENT(INOUT) :: pxtp1c(kbdim,klev,ntrac)      ! tracer concentration in the cloudy part

    !---local data
    REAL(dp) :: zaerodm(kbdim,klev)

    INTEGER :: jn, jk, jtic, jtcl

    REAL(dp), POINTER :: dmelt_p(:,:)

    !---executable procedure

    zaerodm(1:kproma,:) = 0._dp

    !---in-crystal mass:
    DO jn=aero_idx(1),aero_idx(naerospec)
!>>DN [no water tracer in cloud particles]
       IF (jn==id_wat) CYCLE
!<<DN
       jtcl = idt_cd(jn)
       jtic = idt_ic(jn)
       dmelt_p     => dmelt(jn)%ptr
       
        WHERE (pnmelt(1:kproma,:) > zeps .AND. pnactmelt(1:kproma,:) > zeps) 
!>>DNdebug
!          zaerodm(1:kproma,:) = pxtp1c(1:kproma,:,jtic) / pnactmelt(1:kproma,:) * pnmelt(1:kproma,:)
          zaerodm(1:kproma,:) = MAX(0._dp,MIN(pxtp1c(1:kproma,:,jtic),&
               pxtp1c(1:kproma,:,jtic) / pnactmelt(1:kproma,:) * pnmelt(1:kproma,:)))
!<<DNdebug

          !---reduce mass of in-ice crystal aerosol
          pxtp1c(1:kproma,:,jtic) = pxtp1c(1:kproma,:,jtic) - zaerodm(1:kproma,:)

          !---add to mass of in-droplet aerosol
          pxtp1c(1:kproma,:,jtcl) = pxtp1c(1:kproma,:,jtcl) + zaerodm(1:kproma,:)

       END WHERE

       !---vertically integrated diagnostics per species
       DO jk=1,klev
          dmelt_p(1:kproma,krow) = dmelt_p(1:kproma,krow)  + zaerodm(1:kproma,jk)       &
                                   * pdpg(1:kproma,jk)*delta_time/time_step_len
       END DO

       !---3-dimensional diagnostics for sulfate mass
       IF (jn == id_so4)  THEN
          qmelt_ms4cd(1:kproma,:,krow) = qmelt_ms4cd(1:kproma,:,krow)                          &
                                         + zaerodm(1:kproma,:)  *  vphysc%rhoam1(1:kproma,:,krow)         &
                                         *delta_time/time_step_len
       END IF
    END DO

    !---3-dimensional diagnostics (number)
    WHERE (pnmelt(1:kproma,:) > zeps .AND. pnactmelt(1:kproma,:) > zeps) 

       qmelt_incd(1:kproma,:,krow) = qmelt_incd(1:kproma,:,krow)                &
                                     + pnmelt(1:kproma,:) * vphysc%rhoam1(1:kproma,:,krow)      &
                                     *delta_time/time_step_len
        
    END WHERE

  END SUBROUTINE aeroproc_melt
  
  !------------------------------------------------------------------------

  SUBROUTINE aeroproc_evap(kbdim, kproma, klev, krow, pdpg, pxtp1c, pxtp10, pnactevap, pnevap)

    USE mo_kind,            ONLY: dp
    USE mo_tracdef,         ONLY: ntrac
!>>DN
!    USE mo_aero_m7,         ONLY: sizeclass, crdiv, iaits, iaccs, icoas, naermod, nsol
!    USE mo_aero_species,    ONLY: naerospec, aerospec, idt_cd, idt_ic, id_so4
    USE mo_ham,             ONLY: sizeclass, naerocomp, aerocomp, nsol, nclass
    USE mo_ham_m7ctl,       ONLY: crdiv, inucs, iaits, iaccs, icoas
    USE mo_ham_species,     ONLY: idt_cd, idt_ic, id_so4, id_ss, id_du, id_wat
    USE mo_species,         ONLY: naerospec, speclist, aero_idx
    USE mo_ham_streams,     ONLY: rdry
!<<DN
    USE mo_activ,           ONLY: idt_cdnc
    USE mo_time_control,    ONLY: delta_time, time_step_len
    USE mo_vphysc,          ONLY: vphysc
!    USE mo_outctl,          ONLY: lproc_detail

    INTEGER, INTENT(IN)  :: kbdim                            ! geographic block maximum number of locations
    INTEGER, INTENT(IN)  :: kproma                           ! geographic block number of locations on this row
    INTEGER, INTENT(IN)  :: krow                             ! geographic block number
    INTEGER, INTENT(IN)  :: klev                             ! number of vertical levels
    REAL(dp), INTENT(IN) :: pdpg(kbdim,klev)                 ! air mass per unit area
!>>DN
!    REAL(dp), INTENT(IN) :: pnevap(kbdim,klev)               ! cloud droplet evaporation rate
!    REAL(dp), INTENT(IN) :: pnactevap(kbdim,klev)            ! 
    REAL(dp), INTENT(INOUT) :: pnevap(kbdim,klev)               ! cloud droplet evaporation rate
    REAL(dp), INTENT(INOUT) :: pnactevap(kbdim,klev)            ! 
!<<DN
    REAL(dp), INTENT(INOUT) :: pxtp1c(kbdim,klev,ntrac)      ! tracer concentration in the cloudy part
    REAL(dp), INTENT(INOUT) :: pxtp10(kbdim,klev,ntrac)      ! tracer concentration in the cloud-free part

    !---local data
    REAL(dp) :: zvol(kbdim,klev)                             ! volume of dry aerosol in cloud droplet
    REAL(dp) :: zaeromtot(kbdim,klev)                        ! total mass of dry aerosol in cloud droplet
    REAL(dp) :: zaerodm(kbdim,klev)                          ! change of mass due to evaporation
!>>DN
    REAL(dp) :: zaerodn(kbdim,klev,nclass)                          ! change of number due to evaporation
!<<DN

!>>DN
!    REAL(dp) :: crdivm(nsol)                                 ! threshold radii for M7 modes in m
    REAL(dp) :: crdivm(nsol+1)                                 ! threshold radii for M7 modes in m
!<<DN
    REAL(dp) :: zram2cmr_cd                                  ! conversion of radius of average aerosol mass -> count median radius
    REAL(dp) :: zmasschange(kbdim,klev)
    REAL(dp) :: zlog4piover3
    INTEGER :: jn, jl,jk,jt, jtcl, jtnum
    INTEGER :: ispec, imod
!>>DN
    LOGICAL :: lnumupdated_cd(kbdim,klev,nclass),ldnupdate(kbdim,klev) 
    REAL(dp) :: zrdry_incd(kbdim,klev,aero_idx(1):aero_idx(naerospec))               ! radius of released particles corrected for SS and DU
!!<<DN

    REAL(dp), POINTER :: dloss1_p(:,:)
    REAL(dp), POINTER :: dloss2_p(:,:)
    REAL(dp), POINTER :: dloss3_p(:,:)
    REAL(dp), POINTER :: devap_p(:,:)

    !---executable procedure
    zlog4piover3 = LOG(a4piover3)
    zvol(1:kproma,:) = 0._dp
    zaeromtot(1:kproma,:) = 0._dp
    zaerodm(1:kproma,:) = 0._dp
    crdivm(1:nsol) = 0.01_dp*crdiv(1:nsol)
!>>DN
    crdivm(nsol+1)=50.0E-6_dp
    lnumupdated_cd(1:kproma,:,:) = .FALSE.
    ldnupdate(1:kproma,:) = .FALSE.
    zaerodn(1:kproma,:,:) = 0._dp
    pnactevap(1:kproma,:) = MAX(pnevap(1:kproma,:),pnactevap(1:kproma,:))
    pnactevap(1:kproma,:) = MAX(0._dp,pnactevap(1:kproma,:))
!<<DN

    !---conversion factor radius of average mass -> count median radius
!>>DN
!    zram2cmr_cd = 1._dp / ( EXP(1.5_dp*(LOG(sigma_cd))**2) )
    zram2cmr_cd = 1._dp / ( EXP(1.5_dp*(LOG(1.59_dp))**2) )
!<<DN

    !---cutoff limits for 0.5um and Aitken mode minimum (count median radius) converted to radius of average mass
    !zcutoff_0.5 = 0.5_dp * EXP(1.5_dp*(LOG(sigma(icoas)))**2)
    !zcutoff_aitk = (crdiv(iaits)*0.01_dp) * EXP(1.5_dp*(LOG(sigma(iaits)))**2)

    !---total in-cloud droplet aerosol mass
    DO jn = aero_idx(1),aero_idx(naerospec)   
!>>DN [no water tracer in cloud particles]
       IF (jn==id_wat) CYCLE
!<<DN
       jtcl = idt_cd(jn)
       zaeromtot(1:kproma,:) = zaeromtot(1:kproma,:) + MAX(0._dp,pxtp1c(1:kproma,:,jtcl))          
    END DO
       
   !---sum (mass/density) for all species in cloud [ result in m**3/kg(air) ]
!>>DN [moved before change in droplet mass]
   DO jn=aero_idx(1),aero_idx(naerospec) 
      
!>>DN [no water tracer in cloud particles]
      IF (jn==id_wat) CYCLE
!<<DN
      jtcl = idt_cd(jn)
!!>>DN [SO4 in kg(S) kg-1]
!       IF (jn==id_so4) THEN
!          zvol(1:kproma,:) = zvol(1:kproma,:) + pxtp1c(1:kproma,:,jtcl)*astoso4/speclist(jn)%density
!       ELSE
          zvol(1:kproma,:) = zvol(1:kproma,:) + MAX(0._dp,pxtp1c(1:kproma,:,jtcl))/speclist(jn)%density
!       END IF
!!<<DN
   END DO
!<<DN

  !---calculate radius of released particle
   WHERE (pnactevap(1:kproma,:) > zeps_num .AND. zvol(1:kproma,:) > 0._dp)

      !---divide by droplet number concentration [ result in m**3 / drop ]
!      zvol(1:kproma,:) = zvol(1:kproma,:) / pnactevap(1:kproma,:)
     
    !---convert to radius of average mass
      rdry_incd(1:kproma,:,krow) =  zram2cmr_cd * ((zvol(1:kproma,:)/ pnactevap(1:kproma,:)/a4piover3)**(1._dp/3._dp))
    ELSEWHERE
       rdry_incd(1:kproma,:,krow) = 0._dp
    END WHERE
!>>DN [additional diagnostics]
    DO jk=1,klev
!       IF (lproc_detail) THEN
          devap_ncd(1:kproma,krow) = devap_ncd(1:kproma,krow)-pnevap(1:kproma,jk)&
                               *pdpg(1:kproma,jk)*delta_time                                    
!       END IF
    END DO
!<<DN

!>>DN [corrections for too small/too large released particles]
    !use coarse mode sigma, if rdry_incd larger than 0.5 um
    WHERE (rdry_incd(1:kproma,:,krow) > 0.5E-6_dp .AND. pnactevap(1:kproma,:) > zeps_num)  

       rdry_incd(1:kproma,:,krow) =  1._dp / ( EXP(1.5_dp*(LOG(sigma_cd))**2) ) &
                                     * ((zvol(1:kproma,:)/ pnactevap(1:kproma,:)/a4piover3)**(1._dp/3._dp))
    END WHERE
   
    !prescribe rdry=100nm where rdry is too small
    WHERE (rdry_incd(1:kproma,:,krow) < crdivm(iaits))! .AND. &
!         pnevap(1:kproma,:) * (rdry_incd(1:kproma,:,krow)/100.E-9_dp)**3._dp <= 1.E6_dp .AND. &
!         pnactevap(1:kproma,:) > zeps_num .AND. zvol(1:kproma,:) > 0._dp) 
!       pnevap(1:kproma,:)=pnactevap(1:kproma,:) * (rdry_incd(1:kproma,:,krow)/100.E-9_dp)**3._dp
       rdry_incd(1:kproma,:,krow)=100.E-9_dp
       pnevap(1:kproma,:)=zvol(1:kproma,:)/a4piover3/&
            (rdry_incd(1:kproma,:,krow)*EXP(1.5_dp*(LOG(1.59_dp))**2))**3
       pnactevap(1:kproma,:)=pnevap(1:kproma,:)
    END WHERE

   !prescribe rdry=10um here rdry is too large or where previous calculation yielded very high numbers
!    WHERE ((rdry_incd(1:kproma,:,krow) < crdivm(iaits) .AND. &
    WHERE ((rdry_incd(1:kproma,:,krow) == 100.E-9_dp .AND. &
         pnevap(1:kproma,:) > 1.E6_dp) .OR. &
!         pnactevap(1:kproma,:) > zeps_num .AND. zvol(1:kproma,:) > 0._dp) .OR. &
         rdry_incd(1:kproma,:,krow) > 50.E-6_dp) 
!       pnevap(1:kproma,:)=pnactevap(1:kproma,:) * (rdry_incd(1:kproma,:,krow)/10.E-6_dp)**3._dp
       rdry_incd(1:kproma,:,krow)=10.E-6_dp
       pnevap(1:kproma,:)=zvol(1:kproma,:)/a4piover3/&
            (rdry_incd(1:kproma,:,krow)*EXP(1.5_dp*(LOG(sigma_cd))**2))**3
       pnactevap(1:kproma,:)=pnevap(1:kproma,:)
    END WHERE
 
    DO jn = 1,nsol
       ![only calculate zaerodn here; pnevap <=pnactevap]
       WHERE (pnevap(1:kproma,:) > zeps .AND. pnactevap(1:kproma,:) > zeps_num .AND. &
            zaeromtot(1:kproma,:) > 0._dp) 
          zaerodn(1:kproma,:,jn) = MIN(pnevap(1:kproma,:),pnactevap(1:kproma,:))
       END WHERE
    END DO

    DO jn = aero_idx(1),aero_idx(naerospec) 
       zrdry_incd(1:kproma,:,jn)=rdry_incd(1:kproma,:,krow)
       !sea salt and dust: prescribe rdry=rdry_AS or 50nm where rdry is too small
       WHERE ((jn==id_ss .OR. jn==id_du) .AND. &
            rdry_incd(1:kproma,:,krow) < crdivm(iaccs) .AND. &
            pnevap(1:kproma,:) > zeps .AND. pnactevap(1:kproma,:) > zeps_num .AND. &
            zaeromtot(1:kproma,:) > 0._dp)
          zrdry_incd(1:kproma,:,jn)= MAX(rdry(iaccs)%ptr(1:kproma,:,krow),crdivm(iaccs))
          zaerodn(1:kproma,:,iaccs) = MIN(pnevap(1:kproma,:),pnactevap(1:kproma,:)) * &
               (rdry_incd(1:kproma,:,krow)/zrdry_incd(1:kproma,:,jn))**3._dp
       END WHERE
       WHERE ((jn==id_ss .OR. jn==id_du) .AND. &
            rdry_incd(1:kproma,:,krow) < crdivm(iaccs) .AND. &
            pnevap(1:kproma,:) > zeps .AND. pnactevap(1:kproma,:) > zeps_num .AND. & 
            zaeromtot(1:kproma,:) > 0._dp .AND. .NOT. ldnupdate(1:kproma,:))
          zaerodn(1:kproma,:,iaits) = zaerodn(1:kproma,:,iaits)-zaerodn(1:kproma,:,iaccs)
          ldnupdate(1:kproma,:)=.TRUE.
       END WHERE

    END DO
    
    ![only calculate zaerodm here; pnevap <=pnactevap]
    WHERE (pnevap(1:kproma,:) > zeps .AND. pnactevap(1:kproma,:) > zeps_num .AND. &
         zaeromtot(1:kproma,:) > 0._dp) 
       zaerodm(1:kproma,:) = MIN(zaeromtot(1:kproma,:),zaeromtot(1:kproma,:) / &
            pnactevap(1:kproma,:) * pnevap(1:kproma,:))
    END WHERE
!<<DN

    !---attribute number and mass change to correct mode

    !---loop over all aerosol species in all modes
    DO jn=1,naerocomp
       ispec = aerocomp(jn)%spid
       imod = aerocomp(jn)%iclass
       
       !---get tracer index of aerosol number
       jtnum = sizeclass(imod)%idt_no
       
       !---get tracer index to aerosol 
       jt = aerocomp(jn)%idt
       
       DO jk=1,klev
          DO jl=1,kproma
             
             IF ( (imod == iaits .AND. zrdry_incd(jl,jk,ispec) >= crdivm(iaits) .AND. &
                  zrdry_incd(jl,jk,ispec) < crdivm(iaccs)) .OR. &
                  (imod == iaccs .AND. zrdry_incd(jl,jk,ispec) >= crdivm(iaccs) .AND. &
                  zrdry_incd(jl,jk,ispec) < crdivm(icoas)) .OR. &
                  (imod == icoas .AND. zrdry_incd(jl,jk,ispec) >= crdivm(icoas) .AND. &
                  zrdry_incd(jl,jk,ispec) < crdivm(nsol+1)) ) THEN

                   !---increase aerosol number in cloud-free part
!>>DN [update number concentration only once per mode; additional diagnostics]
                IF (.NOT.lnumupdated_cd(jl,jk,imod)) THEN
                   pxtp10(jl,jk,jtnum) = pxtp10(jl,jk,jtnum) + zaerodn(jl,jk,imod)
                   lnumupdated_cd(jl,jk,imod) = .TRUE.
!                   IF (lproc_detail) THEN
                      qevap_incd(jl,jk,krow) = qevap_incd(jl,jk,krow)             &
                                               -zaerodn(jl,jk,imod)*vphysc%rhoam1(jl,jk,krow)    &
                                               *delta_time/time_step_len
                      devap_incd(jl,krow) = devap_incd(jl,krow)-zaerodn(jl,jk,imod)&
                                            *pdpg(jl,jk)*delta_time            
!                   END IF
                END IF
!<<DN

                   !---increase mass in the cloud-free part
!>>DN [index must refer to in-coud species]
!                zmasschange(jl,jk) = pxtp1c(jl,jk,jt)/zaeromtot(jl,jk)*zaerodm(jl,jk)
                IF (zaeromtot(jl,jk)>0._dp.AND.pxtp1c(jl,jk,idt_cd(ispec))>0._dp) THEN
                   zmasschange(jl,jk) = MIN(pxtp1c(jl,jk,idt_cd(ispec)), &
                        pxtp1c(jl,jk,idt_cd(ispec))/zaeromtot(jl,jk)*zaerodm(jl,jk))
                ELSE
                   zmasschange(jl,jk) = 0._dp
                END IF
!<<DN
                pxtp10(jl,jk,jt) = pxtp10(jl,jk,jt) + zmasschange(jl,jk)
                  !---diagnostics

!                IF (lproc_detail) THEN
!>>DNdebug
                   IF (ispec == id_so4) THEN
!                   IF (jn == id_so4) THEN
!<<DNdebug
                      qevap_ms4cd(jl,jk,krow) = qevap_ms4cd(jl,jk,krow)             &
                                                -zmasschange(jl,jk)*vphysc%rhoam1(jl,jk,krow)    &
                                                *delta_time/time_step_len
                   END IF
!                END IF

                   !---TODO (maybe): re-introduce diagnostic of lost SS and DU mass 
             END IF

          END DO
       END DO

!!$       !---aitken mode
!!$       IF (imod == iaits) THEN
!!$          DO jk=1,klev
!!$             DO jl=1,kproma
!!$                
!!$                IF ( rdry_incd(jl,jk,krow) >= crdivm(iaits) .AND. &
!!$                     rdry_incd(jl,jk,krow) < crdivm(iaccs) ) THEN
!!$
!!$                   !---increase aerosol number in cloud-free part
!!$!>>DN [update number concentration only once per mode]
!!$                   IF (.NOT.lnumupdated_cd(imod)) THEN
!!$                      pxtp10(jl,jk,jtnum) = pxtp10(jl,jk,jtnum) + pnevap(jl,jk)
!!$                      lnumupdated_cd(imod) = .TRUE.
!!$                      IF (lproc_detail) THEN
!!$                         qevap_incd(jl,jk,krow) = qevap_incd(jl,jk,krow)             &
!!$                                                  -pnevap(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                  *delta_time/time_step_len
!!$                                                  /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$!<<DN
!!$
!!$                   !---increase mass in the cloud-free part
!!$!>>DN [index must refer to in-coud species]
!!$!                   zmasschange(jl,jk) = pxtp1c(jl,jk,jt)/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$                   zmasschange(jl,jk) = pxtp1c(jl,jk,idt_cd(ispec))/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$!<<DN
!!$                   pxtp10(jl,jk,jt) = pxtp10(jl,jk,jt) + zmasschange(jl,jk)      
!!$
!!$                   !---diagnostics
!!$
!!$                   IF (lproc_detail) THEN
!!$!>>DNdebug
!!$                      IF (ispec == id_so4) THEN
!!$!                      IF (jn == id_so4) THEN
!!$!<<DNdebug
!!$                         qevap_ms4cd(jl,jk,krow) = qevap_ms4cd(jl,jk,krow)             &
!!$                                                   - zmasschange(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                   *delta_time/time_step_len
!!$                                                   /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$
!!$                   !---TODO (maybe): re-introduce diagnostic of lost SS and DU mass 
!!$                END IF
!!$             END DO
!!$          END DO
!!$       END IF
!!$             
!!$       !---accumulation mode
!!$       IF (imod == iaccs) THEN
!!$          DO jk=1,klev
!!$             DO jl=1,kproma
!!$                
!!$                IF ( rdry_incd(jl,jk,krow) >= crdivm(iaccs) .AND. &
!!$                     rdry_incd(jl,jk,krow) < crdivm(icoas) ) THEN
!!$
!!$                   !---increase aerosol number in cloud-free part
!!$!>>DN [update number concentration only once per mode]
!!$                   IF (.NOT.lnumupdated_cd(imod)) THEN
!!$                      pxtp10(jl,jk,jtnum) = pxtp10(jl,jk,jtnum) + pnevap(jl,jk)
!!$                      lnumupdated_cd(imod) = .TRUE.  
!!$                      IF (lproc_detail) THEN
!!$                         qevap_incd(jl,jk,krow) = qevap_incd(jl,jk,krow)             &
!!$                                                  -pnevap(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                  *delta_time/time_step_len
!!$                                                  /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$!<<DN
!!$
!!$                   !---increase mass in the cloud-free part
!!$!>>DN [index must refer to in-coud species]
!!$!                   zmasschange(jl,jk) = pxtp1c(jl,jk,jt)/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$                   zmasschange(jl,jk) = pxtp1c(jl,jk,idt_cd(ispec))/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$!<<DN
!!$                   pxtp10(jl,jk,jt) = pxtp10(jl,jk,jt) + zmasschange(jl,jk)      
!!$
!!$                   !---diagnostics
!!$                   IF (lproc_detail) THEN
!!$!>>DNdebug
!!$                      IF (ispec == id_so4) THEN
!!$!                      IF (jn == id_so4) THEN
!!$!<<DNdebug
!!$                         qevap_ms4cd(jl,jk,krow) = qevap_ms4cd(jl,jk,krow)             &
!!$                                                   - zmasschange(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                   *delta_time/time_step_len
!!$                                                   /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$
!!$                   !---TODO (maybe): re-introduce diagnostic of lost SS and DU mass 
!!$                END IF
!!$             END DO
!!$          END DO
!!$       END IF
!!$             
!!$       !---coarse mode   
!!$       IF (imod == icoas) THEN
!!$          DO jk=1,klev
!!$             DO jl=1,kproma
!!$                
!!$                IF ( rdry_incd(jl,jk,krow) > crdivm(icoas)) THEN
!!$
!!$                   !---increase aerosol number in cloud-free part
!!$!>>DN [update number concentration only once per mode]
!!$                   IF (.NOT.lnumupdated_cd(imod)) THEN
!!$                      pxtp10(jl,jk,jtnum) = pxtp10(jl,jk,jtnum) + pnevap(jl,jk)
!!$                      lnumupdated_cd(imod) = .TRUE.
!!$                      IF (lproc_detail) THEN
!!$                         qevap_incd(jl,jk,krow) = qevap_incd(jl,jk,krow)             &
!!$                                                  -pnevap(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                  *delta_time/time_step_len
!!$                                                  /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$!<<DN
!!$
!!$                   !---increase mass in the cloud-free part
!!$!>>DN [index must refer to in-coud species]
!!$!                   zmasschange(jl,jk) = pxtp1c(jl,jk,jt)/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$                   zmasschange(jl,jk) = pxtp1c(jl,jk,idt_cd(ispec))/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$!<<DN
!!$                   pxtp10(jl,jk,jt) = pxtp10(jl,jk,jt) + zmasschange(jl,jk)      
!!$
!!$                   !---diagnostics
!!$                   IF (lproc_detail) THEN
!!$!>>DNdebug
!!$                      IF (ispec == id_so4) THEN
!!$!                      IF (jn == id_so4) THEN
!!$!<<DNdebug
!!$                         qevap_ms4cd(jl,jk,krow) = qevap_ms4cd(jl,jk,krow)             &
!!$                                                   - zmasschange(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                   *delta_time/time_step_len
!!$                                                   /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$
!!$                   !---TODO (maybe): re-introduce diagnostic of lost SS and DU mass 
!!$                END IF
!!$             END DO
!!$          END DO
!!$       END IF
 
   END DO

    !---TODO: diagnostic of any cases where rdry_incd < aitken mode min, or very large values (> 50um)

!>>DN [diagnostic of any cases where rdry_incd < aitken mode min, or for very large values (> 50um)]
   DO jk=1,klev
      DO jl=1,kproma
         
         IF ( zaeromtot(jl,jk) > 0._dp .AND. &
              pnevap(jl,jk) > zeps .AND. pnactevap(jl,jk)  > zeps_num) THEN
            DO jn=aero_idx(1),aero_idx(naerospec) 
               IF (jn==id_wat) CYCLE
               IF ( zrdry_incd(jl,jk,jn) < crdivm(inucs) ) THEN ! smaller 0.5 nm
                  dloss1_p     => dloss1(jn)%ptr
                  dloss1_p(jl,krow)=dloss1_p(jl,krow)                                &
                       +MIN(pxtp1c(jl,jk,idt_cd(jn)), &
                       pxtp1c(jl,jk,idt_cd(jn))/zaeromtot(jl,jk)*zaerodm(jl,jk))     &
                       *pdpg(jl,jk)*delta_time/time_step_len
                  IF(jn==id_so4) qloss1_ms4cd(jl,jk,krow)=qloss1_ms4cd(jl,jk,krow)   &
                       +MIN(pxtp1c(jl,jk,idt_cd(jn)), &
                       pxtp1c(jl,jk,idt_cd(jn))/zaeromtot(jl,jk)*zaerodm(jl,jk))     &
                       *vphysc%rhoam1(jl,jk,krow)*delta_time/time_step_len
               ELSE IF ( (zrdry_incd(jl,jk,jn) < crdivm(iaits)) .OR. &
                    ((jn==id_ss .OR. jn==id_du) .AND. (zrdry_incd(jl,jk,jn) >= crdivm(iaits) .AND. & ! loss of SS or DU for aitken mode sized particles
                    zrdry_incd(jl,jk,jn) < crdivm(iaccs))) ) THEN ! smaller 5 nm
                  dloss2_p     => dloss2(jn)%ptr
                  dloss2_p(jl,krow)=dloss2_p(jl,krow)                                &
                       +MIN(pxtp1c(jl,jk,idt_cd(jn)), &
                       pxtp1c(jl,jk,idt_cd(jn))/zaeromtot(jl,jk)*zaerodm(jl,jk))     &
                       *pdpg(jl,jk)*delta_time/time_step_len               
                  IF(jn==id_so4) qloss2_ms4cd(jl,jk,krow)=qloss2_ms4cd(jl,jk,krow)   &
                       +MIN(pxtp1c(jl,jk,idt_cd(jn)), &
                       pxtp1c(jl,jk,idt_cd(jn))/zaeromtot(jl,jk)*zaerodm(jl,jk))     &
                       *vphysc%rhoam1(jl,jk,krow)*delta_time/time_step_len
               ELSE IF ( zrdry_incd(jl,jk,jn) >= crdivm(nsol+1) ) THEN ! larger 50 um
                  dloss3_p     => dloss3(jn)%ptr
                  dloss3_p(jl,krow)=dloss3_p(jl,krow)                                &
                       +MIN(pxtp1c(jl,jk,idt_cd(jn)), &
                       pxtp1c(jl,jk,idt_cd(jn))/zaeromtot(jl,jk)*zaerodm(jl,jk))     &
                       *pdpg(jl,jk)*delta_time/time_step_len                             
                  IF(jn==id_so4) qloss3_ms4cd(jl,jk,krow)=qloss3_ms4cd(jl,jk,krow)   &
                       +MIN(pxtp1c(jl,jk,idt_cd(jn)), &
                       pxtp1c(jl,jk,idt_cd(jn))/zaeromtot(jl,jk)*zaerodm(jl,jk))     &
                       *vphysc%rhoam1(jl,jk,krow)*delta_time/time_step_len
               END IF
            END DO
         END IF
         
      END DO
   END DO
!<<DN

!>>DN [in-cloud mass change AFTER the change in cloud-free part]
    !---in-droplet mass:
    DO jn=aero_idx(1),aero_idx(naerospec) 
!>>DN [no water tracer in cloud particles]
      IF (jn==id_wat) CYCLE
!<<DN
      jtcl = idt_cd(jn)
      devap_p     => devap(jn)%ptr
      DO jk=1,klev
         WHERE (pnevap(1:kproma,jk) > zeps .AND. pnactevap(1:kproma,jk)  > zeps_num &
              .AND. zaeromtot(1:kproma,jk) > 0._dp)
            devap_p(1:kproma,krow) = devap_p(1:kproma,krow) + &
                                     MAX(0._dp,MIN(pxtp1c(1:kproma,jk,jtcl), & 
                                     pxtp1c(1:kproma,jk,jtcl)/pnactevap(1:kproma,jk)*pnevap(1:kproma,jk))) &
                                     *pdpg(1:kproma,jk)*delta_time/time_step_len
         END WHERE
      END DO
!>>DNdebug     
!       WHERE (pnevap(1:kproma,:) > zeps_num .AND. pxtp1c(1:kproma,:,idt_cdnc)  > zeps_num) 
       WHERE (pnevap(1:kproma,:) > zeps .AND. pnactevap(1:kproma,:)  > zeps_num &
              .AND. zaeromtot(1:kproma,:) > 0._dp)
!            .AND. rdry_incd(1:kproma,:,krow) >= crdivm(iaits)) 
!<<DNdebug
          !---reduce mass of in-cloud droplet aerosol
          pxtp1c(1:kproma,:,jtcl) = pxtp1c(1:kproma,:,jtcl) - &
                                    MAX(0._dp,MIN(pxtp1c(1:kproma,:,jtcl), & 
                                    pxtp1c(1:kproma,:,jtcl)/pnactevap(1:kproma,:)*pnevap(1:kproma,:)))
!          pxtp1c(1:kproma,:,jtcl) = MAX(pxtp1c(1:kproma,:,jtcl), 0._dp) !DN: leads to mass increase!
       END WHERE
    
   END DO
!<<DN
  END SUBROUTINE aeroproc_evap

  !------------------------------------------------------------------------
 
  SUBROUTINE aeroproc_sub(kbdim, kproma, klev, krow, pdpg, pxtp1c, pxtp10, pnactsub, pnsub)

    USE mo_kind,            ONLY: dp
    USE mo_tracdef,         ONLY: ntrac
!>>DN
!    USE mo_aero_m7,         ONLY: sizeclass, crdiv, iaits, iaccs, icoas, naermod, nsol
!    USE mo_aero_species,    ONLY: naerospec, aerospec, idt_cd, idt_ic, id_so4
    USE mo_ham,             ONLY: sizeclass, naerocomp, aerocomp, nsol, nclass
    USE mo_ham_m7ctl,       ONLY: crdiv, inucs, iaits, iaccs, icoas
    USE mo_ham_species ,    ONLY: idt_cd, idt_ic, id_so4, id_ss, id_du, id_wat
    USE mo_species,         ONLY: naerospec, speclist, aero_idx
    USE mo_ham_streams,     ONLY: rdry
!<<DN
    USE mo_activ,           ONLY: idt_icnc
    USE mo_time_control,    ONLY: delta_time, time_step_len
    USE mo_vphysc,          ONLY: vphysc
!    USE mo_outctl,          ONLY: lproc_detail

    INTEGER, INTENT(IN)  :: kbdim                            ! geographic block maximum number of locations
    INTEGER, INTENT(IN)  :: kproma                           ! geographic block number of locations on this row
    INTEGER, INTENT(IN)  :: krow                             ! geographic block number
    INTEGER, INTENT(IN)  :: klev                             ! number of vertical levels
    REAL(dp), INTENT(IN) :: pdpg(kbdim,klev)                 ! air mass per unit area
!>>DN
!    REAL(dp), INTENT(IN) :: pnsub(kbdim,klev)                ! ice crystal sublimation rate
!    REAL(dp), INTENT(IN) :: pnactsub(kbdim,klev)             ! 
    REAL(dp), INTENT(INOUT) :: pnsub(kbdim,klev)                ! ice crystal sublimation rate
    REAL(dp), INTENT(INOUT) :: pnactsub(kbdim,klev)             ! 
!<<DN
    REAL(dp), INTENT(INOUT) :: pxtp1c(kbdim,klev,ntrac)      ! tracer concentration in the cloudy part
    REAL(dp), INTENT(INOUT) :: pxtp10(kbdim,klev,ntrac)      ! tracer concentration in the cloud-free part
     
    !---local data
    REAL(dp) :: zvol(kbdim,klev)
    REAL(dp) :: zaeromtot(kbdim,klev)
    REAL(dp) :: zaerodm(kbdim,klev)
    REAL(dp) :: zmasschange(kbdim,klev)
!>>DN
    REAL(dp) :: zaerodn(kbdim,klev,nclass)                          ! change of number due to evaporation
!<<DN

!>>DN    
!    REAL(dp):: crdivm(nsol)                                  !threshold radii for M7 modes in m
    REAL(dp):: crdivm(nsol+1)                                  !threshold radii for M7 modes in m
!<<DN
    REAL(dp) :: zram2cmr_ic                                  ! conversion of radius of average aerosol mass -> count median radius
    
    INTEGER :: jn, jl, jk, jt, jtic, jtnum
    INTEGER :: ispec, imod
    REAL(dp), POINTER :: dloss1_p(:,:)
    REAL(dp), POINTER :: dloss2_p(:,:)
    REAL(dp), POINTER :: dloss3_p(:,:)
    REAL(dp), POINTER :: dsubl_p(:,:)
!>>DN
    LOGICAL :: lnumupdated_ic(kbdim,klev,nclass),ldnupdate(kbdim,klev) 
    REAL(dp) :: zrdry_inic(kbdim,klev,aero_idx(1):aero_idx(naerospec))               ! radius of released particles corrected for SS and DU
!<<DN

    !---executable procedure
    zvol(1:kproma,:) = 0._dp
    zaeromtot(1:kproma,:) = 0._dp
    zaerodm(1:kproma,:) = 0._dp
    crdivm(1:nsol) = 0.01_dp*crdiv(1:nsol)
!>>DN
    crdivm(nsol+1)=50.0E-6_dp
    lnumupdated_ic(1:kproma,:,:) = .FALSE.
    ldnupdate(1:kproma,:) = .FALSE.
    zaerodn(1:kproma,:,:) = 0._dp
    pnactsub(1:kproma,:) = MAX(pnsub(1:kproma,:),pnactsub(1:kproma,:))
    pnactsub(1:kproma,:) = MAX(0._dp,pnactsub(1:kproma,:))
    pnactsub(1:kproma,:) = MAX(10._dp/(vphysc%rhoam1(1:kproma,:,krow)+zeps),pnactsub(1:kproma,:))
!<<DN

    !---conversion factor radius of average mass -> count median radius
!>>DN
!    zram2cmr_ic = 1._dp / ( EXP(1.5_dp*(LOG(sigma_ic))**2) )
    zram2cmr_ic = 1._dp / ( EXP(1.5_dp*(LOG(1.59_dp))**2) )
!<<DN

    !---total in-ice crystal aerosol mass
    DO jn = aero_idx(1),aero_idx(naerospec) 
!>>DN [no water tracer in cloud particles]
       IF (jn==id_wat) CYCLE
       jtic = idt_ic(jn)
       zaeromtot(1:kproma,:) = zaeromtot(1:kproma,:) + MAX(0._dp,pxtp1c(1:kproma,:,jtic))          
    END DO
       
   !---sum (mass/density) for all species in cloud [ result in m**3/kg(air) ]
!>>DN [moved before change in droplet mass]
   DO jn=aero_idx(1),aero_idx(naerospec) 
!>>DN [no water tracer in cloud particles]
      IF (jn==id_wat) CYCLE
!<<DN
      jtic = idt_ic(jn)
!!>>DN [SO4 in kg(S) kg-1]
!       IF (jn==id_so4) THEN
!          zvol(1:kproma,:) = zvol(1:kproma,:) + pxtp1c(1:kproma,:,jtic)*astoso4/speclist(jn)%density
!       ELSE
          zvol(1:kproma,:) = zvol(1:kproma,:) + MAX(0._dp,pxtp1c(1:kproma,:,jtic))/speclist(jn)%density
!       END IF
!!<<DN
    END DO
!<<DN

    !---calculate radius of released particle
    WHERE (pnactsub(1:kproma,:) > zeps_num .AND. zvol(1:kproma,:) > 0._dp) 

       !---divide by droplet number concentration [ result in m**3 / drop ]
!       zvol(1:kproma,:) = zvol(1:kproma,:) / pnactsub(1:kproma,:)
       !---convert to radius of average mass
       rdry_inic(1:kproma,:,krow) = zram2cmr_ic * ((zvol(1:kproma,:)/pnactsub(1:kproma,:)/a4piover3)**(1._dp/3._dp))

    ELSEWHERE
       rdry_inic(1:kproma,:,krow) = 0._dp
    END WHERE

!>>DN [corrections for too small/too large released particles]
    !use coarse mode sigma, if rdry_inic larger than 0.5 um
    WHERE (rdry_inic(1:kproma,:,krow) > 0.5E-6_dp .AND. pnactsub(1:kproma,:) > zeps_num)  
       rdry_inic(1:kproma,:,krow) =  1._dp / ( EXP(1.5_dp*(LOG(sigma_ic))**2) ) &
                                     * ((zvol(1:kproma,:)/pnactsub(1:kproma,:)/a4piover3)**(1._dp/3._dp))
    END WHERE

!>>DN [additional diagnostics]
    DO jk=1,klev
!       IF (lproc_detail) THEN
          dsub_nic(1:kproma,krow) = dsub_nic(1:kproma,krow)-pnsub(1:kproma,jk)&
                               *pdpg(1:kproma,jk)*delta_time                                     
!       END IF
    END DO
!<<DN

    !prescribe rdry=100nm where rdry is too small
    WHERE (rdry_inic(1:kproma,:,krow) < crdivm(iaits))! .AND. &
!         pnsub(1:kproma,:) * (rdry_inic(1:kproma,:,krow)/100.E-9_dp)**3._dp <= 1.E6_dp .AND. &
!         pnactsub(1:kproma,:) > zeps_num .AND. zvol(1:kproma,:) > 0._dp)
       rdry_inic(1:kproma,:,krow)=100.E-9_dp
       pnsub(1:kproma,:)=zvol(1:kproma,:)/a4piover3/&
            (rdry_inic(1:kproma,:,krow)*EXP(1.5_dp*(LOG(1.59_dp))**2))**3
       pnactsub(1:kproma,:)=pnsub(1:kproma,:)
    END WHERE
   !prescribe rdry=10um where rdry is too large or where previous calculation yielded very high numbers  
!    WHERE ((rdry_inic(1:kproma,:,krow) < crdivm(iaits) .AND. &
!         pnactsub(1:kproma,:) > zeps_num .AND. zvol(1:kproma,:) > 0._dp) .OR. &
    WHERE ((rdry_inic(1:kproma,:,krow) == 100.E-9_dp .AND. &
         pnsub(1:kproma,:) > 1.E6_dp) .OR. &
         rdry_inic(1:kproma,:,krow) > 50.E-6_dp) 
       rdry_inic(1:kproma,:,krow)=10.E-6_dp      
       pnsub(1:kproma,:)=zvol(1:kproma,:)/a4piover3/&
            (rdry_inic(1:kproma,:,krow)*EXP(1.5_dp*(LOG(sigma_ic))**2))**3
       pnactsub(1:kproma,:)=pnsub(1:kproma,:)
    END WHERE

    DO jn = 1,nsol
       ![only calculate zaerodn here; pnsub <=pnactsub]
       WHERE (pnsub(1:kproma,:) > zeps .AND. pnactsub(1:kproma,:) > zeps_num .AND. &
            zaeromtot(1:kproma,:) > 0._dp) 
          zaerodn(1:kproma,:,jn) = MIN(pnsub(1:kproma,:),pnactsub(1:kproma,:))
       END WHERE
    END DO

    DO jn = aero_idx(1),aero_idx(naerospec) 
       zrdry_inic(1:kproma,:,jn)=rdry_inic(1:kproma,:,krow)
       !sea salt and dust: prescribe rdry=rdry_AS or 50nm where rdry is too small
       WHERE ((jn==id_ss .OR. jn==id_du) .AND. &
            rdry_inic(1:kproma,:,krow) < crdivm(iaccs) .AND. &
            pnsub(1:kproma,:) > zeps .AND. pnactsub(1:kproma,:) > zeps_num .AND. &
            zaeromtot(1:kproma,:) > 0._dp)
          zrdry_inic(1:kproma,:,jn)= MAX(rdry(iaccs)%ptr(1:kproma,:,krow),crdivm(iaccs))
          zaerodn(1:kproma,:,iaccs) = MIN(pnsub(1:kproma,:),pnactsub(1:kproma,:)) * &
               (rdry_inic(1:kproma,:,krow)/zrdry_inic(1:kproma,:,jn))**3._dp
       END WHERE
       WHERE ((jn==id_ss .OR. jn==id_du) .AND. &
            rdry_inic(1:kproma,:,krow) < crdivm(iaccs) .AND. &
            pnsub(1:kproma,:) > zeps .AND. pnactsub(1:kproma,:) > zeps_num .AND. &
            zaeromtot(1:kproma,:) > 0._dp .AND. .NOT. ldnupdate(1:kproma,:))
          zaerodn(1:kproma,:,iaits) = zaerodn(1:kproma,:,iaits)-zaerodn(1:kproma,:,iaccs)
          ldnupdate(1:kproma,:)=.TRUE.
       END WHERE
    END DO
    
   ![only calculate zaerodm here; pnsub <=pnactsub]
    WHERE (pnsub(1:kproma,:) > zeps .AND. pnactsub(1:kproma,:) > zeps_num .AND. &
         zaeromtot(1:kproma,:) > 0._dp) 
       zaerodm(1:kproma,:) = MIN(zaeromtot(1:kproma,:),zaeromtot(1:kproma,:) / &
            pnactsub(1:kproma,:) * pnsub(1:kproma,:)) 
    END WHERE
!<<DN

    !---attribute number and mass change to correct mode

    !---loop over all aerosol species in all modes
    DO jn=1,naerocomp
       ispec = aerocomp(jn)%spid
       imod = aerocomp(jn)%iclass
       
       !---get tracer index of aerosol number
       jtnum = sizeclass(imod)%idt_no
       
       !---get tracer index to aerosol 
       jt = aerocomp(jn)%idt

       DO jk=1,klev
          DO jl=1,kproma
             
             IF ( (imod == iaits .AND. zrdry_inic(jl,jk,ispec) >= crdivm(iaits) .AND. &
                  zrdry_inic(jl,jk,ispec) < crdivm(iaccs)) .OR. &
                  (imod == iaccs .AND. zrdry_inic(jl,jk,ispec) >= crdivm(iaccs) .AND. &
                  zrdry_inic(jl,jk,ispec) < crdivm(icoas)) .OR. &
                  (imod == icoas .AND. zrdry_inic(jl,jk,ispec) >= crdivm(icoas) .AND. &
                  zrdry_inic(jl,jk,ispec) < crdivm(nsol+1)) ) THEN
                
                   !---increase aerosol number in cloud-free part
!>>DN [update number concentration only once per mode; additional diagnostics]
                IF (.NOT.lnumupdated_ic(jl,jk,imod)) THEN
                   pxtp10(jl,jk,jtnum) = pxtp10(jl,jk,jtnum) + zaerodn(jl,jk,imod)
                   lnumupdated_ic(jl,jk,imod) = .TRUE.
!                   IF (lproc_detail) THEN
                      qsub_inic(jl,jk,krow) = qsub_inic(jl,jk,krow)             &
                                              - zaerodn(jl,jk,imod)*vphysc%rhoam1(jl,jk,krow)    &
                                              *delta_time/time_step_len
                      dsub_inic(jl,krow) = dsub_inic(jl,krow)-zaerodn(jl,jk,imod)&
                                          *pdpg(jl,jk)*delta_time           
!                   END IF
                END IF

                   !---increase mass in the cloud-free part
!>>DN [index must refer to in-coud species]
!                   zmasschange(jl,jk) = pxtp1c(jl,jk,jt)/zaeromtot(jl,jk)*zaerodm(jl,jk)
                IF (zaeromtot(jl,jk)>0._dp .AND. pxtp1c(jl,jk,idt_ic(ispec))>0._dp) THEN
                   zmasschange(jl,jk) = MIN(pxtp1c(jl,jk,idt_ic(ispec)), &
                        pxtp1c(jl,jk,idt_ic(ispec))/zaeromtot(jl,jk)*zaerodm(jl,jk))
                ELSE
                   zmasschange(jl,jk) = 0._dp
                END IF
!!<<DN
                pxtp10(jl,jk,jt) = pxtp10(jl,jk,jt) + zmasschange(jl,jk)      

                   !---diagnostics
!                IF (lproc_detail) THEN
!>>DNdebug
                   IF (ispec == id_so4) THEN
!                   IF (jn == id_so4) THEN
!<<DNdebug
                      qsub_ms4ic(jl,jk,krow) = qsub_ms4ic(jl,jk,krow)             &
                                              - zmasschange(jl,jk)*vphysc%rhoam1(jl,jk,krow)    &
                                              *delta_time/time_step_len
                      END IF
                   END IF

                   !---TODO (maybe): re-introduce diagnostic of lost SS and DU mass 
!                END IF

             END DO
          END DO

!!$       !---aitken mode
!!$       IF (imod == iaits) THEN
!!$          DO jk=1,klev
!!$             DO jl=1,kproma
!!$                
!!$                IF ( rdry_inic(jl,jk,krow) >= crdivm(iaits) .AND. &
!!$                     rdry_inic(jl,jk,krow) < crdivm(iaccs) ) THEN
!!$
!!$                   !---increase aerosol number in cloud-free part
!!$!>>DN [update number concentration only once per mode]
!!$                   IF (.NOT.lnumupdated_ic(imod)) THEN
!!$                      pxtp10(jl,jk,jtnum) = pxtp10(jl,jk,jtnum) + pnsub(jl,jk)
!!$                      lnumupdated_ic(imod) = .TRUE.
!!$                      IF (lproc_detail) THEN
!!$                         qsub_inic(jl,jk,krow) = qsub_inic(jl,jk,krow)             &
!!$                                                 -pnsub(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                 *delta_time/time_step_len
!!$                                                 /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$
!!$                   !---increase mass in the cloud-free part
!!$!>>DN [index must refer to in-coud species]
!!$!                   zmasschange(jl,jk) = pxtp1c(jl,jk,jt)/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$                   zmasschange(jl,jk) = pxtp1c(jl,jk,idt_ic(ispec))/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$!<<DN
!!$                   pxtp10(jl,jk,jt) = pxtp10(jl,jk,jt) + zmasschange(jl,jk)      
!!$
!!$                   !---diagnostics
!!$                   IF (lproc_detail) THEN
!!$!>>DNdebug
!!$                      IF (ispec == id_so4) THEN
!!$!                      IF (jn == id_so4) THEN
!!$!<<DNdebug
!!$                         qsub_ms4ic(jl,jk,krow) = qsub_ms4ic(jl,jk,krow)             &
!!$                                                  - zmasschange(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                  *delta_time/time_step_len
!!$                                                  /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$
!!$                   !---TODO (maybe): re-introduce diagnostic of lost SS and DU mass 
!!$                END IF
!!$             END DO
!!$          END DO
!!$       END IF
!!$             
!!$       !---accumulation mode
!!$       IF (imod == iaccs) THEN
!!$          DO jk = 1,klev
!!$             DO jl=1,kproma
!!$                
!!$                IF ( rdry_inic(jl,jk,krow) >= crdivm(iaccs) .AND. &
!!$                     rdry_inic(jl,jk,krow) < crdivm(icoas) ) THEN
!!$
!!$                   !---increase aerosol number in cloud-free part
!!$!>>DN [update number concentration only once per mode]
!!$                   IF (.NOT.lnumupdated_ic(imod)) THEN
!!$                      pxtp10(jl,jk,jtnum) = pxtp10(jl,jk,jtnum) + pnsub(jl,jk)
!!$                      lnumupdated_ic(imod) = .TRUE.
!!$                      IF (lproc_detail) THEN
!!$                         qsub_inic(jl,jk,krow) = qsub_inic(jl,jk,krow)             &
!!$                                                 -pnsub(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                 *delta_time/time_step_len
!!$                                                 /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$
!!$                   !---increase mass in the cloud-free part
!!$!>>DN [index must refer to in-coud species]
!!$!                   zmasschange(jl,jk) = pxtp1c(jl,jk,jt)/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$                   zmasschange(jl,jk) = pxtp1c(jl,jk,idt_ic(ispec))/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$!<<DN
!!$                   pxtp10(jl,jk,jt) = pxtp10(jl,jk,jt) + zmasschange(jl,jk)      
!!$
!!$                   !---diagnostics
!!$                   IF (lproc_detail) THEN
!!$!>>DNdebug
!!$                      IF (ispec == id_so4) THEN
!!$!                      IF (jn == id_so4) THEN
!!$!<<DNdebug
!!$                         qsub_ms4ic(jl,jk,krow) = qsub_ms4ic(jl,jk,krow)             &
!!$                                                  - zmasschange(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                  *delta_time/time_step_len
!!$                                                  /time_step_len
!!$!<<DNdebug
!!$
!!$                      END IF
!!$                   END IF
!!$
!!$                   !---TODO (maybe): re-introduce diagnostic of lost SS and DU mass 
!!$                END IF
!!$             END DO
!!$          END DO
!!$       END IF
!!$             
!!$       !---coarse mode   
!!$       IF (imod == icoas) THEN
!!$          DO jk = 1,klev
!!$             DO jl=1,kproma
!!$                
!!$                IF ( rdry_inic(jl,jk,krow) > crdivm(icoas)) THEN
!!$
!!$                   !---increase aerosol number in cloud-free part
!!$!>>DN [update number concentration only once per mode]
!!$                   IF (.NOT.lnumupdated_ic(imod)) THEN
!!$                      pxtp10(jl,jk,jtnum) = pxtp10(jl,jk,jtnum) + pnsub(jl,jk)
!!$                      lnumupdated_ic(imod) = .TRUE.
!!$                      IF (lproc_detail) THEN
!!$                         qsub_inic(jl,jk,krow) = qsub_inic(jl,jk,krow)             &
!!$                                                 -pnsub(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                 *delta_time/time_step_len
!!$                                                 /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$
!!$                   !---increase mass in the cloud-free part
!!$!>>DN [index must refer to in-coud species]
!!$!                   zmasschange(jl,jk) = pxtp1c(jl,jk,jt)/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$                   zmasschange(jl,jk) = pxtp1c(jl,jk,idt_ic(ispec))/zaeromtot(jl,jk)*zaerodm(jl,jk)
!!$!<<DN
!!$                   pxtp10(jl,jk,jt) = pxtp10(jl,jk,jt) + zmasschange(jl,jk)      
!!$
!!$                   !---diagnostics
!!$                   IF (lproc_detail) THEN
!!$!>>DNdebug
!!$                      IF (ispec == id_so4) THEN
!!$!                      IF (jn == id_so4) THEN
!!$!<<DNdebug
!!$                         qsub_ms4ic(jl,jk,krow) = qsub_ms4ic(jl,jk,krow)             &
!!$                                                  - zmasschange(jl,jk)*rhoam1(jl,jk,krow)    &
!!$!>>DNdebug
!!$!                                                  *delta_time/time_step_len
!!$                                                  /time_step_len
!!$!<<DNdebug
!!$                      END IF
!!$                   END IF
!!$
!!$                   !---TODO (maybe): re-introduce diagnostic of lost SS and DU mass 
!!$                END IF
!!$             END DO
!!$          END DO
!!$       END IF
 
   END DO

   !---TODO: diagnostic of any cases where rdry_inic < aitken mode min, or very large values (> 50um)

!>>DN
   !diagnostic of any cases where rdry_inic < aitken mode min, or for very large values (> 50um)
   DO jk=1,klev
      DO jl=1,kproma
         
         IF ( zaeromtot(jl,jk) > 0._dp .AND. &
              pnsub(jl,jk) > zeps .AND. pnactsub(jl,jk)  > zeps_num) THEN
            DO jn=aero_idx(1),aero_idx(naerospec) 
               IF (jn==id_wat) CYCLE
               IF ( zrdry_inic(jl,jk,jn) < crdivm(inucs) ) THEN ! smaller 0.5 nm
                  dloss1_p     => dloss1(jn)%ptr
                  dloss1_p(jl,krow)=dloss1_p(jl,krow)                                &
                       +MIN(pxtp1c(jl,jk,idt_ic(jn)), &
                       pxtp1c(jl,jk,idt_ic(jn))/zaeromtot(jl,jk)*zaerodm(jl,jk))     &
                       *pdpg(jl,jk)*delta_time/time_step_len
               ELSE IF ( (zrdry_inic(jl,jk,jn) < crdivm(iaits)) .OR. &
                    ((jn==id_ss .OR. jn==id_du) .AND. (zrdry_inic(jl,jk,jn) >= crdivm(iaits) .AND. & ! loss of SS or DU for aitken mode sized particles
                    zrdry_inic(jl,jk,jn) < crdivm(iaccs))) ) THEN ! smaller 5 nm
                  dloss2_p     => dloss2(jn)%ptr
                  dloss2_p(jl,krow)=dloss2_p(jl,krow)                                &
                       +MIN(pxtp1c(jl,jk,idt_ic(jn)), &
                       pxtp1c(jl,jk,idt_ic(jn))/zaeromtot(jl,jk)*zaerodm(jl,jk))     &
                       *pdpg(jl,jk)*delta_time/time_step_len               
               ELSE IF ( zrdry_inic(jl,jk,jn) >= crdivm(nsol+1) ) THEN ! larger 50 um
                  dloss3_p     => dloss3(jn)%ptr
                  dloss3_p(jl,krow)=dloss3_p(jl,krow)                                &
                       +MIN(pxtp1c(jl,jk,idt_ic(jn)), &
                       pxtp1c(jl,jk,idt_ic(jn))/zaeromtot(jl,jk)*zaerodm(jl,jk))     &
                       *pdpg(jl,jk)*delta_time/time_step_len                             
               END IF
            END DO
         END IF
         
      END DO
   END DO
!<<DN

!>>DN [in-cloud mass change AFTER the change in cloud-free part; pnsub <=pnactsub]
    !---in-crystal mass:
    DO jn=aero_idx(1),aero_idx(naerospec)
!>>DN [no water tracer in cloud particles]
       IF (jn==id_wat) CYCLE
!<<DN
       jtic = idt_ic(jn)
       dsubl_p     => dsubl(jn)%ptr
       DO jk=1,klev
          WHERE (pnsub(1:kproma,jk) > zeps .AND. pnactsub(1:kproma,jk) > zeps_num &
               .AND. zaeromtot(1:kproma,jk) > 0._dp)
             dsubl_p(1:kproma,krow) = dsubl_p(1:kproma,krow) + &
                                      MAX(0._dp,MIN(pxtp1c(1:kproma,jk,jtic), & 
                                      pxtp1c(1:kproma,jk,jtic)*pnsub(1:kproma,jk)/ pnactsub(1:kproma,jk))) & 
                                      *pdpg(1:kproma,jk)*delta_time/time_step_len
          END WHERE
       END DO
!>>DNdebug     
!       WHERE (pnsub(1:kproma,:) > zeps_num .AND. pxtp1c(1:kproma,:,idt_icnc) > zeps_num) 
       WHERE (pnsub(1:kproma,:) > zeps .AND. pnactsub(1:kproma,:) > zeps_num &
            .AND. zaeromtot(1:kproma,:) > 0._dp)
!            .AND. rdry_inic(1:kproma,:,krow) >= crdivm(iaits))
!<<DNdebug
          !---reduce mass of in-ice crystal aerosol
          pxtp1c(1:kproma,:,jtic) = pxtp1c(1:kproma,:,jtic) - &
                                    MAX(0._dp,MIN(pxtp1c(1:kproma,:,jtic),& 
                                    pxtp1c(1:kproma,:,jtic)*pnsub(1:kproma,:)/ pnactsub(1:kproma,:)))
          !###DEBUG
!          pxtp1c(1:kproma,:,jtic) = MAX(pxtp1c(1:kproma,:,jtic), 0._dp)!DN: leads to mass increase!      
       END WHERE
   END DO
!<<DN   

  END SUBROUTINE aeroproc_sub

  !------------------------------------------------------------------------

  SUBROUTINE aeroproc_coll(kbdim, kproma, klev, krow, pdpg, pxtp1c)

    USE mo_kind,            ONLY: dp
    USE mo_tracdef,         ONLY: ntrac
    USE mo_time_control,    ONLY: time_step_len, delta_time
    USE mo_ham,             ONLY: sizeclass, naerocomp, aerocomp, nclass
    USE mo_activ,           ONLY: idt_cdnc, idt_icnc
    USE mo_ham_species,     ONLY: idt_cd, idt_ic, id_so4
    USE mo_vphysc,          ONLY: vphysc
  USE mo_physical_constants,     ONLY: tmelt !davidn
  USE mo_echam_cloud_params,         ONLY: cthomi !davidn
   USE mo_ham_m7ctl,       ONLY: icoas, iaits, iaccs !davidn
    USE mo_ham_species,     ONLY: id_bc, id_oc, id_ss, id_du !davidn
 
    INTEGER, INTENT(IN)  :: kbdim                            ! geographic block maximum number of locations
    INTEGER, INTENT(IN)  :: kproma                           ! geographic block number of locations on this row
    INTEGER, INTENT(IN)  :: krow                             ! geographic block number
    INTEGER, INTENT(IN)  :: klev                             ! number of vertical levels
    REAL(dp), INTENT(IN) :: pdpg(kbdim,klev)                 ! air mass per unit area
    REAL(dp), INTENT(INOUT) :: pxtp1c(kbdim,klev,ntrac)      ! tracer concentration in the cloudy part
  
    REAL(dp):: zaercollcd(kbdim,klev), zhelp1cd(kbdim,klev,nclass), zhelp2cd(kbdim,klev)
    REAL(dp):: zaercollic(kbdim,klev), zhelp1ic(kbdim,klev,nclass), zhelp2ic(kbdim,klev)

    LOGICAL :: lcoll_maskcd(kbdim,klev)
    LOGICAL :: lcoll_maskic(kbdim,klev)
    
    INTEGER :: jn, jl, jk, jt, jtnum, jtcl, jtic
    INTEGER :: ispec, imod

    !  REAL, PARAMETER :: collkernel(nclass) &
    !       =(/1.E-11,1.E-12,1.E-13,1.E-14,1.E-12,1.E-13,1.E-14/)  
    !m3/s    !estimates from Wang&Lin
    !  REAL, PARAMETER :: collkernel(nclass)=(/3.E-12,3.E-12,4.E-13,2.E-13,3.E-12,4.E-13,2.E-13/)   
    !m3/s  !for evaporating droplet, Young
!>>DN
!    REAL(dp), PARAMETER :: collkernel(nclass) = (/2.5E-12_dp,2.5E-12_dp,2.E-14_dp,0._dp,2.5E-12_dp,2.E-14_dp,0._dp/)   
!no parameter because nclass is non-constant  
    REAL(dp) :: collkernel_cd(nclass) 
!<<DN
    !m3/s  !for growing droplet, Young (1974)
!davidn
!    REAL(dp), PARAMETER :: collkernel_cd(nclass)=(/3.E-12_dp,3.E-12_dp,4.E-12_dp,2.E-12_dp,3.E-12_dp,4.E-12_dp,2.E-12_dp/)   !(CS,AS,CI,AI)*10
       !m3/s  !for evaporating droplet, Young
!davidn
!>>DN [collision kernel for aerosol-crystal collisions]
!no parameter because nclass is non-constant  
    REAL(dp) :: collkernel_ic(nclass) 
    !m3/s  !for 15 mum growing droplet, Young (1974)
!<<DN

!>>DN [initialize here because nclass is non-constant]
    collkernel_cd(1:nclass) = (/2.5E-12_dp,2.5E-12_dp,2.E-14_dp,0._dp,2.5E-12_dp,2.E-14_dp,0._dp/) 
    collkernel_ic(1:nclass) = (/5.E-11_dp,5.E-11_dp,2.E-12_dp,2.E-13_dp,5.E-11_dp,2.E-12_dp,2.E-13_dp/)  
!<<DN

    !calculate collision rates
    !---aerosol number
    DO imod=1,nclass
       jtnum = sizeclass(imod)%idt_no

!>>DN [add collision kernel for aerosol-crystal collisions]
!       zhelp1cd(1:kproma,:,imod) = (1._dp - EXP(-collkernel(imod)*pxtp1c(1:kproma,:,idt_cdnc)    &
!                                                *rhoam1(1:kproma,:,krow)*time_step_len) )
!
!       zhelp1ic(1:kproma,:,imod) = (1._dp - EXP(-collkernel(imod)*pxtp1c(1:kproma,:,idt_icnc)    &
!                                                *rhoam1(1:kproma,:,krow)*time_step_len) )

       zhelp1cd(1:kproma,:,imod) = (1._dp - EXP(-collkernel_cd(imod)*pxtp1c(1:kproma,:,idt_cdnc)    &
                                                *vphysc%rhoam1(1:kproma,:,krow)*time_step_len) )

       zhelp1ic(1:kproma,:,imod) = (1._dp - EXP(-collkernel_ic(imod)*pxtp1c(1:kproma,:,idt_icnc)    &
                                                *vphysc%rhoam1(1:kproma,:,krow)*time_step_len) )
!<<DN
!>>DNdebug
       zhelp1cd(1:kproma,:,imod) = MAX(0._dp,MIN(1._dp,zhelp1cd(1:kproma,:,imod)))
       zhelp1ic(1:kproma,:,imod) = MAX(0._dp,MIN(1._dp,zhelp1ic(1:kproma,:,imod)))
!<<DNdebug

!>>DN [collision scavenging of cloud droplets before collision scavenging of ice crystals]
       !cloud droplets
       lcoll_maskcd(1:kproma,:) = pxtp1c(1:kproma,:,idt_cdnc) > zeps .AND. pxtp1c(1:kproma,:,jtnum) > 0._dp
       zhelp2cd(1:kproma,:) = pxtp1c(1:kproma,:,jtnum) * zhelp1cd(1:kproma,:,imod)
       zaercollcd(1:kproma,:) = MERGE( MIN(zhelp2cd(1:kproma,:), pxtp1c(1:kproma,:,jtnum)), &
                                           0._dp, lcoll_maskcd(1:kproma,:) )
       !---loss of interstitial number
       pxtp1c(1:kproma,:,jtnum) = pxtp1c(1:kproma,:,jtnum) - zaercollcd(1:kproma,:)
       pxtp1c(1:kproma,:,jtnum) = MAX(pxtp1c(1:kproma,:,jtnum),0._dp)
                                         
       !ice crystals
       lcoll_maskic(1:kproma,:) = pxtp1c(1:kproma,:,idt_icnc) > zeps .AND. pxtp1c(1:kproma,:,jtnum) > 0._dp
       zhelp2ic(1:kproma,:) = pxtp1c(1:kproma,:,jtnum) * zhelp1ic(1:kproma,:,imod)
       zaercollic(1:kproma,:) = MERGE( MIN(zhelp2ic(1:kproma,:), pxtp1c(1:kproma,:,jtnum)), &
                                           0._dp, lcoll_maskic(1:kproma,:) )
       !---loss of interstitial number
       pxtp1c(1:kproma,:,jtnum) = pxtp1c(1:kproma,:,jtnum) - zaercollic(1:kproma,:)
       pxtp1c(1:kproma,:,jtnum) = MAX(pxtp1c(1:kproma,:,jtnum),0._dp)

       !---3-D diagnostics for aerosol number
       qcoll_incd(1:kproma,:,krow) = qcoll_incd(1:kproma,:,krow)  +   &
                                     zaercollcd(1:kproma,:)*vphysc%rhoam1(1:kproma,:,krow)*delta_time/time_step_len
       qcoll_inic(1:kproma,:,krow) = qcoll_inic(1:kproma,:,krow)  +   &
                                     zaercollic(1:kproma,:)*vphysc%rhoam1(1:kproma,:,krow)*delta_time/time_step_len
    END DO
    
    !---aerosol mass
    DO jn=1,naerocomp
       
       ispec = aerocomp(jn)%spid
       imod = aerocomp(jn)%iclass
       
       jtcl = idt_cd(ispec)
       jtic = idt_ic(ispec)

!>>DNdebug
!       jt = aeromode(imod)%idt
       jt = aerocomp(jn)%idt
!<<DNdebug

       !cloud droplets
       lcoll_maskcd(1:kproma,:) = pxtp1c(1:kproma,:,idt_cdnc) > zeps .AND. pxtp1c(1:kproma,:,jt) > 0._dp
       zhelp2cd(1:kproma,:) = pxtp1c(1:kproma,:,jt) * zhelp1cd(1:kproma,:,imod)
       zaercollcd(1:kproma,:) = MERGE( MIN(zhelp2cd(1:kproma,:), pxtp1c(1:kproma,:,jt)), &
                                       0._dp, lcoll_maskcd(1:kproma,:) )
       !---loss of interstitial mass
       pxtp1c(1:kproma,:,jt) = pxtp1c(1:kproma,:,jt) - zaercollcd(1:kproma,:)
       !---gain of in-droplet mass 
       pxtp1c(1:kproma,:,jtcl) = pxtp1c(1:kproma,:,jtcl) + zaercollcd(1:kproma,:) 
 
       !ice crystals
       lcoll_maskic(1:kproma,:) = pxtp1c(1:kproma,:,idt_icnc) > zeps .AND. pxtp1c(1:kproma,:,jt) > 0._dp
       zhelp2ic(1:kproma,:) = pxtp1c(1:kproma,:,jt) * zhelp1ic(1:kproma,:,imod)
       zaercollic(1:kproma,:) = MERGE( MIN(zhelp2ic(1:kproma,:), pxtp1c(1:kproma,:,jt)), &
                                       0._dp, lcoll_maskic(1:kproma,:) )
      !---loss of interstitial mass 
       pxtp1c(1:kproma,:,jt) = pxtp1c(1:kproma,:,jt) - zaercollic(1:kproma,:)
       !---gain of in-crystal mass 
       pxtp1c(1:kproma,:,jtic) = pxtp1c(1:kproma,:,jtic) + zaercollic(1:kproma,:) 

       !---diagnostics
       !   3-D diagnostic for sulfate mass:
       IF (ispec == id_so4) THEN
          qcoll_ms4cd(1:kproma,:,krow) = qcoll_ms4cd(1:kproma,:,krow) +   &
                                         zaercollcd(1:kproma,:)*vphysc%rhoam1(1:kproma,:,krow)*delta_time/time_step_len
          qcoll_ms4ic(1:kproma,:,krow) = qcoll_ms4ic(1:kproma,:,krow) +   &
                                         zaercollic(1:kproma,:)*vphysc%rhoam1(1:kproma,:,krow)*delta_time/time_step_len
       END IF
                   
       !---2-D diagnostics for each species
       DO jk=1,klev
          dcollcd(ispec)%ptr(1:kproma,krow) = dcollcd(ispec)%ptr(1:kproma,krow)    +           &
                                              zaercollcd(1:kproma,jk)*pdpg(1:kproma,jk)*delta_time/time_step_len
          dcollic(ispec)%ptr(1:kproma,krow) = dcollic(ispec)%ptr(1:kproma,krow)    +           &
                                              zaercollic(1:kproma,jk)*pdpg(1:kproma,jk)*delta_time/time_step_len
       END DO

    END DO
!<<DN

  END SUBROUTINE aeroproc_coll

  !------------------------------------------------------------------------

END MODULE mo_ham_m7_processing
