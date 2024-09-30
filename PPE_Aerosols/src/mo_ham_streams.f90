!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_streams.f90
!!
!! \brief
!! Generalized definition of diagnostic output for the aerosol model. Only one
!! diagnostic stream (ham) left. Other output in wetdep, drydep or sedi streams
!!
!! \author Martin G. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin G. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# Martin G. Schultz (FZ Juelich) - original code (2009-09-30)
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!!   -# Philip Stier (University of Oxford) - additional radiation diagnostics (2014-04-04)
!!
!! \limitations
!! None
!!
!! \details
!! This code is based on mo_aero_mem by P. Stier, D. O'Donnell and others. Generalized 
!! by Martin Schultz.
!! *NOTE* All d_emi and emi stuff removed and re-implemented in mo_hammoz_emissions.f90
!!
!! \bibliographic_references
!! None 
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ham_streams

  USE mo_kind,          ONLY: dp
  USE mo_submodel_diag, ONLY: vmem2d, vmem3d
  USE mo_ham_rad_data,  ONLY: Nwv_tot, &
!>>NAJS: optical properties of dry aerosol
                              Niwv_sw_dry, iraddry, &
!<<NAJS
!>>NAJS: lidar backscatter
                              Niwv_sw_beta
!<<NAJS
  USE mo_ham,           ONLY: nmaxclass
  USE mo_ham_ccnctl,    ONLY: nsat

  IMPLICIT NONE

!!  TYPE, PUBLIC :: vmem3d
!!     REAL(dp), POINTER  :: ptr(:,:,:)
!!  END TYPE vmem3d

  PRIVATE

  !--- Service routines

  PUBLIC :: new_stream_ham       ! construct the ham streams
  PUBLIC :: new_stream_ham_rad   ! construct the ham_rad stream

  !-- pointers for diagnostic quantities
  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: rdry(:)          ! used in cloud_cdnc_icnc
  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: rwet(:)          ! used in cuasc and mo_drydep
  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: densaer(:)       ! used in mo_drydep

  !-- some more ...
  REAL(dp), PUBLIC, POINTER :: relhum(:,:,:)             ! ### replace by ECHAM variable ??
  REAL(dp), PUBLIC, POINTER :: ipr(:,:,:)
  REAL(dp), PUBLIC, POINTER :: daylength(:,:)

  !-- chemical production and loss terms
  REAL(dp), PUBLIC, POINTER :: d_prod_ms4as(:,:)
  REAL(dp), PUBLIC, POINTER :: d_prod_ms4cs(:,:)
   
  REAL(dp), PUBLIC, POINTER :: d_prod_so2_dms_oh(:,:)
  REAL(dp), PUBLIC, POINTER :: d_prod_so4_dms_oh(:,:)
  REAL(dp), PUBLIC, POINTER :: d_prod_so2_dms_no3(:,:)
  REAL(dp), PUBLIC, POINTER :: d_prod_so4_so2_oh(:,:)

  REAL(dp), PUBLIC, POINTER :: d_prod_h2so4(:,:,:)
  REAL(dp), PUBLIC, POINTER :: d_prod_so4_liq(:,:,:)

  REAL(dp), PUBLIC, POINTER :: d_cond_so4(:,:)

  REAL(dp), PUBLIC, POINTER :: d_nuc_so4(:,:)

!>>SF
  !-- diagnostics for cloud activation
  REAL(dp), PUBLIC, POINTER :: nbcsol_diag(:,:,:)
  REAL(dp), PUBLIC, POINTER :: nbcinsol_diag(:,:,:)
  REAL(dp), PUBLIC, POINTER :: nbcsol_acc(:,:,:)
  REAL(dp), PUBLIC, POINTER :: nbcsol_ait(:,:,:)
  REAL(dp), PUBLIC, POINTER :: nduinsol_diag(:,:,:)
  REAL(dp), PUBLIC, POINTER :: nbcsol_strat(:,:,:)
  REAL(dp), PUBLIC, POINTER :: nbcsol_cv(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ndusol_strat(:,:,:)
  REAL(dp), PUBLIC, POINTER :: ndusol_cv(:,:,:)
  REAL(dp), PUBLIC, POINTER :: nbcinsol(:,:,:)
  REAL(dp), PUBLIC, POINTER :: nduinsolai(:,:,:)
  REAL(dp), PUBLIC, POINTER :: nduinsolci(:,:,:)
  REAL(dp), PUBLIC, POINTER :: naerinsol(:,:,:)

  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: nact_strat(:)
  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: nact_conv(:)
  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: frac(:)
  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: a(:)
  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: b(:)
  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: sc(:)
  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: rc_strat(:,:)
  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: rc_conv(:,:)
!<<SF

  !-- mo_ham_ccn diagnostics
  TYPE (vmem2d), PUBLIC :: ccn_2d(nsat)
  TYPE (vmem2d), PUBLIC :: ccn_burden(nsat)
  TYPE (vmem3d), PUBLIC :: ccn_3d(nsat)
  REAL(dp), PUBLIC, POINTER :: cn_2d(:,:)
  REAL(dp), PUBLIC, POINTER :: cn_burden(:,:)
  REAL(dp), PUBLIC, POINTER :: cn_3d(:,:,:)

  !-- ham radiation diagnostics (former mo_aero_rad_mem)
  TYPE (vmem3d), PUBLIC :: tau_mode(nmaxclass,Nwv_tot)
  TYPE (vmem3d), PUBLIC :: abs_mode(nmaxclass,Nwv_tot)
  TYPE (vmem3d), PUBLIC :: sigma_mode(nmaxclass,Nwv_tot)
  TYPE (vmem3d), PUBLIC :: omega_mode(nmaxclass,Nwv_tot)
  TYPE (vmem3d), PUBLIC :: asym_mode(nmaxclass,Nwv_tot)
  TYPE (vmem3d), PUBLIC :: nr_mode(nmaxclass,Nwv_tot)
  TYPE (vmem3d), PUBLIC :: ni_mode(nmaxclass,Nwv_tot)

!  TYPE (vmem2d), PUBLIC :: tau_2d_mode(nmaxclass,Nwv_tot)
  TYPE (vmem2d), PUBLIC :: sigma_2d_mode(nmaxclass,Nwv_tot)
  TYPE (vmem2d), PUBLIC :: omega_2d_mode(nmaxclass,Nwv_tot)
  TYPE (vmem2d), PUBLIC :: asym_2d_mode(nmaxclass,Nwv_tot)
  TYPE (vmem2d), PUBLIC :: nr_2d_mode(nmaxclass,Nwv_tot)
  TYPE (vmem2d), PUBLIC :: ni_2d_mode(nmaxclass,Nwv_tot)
  TYPE (vmem2d), PUBLIC :: tau_2d(Nwv_tot)
  TYPE (vmem2d), PUBLIC :: abs_2d(Nwv_tot)
  TYPE (vmem2d), PUBLIC :: ant_2d(Nwv_tot)
!>>NAJS: optical properties of dry aerosol
  TYPE (vmem2d), PUBLIC :: tau_dry(Niwv_sw_dry)
  TYPE (vmem2d), PUBLIC :: abs_dry(Niwv_sw_dry)
  TYPE (vmem3d), PUBLIC :: ext_dry(Niwv_sw_dry)
!<<NAJS
!>>DN: AeroCom
  TYPE (vmem3d), PUBLIC :: abs_lev_dry(Niwv_sw_dry)
!<<DN

  TYPE (vmem2d), ALLOCATABLE, PUBLIC :: tau_comp(:,:)
  TYPE (vmem2d), ALLOCATABLE, PUBLIC :: abs_comp(:,:)

  REAL(dp), PUBLIC, POINTER :: ang(:,:)
  REAL(dp), PUBLIC, POINTER :: ai(:,:)

!>>NAJS: lidar backscatter
  TYPE (vmem3d), PUBLIC :: attbeta_tot_space(Niwv_sw_beta) ! 3D attenuated aerosol and gas backscatter from space
  TYPE (vmem3d), PUBLIC :: attbeta_tot_earth(Niwv_sw_beta) ! 3D attenuated aerosol and gas backscatter from earth
  TYPE (vmem3d), PUBLIC :: beta_tot(Niwv_sw_beta) ! 3D aerosol and gas backscatter
  TYPE (vmem3d), PUBLIC :: alfa_aer(Niwv_sw_beta) ! 3D aerosol extinction
  TYPE (vmem3d), PUBLIC :: alfa_gas(Niwv_sw_beta) ! 3D gaseous extinction
  TYPE (vmem3d), PUBLIC :: lidar_ratio_aer(Niwv_sw_beta) ! 3D LIDAR ratio for aerosol
!<<NAJS


!----------------------------------------------------------------------------------------------------------------

  CONTAINS

!
!! brief: define diagnostics for the HAM stream
!

  SUBROUTINE new_stream_ham (nclass)

    USE mo_linked_list,   ONLY: t_stream
    USE mo_filename,      ONLY: trac_filetype
    USE mo_linked_list,   ONLY: SURFACE, HYBRID
    USE mo_memory_base,   ONLY: AUTO, new_stream, default_stream_setting,     &
                                add_stream_reference, add_stream_element
!!    USE mo_submodel_diag, ONLY: t_diag_list,                 &
!!                                new_diag_list, new_diag,     & 
!!                                BYTRACER, BYSPECIES, BYMODE
    USE mo_tracer,        ONLY: get_ntrac
    USE mo_species,       ONLY: get_nspec
    USE mo_ham,           ONLY: sizeclass, nccndiag
!>>SF
    USE mo_param_switches, ONLY: ncd_activ, nactivpdf
!<<SF
    USE mo_activ,         ONLY: nw
    USE mo_ham_ccnctl,    ONLY: zsat


    IMPLICIT NONE

    !-- parameters
    INTEGER,          INTENT(in)      :: nclass         ! number of aerosol modes (or bins)
    TYPE (t_stream),  POINTER         :: stream_ham  ! stream reference for internal use

    !-- local variables
    INTEGER           :: jt, jclass, jw, jsat
    CHARACTER(len=10) :: cbin, csat

    !--- executable procedure ---------


    !--- Allocate memory for aerosol properties
    IF (.NOT. ALLOCATED(rdry)) ALLOCATE(rdry(nclass))
    IF (.NOT. ALLOCATED(rwet)) ALLOCATE(rwet(nclass))
    IF (.NOT. ALLOCATED(densaer)) ALLOCATE(densaer(nclass))
!>>SF
    IF (.NOT. ALLOCATED(nact_strat)) ALLOCATE(nact_strat(nclass))
    IF (.NOT. ALLOCATED(nact_conv)) ALLOCATE(nact_conv(nclass))
    IF (.NOT. ALLOCATED(frac)) ALLOCATE(frac(nclass))

    IF (ncd_activ == 2 .OR. nccndiag > 0) THEN
       IF (.NOT. ALLOCATED(a)) ALLOCATE(a(nclass))
       IF (.NOT. ALLOCATED(b)) ALLOCATE(b(nclass))
    ENDIF

    IF (ncd_activ == 2) THEN
       IF (.NOT. ALLOCATED(sc)) ALLOCATE(sc(nclass))
       IF (.NOT. ALLOCATED(rc_strat)) ALLOCATE(rc_strat(nclass,nw))
       IF (.NOT. ALLOCATED(rc_conv)) ALLOCATE(rc_conv(nclass,nw))
    ENDIF
!<<SF

    !--- Create the ham streams
    !>>dod deleted laero_out
    CALL new_stream(stream_ham, 'ham', filetype=trac_filetype)
    CALL default_stream_setting (stream_ham,              &
                                 lrerun    = .FALSE.,      &
                                 laccu     = .FALSE. ,     &
                                 lpost     =  .TRUE.,      &
                                 leveltype =  SURFACE,     &
                                 table     =  199,         &
                                 code      =  AUTO         )

    !>>dod redmine #113 add lpost=T
    CALL add_stream_reference (stream_ham, 'geosp'   ,'g3b',    lpost=.TRUE.)
    CALL add_stream_reference (stream_ham, 'lsp'     ,'sp',     lpost=.TRUE.)
    CALL add_stream_reference (stream_ham, 'aps'     ,'g3b',    lpost=.TRUE.)
    CALL add_stream_reference (stream_ham, 'gboxarea','geoloc', lpost=.TRUE.)
    CALL add_stream_reference (stream_ham, 'aclcv'   ,'g3b',    lpost=.TRUE.)
    !<<dod

    !---------------------------------------------------------------------------------
    !--- Add non accumulated auxiliary fields to ham stream
    CALL default_stream_setting (stream_ham, lpost =.FALSE., lrerun=.TRUE.)

    !--- 4.1) Auxiliary fields for dust emissions:
    CALL add_stream_element (stream_ham, 'daylength',     daylength,       &
                             longname='relative daylength', units='1')

    !--- 4.3) Hydrological parameters:

    CALL default_stream_setting (stream_ham, leveltype = HYBRID, lpost=.TRUE.)

    CALL add_stream_element (stream_ham, 'relhum',        relhum,          &
                             longname='ambient relative humidity', units='%')

    !--- 4.4) Aerosol Properties:

    CALL default_stream_setting (stream_ham, leveltype = HYBRID, lpost=.TRUE.)

    DO jt = 1, nclass
       CALL add_stream_element (stream_ham, 'rdry_'//TRIM(sizeclass(jt)%shortname), rdry(jt)%ptr,     &
                                longname='dry number median radius - '//TRIM(sizeclass(jt)%shortname), &
                                units='m' )
       CALL add_stream_element (stream_ham, 'rwet_'//TRIM(sizeclass(jt)%shortname), rwet(jt)%ptr,     &
                                longname='wet number median radius - '//TRIM(sizeclass(jt)%shortname), &
                                units='m', lpost=sizeclass(jt)%lsoluble )
       !<<dod
       CALL add_stream_element (stream_ham, 'densaer_'//TRIM(sizeclass(jt)%shortname), densaer(jt)%ptr, &
                                longname='wet density'//TRIM(sizeclass(jt)%shortname),                 &
                                units='kg m-3' )

    ENDDO

    !--- Ionisation rate (non accumulated)
    CALL add_stream_element (stream_ham, 'ipr', ipr,                                   &
                             longname='ion pair production rate', units='m-3 s-1',      &
                             laccu=.FALSE., lrerun=.FALSE., leveltype=HYBRID  )


    !---------------------------------------------------------------------------------
    !--- Accumulated diagnostics elements

    CALL default_stream_setting (stream_ham, units     = 'kg m-2 s-1',&
                                       lrerun    = .TRUE. ,     &
                                       laccu     = .TRUE. ,     &
                                       lpost     = .TRUE. ,     &
                                       leveltype = SURFACE,     &
                                       table     = 199,         &
                                       code      = AUTO         )

    !--- 2.2) Mass diagnostics:

    CALL add_stream_element (stream_ham, 'D_PROD_MS4AS',      d_prod_ms4as,                             &
                             longname='sulfate production liquid phase acc',    units='kg(SO4) m-2 s-1' )

    CALL add_stream_element (stream_ham, 'D_PROD_MS4CS',      d_prod_ms4cs,                             &
                             longname='sulfate production liquid phase coarse', units='kg(SO4) m-2 s-1' )

    CALL add_stream_element (stream_ham, 'D_PROD_SO2_DMS_OH', d_prod_so2_dms_oh,                        &
                             longname='sulfur production gas phase via DMS+OH', units='kg(SO2) m-2 s-1' )
    
    CALL add_stream_element (stream_ham, 'D_PROD_SO4_DMS_OH', d_prod_so4_dms_oh,                        &
                             longname='sulfate production gas phase via DMS+OH',units='kg(SO4) m-2 s-1' )

    CALL add_stream_element (stream_ham, 'D_PROD_SO2_DMS_NO3',d_prod_so2_dms_no3,                       &
                             longname='sulfur production gas phase via DMS+NO3',units='kg(SO2) m-2 s-1' )

    CALL add_stream_element (stream_ham, 'D_PROD_SO4_SO2_OH', d_prod_so4_so2_oh,                        &
                             longname='sulfate production gas phase via SO2+OH',units='kg(SO4) m-2 s-1' )

    CALL add_stream_element (stream_ham, 'D_COND_SO4',        d_cond_so4,                               &
                             longname='condensation of sulfate on aerosol',     units='kg(SO4) m-2 s-1' )

    CALL add_stream_element (stream_ham, 'D_NUC_SO4',         d_nuc_so4,                                &
                             longname='nucleation of sulfate',                  units='kg(SO4) m-2 s-1' )

    !--- 3D fields:

    CALL default_stream_setting (stream_ham, leveltype = HYBRID)

    CALL add_stream_element (stream_ham, 'D_PROD_SO4_LIQ',    d_prod_so4_liq,                           &
                             longname='sulfate production liquid phase',        units='kg(SO4) m-3 s-1' )

    CALL add_stream_element (stream_ham, 'D_PROD_H2SO4',    d_prod_h2so4,                           &
                             longname='sulfate production gas phase',           units='kg(SO4) m-3 s-1' )

!>>SF moved this aerosol-dependent stuff out of the activ stream 
!     (as for other quantities above, this is in fact still HAM-dep more than simply aerosol-dep)
    CALL add_stream_element (stream_ham,   'NBCSOL_DIAG',   nbcsol_diag,                         &
                             longname='Number of internally mixed BC particles',units='m-3' )

    CALL add_stream_element (stream_ham,   'NBCINSOL_DIAG',   nbcinsol_diag,                     &
                             longname='Number of externally mixed BC particles',units='m-3' )

    CALL add_stream_element (stream_ham,   'NBCSOL_ACC',   nbcsol_acc,                           &
                            longname='Number of internally mixed accu mode BC',units='m-3' )
  
    CALL add_stream_element (stream_ham,   'NBCSOL_AIT',   nbcsol_ait,                         &
                             longname='Number internally mixed Aitken mode BC',units='m-3'  )

    CALL add_stream_element (stream_ham,   'NDUINSOL_DIAG',   nduinsol_diag,                     &
                             longname='Number of externally mixed dust particles',units='m-3')

    CALL default_stream_setting (stream_ham, laccu=.FALSE.)

    CALL add_stream_element (stream_ham,   'NBCSOL_STRAT', nbcsol_strat,            &
         longname='int mixed BC aerosols from strat. activ. part.', units='m-3')

    CALL add_stream_element (stream_ham,   'NDUSOL_STRAT', ndusol_strat,            &
         longname='int mixed dust aerosols from strat. activ. part.', units='m-3')

    CALL add_stream_element (stream_ham,   'NBCINSOL', nbcinsol,                    &
         longname='ext mixed BC aerosols for frz in stream_ham', units='m-3')

    CALL add_stream_element (stream_ham,   'NDUINSOLAI', nduinsolai,                &
         longname='ext mixed accu. dust AP for frz in stream_ham', units='m-3')

    CALL add_stream_element (stream_ham,   'NDUINSOLCI', nduinsolci,                &
         longname='ext mixed coarse dust AP for frz in stream_ham', units='m-3')

    CALL add_stream_element (stream_ham,   'NAERINSOL', naerinsol,                  &
         longname='total number of insoluble AP for frz in stream_ham', units='m-3')

    CALL default_stream_setting (stream_ham, lrerun=.FALSE.)
  
    DO jclass=1, nclass
       !>>dod #377
       IF (sizeclass(jclass)%lactivation) THEN
          CALL add_stream_element (stream_ham, 'NACT_STRAT_'//TRIM(sizeclass(jclass)%shortname),nact_strat(jclass)%ptr, &
                                longname='number of activated particles stratiform', units='m-3', lrerun=.TRUE.)
  
          CALL add_stream_element (stream_ham, 'NACT_CONV_'//TRIM(sizeclass(jclass)%shortname), nact_conv(jclass)%ptr, &
                                longname='number of activated particles convective', units='m-3', lrerun=.TRUE.)
  
          CALL add_stream_element (stream_ham, 'FRAC_'//TRIM(sizeclass(jclass)%shortname),      frac(jclass)%ptr, &
                                longname='fraction of activated particles',units='%', lrerun=.TRUE.)
       !<<SF #384 
  
          IF (ncd_activ == 2 .OR. nccndiag > 0) THEN
  
             CALL add_stream_element (stream_ham, 'A_'//TRIM(sizeclass(jclass)%shortname),         a(jclass)%ptr, &
                                   longname='curvature parameter',           units='m', lrerun=.TRUE.)

             CALL add_stream_element (stream_ham, 'B_'//TRIM(sizeclass(jclass)%shortname),         b(jclass)%ptr, &
                                   longname='hygroscopicity parameter',      units='m3', lrerun=.TRUE.)
       END IF

       IF (ncd_activ == 2) THEN
  
          CALL add_stream_element (stream_ham, 'SC_'//TRIM(sizeclass(jclass)%shortname),        sc(jclass)%ptr, &
                                   longname='critical supersaturation',      units='1')

             IF (nactivpdf == 0) THEN
                CALL add_stream_element (stream_ham, 'RC_STRAT_'//TRIM(sizeclass(jclass)%shortname), &
                                      rc_strat(jclass,1)%ptr,                                     &
                                      longname='critical radius stratiform',                      &
                                      units='m', lrerun=.TRUE.)
   
                CALL add_stream_element (stream_ham, 'RC_CONV_'//TRIM(sizeclass(jclass)%shortname), &
                                      rc_conv(jclass,1)%ptr,                                     &
                                      longname='critical radius convective',                     &
                                      units='m', lrerun=.TRUE.)
             ELSE IF (nactivpdf < 0) THEN
                DO jw=1,nw
                   WRITE (cbin, "(I2.2)") jw
                   CALL add_stream_element (stream_ham,                                                      &
                                        'RC_STRAT_'//TRIM(sizeclass(jclass)%shortname)//'_'//TRIM(cbin), &
                                        rc_strat(jclass,jw)%ptr,                                         &
                                        longname='critical radius stratiform, vertical velocity bin '//TRIM(cbin), &
                                        units='m', lrerun=.TRUE.)
   
                   CALL add_stream_element (stream_ham,                                                     &
                                        'RC_CONV_'//TRIM(sizeclass(jclass)%shortname)//'_'//TRIM(cbin), &
                                        rc_conv(jclass,jw)%ptr,                                         &
                                        longname='critical radius convective, vertical velocity bin '//TRIM(cbin), &
                                        units='m', lrerun=.TRUE.)
                ENDDO
             END IF
  
          END IF
       END IF
       !<<dod
    ENDDO
   
!<<SF

    !-- Stream entries for mo_ham_ccn (see also #586)
    IF (nccndiag > 0) THEN

      ! Main CCN and CN diagnistics at surface if nccndiag in {1,3,5};
      ! or 3D if nccndiag in {2,4,6})

      CALL default_stream_setting (stream_ham, lrerun = .TRUE., laccu = .FALSE.)
      SELECT CASE (nccndiag)
         CASE (1,3,5)
           CALL default_stream_setting (stream_ham, leveltype = SURFACE)
         CASE (2,4,6)
           CALL default_stream_setting (stream_ham, leveltype = HYBRID)
      END SELECT

      ! CCN diagnostics for each supersaturation
      ! (surface if nccndiag in {1,3,5}; 3D if in {2,4,6})
      DO jsat=1,nsat

        WRITE(csat,'(F7.3)') zsat(jsat)*100.0

        SELECT CASE (nccndiag)
           CASE (1,3,5,7)
             CALL add_stream_element (stream_ham, 'CCN_'//TRIM(ADJUSTL(csat)),              &
                                      ccn_2d(jsat)%ptr, units='m-3',                        &
                                      longname='Cloud Condensation Nuclei at S='//csat//'%' )
           CASE (2,4,6)
             CALL add_stream_element (stream_ham, 'CCN_'//TRIM(ADJUSTL(csat)),              &
                                      ccn_3d(jsat)%ptr, units='m-3',                        &
                                      longname='Cloud Condensation Nuclei at S='//csat//'%' )
        END SELECT

      ENDDO ! jsat

      ! CN diagnostics if nccndiag in {3,4,5,6}
      ! (surface if in {3,5}; 3D if in {4,6})
      SELECT CASE (nccndiag)
         CASE (3,5)
           CALL add_stream_element (stream_ham, 'CN', cn_2d, units='m-3', &
                                    longname='Condensation Nuclei'        )
         CASE (4,6)
           CALL add_stream_element (stream_ham, 'CN', cn_3d, units='m-3', &
                                    longname='Condensation Nuclei'        )
      END SELECT

      ! CCN and CN column burden diagnostics if nccndiag in {5,6}
      SELECT CASE (nccndiag)
         CASE (5,6)
            CALL default_stream_setting (stream_ham, leveltype = SURFACE)
            DO jsat=1,nsat
               WRITE(csat,'(F7.3)') zsat(jsat)*100.0
               CALL add_stream_element (stream_ham, 'CCN_BURDEN_'//TRIM(ADJUSTL(csat)),     &
                                        ccn_burden(jsat)%ptr, units='m-2',                  &
                                        longname='Cloud Condensation Nuclei burden at S='//csat//'%' )
            ENDDO
   	  CASE (7)
            CALL add_stream_element (stream_ham, 'CN_BURDEN', cn_burden, units='m-2', &
                                     longname='Condensation Nuclei Burden'            )
      END SELECT

    END IF ! nccndiag
    !-- end entries for mo_ham_ccn

  END SUBROUTINE new_stream_ham

!--------------------------------------------------------------------------------------------------------------------

  SUBROUTINE new_stream_ham_rad

  ! *new_stream_ham_rad* constructs diagnostics stream rad
  !                      and defines stream elemenst for
  !                      HAM aerosol optical properties
  !
  ! Author:
  ! -------
  ! Philip Stier, MPI-Met, Hamburg          05/2003
  !
  ! Interface:
  ! ----------
  ! *new_stream_ham_rad* is called from *init_subm_memory*
  !                             in *mo_submodel_interface*

  USE mo_memory_base,   ONLY: new_stream, add_stream_element, &
                              default_stream_setting, add_stream_reference, AUTO
  USE mo_linked_list,   ONLY: t_stream, SURFACE, HYBRID
  USE mo_filename,      ONLY: trac_filetype
  USE mo_species,       ONLY: speclist
  USE mo_ham_rad_data,  ONLY: lambda, nradang
  USE mo_ham_rad_data,  ONLY: Nwv_tot, nraddiagwv, iradbeta, Niwv_sw_beta
  USE mo_exception,     ONLY: finish
  USE mo_ham,           ONLY: nclass, nraddiag, sizeclass, nrad, subm_naerospec, subm_aerospec, &
!>>NAJS: optical properties of dry aerosol
                              nraddry,                                                          &
!<<NAJS
!>>NAJS: lidar backscatter
                              nradbeta
!<<NAJS


  IMPLICIT NONE

  INTEGER :: jclass, jn, jwv

  CHARACTER(len=30) :: cwv, cwv2
  CHARACTER(len=18) :: comp_str

  TYPE (t_stream),  POINTER :: rad

  !---executable procedure

  !---allocate arrays for per-compound diagnostics (only diagnosed wvs are actually allocated as stream

  ALLOCATE(tau_comp(subm_naerospec,Nwv_tot))
  ALLOCATE(abs_comp(subm_naerospec,Nwv_tot))

  !--- Create new stream:

  CALL new_stream (rad ,'rad',filetype=trac_filetype)

  !--- 1) Add standard fields for post-processing:

  !>>dod redmine #113 add lpost=T
  CALL add_stream_reference (rad, 'geosp'   ,'g3b',    lpost=.TRUE.)
  CALL add_stream_reference (rad, 'lsp'     ,'sp',     lpost=.TRUE.)
  CALL add_stream_reference (rad, 'aps'     ,'g3b',    lpost=.TRUE.)
  CALL add_stream_reference (rad, 'gboxarea','geoloc', lpost=.TRUE.)
  !<<dod

  CALL default_stream_setting (rad, lrerun   = .TRUE.,   &
                                    laccu    = .FALSE.,  &
                                    lpost    = .TRUE.,   &
                                    table    = 200,      &
                                    code     = AUTO      )
  !<<dod
  !--- Wavelength loop, only requested wavelengths are allocated:

  DO jwv=1, Nwv_tot

     IF (nraddiagwv(jwv)>0) THEN

        WRITE(cwv,'(I6)') INT(lambda(jwv)*1.E9_dp)
        cwv=TRIM(ADJUSTL(cwv))//'nm'

        !--- 2.1) Diagnostics for each mode:

        DO jclass=1, nclass

           IF( nrad(jclass)>0 ) THEN

              CALL default_stream_setting (rad, leveltype = HYBRID )

              CALL add_stream_element(rad, 'TAU_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv,  &
                                      tau_mode(jclass,jwv)%ptr, units='1',     &
                                      longname='Optical thickness '//TRIM(sizeclass(jclass)%shortname)//' '//cwv )

              IF (nraddiagwv(jwv)>1) THEN

                 CALL add_stream_element(rad, 'ABS_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv,   &
                                         abs_mode(jclass,jwv)%ptr, units='1',     &
                                         longname='Absorption Optical thickness '   &
                                         //TRIM(sizeclass(jclass)%shortname)//' '//cwv)

                 CALL default_stream_setting (rad, leveltype = SURFACE)

                 IF (nraddiag>0) THEN

                    ! THANOS
!                    IF (jwv <= Nwv_sw+Nwv_sw_opt) THEN
!                      CALL add_stream_element(rad,'TAU_2D_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv, &
!                                     tau_2d_mode(jclass,jwv)%ptr, units='1',       &
!                                     longname='Optical thickness'//TRIM(sizeclass(jclass)%shortname)//' '//cwv )
!                    END IF

                    CALL add_stream_element(rad, 'SIGMA_2D_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv, &
                                            sigma_2d_mode(jclass,jwv)%ptr,units='m-2', &
                                            longname='Extinction cross section per particle 2D '//cwv, &
                                            contnorest=.TRUE.)

!                    CALL add_stream_element(rad, 'TAU_2D_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv, &
!                                            tau_2d_mode(jclass,jwv)%ptr,units='m-2', &
!                                            longname='Mode optical thickness 2D '//cwv, &
!                                            contnorest=.TRUE.)

                    CALL add_stream_element(rad, 'OMEGA_2D_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv, &
                                            omega_2d_mode(jclass,jwv)%ptr,units='1',   &
                                            longname='Single scattering albedo 2D '//cwv,  &
                                            contnorest=.TRUE.)

                    CALL add_stream_element(rad, 'ASYM_2D_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv, &
                                           asym_2d_mode(jclass,jwv)%ptr, units='1',   &
                                            longname='Asymetry Factor 2D '//cwv,    &
                                            contnorest=.TRUE. )

                    CALL add_stream_element(rad, 'NR_2D_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv, &
                                            nr_2d_mode(jclass,jwv)%ptr,   units='1',   &
                                            longname='Refractive Index - real part 2D '//cwv,   &
                                            contnorest=.TRUE.  )

                    CALL add_stream_element(rad, 'NI_2D_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv,  &
                                            ni_2d_mode(jclass,jwv)%ptr,   units='1',   &
                                            longname='Refractive Index - imaginary part 2D '//cwv, &
                                            contnorest=.TRUE.    )

                 END IF

                 IF (nraddiag==2) THEN

                    CALL default_stream_setting (rad, leveltype = HYBRID )

                    CALL add_stream_element(rad, 'SIGMA_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv,  &
                                            sigma_mode(jclass,jwv)%ptr,   units='cm+2 part-1', &
                                            longname='Extinction cross section per particle '//cwv  )

                    CALL add_stream_element(rad, 'OMEGA_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv,  &
                                            omega_mode(jclass,jwv)%ptr,   units='1',           &
                                            longname='Single scattering albedo '//cwv   )

                    CALL add_stream_element(rad, 'ASYM_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv,  &
                                            asym_mode(jclass,jwv)%ptr,    units='1',           &
                                            longname='Asymetry Factor '//cwv           )

                    CALL add_stream_element(rad, 'NR_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv,  &
                                            nr_mode(jclass,jwv)%ptr,      units='1',           &
                                            longname='Refractive Index - real part '//cwv       )

                    CALL add_stream_element(rad, 'NI_MODE_'//TRIM(sizeclass(jclass)%shortname)//'_'//cwv,  &
                                            ni_mode(jclass,jwv)%ptr,      units='1',           &
                                            longname='Refractive Index - imaginary part '//cwv   )

                 END IF ! nraddiag

              END IF !nraddiagwv>1

           END IF
        ENDDO

        !--- 2.3) 2D fields:

        CALL default_stream_setting (rad, leveltype = SURFACE )

        IF (nraddiagwv(jwv)>0) THEN

           CALL add_stream_element(rad, 'TAU_2D'//'_'//cwv,  tau_2d(jwv)%ptr,     units='1',  &
                                   longname='Optical thickness - total '//cwv    )

        END IF

        IF (nraddiagwv(jwv)>1) THEN

           !--- Diagnostics for each compound:

           DO jn=1, subm_naerospec
              comp_str = TRIM(ADJUSTL(speclist(subm_aerospec(jn))%shortname))

              CALL add_stream_element(rad, 'TAU_COMP_'//TRIM(ADJUSTL(comp_str))//'_'//cwv,  &
                                      tau_comp(jn,jwv)%ptr,     units='1', &
                                      longname='Optical thickness '//TRIM(ADJUSTL(comp_str))//' '//cwv)

              CALL add_stream_element(rad, 'ABS_COMP_'//TRIM(ADJUSTL(comp_str))//'_'//cwv,  &
                                      abs_comp(jn,jwv)%ptr,     units='1', &
                                      longname='Absorption optical thickness '//TRIM(ADJUSTL(comp_str))//' '//cwv )

           ENDDO

           !--- 2.3) 2D total fields:

           CALL add_stream_element(rad, 'ABS_2D'//'_'//cwv,     abs_2d(jwv)%ptr,     units='1',  &
                                   longname='Absorption optical thickness - total '//cwv              )


        END IF !nraddiagwv>1

     END IF !nraddiagwv>0
  ENDDO !jwv

!>>NAJS: optical properties of dry aerosol
  IF (nraddry.GT.0) THEN 
    DO jwv=1,Niwv_sw_dry 
      IF (iraddry(jwv)/=0) THEN

        WRITE(cwv,'(I6)') INT(lambda(iraddry(jwv))*1.E9_dp)
        cwv=TRIM(ADJUSTL(cwv))//'nm'

        CALL default_stream_setting (rad, leveltype = SURFACE )
        CALL add_stream_element(rad, 'TAU_DRY_2D'//'_'//cwv,     tau_dry(jwv)%ptr,     units='1',  &
                                longname='Dry optical thickness '//cwv                            )

        IF (nraddry.GT.1) THEN
          CALL add_stream_element(rad, 'ABS_DRY_2D'//'_'//cwv,     abs_dry(jwv)%ptr,     units='1',  &
                                  longname='Dry absorption optical thickness '//cwv               )
        ENDIF

        IF (nraddry.GT.2) THEN
          CALL default_stream_setting (rad, leveltype = HYBRID )
          CALL add_stream_element(rad, 'ALFA_AER_DRY'//'_'//cwv,     ext_dry(jwv)%ptr,     units='1/m',  &
                                  longname='Dry extinction '//cwv                              )
!>>DN: AeroCom
          CALL add_stream_element(rad, 'ABS_AER_DRY'//'_'//cwv,     abs_lev_dry(jwv)%ptr,     units='1/m',  &
                                  longname='Dry absorption '//cwv                              )
!<<DN
        ENDIF

      ENDIF
    ENDDO ! jwv
  ENDIF
!<<NAJS


  CALL default_stream_setting (rad, leveltype = SURFACE )

  IF (nradang(1)/=0 .AND. nradang(2)/=0) THEN

     IF (nraddiagwv(nradang(1))>0 .AND. nraddiagwv(nradang(2))>0 ) THEN

        WRITE(cwv,'(I6)') INT(lambda(nradang(1))*1.E9_dp)
        cwv=TRIM(ADJUSTL(cwv))//'nm'

        WRITE(cwv2,'(I6)') INT(lambda(nradang(2))*1.E9_dp)
        cwv2=TRIM(ADJUSTL(cwv2))//'nm'

        CALL add_stream_element(rad, 'ANG_'//TRIM(ADJUSTL(cwv))//'_'//TRIM(ADJUSTL(cwv2)), ang, units='1', &
                               longname='Angstroem parameter between '//TRIM(ADJUSTL(cwv))//' and '//TRIM(ADJUSTL(cwv2)) )

        CALL add_stream_element(rad, 'AI_'//TRIM(ADJUSTL(cwv))//'_'//TRIM(ADJUSTL(cwv2)), ai, units='1', &
                               longname='Aerosol index corresponding to '//TRIM(ADJUSTL(cwv))//' and '//TRIM(ADJUSTL(cwv2)) )        

     ELSE
        CALL finish('construct_stream_ham_rad:','Angstroem parameter between '//TRIM(ADJUSTL(cwv))  &
                    //' and '//TRIM(ADJUSTL(cwv2))//' requested but inconsistent nraddiagwv')
     END IF

  END IF

!>>NAJS: lidar backscatter
  !--- 3) Aerosol and gaseous backscatter:
  IF (nradbeta.GT.0) THEN
    CALL default_stream_setting (rad, leveltype = HYBRID )
    DO jwv=1,Niwv_sw_beta 
      IF (iradbeta(jwv)/=0) THEN
  
        WRITE(cwv,'(I6)') INT(lambda(iradbeta(jwv))*1.E9_dp)
        cwv=TRIM(ADJUSTL(cwv))//'nm'
    
        CALL add_stream_element(rad, 'ALFA_AER'//'_'//cwv, alfa_aer(jwv)%ptr, units='1/m', &
                                longname='Aerosol extinction at '//TRIM(ADJUSTL(cwv)) )
        IF (nradbeta.GT.1) THEN
          CALL add_stream_element(rad, 'BETA_TOT'//'_'//cwv, beta_tot(jwv)%ptr, units='1/m/sr', &
                                  longname='Aerosol & gas backscatter at '//TRIM(ADJUSTL(cwv)) )
          CALL add_stream_element(rad, 'ALFA_GAS'//'_'//cwv, alfa_gas(jwv)%ptr, units='1/m', &
                                  longname='Gaseous extinction at '//TRIM(ADJUSTL(cwv)) )
          CALL add_stream_element(rad, 'LIDAR_RATIO_AER'//'_'//cwv, lidar_ratio_aer(jwv)%ptr, units='sr', &
                                  longname='Aerosol LIDAR ratio at '//TRIM(ADJUSTL(cwv)) )
        ENDIF ! nradbeta.GT.1
        IF (nradbeta.GT.2) THEN
          CALL add_stream_element(rad, 'ATTBETA_TOT_SPACE'//'_'//cwv, attbeta_tot_space(jwv)%ptr, units='1/m/sr', &
                                  longname='Attenuated aerosol & gas backscatter from space at '//TRIM(ADJUSTL(cwv)) )
          CALL add_stream_element(rad, 'ATTBETA_TOT_EARTH'//'_'//cwv, attbeta_tot_earth(jwv)%ptr, units='1/m/sr', &
                                  longname='Attenuated aerosol & gas backscatter from Earth at '//TRIM(ADJUSTL(cwv)) )
        ENDIF ! nradbeta.GT.2
 
      ENDIF ! iradbeta(jwv).GT.0
    ENDDO ! jwv
  ENDIF ! nradbeta.GT.0
!<<NAJS

  END SUBROUTINE new_stream_ham_rad



END MODULE mo_ham_streams
