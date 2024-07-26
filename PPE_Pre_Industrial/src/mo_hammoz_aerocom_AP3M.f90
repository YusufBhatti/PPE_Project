!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_AP3M.f90
!!
!! \brief
!! Module for AeroCom Phase 3 experiments monthly (AP3M) diagnostics
!!
!! \author D. Neubauer (ETH Zurich)
!!  adapted from AerChemMIP diagnostics from S. Ferrachat
!!
!! D. Neubauer, david.neubauer@env.ethz.ch
!!
!! \revision_history
!!   -# D. Neubauer (ETH Zurich) - original code (2019-03-12)
!!
!! \limitations
!! None
!!
!! \details
!! None
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
MODULE mo_hammoz_aerocom_AP3M

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  
  IMPLICIT NONE  

  PRIVATE

  PUBLIC :: construct_AP3M_stream
  PUBLIC :: init_AP3M
  PUBLIC :: update_AP3M_diags

  !-- Define data structures to handle the PMxx calcs
  INTEGER, PARAMETER :: npm = 3 ! number of PMxx calculations

  REAL(dp), PARAMETER :: dcrit_pm(npm) = (/ 1.e-6_dp, 2.5e-6_dp, 10.e-6_dp /) ! critical diams
                                                                        ! for PMxx calculations
  TYPE t_class
      INTEGER :: iclass                      !< aerosol class index (== M7 mode)
      REAL(dp), ALLOCATABLE :: radius(:,:,:) !< class radius (may be dry or wet) [m]
      REAL(dp), ALLOCATABLE :: frac(:,:,:)   !< distribution fraction above critical diam [1]
  END TYPE t_class

  TYPE t_pm
      REAL(dp) :: rcrit                       ! critical radius [m]
      INTEGER :: nclass                       ! number of relevant classes
      TYPE(t_class), ALLOCATABLE :: lclass(:) ! class-related infos
      LOGICAL :: ldry                         ! dry vs wet
      REAL(dp), POINTER :: ptr(:,:,:)       ! pointer to corresponding stream var for final diags
  END TYPE t_pm

  TYPE(t_pm) :: pm_info(npm) ! holds data for the PMxx calculations

  PUBLIC :: t_pm, t_class ! needed for PM calcs in other modules

  TYPE (t_stream), PUBLIC, POINTER :: acp3m
  
  REAL(dp), PUBLIC, POINTER :: co(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rn222(:,:,:)
  REAL(dp), PUBLIC, POINTER :: pb210(:,:,:)
  REAL(dp), PUBLIC, POINTER :: loadlt1du(:,:)
  REAL(dp), PUBLIC, POINTER :: loadlt25du(:,:)
  REAL(dp), PUBLIC, POINTER :: loadlt10du(:,:)
  REAL(dp), PUBLIC, POINTER :: loadlt1ss(:,:)
  REAL(dp), PUBLIC, POINTER :: loadlt25ss(:,:)
  REAL(dp), PUBLIC, POINTER :: mmr_class(:,:,:,:)
  REAL(dp), PUBLIC, POINTER :: mmr_classss(:,:,:,:)
  REAL(dp), PUBLIC, POINTER :: mmr_classdu(:,:,:,:)
  REAL(dp), PUBLIC, POINTER :: mmrpm10(:,:,:)
  REAL(dp), PUBLIC, POINTER :: mmrpm2p5(:,:,:)
  REAL(dp), PUBLIC, POINTER :: mmrpm1(:,:,:)
  REAL(dp), PUBLIC, POINTER :: od550csaer(:,:)
  REAL(dp), PUBLIC, POINTER :: od550csaerfreq(:,:)
  REAL(dp), PUBLIC, POINTER :: od550lt1aer(:,:)
  REAL(dp), PUBLIC, POINTER :: abs550lt1aer(:,:)
  REAL(dp), PUBLIC, POINTER :: nh50(:,:,:)
  REAL(dp), PUBLIC, POINTER :: aoanh(:,:,:)
  REAL(dp), PUBLIC, POINTER :: rldscs(:,:)
  REAL(dp), PUBLIC, POINTER :: rluscs(:,:)

  INTEGER :: idt_pb210  ! tracer index for Lead 210
  REAL(dp), PARAMETER :: co_decay = -1._dp / 50._dp / 86400._dp !< decay constant for CO corresponding to
                                                              !! a lifetime of 50 days, expressed in [s-1]
  REAL(dp), PARAMETER :: rn222_decay = -2.11E-6_dp !< decay constant for Radon 2222, expressed in [s-1]
  !-- To handle age of northern hemisphere air:
  INTEGER :: idt_aoanh !< tracer index for age of northern hemisphere air

  REAL(dp), PARAMETER :: sec_per_year = 3600._dp * 24._dp * 365.25 !< seconds to year conversion
                                                                   !! [s yr-1]
  REAL(dp), PARAMETER :: lat_emi_min = 30._dp
  REAL(dp), PARAMETER :: lat_emi_max = 50._dp

  LOGICAL, ALLOCATABLE  :: emi_mask(:,:) !< mask for defining the emission region of 
                                         !! the northern hemisphere age of air tracer

   !-- To handle the artificial tracer with 50 day lifetime                                       
  INTEGER :: idt_nh50  !< tracer index for nh50 tracer
  
  REAL(dp), PARAMETER :: c_0 = 100.e-9_dp !< Initial concentration for artificial tracer nh50
                                          !! [mol mol-1]
  REAL(dp), PARAMETER :: c_decay = -1._dp / 50._dp / 86400._dp !< decay constant corresponding to
                                                              !! a lifetime of 50 days, expressed in [s-1]

  CONTAINS

  !------------------------------------------------
  SUBROUTINE init_AP3M
    
    USE mo_ham_m7ctl,     ONLY: icoas, icoai
    USE mo_decomposition, ONLY: dc => local_decomposition, global_decomposition
    USE mo_gaussgrid, ONLY: philat
    USE mo_transpose, ONLY: scatter_gp
    USE mo_mpi,       ONLY: p_parallel_io
    USE mo_tracer,    ONLY: new_tracer, get_tracer
    USE mo_tracdef,   ONLY: RESTART, CONSTANT, INITIAL, ON, OFF, ntrac
    
    INTEGER :: ipm, jc, ierr, jlon

     !SFnote: the following 2 vars are REALs instead of LOGICALs 
     !        because scatter_gp does not handle logicals
     REAL(dp), POINTER     :: zmask_glo(:,:) !< tracer emission mask (global)
     REAL(dp), ALLOCATABLE :: zmask_loc(:,:) !< tracer emission mask (local)

     LOGICAL, ALLOCATABLE :: lo_glo(:,:)

    !-- Populate pm_info data struc:
    DO ipm=1,npm
       pm_info(ipm) %rcrit = dcrit_pm(ipm) * 0.5_dp ! diam -> radius
       
       pm_info(ipm) %nclass = 2
       ALLOCATE(pm_info(ipm) %lclass( pm_info(ipm) %nclass ))
       pm_info(ipm) %lclass(1) %iclass = icoas
       pm_info(ipm) %lclass(2) %iclass = icoai
       
       IF (dcrit_pm(ipm) == 10.e-6_dp) THEN
          pm_info(ipm) %ldry = .TRUE.
       ELSE
          pm_info(ipm) %ldry = .TRUE.
       ENDIF
       
       DO jc=1,pm_info(ipm) %nclass
          ALLOCATE(pm_info(ipm) %lclass(jc) %radius(dc%nproma, dc%nlev, dc%ngpblks))
          ALLOCATE(pm_info(ipm) %lclass(jc) %frac(dc%nproma, dc%nlev, dc%ngpblks))
       ENDDO
       
    ENDDO
         
    !-- Define a tracer for PB210
     CALL new_tracer('PB210', 'acp3m', ierr=ierr &
          ,units='kg kg-1' &
          ,nwrite=OFF &
          ,ninit=RESTART+CONSTANT &
          ,longname='Lead 210' &
          ,idx = idt_pb210 &
          )
     
     !-- Define an emission mask for the above tracer
     !   (to emit the tracer only in a given lat band)
     IF (p_parallel_io) THEN
         ALLOCATE (zmask_glo(dc%nlon,dc%nlat))
         ALLOCATE (lo_glo(dc%nlon,dc%nlat))

         DO jlon=1,dc%nlon
             lo_glo(jlon,1:dc%nlat) = (&
                       philat(1:dc%nlat) >= lat_emi_min &
                 .AND. philat(1:dc%nlat) <= lat_emi_max)
         ENDDO

         zmask_glo(:,:) = MERGE(1._dp, 0._dp, lo_glo(:,:))
     ENDIF

     ALLOCATE(zmask_loc(dc%nproma,dc% ngpblks))
     ALLOCATE(emi_mask(dc%nproma,dc% ngpblks))

     CALL scatter_gp(zmask_glo,zmask_loc,global_decomposition)

     emi_mask(:,:) = (zmask_loc(:,:) == 1._dp)

     IF (p_parallel_io) THEN
         DEALLOCATE(zmask_glo, zmask_loc, lo_glo)
     END IF

     !-- Define a passive tracer to track the ideal age of northern hemisphere air
     CALL new_tracer('aoanh', 'AP3', ierr=ierr &
         ,units='yr' &
         ,nwrite=OFF &
         ,ninit=RESTART+CONSTANT &
         ,longname='mean age of stratospheric air' &
         ,idx = idt_aoanh &
         )

     !-- Define a passive tracer with a fixed decay
     CALL new_tracer('nh50', 'AP3', ierr=ierr &
         ,units='mol mol-1' &
         ,nwrite=OFF &
         ,ninit=RESTART+CONSTANT &
         ,longname='Artificial tracer with 50 day lifetime' &
         ,idx = idt_nh50 &
         )

  END SUBROUTINE init_AP3M
  !------------------------------------------------

  SUBROUTINE construct_AP3M_stream

    USE mo_hammoz_aerocom_data, ONLY: AEROCOM_FILETYPE 
    USE mo_memory_base,         ONLY: new_stream, add_stream_element, AUTO,  &
                                    default_stream_setting, add_stream_reference
    USE mo_time_event,          ONLY: io_time_event
    USE mo_ham,                 ONLY: nclass
    USE mo_exception,           ONLY: message, em_error

    TYPE(io_time_event) :: put_interval
    INTEGER :: ipm
    
    !-- set output interval
    put_interval%counter      = 1
    put_interval%unit         = 'months'
    put_interval%adjustment   = 'last'
    put_interval%offset       = 0
    
    !-- Create new stream:
    CALL new_stream (acp3m ,'acp3m', &
                     filetype = AEROCOM_FILETYPE, &
                     lrerun = .TRUE., &
                     interval = put_interval, &
                     lpost = .TRUE.)
  
    !-- Add standard fields for post-processing:
    CALL default_stream_setting (acp3m, &
                                 table = 199, &
                                 code = AUTO )

    !-- Basic vars (may sometimes be necessary)
    CALL add_stream_reference (acp3m, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (acp3m, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (acp3m, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (acp3m, 'gboxarea','geoloc',lpost=.TRUE.)
  
    CALL add_stream_element (acp3m,'CO',co, &
         longname='Carbon monoxide - transport tracer', &
         units='kg kg-1',   &
         laccu = .TRUE.,&
         lpost = .TRUE., &
         lrerun = .TRUE. )
    
    CALL add_stream_element (acp3m,'RN222',rn222, &
         longname='Radon 222', &
         units='kg kg-1',   &
         laccu = .TRUE.,&
         lpost = .TRUE., &
         lrerun = .TRUE. )
    
    CALL add_stream_element (acp3m,'PB210',pb210, &
         longname='Lead 210 - deposition tracer', &
         units='kg kg-1',   &
         laccu = .TRUE.,&
         lpost = .TRUE., &
         lrerun = .TRUE. )
    
    CALL add_stream_element (acp3m, 'loadlt1du', loadlt1du, &
        longname = 'Column pm1 dust mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

    CALL add_stream_element (acp3m, 'loadlt25du', loadlt25du, &
        longname = 'Column pm2p5 dust mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

    CALL add_stream_element (acp3m, 'loadlt10du', loadlt10du, &
        longname = 'Column pm10 dust mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

    CALL add_stream_element (acp3m, 'loadlt1ss', loadlt1ss, &
        longname = 'Column pm1 seasalt mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

    CALL add_stream_element (acp3m, 'loadlt25ss', loadlt25ss, &
        longname = 'Column pm2p5 seasalt mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .FALSE. )

    CALL add_stream_element (acp3m, 'od550csaer', od550csaer, &
        longname = 'ambient aerosol optical thickness at 550 nm', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3m, 'od550csaerfreq', od550csaerfreq, &
        longname = 'frequency of ambient aerosol optical thickness at 550 nm', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3m, 'od550lt1aer', od550lt1aer, &
        longname = 'ambient fine mode aerosol optical thickness at 550 nm', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3m, 'abs550lt1aer', abs550lt1aer, &
        longname = 'ambient fine mode aerosol absorption optical thickness at 550 nm', &
        units = '1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3m, 'mmrpm10', mmrpm10, &
        longname = 'PM10 mass mixing ratio', &
        units = 'kg kg-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3m, 'mmrpm2p5', mmrpm2p5, &
        longname = 'PM2.5 mass mixing ratio', &
        units = 'kg kg-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3m, 'mmrpm1', mmrpm1, &
        longname = 'PM1.0 mass mixing ratio', &
        units = 'kg kg-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3m, 'nh50', nh50, &
        longname = 'Artificial tracer with 50 day lifetime', &
        units = 'mol mol-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3m, 'aoanh', aoanh, &
        longname = 'Tracer age of air Northern Hemisphere', &
        units = 'yr', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3m, 'rldscs', rldscs, &
        longname = 'Surface downwelling clear-sky longwave radiation', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acp3m, 'rluscs', rluscs, &
        longname = 'Surface upwelling clear-sky longwave radiation', &
        units = 'W m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    !-- Save per-class mmr for PM calcs in other modules
    CALL add_stream_element (acp3m, 'mmr_class', mmr_class, &
        longname = 'mass mixing ratio per aerosol class', &
        units = 'kg kg-1', &
        ktrac = nclass, &
        laccu = .FALSE., &
        lpost = .FALSE., &
        lrerun = .FALSE. )

    CALL add_stream_element (acp3m, 'mmr_classdu', mmr_classdu, &
        longname = 'Dust mass mixing ratio per aerosol class', &
        units = 'kg kg-1', &
        ktrac = nclass, &
        laccu = .FALSE., &
        lpost = .FALSE., &
        lrerun = .FALSE. )

    CALL add_stream_element (acp3m, 'mmr_classss', mmr_classss, &
        longname = 'Seasalt mixing ratio per aerosol class', &
        units = 'kg kg-1', &
        ktrac = nclass, &
        laccu = .FALSE., &
        lpost = .FALSE., &
        lrerun = .FALSE. )

    !-- For the PMxx calcs
    DO ipm=1,npm

        IF (dcrit_pm(ipm) == 1.e-6_dp) THEN
            pm_info(ipm) %ptr => mmrpm1
        ELSEIF (dcrit_pm(ipm) == 2.5e-6_dp) THEN
            pm_info(ipm) %ptr => mmrpm2p5
        ELSEIF (dcrit_pm(ipm) == 10.e-6_dp) THEN
            pm_info(ipm) %ptr => mmrpm10
        ELSE
            CALL message('construct_AP3M_stream','Unknown critical diameter',level=em_error)
         ENDIF

      END DO

  END SUBROUTINE construct_AP3M_stream

  SUBROUTINE update_AP3M_diags(kproma, kbdim, klev, krow)

    USE mo_time_control, ONLY: delta_time
    USE mo_memory_g1a,   ONLY: xtm1
    USE mo_scan_buffer,  ONLY: xtte
    USE mo_ham_m7_trac,  ONLY: idt_co, idt_rn222
    USE mo_hammoz_aerocom_HEaci, ONLY: dpg, clt_inst
    USE mo_ham,          ONLY: nclass, naerocomp, aerocomp, nrad
    USE mo_ham_rad_data, ONLY: Nwv_sw
    USE mo_ham_species,  ONLY: id_du, id_ss
    USE mo_ham_streams,  ONLY: abs_2d, tau_2d, tau_comp, rwet, rdry, tau_mode, abs_mode
    USE mo_ham_tools,    ONLY: ham_m7_logtail
    USE mo_memory_cfdiag,ONLY: irlucs, irldcs

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow

    INTEGER :: iclass,it, ipm, jc, jn, ispec, jt, jk
    REAL(dp) :: zmmr_instdu(kbdim,klev),zmmr_instss(kbdim,klev),zmmr_inst(kbdim,klev)
    REAL(dp) :: zfrac(kbdim, klev, nclass) !< fraction of distribution relevant for mmr sum
                                           !! in the PMxx calcs
    REAL(dp) :: zrcrit(kbdim,klev)         !< critical radius (ham_m7_logtail
                                           !! requires a 2D var...)
    REAL(dp) :: odlt1_inst(kbdim) !< instantaneous version of od550lt1aer
    REAL(dp) :: abslt1_inst(kbdim) !< instantaneous version of abs550lt1aer

    !-- aoanh
    xtm1(1:kproma,:,idt_aoanh,krow) = xtm1(1:kproma,:,idt_aoanh,krow) &
        + delta_time / sec_per_year
    xtm1(1:kproma,klev,idt_aoanh,krow) = MERGE(0._dp, xtm1(1:kproma,klev,idt_aoanh,krow), &
        emi_mask(1:kproma,krow))
    xtte(1:kproma,1:klev,idt_aoanh,krow) = MERGE(0._dp, xtte(1:kproma,1:klev,idt_aoanh,krow), &
        xtm1(1:kproma,1:klev,idt_aoanh,krow)<0._dp)
    xtm1(1:kproma,1:klev,idt_aoanh,krow) = MERGE(0._dp, xtm1(1:kproma,1:klev,idt_aoanh,krow), &
        xtm1(1:kproma,1:klev,idt_aoanh,krow)<0._dp)
    xtte(1:kproma,klev,idt_aoanh,krow) = 0._dp

    aoanh(1:kproma,:,krow) = aoanh(1:kproma,:,krow) &
        + xtm1(1:kproma,:,idt_aoanh,krow) * delta_time

    !-- nh50
    xtm1(1:kproma,:,idt_nh50,krow) = xtm1(1:kproma,:,idt_nh50,krow) &
        * EXP(delta_time * c_decay)
    xtm1(1:kproma,klev,idt_nh50,krow) = MERGE(c_0, xtm1(1:kproma,klev,idt_nh50,krow), &
        emi_mask(1:kproma,krow))
    xtte(1:kproma,1:klev,idt_nh50,krow) = MERGE(0._dp, xtte(1:kproma,1:klev,idt_nh50,krow), &
         xtm1(1:kproma,1:klev,idt_nh50,krow)<0._dp)
    xtm1(1:kproma,1:klev,idt_nh50,krow) = MERGE(0._dp, xtm1(1:kproma,1:klev,idt_nh50,krow), &
         xtm1(1:kproma,1:klev,idt_nh50,krow)<0._dp)
    xtm1(1:kproma,1:klev,idt_nh50,krow) = MERGE(0._dp, xtm1(1:kproma,1:klev,idt_nh50,krow), &
         xtm1(1:kproma,1:klev,idt_nh50,krow)>1._dp)
    xtte(1:kproma,klev,idt_nh50,krow) = 0._dp

    nh50(1:kproma,:,krow) = nh50(1:kproma,:,krow) &
        + xtm1(1:kproma,:,idt_nh50,krow) * delta_time

    !-- Carbon Monoxide
    xtm1(1:kproma,:,idt_co,krow) = xtm1(1:kproma,:,idt_co,krow) &
         * EXP(delta_time * co_decay)
    xtte(1:kproma,1:klev,idt_co,krow) = MERGE(0._dp, xtte(1:kproma,1:klev,idt_co,krow), &
         xtm1(1:kproma,1:klev,idt_co,krow)<0._dp)
    xtm1(1:kproma,1:klev,idt_co,krow) = MERGE(0._dp, xtm1(1:kproma,1:klev,idt_co,krow), &
         xtm1(1:kproma,1:klev,idt_co,krow)<0._dp)

    co(1:kproma,:,krow) = co(1:kproma,:,krow) &
         + xtm1(1:kproma,:,idt_co,krow) * delta_time

    !-- Lead 210
    xtm1(1:kproma,:,idt_pb210,krow) = xtm1(1:kproma,:,idt_pb210,krow) + &
         xtm1(1:kproma,:,idt_rn222,krow) * (1._dp - EXP(delta_time * rn222_decay))
    xtte(1:kproma,1:klev,idt_pb210,krow) = MERGE(0._dp, xtte(1:kproma,1:klev,idt_pb210,krow), &
         xtm1(1:kproma,1:klev,idt_pb210,krow)<0._dp)
    xtm1(1:kproma,1:klev,idt_pb210,krow) = MERGE(0._dp, xtm1(1:kproma,1:klev,idt_pb210,krow), &
         xtm1(1:kproma,1:klev,idt_pb210,krow)<0._dp)

    pb210(1:kproma,:,krow) = pb210(1:kproma,:,krow) &
         + xtm1(1:kproma,:,idt_pb210,krow) * delta_time

    !-- Radon 222
    xtm1(1:kproma,:,idt_rn222,krow) = xtm1(1:kproma,:,idt_rn222,krow) &
         * EXP(delta_time * rn222_decay)
    xtte(1:kproma,1:klev,idt_rn222,krow) = MERGE(0._dp, xtte(1:kproma,1:klev,idt_rn222,krow), &
         xtm1(1:kproma,1:klev,idt_rn222,krow)<0._dp)
    xtm1(1:kproma,1:klev,idt_rn222,krow) = MERGE(0._dp, xtm1(1:kproma,1:klev,idt_rn222,krow), &
         xtm1(1:kproma,1:klev,idt_rn222,krow)<0._dp)

    rn222(1:kproma,:,krow) = rn222(1:kproma,:,krow) &
         + xtm1(1:kproma,:,idt_rn222,krow) * delta_time
    
    !-- mmrpm1, mmrpm2.5, mmrpm10

    !--- Compute mass mixing ratio for each mode 
    mmr_class(1:kproma,:,:,krow) = 0._dp
    mmr_classss(1:kproma,:,:,krow) = 0._dp
    mmr_classdu(1:kproma,:,:,krow) = 0._dp
    DO jn = 1,naerocomp
       iclass = aerocomp(jn)%iclass
       ispec = aerocomp(jn)%spid
       jt = aerocomp(jn)%idt
       mmr_class(1:kproma,:,iclass,krow) = mmr_class(1:kproma,:,iclass,krow) &
           + xtm1(1:kproma,:,jt,krow)
       IF (ispec == id_ss) THEN
          mmr_classss(1:kproma,:,iclass,krow) = mmr_classss(1:kproma,:,iclass,krow) &
               + xtm1(1:kproma,:,jt,krow)
       ELSE IF (ispec == id_du) THEN
          mmr_classdu(1:kproma,:,iclass,krow) = mmr_classdu(1:kproma,:,iclass,krow) &
               + xtm1(1:kproma,:,jt,krow)
       END IF
       zfrac(1:kproma,:,iclass) = 1._dp ! will be updated for relevant classes
    END DO
 
    !---
    DO ipm=1,npm

        DO jc=1,pm_info(ipm) %nclass

            !--- Set radii distribs
            iclass = pm_info(ipm) %lclass(jc) %iclass
            pm_info(ipm) %lclass(jc) %radius(:,:,krow) = &
                MERGE(rdry(iclass)%ptr(:,:,krow), &
                      rwet(iclass)%ptr(:,:,krow), &
                      pm_info(ipm) %ldry)

            !--- Compute fraction of distrib above threshold
            zrcrit(1:kproma,:) = pm_info(ipm) %rcrit
            CALL ham_m7_logtail(kproma, kbdim, klev, krow, iclass, .FALSE., &
                pm_info(ipm) %lclass(jc) %radius(:,:,krow), &
                zrcrit, &
                pm_info(ipm) %lclass(jc) %frac(:,:,krow) )

            !--- Fraction below threshold
            zfrac(1:kproma,:,iclass) = (1._dp - pm_info(ipm) %lclass(jc) %frac(1:kproma,:,krow))
        ENDDO

        !-- Sum up weighted mmr's to get the instantaneous PMxx mmr
        zmmr_inst(1:kproma,:) = 0._dp
        zmmr_instss(1:kproma,:) = 0._dp
        zmmr_instdu(1:kproma,:) = 0._dp
        DO jc=1,nclass
           zmmr_inst(1:kproma,:) = zmmr_inst(1:kproma,:) + &
                mmr_class(1:kproma,:,jc,krow) * zfrac(1:kproma,:,jc)
           zmmr_instss(1:kproma,:) = zmmr_instss(1:kproma,:) + &
                mmr_classss(1:kproma,:,jc,krow) * zfrac(1:kproma,:,jc)
           zmmr_instdu(1:kproma,:) = zmmr_instdu(1:kproma,:) + &
                mmr_classdu(1:kproma,:,jc,krow) * zfrac(1:kproma,:,jc)
        ENDDO

          !-- Diagnose accumulated PMxx
        pm_info(ipm) %ptr(1:kproma,:,krow) = pm_info(ipm) %ptr(1:kproma,:,krow) + &
             zmmr_inst(1:kproma,:) * delta_time 
        IF (dcrit_pm(ipm) == 1.0e-6_dp) THEN
           DO jk=1,klev
              loadlt1ss(1:kproma,krow) = loadlt1ss(1:kproma,krow) + & 
                   zmmr_instss(1:kproma,jk) * dpg(1:kproma,jk,krow) * delta_time 
              loadlt1du(1:kproma,krow) = loadlt1du(1:kproma,krow) + & 
                   zmmr_instdu(1:kproma,jk) * dpg(1:kproma,jk,krow) * delta_time
           END DO      
        ELSE IF (dcrit_pm(ipm) == 2.5e-6_dp) THEN
           DO jk=1,klev
              loadlt25ss(1:kproma,krow) = loadlt25ss(1:kproma,krow) + & 
                   zmmr_instss(1:kproma,jk) * dpg(1:kproma,jk,krow) * delta_time 
              loadlt25du(1:kproma,krow) = loadlt25du(1:kproma,krow) + & 
                   zmmr_instdu(1:kproma,jk) * dpg(1:kproma,jk,krow) * delta_time
           END DO
        ELSE IF (dcrit_pm(ipm) == 10.0e-6_dp) THEN
           DO jk=1,klev
              loadlt10du(1:kproma,krow) = loadlt10du(1:kproma,krow) + & 
                   zmmr_instdu(1:kproma,jk) * dpg(1:kproma,jk,krow) * delta_time
           END DO
        END IF
    ENDDO
    
    !-- od550csaer
    od550csaer(1:kproma,krow) = od550csaer(1:kproma,krow) + &
         MERGE(tau_2d(Nwv_sw+1)%ptr(1:kproma,krow) * &
         delta_time,0._dp, clt_inst(1:kproma,krow)/=1._dp)

    !-- divide od550csaer by od550csaerfreq to get od550csaer in (partly) clear-sky conditions only
    od550csaerfreq(1:kproma,krow) = od550csaerfreq(1:kproma,krow) + &
         MERGE(delta_time,0._dp,clt_inst(1:kproma,krow)/=1._dp)

    !-- od550lt1aer

    !--- Initialize fraction of numb. conc. for each mode 
    DO jn = 1,naerocomp
       iclass = aerocomp(jn)%iclass
       zfrac(1:kproma,:,iclass) = 1._dp
    END DO

    !--- Updates zfrac for relevant classes
    DO jc=1,pm_info(1) %nclass

        !--- Set radii distribs
        iclass = pm_info(1) %lclass(jc) %iclass
        pm_info(1) %lclass(jc) %radius(:,:,krow) = &
            MERGE(rwet(iclass)%ptr(:,:,krow), &
                  rdry(iclass)%ptr(:,:,krow), &
                  .TRUE.)

        !--- Compute fraction of distrib above threshold
        zrcrit(1:kproma,:) = pm_info(1) %rcrit
        CALL ham_m7_logtail(kproma, kbdim, klev, krow, iclass, .TRUE., &
            pm_info(1) %lclass(jc) %radius(:,:,krow), &
            zrcrit, &
            pm_info(1) %lclass(jc) %frac(:,:,krow) )

        !--- Fraction below threshold
        zfrac(1:kproma,:,iclass) = (1._dp - pm_info(1) %lclass(jc) %frac(1:kproma,:,krow))
    ENDDO

    !-- Sum up weighted optical thicknesses per class / per level to get the instantaneous od
    ! for part. smaller than 1 um
    odlt1_inst(1:kproma) = 0._dp
    abslt1_inst(1:kproma) = 0._dp
    DO jc=1,nclass
        IF( nrad(jc)>0 ) THEN
            DO jk=1,klev
                odlt1_inst(1:kproma) = odlt1_inst(1:kproma) + &
                    tau_mode(jc,Nwv_sw+1)%ptr(1:kproma,jk,krow) * zfrac(1:kproma,jk,jc)
                abslt1_inst(1:kproma) = abslt1_inst(1:kproma) + &
                    abs_mode(jc,Nwv_sw+1)%ptr(1:kproma,jk,krow) * zfrac(1:kproma,jk,jc)
            ENDDO
        ENDIF
    ENDDO

    !-- Diagnose accumulated odlt1 
    od550lt1aer(1:kproma,krow) = od550lt1aer(1:kproma,krow) + &
        odlt1_inst(1:kproma) * delta_time
    
    !-- abs550aer
    abs550lt1aer(1:kproma,krow) = abs550lt1aer(1:kproma,krow) &
         + abslt1_inst(1:kproma) * delta_time
    
    !-- rluscs,rldscs
    rluscs(1:kproma,krow) = rluscs(1:kproma,krow) + delta_time*irlucs(1:kproma,klev+1,krow)
    rldscs(1:kproma,krow) = rldscs(1:kproma,krow) + delta_time*irldcs(1:kproma,klev+1,krow)

  END SUBROUTINE update_AP3M_diags

END MODULE mo_hammoz_aerocom_AP3M
