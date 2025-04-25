!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_cosp_simulator
 !
 ! Purpose:
 ! --------
 ! Model output for comparisons with CloudSat, CALIPSO, and ISCCP satellite data
 !
 ! The COSP simulator is (c) copyright by the British Crown / Met Office 2008
 ! See http://cfmip.metoffice.com/cosp/cosp.v0.3/  and
 ! Refer to Met_Office_licencse.text for details
 !
 ! Version v0.1.beta CICCS implemented by Christine Nam, July    2008
 ! Version v0.3.beta COSP - Lidar implemented by C.Nam, October 2008
 ! Version v0.3.beta COSP - Radar implemented by C.Nam, January 2009
 ! Version v1.0 COSP - Lidar updated by C.Nam, April 2009
 ! Version v1.1 COSP - Lidar updated by C.Nam, May 2009
 ! Version v1.1 COSP - Radar updated by C.Nam, October 2009
 ! Version v1.2.1 COSP - Lidar & Radar updat implemented by C.Nam, Febuary 2010
 ! Zak Kipling: bugfix for allowing COSP output on 40 levels. (#474)
 !
 !  m300111:
 ! Initial implementation into ECHAM6,  update to v1.3, incl. ISCCP simulator, 
 !  time control, openMP compatible, some bug fixes, de-implemented radar,
 !  threw away cosp types & convective cloud 
 ! 
 ! todo: further evaluation, sensitivity microphysics assumptions 
 !  
 !Uses ECHAM modules
  USE mo_kind,          ONLY: wp
  USE mo_linked_list,   ONLY: t_stream,  HYBRID, SURFACE, BELOWSUR, &
                              ABOVESUR10, COSP_LEVELS, NETCDF, GRIB
  USE mo_memory_base,   ONLY: new_stream, add_stream_element, &
                              default_stream_setting, add_stream_reference, AUTO
  USE mo_cosp_metoff_cosp

  USE mo_tr_omp_decomposition, ONLY: omp_get_thread_num

!  USE mo_radiation,  ONLY: emsfc_lw => cemiss  

  IMPLICIT NONE

  INTEGER, SAVE      :: cospo_index ! event index of the COSP stream

  TYPE, PUBLIC :: vmem3d
     REAL(wp), POINTER :: ptr(:,:,:)
  END TYPE vmem3d
  
  TYPE, PUBLIC :: vmem2d
     REAL(wp), POINTER :: ptr2(:,:)
  END TYPE vmem2d

  PUBLIC:: cosp_initialize
  PUBLIC:: construct_stream_cosp
  PUBLIC:: call_cospsimulator

  LOGICAL :: locosp = .FALSE.     ! Doing satellite output?
  LOGICAL, SAVE :: Llidar_sim,Lisccp_sim,Lstats,Llidar_cfad
  LOGICAL :: l_fixed_reff ! use fixed effective radii for testing only
  LOGICAL :: extra_output
  LOGICAL :: use_netcdf
  INTEGER :: offl2dout

  INTEGER :: Ncolumns          ! number of sub-columns used for each profile (namelist input)
  INTEGER, PARAMETER :: isccp_overlap = 3       ! Overlap Type: 1=max, 2=rand, 3=max/rand
  LOGICAL, PARAMETER :: use_reff = .TRUE. ! .true.: use effective radii from model  
                                          ! .false.: defaults

  INTEGER, PARAMETER :: lidar_ice_type = 0 ! (0=ice-spheres ; 1=ice-non-spherical)

  INTEGER, PARAMETER :: isccp_top_height = 1 ! 1 = adjust cloud top pressur using IR brightness 
                                             ! temperature and visible optical depth
  INTEGER, PARAMETER :: isccp_top_height_direction = 2 ! 2= find the lowest pressure level with interpolated  
                                             ! temperature equal to the radiance determined cloud-top temperature
  REAL(wp), PARAMETER :: isccp_emsfc_lw = 0.996_wp  ! Surface emissivity at 10.5 micron (fraction)
                                           ! ... pretty close, check sensitivity ...!
    
  TYPE (vmem3d), PRIVATE, POINTER :: cosp_cfadsr(:)
                                             ! COSP CFADSR for the 15 srbins of lidar scattering ratio
  TYPE (vmem2d), PRIVATE, POINTER :: isccp_cldtypes(:)   
                                             ! ISCCP cloud fraction for the 49 ISCCP cloud types
  TYPE (vmem2d), PRIVATE, POINTER :: cosp_parasol_refl(:)   
                                             ! PARASOL reflectance for 5 bins of solar zenith angle

  REAL(wp), PRIVATE, POINTER :: cosp_isccp_totalcloud(:,:)    ! ISCCP total cloud cover 
  REAL(wp), PRIVATE, POINTER :: cosp_isccp_meanptop(:,:)      ! mean ISCCP cloud top pressure
  REAL(wp), PRIVATE, POINTER :: cosp_isccp_meanalbedocld(:,:) ! mean ISCCP cloud albedo
  REAL(wp), PRIVATE, POINTER :: cosp_isccp_meantaucld(:,:)    ! mean ISCCP cloud optical thickness

  REAL(wp), PRIVATE, POINTER :: cosp_calipso_cloud_fraction(:,:,:) ! COSP Cloud frequency of occurance as seen by CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_lowcloud(:,:)           ! COSP Low-level cloud cover from CALIPSO  
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_midcloud(:,:)           ! COSP Mid-level cloud cover from CALIPSO  
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_highcloud(:,:)          ! COSP High-level cloud cover from CALIPSO
  REAL(wp), PRIVATE, POINTER :: cosp_lidar_totalcloud(:,:)         ! COSP Total cloud cover from CALIPSO

  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cosp_reffl     ! Liquid water droplet effective radius [um]
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: cosp_reffi     ! Ice crystal effective radius [um]

  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: xi_cosp     
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: xl_cosp     
  REAL(wp), PUBLIC, POINTER, DIMENSION(:,:,:) :: tm1_cosp     


  REAL(wp), PUBLIC,  POINTER, DIMENSION(:,:) :: cosp_sunlit    ! sunlit grid points
  REAL(wp), PUBLIC,  POINTER, DIMENSION(:,:) :: cosp_sunlit_av ! average sunlit fraction (needed for correctly 
                                                               !  averaging output)

  REAL(wp), PUBLIC,  POINTER, DIMENSION(:,:,:) :: N_miss_lidar   ! lidar simulator: number of missing points
  REAL(wp), PUBLIC,  POINTER, DIMENSION(:,:) :: N_miss_lidar_mid ! lidar simulator: number of missing points 
                                                                 !  mid level clouds
  REAL(wp), PUBLIC,  POINTER, DIMENSION(:,:) :: N_miss_lidar_low ! lidar simulator: number of missing points 
                                                                 !  low level clouds

  REAL(wp), PUBLIC,  POINTER, DIMENSION(:,:,:) :: cisccp_cldtau3d  ! Cloud optical thickness
  REAL(wp), PUBLIC,  POINTER, DIMENSION(:,:,:) :: cisccp_cldemi3d  ! Cloud emissivity @ 10.5 µm

  REAL(wp), PUBLIC,  POINTER, DIMENSION(:,:,:) :: cosp_f3d ! 3-d cloud fraction as seen in radiation
  
  REAL(wp), SAVE :: od_cosp ! fraction of day calls to satellite simulator 
                            ! (currently equal to radiation time step). 

  TYPE (t_stream), PUBLIC, POINTER :: scosp


CONTAINS
  
!--------------------------------------------------------------------------------------------------------------------

  SUBROUTINE cosp_initialize

   USE mo_mpi,          ONLY: p_parallel, p_parallel_io, p_bcast, p_io
   USE mo_namelist,     ONLY: open_nml, position_nml, POSITIONED
   USE mo_time_control, ONLY: trigrad
   USE mo_exception,    ONLY: finish

   IMPLICIT NONE 

   INTEGER :: ierr, inml, iunit

   NAMELIST /cospctl/ locosp,  Llidar_sim, Llidar_cfad, Lisccp_sim, &
                       l_fixed_reff, Ncolumns,extra_output, use_netcdf, &
                       offl2dout

    Llidar_sim = .true.
    Lisccp_sim = .true.
    Llidar_cfad = .false.
    l_fixed_reff = .false.
    extra_output = .false.
    use_netcdf  = .true.
    Ncolumns = 12
    offl2dout=-1
    
      !--- Local variables:
     IF (p_parallel_io) THEN
       inml = open_nml ('namelist.echam')
       iunit = position_nml ('COSPCTL', inml, status=ierr)
       SELECT CASE (ierr)
       CASE (POSITIONED)
         READ (iunit, cospctl)          
       END SELECT
     END IF

    IF (offl2dout .GT. 0) THEN
      extra_output = .true.
    END IF

      !--- 2) Broadcast over processors:
     IF (p_parallel) THEN
         CALL p_bcast (locosp, p_io)
         CALL p_bcast (Llidar_sim, p_io)
         CALL p_bcast (Llidar_cfad, p_io)
         CALL p_bcast (Lisccp_sim, p_io)
         CALL p_bcast (l_fixed_reff, p_io)
         CALL p_bcast (Ncolumns, p_io)
         CALL p_bcast (extra_output, p_io)
         CALL p_bcast (offl2dout, p_io)
      END IF ! p_parallel

    Lstats = .FALSE.
    IF ((Llidar_sim).OR.(Lisccp_sim)) Lstats = .TRUE.
   
   IF (.NOT. locosp) RETURN

    ! timestep for satellite simulator call
   IF ( trigrad%unit .EQ. 'hours' ) THEN
     od_cosp =  3600._wp * REAL(trigrad%counter, wp) / 86400._wp
   ELSE
     CALL finish('mo_cosp_simulator ','expecting different trigrad unit')
   END IF
 
END SUBROUTINE cosp_initialize

!-----------------------------------------------------------------------------------------------------------------

  SUBROUTINE construct_stream_cosp

    USE mo_time_event,  ONLY: io_time_event
    USE mo_exception,   ONLY: finish

    IMPLICIT NONE

    INTEGER :: i
    CHARACTER*40 :: name
    !
    ! Local
    ! 

!--- 1) Construct the cosp stream: ------------------------------------------------------------------------------

! do not change output frequency
    IF (use_netcdf) THEN
     IF (offl2dout .LT. 0 ) THEN !default
      CALL new_stream (scosp,'cosp',filetype=NETCDF, interval=io_time_event(1,'days','last',0) )
     ELSE
      CALL new_stream (scosp,'cosp',filetype=NETCDF, interval=io_time_event(offl2dout,'minutes','last',0) )
     END IF 
    ELSE
      CALL finish('mo_cosp_simulator ','GRIB OUTPUT NOT YET IMPLEMENTED')
     CALL new_stream (scosp,'cosp',filetype=GRIB, interval=io_time_event(1,'days','last',0) )
    END IF
    cospo_index  = scosp%post_idx

    !--- Add standard fields for post-processing:

    CALL add_stream_reference (scosp, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (scosp, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (scosp, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (scosp, 'gboxarea','geoloc',lpost=.TRUE.)

    IF (offl2dout .GT. 0 ) THEN
        CALL add_stream_reference (scosp, 'u10'   ,'g3b'   ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'v10'   ,'g3b'   ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'slm'   ,'g3b'   ,lpost=.TRUE.)

        CALL add_stream_reference (scosp, 'aclc'   ,'g3b'   ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'ao3'    ,'g3b'   ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'xl'     ,'gl'    ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'xi'     ,'gl'    ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'q'      ,'gl'    ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'tm1'    ,'g1a'   ,lpost=.TRUE.)
        CALL add_stream_reference (scosp, 'relhum' ,'g3b'   ,lpost=.TRUE.)

    CALL add_stream_element (scosp, 'xi_cosp', xi_cosp,   laccu     = .FALSE.,&
         leveltype = HYBRID, lpost=.TRUE., &
         longname='xi_cosp', units='')

    CALL add_stream_element (scosp, 'xl_cosp', xl_cosp,   laccu     = .FALSE.,&
         leveltype = HYBRID, lpost=.TRUE., &
         longname='xl_cosp', units='')

    CALL add_stream_element (scosp, 'tm1_cosp', tm1_cosp,   laccu     = .FALSE.,&
         leveltype = HYBRID, lpost=.TRUE., &
         longname='tm1_cosp', units='', lrerun=.TRUE. )

    END IF

!--- 2) Add stream elements: ------------------------------------------------------------------------------------

    CALL default_stream_setting (scosp,  units     = '',           &
                                         lrerun    = .FALSE. ,     &
                                         lpost     = .TRUE. ,     &
                                         laccu     = .FALSE. ,    &
                                         reset     = 1.e-50_wp,   &
                                         leveltype = SURFACE ,    &
                                         contnorest = .TRUE. )

!! CNAM: Check the units below 
!! lpost set to TRUE for debugging 
    CALL add_stream_element (scosp, 'reffl', cosp_reffl,   laccu     = .FALSE.,&
         leveltype = HYBRID, lpost=extra_output, &
         longname='Liquid water droplet effective radius', units='um')
    CALL add_stream_element (scosp, 'reffi', cosp_reffi,   laccu     = .FALSE.,&
         leveltype = HYBRID, lpost=extra_output, &
         longname='Ice crystal effective radius', units='um')
    CALL add_stream_element (scosp, 'cosp_sunlit', cosp_sunlit,  &
                              longname='COSP sunlit fraction')
    CALL add_stream_element (scosp, 'cosp_sunlit_av', cosp_sunlit_av,  &
                              longname='COSP average sunlit fraction')
    CALL add_stream_element (scosp, 'cosp_f3d', cosp_f3d, &
                              leveltype = HYBRID, &
                              longname='COSP cloud fraction',lpost=extra_output)

   IF (Llidar_sim) THEN
      ! 2d mean quantities
      !! CNAM: Don't know how to add 1D stream
      !!    CALL add_stream_element (scosp, 'srbval', cosp_srbval,  &
      !!                            longname='Values of lowest edges of cfad_sr', units='1')
      CALL add_stream_element (scosp, 'calipso_cloud_fraction', cosp_calipso_cloud_fraction,  &
                           leveltype = COSP_LEVELS, &
                           longname='CALIPSO Lidar Cloud Fraction (532nm)', units='1')
      CALL add_stream_element (scosp, 'lidar_lowcloud', cosp_lidar_lowcloud,  &
                           longname='CALIPSO Low-level cloud fraction', units='1')
      CALL add_stream_element (scosp, 'lidar_midcloud', cosp_lidar_midcloud,  &
                           longname='CALIPSO Mid-level cloud fraction', units='1')
      CALL add_stream_element (scosp, 'lidar_highcloud', cosp_lidar_highcloud,  &
                           longname='CALIPSO High-level cloud fraction', units='1')
      CALL add_stream_element (scosp, 'lidar_totalcloud', cosp_lidar_totalcloud,  &
                           longname='CALIPSO Total cloud fraction', units='1')


    IF (Llidar_cfad) THEN
      ALLOCATE(cosp_cfadsr(15))
          DO i=1,15
                WRITE(name, '(a7,i2.2)') 'cfad_sr',i
             CALL add_stream_element (scosp, name, cosp_cfadsr(i)%ptr,  &
                  leveltype = COSP_LEVELS, longname='COSP CFAD Lidar Scattering Ratio 532nm')
          ENDDO
    END IF 

      CALL add_stream_element (scosp, 'N_miss_lidar_mid', N_miss_lidar_mid,  &
                            longname='COSP lidar number of missing values mid level cloud')
      CALL add_stream_element (scosp, 'N_miss_lidar_low', N_miss_lidar_low,  &
                            longname='COSP lidar number of missing values low level cloud')
      CALL add_stream_element (scosp, 'N_miss_lidar', N_miss_lidar,  &
                            leveltype = COSP_LEVELS, &
                            longname='COSP lidar number of missing values')

      ALLOCATE(cosp_parasol_refl(PARASOL_NREFL))

      DO i=1, PARASOL_NREFL
         WRITE(name, '(a11,i1)') 'parasolRefl',i
          CALL add_stream_element (scosp, name, cosp_parasol_refl(i)%ptr2, &
                               longname='Parasol like mono-directional reflectance', units='1')
      ENDDO

   ENDIF! Llidar_sim
 
   IF (Lisccp_sim) THEN

      CALL add_stream_element (scosp, 'cltisccp', cosp_isccp_totalcloud,  &
                           longname='ISCCP total cloud area fraction', units='1')
      CALL add_stream_element (scosp, 'pctisccp', cosp_isccp_meanptop,  &
                           longname='ISCCP cloud top pressure', units='Pa')
      CALL add_stream_element (scosp, 'albisccp', cosp_isccp_meanalbedocld,  &
                           longname='ISCCP mean cloud albedo', units='1')
      CALL add_stream_element (scosp, 'tauisccp', cosp_isccp_meantaucld,  &
                           longname='ISCCP mean cloud tau', units='1')
      CALL add_stream_element (scosp, 'cisccp_tau3d', cisccp_cldtau3d, &
                           leveltype = HYBRID, lpost=extra_output,  &
                           longname='COSP-ISCCP cloud optical thickness')
      CALL add_stream_element (scosp, 'cisccp_emi3d', cisccp_cldemi3d, &
                           leveltype = HYBRID, lpost=extra_output, &
                           longname='COSP-ISCCP cloud optical thickness')

      ALLOCATE(isccp_cldtypes(49))
      DO i=1,49
            WRITE(name, '(a7,i2.2)') 'cldtype',i
         CALL add_stream_element (scosp, name, isccp_cldtypes(i)%ptr2,  &
                longname='Fraction covered by ISCCP cloud type')
      ENDDO

   END IF!Lisccp_sim

  END SUBROUTINE construct_stream_cosp

!--------------------------------------------------------------------------------------------------------------------

  SUBROUTINE call_cospsimulator(             &
               kproma,          klev,                    jrow,                &
               p,              ph,                      pgeo,                &
               pgeospm,        pslm,                    ptm1,                &
               qm,                                      xlm1,                &
               xim1,           zfrl,                    zfrw,                &
               zfri,           tslm1,                   tsw,                 &
               tsi                                                           )


    USE mo_physical_constants,ONLY: grav        ! Grativation acceleration
    USE mo_cosp_constants,    ONLY: R_UNDEF
    USE mo_time_control,      ONLY: l_trigrad,  l_putdata
    USE mo_exception,         ONLY: message, message_text
    USE mo_random_numbers,    ONLY: set_seed_random

    IMPLICIT NONE

    !
    ! Input to COSP-Simulator set here
    !
    INTEGER                             :: kproma, klev, jrow ! Dimensions
    REAL(wp), DIMENSION(kproma, klev)    :: height       ! Height of model levels [m]
    REAL(wp), DIMENSION(kproma, klev+1)  :: height_half  ! Height @ layer interfaces  [m]

    REAL(wp), DIMENSION(kproma, klev)    :: p            ! Pressure @ full levels [Pa]
    REAL(wp), DIMENSION(kproma, klev+1)  :: ph           ! Pressure @ layer interfaces [Pa]
    REAL(wp), DIMENSION(kproma, klev)    :: pgeo         ! geopotential height above surface
    REAL(wp), DIMENSION(kproma)  :: pgeospm ! surface geopotential height
    REAL(wp), DIMENSION(kproma)  :: pgeospm_l ! surface geopotential height, 0 over water


    REAL(wp), DIMENSION(kproma)  :: pslm

    REAL(wp), DIMENSION(kproma)  :: zfrl, zfrw, zfri,  tslm1, tsw, tsi

 
    REAL(wp), DIMENSION(kproma, klev)    :: ptm1         ! Temperature at model level [K]
    REAL(wp), DIMENSION(kproma, klev)    :: qm           ! Water vapour mixing ratio [kg kg-1]

    !   REAL(wp), DIMENSION(kproma, 1)   :: psfc         ! Surface pressure [Pa]

    REAL(wp), DIMENSION(kproma, klev)    :: xlm1         ! Mixing ratio large scale cloud liquid [kg kg-1]
    REAL(wp), DIMENSION(kproma, klev)    :: xim1         ! Mixing ratio large scale cloud ice [kg kg-1]

    REAL(wp), DIMENSION(:,:,:), POINTER :: cs
    REAL(wp), DIMENSION(:,:), POINTER  :: ct

! gbx gridbox

    REAL(wp) :: p_l(kproma,klev) 
    REAL(wp) :: ph_l(kproma,klev)
    REAL(wp) :: T(kproma,klev) 
    REAL(wp) :: sh(kproma,klev)
    REAL(wp) :: dtau_s(kproma,klev)
    REAL(wp) :: dem_s(kproma,klev)
    REAL(wp) :: tca(kproma,klev)
    REAL(wp) :: sunlit(kproma)
    REAL(wp) :: skt(kproma)
    REAL(wp) :: Reff(kproma,klev,N_hydro)
    REAL(wp) :: zlev_half(kproma,klev)
    REAL(wp) :: zlev(kproma,klev)
    REAL(wp) :: psfc(kproma)
    REAL(wp) :: land(kproma)
    REAL(wp) :: mr_hydro(kproma,klev,N_hydro)

! sglidar output from lidar
    REAL(wp) :: beta_mol(kproma,klev)
    REAL(wp) :: beta_tot(kproma,Ncolumns,klev)
    REAL(wp) :: tau_tot(kproma,Ncolumns,klev)
    REAL(wp) :: refl(kproma,Ncolumns,PARASOL_NREFL)

! stlidar lidar stats
    REAL(wp) :: srbval(SR_BINS)
    REAL(wp) :: cfad_sr(kproma,SR_BINS,Nlr) 
    REAL(wp) :: lidarcld(kproma,Nlr)
    REAL(wp) :: cldlayer(kproma,LIDAR_NCAT)
    REAL(wp) :: parasolrefl(kproma,PARASOL_NREFL)

! isccp
    REAL(wp) :: fq_isccp(kproma,7,7), totalcldarea(kproma)
    REAL(wp) :: meanptop(kproma), meantaucld(kproma)
    REAL(wp) :: meantb(kproma), meantbclr(kproma)
    REAL(wp) :: boxtau(kproma,Ncolumns), boxptop(kproma,Ncolumns)
    REAL(wp) :: meanalbedocld(kproma)

    INTEGER:: k,j,i, itau, ipres, jk, jkb, jl

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    IF (l_putdata(cospo_index) .AND. jrow .eq.1 ) THEN
       WRITE(message_text,*) 'output step', jrow
       CALL message('cosp ',message_text) 
    END IF

    !only call at radiation step
    IF ( .NOT.  l_trigrad ) THEN
       RETURN
    END IF



    IF (jrow .eq.1)THEN
       WRITE(message_text,*) 'called cosp',  od_cosp
       CALL message('cosp: ',message_text) 
    END IF

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!gbx

    zlev      = 0.0_wp
    zlev_half = 0.0_wp
    p_l       = 0.0_wp
    ph_l      = 0.0_wp
    T         = 0.0_wp
    sh        = 0.0_wp
    dtau_s    = 0.0_wp
    dem_s     = 0.0_wp
    tca       = 0.0_wp
    Reff      = 0.0_wp
    psfc      = 0.0_wp
    land      = 0.0_wp
    sunlit    = 0.0_wp
    skt       = 0.0_wp
    mr_hydro  = 0.0_wp

!sglidar
    srbval    = 0.0_wp
    cfad_sr   = 0.0_wp
    lidarcld  = 0.0_wp
    cldlayer  = 0.0_wp
    parasolrefl  = 0.0_wp
    beta_mol   = 0.0_wp
    beta_tot   = 0.0_wp
    tau_tot    = 0.0_wp
    refl       = 0.0_wp ! parasol

!stlidar
!    srbval    = 0.0_wp
!    cfad_sr   = 0.0_wp
!    lidarcld  = 0.0_wp
!    cldlayer  = 0.0_wp
!    parasolrefl  = 0.0_wp

!isccp
    fq_isccp     = 0.0_wp
    totalcldarea = 0.0_wp
    meanptop     = 0.0_wp
    meantaucld   = 0.0_wp
    meantb       = 0.0_wp
    meantbclr    = 0.0_wp
    boxtau       = 0.0_wp
    boxptop      = 0.0_wp
    meanalbedocld= 0.0_wp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Height of model levels [m]: ECHAM5 Top of atmosphere is 1, surface is klev (i.e 19)
   DO jl = 1,kproma
     IF (pslm(jl) .GT. 0._wp) THEN
        pgeospm_l(jl) = pgeospm( jl)
     ELSE
        pgeospm_l(jl) =  0._wp
     END IF
   END DO

   DO jk = 1,klev
      DO jl = 1,kproma
          height(jl,jk) = (pgeo(jl,jk)  +  pgeospm_l( jl) )/grav
      END DO
  END DO
    !cms++

    DO jk = 2,klev
       DO jl = 1,kproma
          height_half(jl,jk) = 0.5_wp*(height(jl,jk)+height(jl,jk-1))
       END DO
    END DO
    DO jl = 1,kproma
       height_half(jl,1) =  height(jl,1) + ( height(jl,1) - height_half(jl,2) ) 
       height_half(jl,klev+1) =   pgeospm_l( jl ) /grav 
    END DO

    !cms--

    ! Relative Humidity of model levels [%]
    ! [in echam src directory, physc.f90 has variable relhum]

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Create gbx structure
    ! To see what the variables mean:
    ! CONSTRUCT_COSP_GRIDBOX(radar_freq,surface_radar, use_mie_tables, use_gas_abs, do_ray,melt_lay, k2, &  (no)
    !               Npoints,Nlevels,Ncolumns,Nhydro,
    !  Nprmts_max_hydro, (no)
    !  Naero,
    !  Nprmts_max_aero,  (no)
    !  lidar_ice_type, &
    !   isccp_top_height,isccp_top_height_direction,isccp_overlap,isccp_emsfc_lw, &
    !               use_reff,y)
    ! *** Following values from cosp_input_nl.txt 
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
    ! npoints: kproma
    ! nlevels: klev

 !   CALL construct_cosp_gridbox(kproma,klev,Ncolumns,N_HYDRO,1,0,    &
 !        isccp_top_height,isccp_top_height_direction,overlap,emsfc_lw,                             &
 !        .true.,gbx)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Code to populate input structure
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    zlev(1:kproma,1:klev) = height(1:kproma,klev:1:-1)
    zlev_half(1:kproma,1:klev) = height_half(1:kproma,klev+1:2:-1)        
    p_l(1:kproma,1:klev) = p(1:kproma,klev:1:-1)
    ph_l(1:kproma,1:klev) = ph(1:kproma,klev+1:2:-1)
    T(1:kproma,1:klev) = ptm1(1:kproma,klev:1:-1)
    tca(1:kproma,1:klev) = cosp_f3d(1:kproma,klev:1:-1,jrow) !tca(1:kproma,klev:1:-1)
    psfc(1:kproma) = ph(1:kproma,klev+1) 
    land(1:kproma) = pslm(1:kproma)        

    IF ( l_fixed_reff ) THEN
       ! Constant values for cloud droplet effective radius of liquid and ice in the test
       Reff(1:kproma,1:klev,I_LSCLIQ)=10.0e-6_wp
       Reff(1:kproma,1:klev,I_LSCICE)= 40.0e-6_wp
    ELSE
       Reff(1:kproma,1:klev,I_LSCLIQ)=cosp_reffl(1:kproma,klev:1:-1,jrow)*1.e-6_wp
       Reff(1:kproma,1:klev,I_LSCICE)=cosp_reffi(1:kproma,klev:1:-1,jrow)*1.e-6_wp 
    END IF
    
    sunlit(1:kproma)   = cosp_sunlit(1:kproma,jrow)  

    ! ISCCP simulator
    IF (Lisccp_sim) THEN
       skt(1:kproma) = zfrl(1:kproma)*tslm1(1:kproma)              &
             +zfri(1:kproma)*tsi(1:kproma)                &
             +zfrw(1:kproma)*tsw(1:kproma)

           ! old ptm1(1:kproma,klev)
    
       sh(1:kproma,1:klev) = qm(1:kproma,klev:1:-1)   

       DO jk=1,klev
          jkb = klev+1-jk
          DO jl=1,kproma
             IF ( cosp_f3d(jl,jkb,jrow) .GT. 0._wp ) THEN
                !cms: cisccp_cldtau3d is already in-cloud, (and the all-sky comment
                ! in mo_srtm.f90 close to line 1111 is not correct)
                dtau_s(jl,jk) = cisccp_cldtau3d(jl,jkb,jrow)
                dem_s(jl,jk) = cisccp_cldemi3d(jl,jkb,jrow)
             ELSE
                dtau_s(jl,jk) = 0._wp
                dem_s(jl,jk) = 0._wp
             END IF
          END DO
       END DO

    END IF
  
! hydrometeors
    DO jk=1,klev
       jkb = klev+1-jk
       DO jl=1,kproma
!!          IF ( cosp_f3d(jl,jkb,jrow).GT.0._wp ) THEN
             mr_hydro(jl,jk,I_LSCLIQ) = xlm1(jl,jkb) 
!!! / cosp_f3d(jl,jkb,jrow)
             mr_hydro(jl,jk,I_LSCICE) = xim1(jl,jkb) 
!!/ & 
!!                  cosp_f3d(jl,jkb,jrow)
!!          ELSE
!!             mr_hydro(jl,jk,I_LSCLIQ) = 0._wp
!!             mr_hydro(jl,jk,I_LSCICE) = 0._wp
!!          END IF
       END DO
    END DO


! extra diagnostic for offl test only
   IF  ( offl2dout .GT. 0) THEN
     DO jk=1,klev
       DO jl=1,kproma
          tm1_cosp(jl,jk,jrow) =ptm1 (jl,jk)
          xi_cosp(jl,jk,jrow) = xim1 (jl,jk)
          xl_cosp(jl,jk,jrow) = xlm1 (jl,jk)
       END DO
      END DO
   END IF
    
   !
   ! Initialize global random number generator used inside COSP
   !   This is not thread or parallel safe but affects only COSP diagnostics
   !
   CALL set_seed_random(INT((T(1,1:klev) - INT(T(1,1:klev)))*10000000))
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Call simulator
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL cosp(kproma,klev, Ncolumns,                        &
              Lisccp_sim,  Llidar_sim, Llidar_cfad, Lstats,        &
!gbx
              use_reff, lidar_ice_type,                    &
              isccp_top_height, isccp_top_height_direction, &
              isccp_overlap, isccp_emsfc_lw,               &
              ph_l, p_l, T, Reff, tca,                     &             
              psfc, mr_hydro, sh, dem_s, dtau_s,           &
              sunlit, skt, zlev, zlev_half, land,          &
!sglidar
              beta_mol, beta_tot, refl, tau_tot,           &
!stlidar
              srbval, cfad_sr, lidarcld,                   &
              cldlayer, parasolrefl,                       &
!isccp
              fq_isccp,  totalcldarea,                     &
              meanptop,  meantaucld,                       &
              meantb, meantbclr, boxtau,  boxptop,         &
              meanalbedocld                                )
   

    DO j=1,kproma
       cosp_sunlit_av(j,jrow) = cosp_sunlit_av(j,jrow)+   cosp_sunlit(j,jrow) *  od_cosp
    END DO

    IF (Llidar_sim) THEN
       ! For 'cosp_calipso_cloud_fraction'
       DO k=1,Nlr
          DO j=1,kproma
             IF (lidarcld(j,k) .NE.  R_UNDEF  ) THEN
                cosp_calipso_cloud_fraction(j,k,jrow)=cosp_calipso_cloud_fraction(j,k,jrow) + lidarcld(j,k) *  od_cosp
             ELSE
                N_miss_lidar(j,k,jrow) =  N_miss_lidar(j,k,jrow) +  1._wp *  od_cosp
             END IF
          END DO
       ENDDO

       ! For 'cosp_lidar_parasol_refl'
       DO k=1,PARASOL_NREFL
           ct => cosp_parasol_refl(k)%ptr2
          DO j=1,kproma
             ct(j,jrow) = ct(j,jrow) + parasolrefl(j,k) *  od_cosp
          END DO
       END DO

       DO j=1,kproma
          IF ( cldlayer(j,1) .NE. R_UNDEF  ) THEN
             cosp_lidar_lowcloud(j,jrow)= cosp_lidar_lowcloud(j,jrow) &
                  + cldlayer(j,1) *  od_cosp
          ELSE
             N_miss_lidar_low(j,jrow)= N_miss_lidar_low(j,jrow) + 1._wp *  od_cosp
          END IF
       END DO

       !For 'cosp_lidar_midcloud'
       DO j=1,kproma
          IF ( cldlayer(j,2) .NE. R_UNDEF  ) THEN
             cosp_lidar_midcloud(j,jrow)= cosp_lidar_midcloud(j,jrow) + &
                  cldlayer(j,2) *  od_cosp 
          ELSE
             N_miss_lidar_mid(j,jrow)= N_miss_lidar_mid(j,jrow) + 1._wp *  od_cosp
          END IF
       END DO

       !For 'cosp_lidar_highcloud'
       DO j=1,kproma
          cosp_lidar_highcloud(j,jrow)= cosp_lidar_highcloud(j,jrow) + cldlayer(j,3) *  od_cosp
       ENDDO

       !For 'cosp_lidar_totalcloud'
       DO j=1,kproma
          cosp_lidar_totalcloud(j,jrow)= cosp_lidar_totalcloud(j,jrow) + cldlayer(j,4) *  od_cosp
       ENDDO

       !For CFAD of Lidar scattering ratio
     IF (Llidar_cfad) THEN
       DO i=1,15
          DO k=1,Nlr
             DO j=1,kproma
                cs => cosp_cfadsr(i)%ptr
                cs(j,k,jrow)=cs(j,k,jrow) + cfad_sr(j,i,k) *  od_cosp
             END DO
          ENDDO
       ENDDO
     END IF
    END IF !Llidar_sim
   

    IF (Lisccp_sim) THEN

       DO j=1,kproma
          cosp_isccp_totalcloud(j,jrow) = cosp_isccp_totalcloud(j,jrow)+ totalcldarea(j)  *  od_cosp
       ENDDO

       ! weight by totalcloud for time mean
       DO j=1,kproma
          cosp_isccp_meanptop(j,jrow) = cosp_isccp_meanptop(j,jrow)+  totalcldarea(j) * meanptop(j)  *  od_cosp
       ENDDO

       ! weight by totalcloud for time mean
       DO j=1,kproma
          cosp_isccp_meanalbedocld(j,jrow) = cosp_isccp_meanalbedocld(j,jrow) &
               + totalcldarea(j) * meanalbedocld(j)  *  od_cosp
       ENDDO

       ! weight by totalcloud for time mean
       DO j=1,kproma
          cosp_isccp_meantaucld(j,jrow) = cosp_isccp_meantaucld(j,jrow)   &
               +   totalcldarea(j) * meantaucld(j)  *  od_cosp
       ENDDO

       i=0
       DO itau=1,7
          DO ipres=1,7
             i=i+1
             ct => isccp_cldtypes(i)%ptr2
             ct(1:kproma,jrow) =  ct(1:kproma,jrow) + fq_isccp(1:kproma,itau,ipres)  *  od_cosp
          ENDDO
       ENDDO

    END IF

  END SUBROUTINE call_cospsimulator
!--------------------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------------------------

END MODULE mo_cosp_simulator
