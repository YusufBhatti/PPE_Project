!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham_dust.f90
!!
!! \brief
!! mo_ham_dust contains the MPI BGC dust emission scheme.
!!
!! \author Martin Werner, MPI Biogeochemistry
!!
!! \responsible_coder
!! Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Werner (MPI BGC) - original code (2002-10)
!!   -# P. Stier (MPI-Met) - fitted parameters for emissions size distribution
!!                           introduced scale factor for wind stress threshold (2003)
!!   -# M. Werner (MPI BGC) - changed scheme from 0.5x0.5 degree grid to default
!!                            Echam5-HAM grid (to save computational costs) (2004-10)
!!   -# M. Schultz (FZ Juelich) - merged all bgc_dust modules into one,
!!                                also added setdust routine and ham_dustctl namelist (2009-10-01)
!!   -# C. Siegenthaler (C2SM-ETHZ) - implementation of a regional tuning factor nduscale 
!!                                    Issue 433 (06-2015)
!!   -# S. Ferrachat (ETH Zurich) - Improved support for namelist settings (See #479)
!!                                - Soil type input file merging (#253)
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!!   - Marticorena, B., and G. Bergametti (1995), Modeling the atmospheric dust cycle: 1. 
!!     Design of a soil-derived dust emission scheme, 
!!     J. Geophys. Res., 100(D8), 16,415?C16,430.
!!   - Marticorena, B., G. Bergametti, B. Aumont, Y. Callot, C. N'Doume, and M. Legrand (1997), 
!!     Modeling the atmospheric dust cycle 2. Simulation of Saharan dust sources, 
!!     J. Geophys. Res., 102(D4), 4387?C4404.
!!   - Tegen, I., S. P. Harrison, K. Kohfeld, I. C. Prentice, M. Coe, and M. Heimann (2002), 
!!     Impact of vegetation and preferential source areas on global dust aerosol: Results 
!!     from a model study, 
!!     J. Geophys. Res., 107(D21), 4576, doi:10.1029/2001JD000963. 
!!   - Cheng, T., Peng, Y., Feichter, J., and Tegen, I.: 
!!     An improvement on the dust emission scheme in the global aerosol-climate model ECHAM5-HAM, 
!!     Atmos. Chem. Phys., 8, 1105-1117, 2008. 
!!     http://www.atmos-chem-phys.net/8/1105/2008/acp-8-1105-2008.html
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

MODULE mo_ham_dust

  
  USE mo_kind, ONLY: dp
  USE mo_mpi,  ONLY: p_parallel, p_parallel_io, p_bcast, p_io, p_pe
  USE mo_ham,  ONLY: ndust !SF #479

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setdust,                        &
            bgc_dust_initialize,            &
            bgc_dust_init_diag,             &
            bgc_dust_read_monthly,          &
            bgc_dust_calc_emis,             &
            bgc_dust_trastat,               &
            bgc_dust_cleanup

  PUBLIC :: ndurough

  PUBLIC :: flux_6h, ntrace         ! emission mass flux 

  PUBLIC :: Dmin, Dmax, nbinit, Dstep

  INTEGER :: ibc_regint ! boundary condition indices for regional mask file (!csld #433) 

! *** dust parameter
  REAL(dp), PARAMETER     :: a_rnolds=1331.647_dp          ! Reynolds constant
  REAL(dp), PARAMETER     :: aeff=0.35_dp                  ! efficient fraction

  REAL(dp), PARAMETER     :: b_rnolds=0.38194_dp           ! Reynolds constant
  
  REAL(dp), PARAMETER     :: d_thrsld=0.00000231_dp        ! thresold value
  REAL(dp), PARAMETER     :: Dmin=0.00002_dp               ! minimum particules diameter (cm)
  REAL(dp), PARAMETER     :: Dmax=0.130_dp                 ! maximum particules diameter (cm)
  REAL(dp), PARAMETER     :: Dstep=0.0460517018598807_dp   ! diameter increment (cm)
      
  INTEGER, PARAMETER      :: ntrace=8                      ! number of tracers
  INTEGER, PARAMETER      :: nbin=24                       ! number of bins per tracer
  INTEGER, PARAMETER      :: nbinit=24                     ! number of bins per tracer
                                                           !alaak: needed for salsa interface
  INTEGER, PARAMETER      :: nclass=ntrace*nbin            ! number of particle classes
  INTEGER, PARAMETER      :: nats =17                      ! number of soil types
  INTEGER, PARAMETER      :: nmode=4
  INTEGER, PARAMETER      :: nspe = nmode*3+2
    
  REAL(dp), PARAMETER     :: roa=0.001227_dp               ! air density (g/cm-3)
  REAL(dp), PARAMETER     :: rop=2.65_dp                   ! particle density (g/cm-3)
  
  REAL(dp), PARAMETER     :: umin=21._dp                   ! minimum threshold friction windspeed (cm/s)
 
  REAL(dp), PARAMETER     :: vk=0.4_dp                     ! Von Karman constant: 0.4 (0.35 <--> 0.42 see Stull) 
  
  REAL(dp), PARAMETER     :: w0=0.99_dp                    ! threshold of relative soil humidity 
  
  REAL(dp), PARAMETER     :: x_rnolds=1.561228_dp          ! Reynolds constant
  REAL(dp), PARAMETER     :: xeff=10._dp                   ! efficient fraction
  
  REAL(dp), PARAMETER     :: ZZ=1000._dp                   ! wind measurment height (cm)
   
  REAL(dp), PARAMETER     :: Z0s=0.001_dp                  ! roughness length of surface without obstacles(cm)
                                                            ! (see: thesis of B.Marticorena, p.85) 
  ! Attention: parameters below listed are only for sensitivity tests
  !            and should not be changed in the control run
                                       !
  REAL(dp):: ndurough     = 1.0E-03_dp ! Surface roughness length (cm)
                                       !
                                       !    ndurough = 0 A monthly mean satellite derived (Prigent et al.,
                                       !                 JGR 2005) surface roughness length map is used
                                       !             > 0 The globally constant surface roughness length
                                       !                 ndurough (cm) us used.
                                       ! default 0.001cm
                                       !
  REAL(dp):: nduscale_reg(8) = 8.6E-01_dp ! Regional scale factor for threshold wind friction velocity !csld #433
                                      !  
                                      ! The indices correspond to the following regions
                                      !     1       All the other locations than the following regions
                                      !     2       North america
                                      !     3       South America 
                                      !     4       North Africa 
                                      !     5       South Africa 
                                      !     6       Middle East 
                                      !     7       Asia
                                      !     8       Australia 
  REAL(dp):: r_dust_lai   = 1.0E-10_dp ! Parameter for the threshold lai value (unitless)
                                       ! default 1.E-10
                                       !
  REAL(dp):: r_dust_umin  = 2.1E+01_dp ! Minimum threshold friction windspeed (cm/s)
                                       ! default 21. cm/s
                                       !
  REAL(dp):: r_dust_z0s   = 1.0E-03_dp ! z0s (cm)
                                       ! default 0.001 cm
                                       !
  REAL(dp):: r_dust_scz0  = 1.0E+00_dp ! Scale factor of satellite z0 (unitless)
                                       ! default 1.
                                       !
  REAL(dp):: r_dust_z0min = 1.0E-05_dp ! Parameter for minimum of z0 (cm)
                                       ! default 1.E-05 cm
                                       !
  REAL(dp):: r_dust_sf13  = 1.0E+00_dp ! duscale over Takelimakan desert (unitless)
  REAL(dp):: r_dust_sf14  = 1.0E+00_dp ! duscale over Loess (unitless)
  REAL(dp):: r_dust_sf15  = 1.0E+00_dp ! duscale over Gobi desert (unitless)
  REAL(dp):: r_dust_sf16  = 1.0E+00_dp ! duscale over other mixture soils (unitless)
  REAL(dp):: r_dust_sf17  = 1.0E+00_dp ! duscale over desert and sand land (unitless)
                                       ! default 1.
                                       !
  REAL(dp):: r_dust_af13  = 1.9E-06_dp ! Parameter for the alfa value over Takelimakan desert (unitless)
  REAL(dp):: r_dust_af14  = 1.9E-04_dp ! Parameter for the alfa value over Loess (unitless)
  REAL(dp):: r_dust_af15  = 3.9E-05_dp ! Parameter for the alfa value over Gobi desert  (unitless)
  REAL(dp):: r_dust_af16  = 3.1E-05_dp ! Parameter for the alfa value over other mixture soils (unitless)
  REAL(dp):: r_dust_af17  = 2.8E-06_dp ! Parameter for the alfa value over desert and sand land  (unitless)
                                       ! default values taken from Cheng et al.(2008)
                                       !
  INTEGER :: k_dust_smst  = 1          ! Effect of soil moisture on threshold wind friction velocity
                                       !    0: on
                                       !    1: off
                                       ! default 1
                                       !
  INTEGER :: k_dust_easo  = 1          ! Including the new East-Asia soil type
                                       !    0: on, cheng's implementation, with a bug
                                       !    1: off
                                       !    2: on, bug-removed-version 0
                                       ! default 1
                                       !

     
! *** dust emssion scheme related variables
  REAL(dp), ALLOCATABLE    :: biome(:,:)            ! biome distribution 

  REAL(dp)                 :: dustsum(4)            ! total dust flux (for different paricle sizes)

  !>>dod omp bugfix
  REAL(dp), ALLOCATABLE    :: flux_6h(:,:,:)        ! 6h flux (for tm3 use)
  !$OMP THREADPRIVATE (flux_6h)
  !<<dod

  REAL(dp), ALLOCATABLE    :: flux_ann(:,:)         ! mean annual flux 
  REAL(dp), ALLOCATABLE    :: flux_a1(:,:)          ! for statistics of small particles
  REAL(dp), ALLOCATABLE    :: flux_a2(:,:)          ! for statistics of small particles
  REAL(dp), ALLOCATABLE    :: flux_a10(:,:)         ! for statistics of small particles

  INTEGER, ALLOCATABLE     :: idust(:,:)            ! dust potential source (1: source, 0: no source)
   
  REAL(dp), ALLOCATABLE    :: k_fpar_eff(:,:)       ! KERNEL effective FPAR

  REAL(dp), ALLOCATABLE         :: mat_s1(:,:)    ! soil type #1 area (in relative percentage)
  REAL(dp), ALLOCATABLE, TARGET :: mat_s2(:,:)    ! soil type #2 area (in relative percentage)
  REAL(dp), ALLOCATABLE, TARGET :: mat_s3(:,:)    ! soil type #3 area (in relative percentage)
  REAL(dp), ALLOCATABLE, TARGET :: mat_s4(:,:)    ! soil type #4 area (in relative percentage)
  REAL(dp), ALLOCATABLE, TARGET :: mat_s6(:,:)    ! soil type #6 area (in relative percentage)
  REAL(dp), ALLOCATABLE, TARGET :: mat_s13(:,:)   ! soil type #13 area (in relative percentage)
  REAL(dp), ALLOCATABLE, TARGET :: mat_s14(:,:)   ! soil type #14 area (in relative percentage)
  REAL(dp), ALLOCATABLE, TARGET :: mat_s15(:,:)   ! soil type #15 area (in relative percentage)
  REAL(dp), ALLOCATABLE, TARGET :: mat_s16(:,:)   ! soil type #16 area (in relative percentage)
  REAL(dp), ALLOCATABLE, TARGET :: mat_s17(:,:)   ! soil type #17 area (in relative percentage)
  REAL(dp), ALLOCATABLE         :: mat_psrc(:,:)  ! preferential source area (in relative percentage)
  REAL(dp), ALLOCATABLE         :: mat_msg(:,:)   ! MSG Saharan pref. source area (in relative percent.) !BH #382
 
  REAL(dp), ALLOCATABLE    :: Z01(:,:)              ! surface rough length (see: Thesis of B.Marticorena)
  REAL(dp), ALLOCATABLE    :: Z02(:,:)              ! surface rough length (see: Thesis of B.Marticorena)
  REAL(dp), ALLOCATABLE    :: nduscale_2d(:,:)      ! 2D nduscale !csld #433

  REAL(dp)                 :: srel(nats,nclass)     ! relative surface 
  REAL(dp)                 :: srelV(nats,nclass)
  REAL(dp)                 :: su_srelV(nats,nclass)
   
  REAL(dp)                 :: Uth(nclass)           ! threshold friction velocity

  !>>SF #253
  INTEGER, PARAMETER :: nsoil_types = 9

  TYPE mat_general
       REAL(dp), POINTER :: ptr(:,:)
  END TYPE mat_general

  TYPE(mat_general) :: mat_all(nsoil_types)

  CHARACTER(len=2),PARAMETER :: jtypes(nsoil_types) = (/ &
                                                ' 2',' 3',' 4',' 6','13','14','15','16','17' /) ! Soil types labels
  !<<SF #253

  ! -- diagnostic fields
  REAL(dp), POINTER        :: ustar_acrit(:,:)

  ! --- from former mo_bgc_dust_data:
  ! combimax    value   int.    number of different biome types
  INTEGER, PARAMETER :: combimax=29
  ! active vegetation types
  INTEGER :: active(combimax)
  ! soil specs
  REAL(dp) :: solspe(17,14)

  CONTAINS

  !--- set dust_data

  SUBROUTINE set_dust_data

  ! Description:
  ! ------------
  ! *mo_bgc_dust_data* holds some classification dependent biome and soil data 
  !
  ! Author:
  ! -------
  ! Martin Werner, MPI Biogeochemistry (MPI BGC), Jena, May 2003.
  ! Tiantao Cheng, MPI-M, Hamburg, 2006-2007.
  !      add East-Asia soil type in solspe (type 13-17)
  !
  !   The 28 different biome types output by BIOME4
  !-----------------------------------------------------------------
  !     1       Tropical evergreen broadleaf forest
  !     2       Tropical semi-evergreen broadleaf forest
  !     3       Tropical deciduous broadleaf forest and woodland
  !     4       Temperate deciduous broadleaf forest
  !     5       Temperate evergreen needleleaf forest
  !     6       Warm-temperate evergreen broadleaf and mixed forest
  !     7       Cool mixed forest
  !     8       Cool evergreen needleleaf forest
  !     9       Cool-temperate evergreen needleleaf and mixed forest
  !     10      Cold evergreen needleleaf forest
  !     11      Cold deciduous forest
  !     12      Tropical savanna
  !     13      Tropical xerophytic shrubland
  !     14      Temperate xerophytic shrubland
  !     15      Temperate sclerophyll woodland and shrubland
  !     16      Temperate deciduous broadleaf savanna
  !     17      Temperate evergreen needleleaf open woodland
  !     18      Cold parkland
  !     19      Tropical grassland
  !     20      Temperate grassland
  !     21      Desert
  !     22      Graminoid and forb tundra
  !     23      Low and high shrub tundra
  !     24      Erect dwarf-shrub tundra
  !     25      Prostrate dwarf-shrub tundra
  !     26      Cushion-forb tundra
  !     27      Barren
  !     28      Ice
  !    (29)     Water (this is implied)

  !
  ! active vegetation types:
  ! ------------------------
  ! 
  active = (/ &
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
       1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, & 
       1, 1, 1, 0, 0 /)

!!mgs!!   commented out because shrubs is not used anywhere
!!mgs!!   ! Biomes including shrubs (active=1 only) 
!!mgs!!   ! ------------------------
!!mgs!!   ! 
!!mgs!!   REAL(dp) :: shrub(combimax) = (/ &       
!!mgs!!        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
!!mgs!!        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
!!mgs!!        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &       
!!mgs!!        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
!!mgs!!        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)

  !----------------------------------------------------------------------------
  ! solspe --> SOIL CARACTERISTICS:
  ! -------------------------------
  !
  !  ZOBLER texture classes |
  !
  ! SOLSPE: for 4 populations : values = 3*(Dmed sig p); ratio of fluxes; 
  !                                      residual moisture
  
  !   Populations: Coarse sand, medium/fine sand, Silt, Clay 
  !
  !     soil type 1 : Coarse     
  !     soil type 2 : Medium     
  !     soil type 3 : Fine     
  !     soil type 4 : Coarse Medium      
  !     soil type 5 : Coarse Fine 
  !     soil type 6 : Medium Fine 
  !     soil type 7 : Coarse_dp, Medium_dp, Fine 
  !     soil type 8 : Organic
  !     soil type 9 : Ice
  !     soil type 10 : Potential Lakes (additional)
  !     soil type 11 : Potential Lakes (clay)
  !     soil type 12 : Potential Lakes Australia
  !     soil type 13 : Taklimakan desert
  !     soil type 14 : Asian Loess
  !     soil type 15 : Asian Gobi
  !     soil type 16 : Asian other mixture
  !     soil type 17 : Asian desert
  !----------------------------------------------------------------------------

  solspe = &
       RESHAPE ( (/ &
       0.0707_dp, 2.0_dp, 0.43_dp, 0.0158_dp, 2.0_dp, 0.40_dp, 0.0015_dp, &!1 
       2.0_dp, 0.17_dp, 0.0002_dp, 2.0_dp, 0.00_dp, 2.1e-06_dp, 0.20_dp,  &  
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.37_dp, 0.0015_dp, &!2 
       2.0_dp, 0.33_dp, 0.0002_dp, 2.0_dp, 0.30_dp, 4.0e-06_dp, 0.25_dp,  &  
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.00_dp, 0.0015_dp, &!3 
       2.0_dp, 0.33_dp, 0.0002_dp, 2.0_dp, 0.67_dp, 1.0e-07_dp, 0.50_dp,  &   
       0.0707_dp, 2.0_dp, 0.10_dp, 0.0158_dp, 2.0_dp, 0.50_dp, 0.0015_dp, &!4 
       2.0_dp, 0.20_dp, 0.0002_dp, 2.0_dp, 0.20_dp, 2.7e-06_dp, 0.23_dp,  &   
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.50_dp, 0.0015_dp, &!5 
       2.0_dp, 0.12_dp, 0.0002_dp, 2.0_dp, 0.38_dp, 2.8e-06_dp, 0.25_dp,  &   
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.27_dp, 0.0015_dp, &!6 
       2.0_dp, 0.25_dp, 0.0002_dp, 2.0_dp, 0.48_dp, 1.0e-07_dp, 0.36_dp,  &   
       0.0707_dp, 2.0_dp, 0.23_dp, 0.0158_dp, 2.0_dp, 0.23_dp, 0.0015_dp, &!7 
       2.0_dp, 0.19_dp, 0.0002_dp, 2.0_dp, 0.35_dp, 2.5e-06_dp, 0.25_dp,  &   
       0.0707_dp, 2.0_dp, 0.25_dp, 0.0158_dp, 2.0_dp, 0.25_dp, 0.0015_dp, &!8 
       2.0_dp, 0.25_dp, 0.0002_dp, 2.0_dp, 0.25_dp, 0.0e-00_dp, 0.50_dp,  &   
       0.0707_dp, 2.0_dp, 0.25_dp, 0.0158_dp, 2.0_dp, 0.25_dp, 0.0015_dp, &!9 
       2.0_dp, 0.25_dp, 0.0002_dp, 2.0_dp, 0.25_dp, 0.0e-00_dp, 0.50_dp,  &   
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.00_dp, 0.0015_dp, &!10: 100% silt 
       2.0_dp, 1.00_dp, 0.0002_dp, 2.0_dp, 0.00_dp, 1.0e-05_dp, 0.25_dp,  &
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.00_dp, 0.0015_dp, &!11: 100% clay
       2.0_dp, 0.00_dp, 0.0002_dp, 2.0_dp, 1.00_dp, 1.0e-05_dp, 0.25_dp,  & 
       0.0707_dp, 2.0_dp, 0.00_dp, 0.0158_dp, 2.0_dp, 0.00_dp, 0.0027_dp, &!12: 100% silt 
       2.0_dp, 1.00_dp, 0.0002_dp, 2.0_dp, 0.00_dp, 1.0e-05_dp, 0.25_dp,  & 
       0.0442_dp, 1.5_dp, 0.03_dp, 0.0084_dp, 1.5_dp, 0.85_dp, 0.0015_dp, &!13: taklimakan
       2.0_dp, 0.11_dp, 0.0002_dp, 2.0_dp, 0.02_dp, 1.9e-06_dp, 0.12_dp,  & 
       0.0450_dp, 1.5_dp, 0.00_dp, 0.0070_dp, 1.5_dp, 0.33_dp, 0.0015_dp, &!14: loess 
       2.0_dp, 0.50_dp, 0.0002_dp, 2.0_dp, 0.17_dp, 1.9e-04_dp, 0.15_dp,  & 
       0.0457_dp, 1.8_dp, 0.31_dp, 0.0086_dp, 1.5_dp, 0.22_dp, 0.0015_dp, &!15: gobi 
       2.0_dp, 0.34_dp, 0.0002_dp, 2.0_dp, 0.12_dp, 3.9e-05_dp, 0.13_dp,  & 
       0.0293_dp, 1.8_dp, 0.39_dp, 0.0090_dp, 1.5_dp, 0.16_dp, 0.0015_dp, &!16: other mixture soils 
       2.0_dp, 0.35_dp, 0.0002_dp, 2.0_dp, 0.10_dp, 3.1e-05_dp, 0.13_dp,  & 
       0.0305_dp, 1.5_dp, 0.46_dp, 0.0101_dp, 1.5_dp, 0.41_dp, 0.0015_dp, &!17: desert and sand land
       2.0_dp, 0.10_dp, 0.0002_dp, 2.0_dp, 0.03_dp, 2.8e-06_dp, 0.12_dp   &
        /), (/ 17, 14 /), ORDER=(/2,1/))


  END SUBROUTINE set_dust_data
!==================================================================================================== 
!>>csld #433
  SUBROUTINE set_bc_nduscale_reg
  
    ! Description:
    ! ------------
    ! *set_bc_nduscale_reg* define the boundary condition for reading the input file containing the different regions (#433)
    !
    ! Authors:
    ! --------
    ! C. Siegenthaler, C2SM-ETHZ, June 2015 : original source
    !
    ! Interface:
    ! ----------
    ! *set_bc_nduscale_reg* is called from *bgc_dust_initialize*.
    !
    ! Method:
    ! -------
    ! Reads regions defined in external netcdf file
    !
    ! The integer contained in the external netcdf file dust_regions.nc describes the following  8 regions
    !   (the integer correspond also to the index of the region in nduscale_reg)
    !-----------------------------------------------------------------
    !     1       All the other locations than the following regions
    !     2       North america
    !     3       South America 
    !     4       North Africa 
    !     5       South Africa 
    !     6       Middle East 
    !     7       Asia
    !     8       Australia 

    USE mo_boundary_condition,       ONLY: bc_nml, bc_define, p_bcast_bc, BC_BOTTOM,BC_REPLACE
    USE mo_external_field_processor, ONLY: EF_FILE, EF_LONLAT, EF_IGNOREYEAR, EF_NOINTER
    USE mo_decomposition,     ONLY: lc => local_decomposition, global_decomposition


    TYPE(bc_nml) :: bc_regint            

    ! allocate the global 2D tuning array nduscale_2d
    IF (.NOT. ALLOCATED(nduscale_2d))   ALLOCATE (nduscale_2d(lc%nproma, lc%ngpblks))    

    ! Set defaults
    ibc_regint   = -1

    ! define the BC for regions
    bc_regint%bc_domain = BC_BOTTOM
    bc_regint%ef_type = EF_FILE
    bc_regint%ef_template = 'dust_regions.nc' 
    bc_regint%ef_varname = 'regions'
    bc_regint%ef_geometry = EF_LONLAT
    bc_regint%ef_timedef = EF_IGNOREYEAR
    bc_regint%ef_interpolate = EF_NOINTER 
    bc_regint%bc_mode = BC_REPLACE

    IF (p_parallel) THEN
       CALL p_bcast_bc (bc_regint,   p_io)
    ENDIF

    ! assign index for BC
    ibc_regint   = bc_define('regions for different nduscale', bc_regint, 2, .TRUE.)
    
     END SUBROUTINE set_bc_nduscale_reg
     

     SUBROUTINE comp_nduscale_reg(kproma, krow)  
    ! Description:
    ! ------------
    ! *comp_nduscale_reg* computes a location-dependent tuning factor nduscale to resolve Issue #433
    !
    ! Authors:
    ! --------
    ! C. Siegenthaler, C2SM-ETHZ, June 2015 : original source
    !
    ! Interface:
    ! ----------
    ! *comp_nduscale_reg* is called from *bgc_dust_calc_emiss*.
    !
    ! Method:
    ! -------
    ! Get regions defined in external netcdf file through bc_apply, derive 2D nduscale by
    ! assigning the appropriate nduscale for each location
    !
       USE mo_boundary_condition,  ONLY: bc_apply
       
       INTEGER, INTENT(in)          :: kproma, krow
       
       INTEGER        ::  jreg     
       INTEGER        ::  nreg     = 8   ! number of regions
       
       REAL(dp)       ::  zregint(kproma) 
 
       
       ! get the regions integer from the input file
       CALL bc_apply(ibc_regint,   kproma, krow, zregint)
       
       ! assign the right tuning factor for each location depending on in which region it is lying
       DO jreg=1,nreg
          nduscale_2d(1:kproma,krow)=MERGE(nduscale_reg(jreg), nduscale_2d(1:kproma,krow), INT(zregint) .EQ. jreg)
       END DO

  END SUBROUTINE comp_nduscale_reg
!<<csld #433
!====================================================================================================

  SUBROUTINE bgc_read_annual_fields

    ! Description:
    ! ------------
    ! *bgc_read_annual_fields* reads the prescribed annual field of preferential dust source area
    !
    ! Authors:
    ! --------
    ! M. Werner, MPI BGC, November 2002: original source
    ! P. Stier,  MPI Met, January  2004: replaced by fast NetCDF read-in
    ! M. Werner, MPI BGC, October  2004: changed emission cheme from 0.5x0.5 degree grid to 
    !                                    default Echam5-HAM grid (to save computational costs)
    !                                 => annual fields of different biome types are not read
    !                                    in this version any longer but may be added again in future releases
    !
    ! Interface:
    ! ----------
    ! *bgc_read_annual_fields* is called from *call_read_forcing* in file *call_submodels.f90*.
    !
    ! Method:
    ! -------
    ! Reads biome and soil fields specified in the three different namlists 
    ! *biome_fields.nml*, *soil_types.nml*, *pot_sources.nml*
    ! and stores them in the corresponding variables.
    !
    ! This routine is based on *xt_read_emiss_fields* by Philip Stier,  MPI-Met, HH.
    ! 

    USE mo_exception,         ONLY: finish
    USE mo_control,           ONLY: nlon, ngl
    USE mo_read_netcdf77,     ONLY: read_var_nf77_2d, read_var_nf77_3d
    USE mo_decomposition,     ONLY: lc => local_decomposition, global_decomposition
    USE mo_transpose,         ONLY: scatter_gp
    USE mo_exception,         ONLY: message, message_text, em_info, em_param
    USE mo_util_string,       ONLY: separator

    IMPLICIT NONE

    LOGICAL        :: lex                                 ! flag to test if file exists
    INTEGER        :: ierr
    INTEGER        :: i, j
    
    CHARACTER(512)   :: cfile                               ! filename

    !>>Kai Zhang, 2009-02
    REAL(dp)       :: frac_eastasia(lc%nproma)            ! total fraction of east asia soil type (from tiantao) 
    !<<Kai Zhang, 2009-02

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:), zin_degen_time(:,:,:), zin_all(:,:,:,:)
    REAL(dp), POINTER :: gl_data(:,:)

    IF (.NOT. ALLOCATED(biome))    ALLOCATE (biome(nlon, ngl))
    IF (.NOT. ALLOCATED(mat_s1))   ALLOCATE (mat_s1(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(mat_s2))   ALLOCATE (mat_s2(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(mat_s3))   ALLOCATE (mat_s3(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(mat_s4))   ALLOCATE (mat_s4(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(mat_s6))   ALLOCATE (mat_s6(lc%nproma, lc%ngpblks))
!   read asian soil type data (T. Cheng)
    IF (.NOT. ALLOCATED(mat_s13))  ALLOCATE (mat_s13(lc%nproma, lc%ngpblks))    
    IF (.NOT. ALLOCATED(mat_s14))  ALLOCATE (mat_s14(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(mat_s15))  ALLOCATE (mat_s15(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(mat_s16))  ALLOCATE (mat_s16(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(mat_s17))  ALLOCATE (mat_s17(lc%nproma, lc%ngpblks))

    IF (.NOT. ALLOCATED(mat_psrc)) ALLOCATE (mat_psrc(lc%nproma, lc%ngpblks))

    !>>Bernd Heinold, 2014-03 #382
    IF (ndust == 5) THEN
       IF (.NOT. ALLOCATED(mat_msg))  ALLOCATE (mat_msg(lc%nproma, lc%ngpblks))
    ENDIF

    IF (.NOT. ALLOCATED(Z01)) ALLOCATE (Z01(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(Z02)) ALLOCATE (Z02(lc%nproma, lc%ngpblks))
    !<<Bernd Heinold, 2014-03 #382

    IF (p_parallel_io) THEN
      ALLOCATE (zin(lc%nlon,lc%nlat))
      ALLOCATE (zin_degen_time(lc%nlon,lc%nlat,1))    !SF #227
      ALLOCATE (zin_all(nsoil_types,lc%nlon,lc%nlat,1)) !SF #253
      !SFNote: Both zin_degen_time and zin_all have a 'time' degenerate dimension to account for
      !        the new dust_preferential_sources.nc (#227) and the new merged soil types file (#253)
      !        in waiting for a proper handling by the boundary condition scheme (ToDo)
    END IF

    !<<SF #253
    !--- Map soil type matrices onto mat_all
    DO i=1,nsoil_types
       NULLIFY(mat_all(i)%ptr)
    ENDDO
    mat_all(1)%ptr => mat_s2(:,:)
    mat_all(2)%ptr => mat_s3(:,:)
    mat_all(3)%ptr => mat_s4(:,:)
    mat_all(4)%ptr => mat_s6(:,:)
    mat_all(5)%ptr => mat_s13(:,:)
    mat_all(6)%ptr => mat_s14(:,:)
    mat_all(7)%ptr => mat_s15(:,:)
    mat_all(8)%ptr => mat_s16(:,:)
    mat_all(9)%ptr => mat_s17(:,:)
    !<<SF #253

    !--- Set biome type distribution:

    biome(:,:)=21._dp      ! Biome type is set to dessert (= biome type 21) for all grid points

    !--- Set initial soil type distribution:

    mat_s1(:,:)=1._dp      ! Soil type is initialized to coarse texture (= soil type 1) for all grid points

    !--- Read other soil type distributions:
    
    CALL message('',separator)
    CALL message('bgc_read_annual_fields','Reading soil type distributions',level=em_info)

    IF (p_parallel_io) THEN
       cfile='soil_type_all.nc'
       INQUIRE (file=TRIM(cfile), exist=lex)
       IF (lex) THEN
          DO i=1,nsoil_types
             WRITE(message_text,'(a,a,a,a)') 'reading soil type #',ADJUSTL(TRIM(jtypes(i))), &
                                             ' distribution from: ', TRIM(cfile)
             CALL message('',message_text,level=em_param)
             CALL read_var_nf77_3d (cfile, "lon", "lat", "time", "type"//ADJUSTL(TRIM(jtypes(i))), &
                                    zin_all(i,:,:,:), ierr)
          ENDDO
       ELSE
          CALL finish('bgc_read_annual_fields','missing input file '//TRIM(cfile))
       END IF
    END IF
   
    DO i=1,nsoil_types
       NULLIFY(gl_data)
       IF (p_pe == p_io) gl_data => zin_all(i,:,:,1)
       CALL scatter_gp (gl_data, mat_all(i)%ptr(:,:), global_decomposition)
    ENDDO

    !--- Read preferential sources distribution:

    IF (p_parallel_io) THEN
       cfile='dust_preferential_sources.nc' !SF partial resolution of #227
       WRITE(message_text,'(a,a)') 'reading potential sources distribution ', TRIM(cfile)
       CALL message('',message_text,level=em_param)
       INQUIRE (file=TRIM(cfile), exist=lex)
       IF (lex) THEN
          CALL read_var_nf77_3d (cfile, "lon", "lat", "time", "source", zin_degen_time, ierr)
       ELSE
          CALL finish('bgc_read_annual_fields','missing input file '//TRIM(cfile))
       END IF
    END IF

    NULLIFY (gl_data)
    IF (p_pe == p_io) gl_data => zin_degen_time(:,:,1)

    CALL scatter_gp (gl_data, mat_psrc(:,:), global_decomposition)

    !>>Bernd Heinold, 2014-03 #382
       !--- Read Saharan dust source activation frequency map from MSG-SEVIRI dust index
       !    (Schepanski et al., GRL 2007; RSE 2012):
    IF (ndust == 5) THEN

       IF (p_parallel_io) THEN
          cfile='dust_msg_pot_sources.nc'
          WRITE(message_text,'(a,a)') 'reading potential Saharan sources', TRIM(cfile)
          CALL message('',message_text,level=em_param)
          INQUIRE (file=TRIM(cfile), exist=lex)
          IF (lex) THEN
             CALL read_var_nf77_2d (cfile, "lon", "lat", "dsaf", zin, ierr)
          ELSE
             CALL finish('bgc_read_annual_fields','missing input file '//TRIM(cfile))
          END IF
       END IF
      
       NULLIFY (gl_data)
       IF (p_pe == p_io) gl_data => zin(:,:)
       CALL scatter_gp (gl_data, mat_msg(:,:), global_decomposition)

    ENDIF
    !<<Bernd Heinold, 2014-03 #382

    IF (p_parallel_io) THEN
      DEALLOCATE (zin)
      DEALLOCATE (zin_degen_time)
      DEALLOCATE(zin_all)
    ENDIF

    CALL message('',separator)

    !>>Kai Zhang, 2009-02

    !if East Asia soil type is used, set original soil type fraction to zero in this region

    IF(k_dust_easo.eq.2) THEN

       DO j = 1,lc%ngpblks !krow

       frac_eastasia(:) = mat_s13(:,j) + mat_s14(:,j) + mat_s15(:,j) + mat_s16(:,j) + mat_s17(:,j) 

       DO i = 1,lc%nproma  !nproma 
          IF (frac_eastasia(i).gt.0) THEN
             mat_s1  (i,j) = frac_eastasia(i) 
             mat_s2  (i,j) = 0._dp
             mat_s3  (i,j) = 0._dp
             mat_s4  (i,j) = 0._dp
             mat_s6  (i,j) = 0._dp
             mat_psrc(i,j) = 0._dp
          END IF
       END DO 

       END DO 

    END IF 

    !<<Kai Zhang, 2009-02

    !>>Bernd Heinold, 2014-03 #382
     
    ! [ndust = 5]:
    ! Alternative approach for the Tegen et al. (2002) scheme using satellite observed dust-source-activation 
    ! frequencies to prescribe potential dust sources in the Sahara. The DSAF map is based on the IR
    ! dust index product of SEVIRI observations aboard the Meteosat Second Generation (MSG) satellite 
    ! (Schepanski et al., GRL 2007; RSE 2012).
    !
    ! Dust emission is calculated for grid cells, where at least 1 percent of the observations from March 2006 to
    ! February 2010 show dust source activations. The surface roughness in those areas is set to a constant 
    ! value of r_dust_z0s = 0.001 cm. 

    !--- Adapt soil type and roughness settings in the Saharan region

    IF (ndust == 5) THEN

       DO j = 1,lc%ngpblks !krow
         DO i = 1,lc%nproma  !nproma

           IF (mat_msg(i,j).gt.0) THEN
             mat_s1  (i,j) = 0._dp
             mat_s2  (i,j) = 0._dp
             mat_s3  (i,j) = 0._dp
             mat_s4  (i,j) = 0._dp
             mat_s6  (i,j) = 0._dp
             mat_psrc(i,j) = 0._dp

             IF (mat_msg(i,j).ge.0.01) THEN
               mat_psrc(i,j) = 1._dp
               Z01     (i,j) = r_dust_z0s
               Z02     (i,j) = r_dust_z0s
             END IF
           END IF

         END DO
       END DO

    END IF

    !<<Bernd Heinold, 2014-03 #382
    
    CALL set_bc_nduscale_reg !csld #433

  END SUBROUTINE bgc_read_annual_fields


!====================================================================================================

!!mgs!! cleanup of interface
  SUBROUTINE bgc_dust_read_monthly(kdate, ndurough)

    INTEGER, INTENT(in)       :: kdate
    REAL(dp), INTENT(in)      :: ndurough
   
    ! --- read effective monthly FPAR field
    CALL bgc_read_fpar_field(kdate)

    ! --- set constant surface roughness if desired
    IF (ndurough>0.0_dp) CALL bgc_set_constant_surf_rough(ndurough)
 
  END SUBROUTINE bgc_dust_read_monthly

!====================================================================================================

  SUBROUTINE bgc_read_fpar_field(kdate)

    ! Description:
    ! ------------
    ! *bgc_read_fpar_field* reads the effective FPAR field (NDVI satellite data or KERNEL output)
    !
    ! Authors:
    ! --------
    ! M. Werner, MPI BGC, November 2002: original source
    ! P. Stier,  MPI Met, January  2004: replaced by fast NetCDF read-in
    ! M. Werner, MPI BGC, October  2004: changed emission cheme from 0.5x0.5 degree grid to 
    !                                    default Echam5-HAM grid (to save computational costs)
    !
    ! Interface:
    ! ----------
    ! *bgc_read_fpar_field* is called from *call_read_forcing* in file *call_submodels.f90*.
    !
    ! Method:
    ! -------
    ! Reads daily (or monthly) effective FPAR fields specified in the namlist 
    ! *fpar_eff.nml* and stores it in the variable *k_fpar_eff*.
    !

    USE mo_exception,     ONLY: finish
    USE mo_read_netcdf77, ONLY: read_var_hs_nf77_2d
    USE mo_decomposition, ONLY: lc => local_decomposition, global_decomposition
    USE mo_transpose,     ONLY: scatter_gp
    USE mo_exception,     ONLY: message, message_text, em_info
    USE mo_util_string,   ONLY: separator

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kdate
    INTEGER             :: ierr
    LOGICAL             :: lex    ! flag to test if file exists
    CHARACTER(512)      :: cfile

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:)
    REAL(dp), POINTER :: gl_data(:,:)


    IF (.NOT. ALLOCATED(k_fpar_eff)) ALLOCATE (k_fpar_eff(lc%nproma, lc%ngpblks))

    !>>T. Cheng 2006-2007

    !read surface roughness length data 

    IF (p_parallel_io) THEN
      ALLOCATE (zin(lc%nlon,lc%nlat))
    END IF

    CALL message('',separator)
    CALL message('','')

    IF (p_parallel_io) THEN

       !--- Read LAI fields

       cfile='dust_potential_sources.nc' !csld / SF partial resolution of #227

       WRITE(message_text,'(a,a)') 'reading LAI fields ', TRIM(cfile)
       CALL message('bgc_read_fpar_field',message_text,level=em_info)

       INQUIRE (file=TRIM(cfile), exist=lex)

       IF (lex) THEN

          CALL read_var_hs_nf77_2d (cfile, "lon", "lat", "time", &
                                    kdate, "pot_source", zin, ierr  ) !csld / SF partial resolution of #227

       ELSE

          CALL finish('bgc_read_fpar_field','missing input file '//TRIM(cfile))

       END IF

    END IF

    NULLIFY (gl_data)

    IF (p_pe == p_io) gl_data => zin(:,:)

    CALL scatter_gp (gl_data, k_fpar_eff(:,:), global_decomposition)


    !--- Read surface roughness length field:

    IF (p_parallel_io) THEN

       cfile='surface_rough_12m.nc'
       
       WRITE(message_text,'(a,a)') 'reading surface roughness field', TRIM(cfile)
       CALL message('bgc_read_fpar_field',message_text,level=em_info)

       INQUIRE (file=TRIM(cfile), exist=lex)

       IF (lex) THEN
          CALL read_var_hs_nf77_2d (cfile, "lon", "lat", "time",  &
                                    kdate, "surfrough",  zin, ierr  )
       ELSE
          CALL finish('bgc_read_fpar_field','missing input file '//TRIM(cfile))
       END IF

    END IF

    NULLIFY (gl_data)

    IF (p_pe == p_io) gl_data => zin(:,:)

    CALL scatter_gp (gl_data, Z02(:,:), global_decomposition)

    IF (p_parallel_io) THEN
    
       !---Read surface roughness length field on 720x360 degree grid:

       cfile='surface_rough_12m.nc'

       WRITE(message_text,'(a,a)') 'reading surface roughness field', TRIM(cfile)
       CALL message('bgc_read_fpar_field',message_text,level=em_info)

       INQUIRE (file=TRIM(cfile), exist=lex)

       IF (lex) THEN
          CALL read_var_hs_nf77_2d (cfile, "lon", "lat", "time",  &
                                 kdate, "surfrough", zin, ierr )
       ELSE
          CALL finish('bgc_read_fpar_field','missing input file '//TRIM(cfile))
       END IF

    END IF
 
    NULLIFY (gl_data)
 
    IF (p_pe == p_io) gl_data => zin(:,:)

    CALL scatter_gp (gl_data, Z01(:,:), global_decomposition)  

    IF (p_parallel_io) THEN
       DEALLOCATE (zin)
    ENDIF

    !<<T. Cheng 2006-2007

    !>>Kai Zhang, 2009-02

    Z01(:,:) = max (Z01(:,:), r_dust_z0min)
    Z02(:,:) = max (Z02(:,:), r_dust_z0min)
    
    Z01(:,:) = Z01(:,:) * r_dust_scz0
    Z02(:,:) = Z02(:,:) * r_dust_scz0

    !<<Kai Zhang, 2009-02

    CALL message('','')
    CALL message('',separator)

  END SUBROUTINE bgc_read_fpar_field
  

!====================================================================================================


  SUBROUTINE bgc_set_constant_surf_rough(roughness_length)

    ! Description:
    ! ------------
    ! *bgc_set_constant_surf_rough* overwrites the surface roughness length
    ! fields with a given constant surface roughness length (cm).
    !
    ! Authors:
    ! --------
    ! J. Kazil, MPI-M, November 2008
    !
    ! Interface:
    ! ----------
    ! *bgc_set_constant_surf_rough* is called from *call_read_forcing* in file *call_submodels.f90*.
    !
    ! Method:
    ! -------
    ! Fills the surface roughness fields Z01 and Z02 with a given constant.

    IMPLICIT NONE

    REAL(dp):: roughness_length
    
    Z01 = roughness_length
    Z02 = roughness_length
    
  END SUBROUTINE bgc_set_constant_surf_rough
  

!====================================================================================================


  SUBROUTINE setdust

    ! Description:
    ! ------------
    ! *setdust* gets pre-configured defaults and reads the ham_dustctl namelist for potential overwritting
    !
    ! author: m.schultz@fz-juelich.de

    USE mo_namelist,            ONLY: open_nml, position_nml, POSITIONED, MISSING
    USE mo_exception,           ONLY: finish, message, message_text, em_info, em_warn, em_param
    USE mo_submodel,            ONLY: print_value
    USE mo_util_string,         ONLY: separator

    INCLUDE 'ham_dustctl.inc'

    !--- Local variables

    INTEGER :: ierr, inml, iunit

    !--- Set defaults based on context !SF #479
    CALL get_dust_namelist_defaults
    
    !--- read namelist
    CALL message('',separator)
    CALL message('setdust', 'Reading namelist ham_dustctl...', level=em_info)

    IF (p_parallel_io) THEN
       inml = open_nml('namelist.echam')
       iunit = position_nml ('HAM_DUSTCTL', inml, status=ierr)
       SELECT CASE (ierr)
       CASE (POSITIONED)
          READ (iunit, ham_dustctl)
       CASE (MISSING)
          WRITE(message_text,'(a,i0,a)') 'Namelist ham_dustctl not found. Will use default values'
          CALL message('setdust', message_text, level=em_warn)
       CASE DEFAULT      ! LENGTH_ERROR or READ_ERROR
          WRITE(message_text,'(a,i0)') 'Namelist ham_dustctl not correctly read! ierr = ', ierr
          CALL finish('setdust', message_text)
       END SELECT
    ENDIF

    !--- 2) Broadcast over processors:
    IF (p_parallel) THEN
       CALL p_bcast (ndurough,  p_io)
       CALL p_bcast (nduscale_reg, p_io)    !csld #433
       CALL p_bcast (ndust,   p_io)
       CALL p_bcast (r_dust_lai,p_io)
       CALL p_bcast (r_dust_umin,p_io)
       CALL p_bcast (r_dust_z0s,p_io)
       CALL p_bcast (r_dust_scz0,p_io)
       CALL p_bcast (r_dust_z0min,p_io)
       CALL p_bcast (r_dust_sf13,p_io)
       CALL p_bcast (r_dust_sf14,p_io)
       CALL p_bcast (r_dust_sf15,p_io)
       CALL p_bcast (r_dust_sf16,p_io)
       CALL p_bcast (r_dust_sf17,p_io)
       CALL p_bcast (r_dust_af13,p_io)
       CALL p_bcast (r_dust_af14,p_io)
       CALL p_bcast (r_dust_af15,p_io)
       CALL p_bcast (r_dust_af16,p_io)
       CALL p_bcast (r_dust_af17,p_io)
       CALL p_bcast (k_dust_smst, p_io)
       CALL p_bcast (k_dust_easo, p_io)
    END IF

    !--- Report parameter settings
    CALL message('','',level=em_param)
    CALL message('setdust','Parameter settings for BGC dust emission model',level=em_info)
    CALL message('','',level=em_param)
    IF (ndurough==0.0_dp)  THEN
      CALL message('','Using satellite-derived surface roughness length map (ndurough = 0.)', level=em_param)
      CALL message('','(Prigent et al., JGR 2005)', level=em_param)
    ELSE
      CALL print_value(' Globally constant surface roughness length (ndurough)  = ', ndurough)
    ENDIF
    CALL message('', '', level=em_param)
!>>csld #433
    CALL message('', 'Regional threshold wind friction velocity:', level=em_param)
    CALL print_value('nduscale_reg(2) (North America) = ', nduscale_reg(2)) 
    CALL print_value('nduscale_reg(3) (South America) = ', nduscale_reg(3)) 
    CALL print_value('nduscale_reg(4) (North Africa)  = ', nduscale_reg(4)) 
    CALL print_value('nduscale_reg(5) (South Africa)  = ', nduscale_reg(5)) 
    CALL print_value('nduscale_reg(6) (Middle East)   = ', nduscale_reg(6))
    CALL print_value('nduscale_reg(7) (Asia)          = ', nduscale_reg(7)) 
    CALL print_value('nduscale_reg(8) (Australia)     = ', nduscale_reg(8))
    CALL print_value('nduscale_reg(1) (elsewhere)     = ', nduscale_reg(1)) 
!<<csld #433
    CALL message('', '', level=em_param)
    CALL print_value(' r_dust_lai      = ', r_dust_lai)
    CALL print_value(' r_dust_umin     = ', r_dust_umin)
    CALL print_value(' r_dust_z0s      = ', r_dust_z0s)
    CALL print_value(' r_dust_scz0     = ', r_dust_scz0)
    CALL print_value(' r_dust_z0min    = ', r_dust_z0min)
    CALL print_value(' r_dust_sf13     = ', r_dust_sf13)
    CALL print_value(' r_dust_sf14     = ', r_dust_sf14)
    CALL print_value(' r_dust_sf15     = ', r_dust_sf15)
    CALL print_value(' r_dust_sf16     = ', r_dust_sf16)
    CALL print_value(' r_dust_sf17     = ', r_dust_sf17)
    CALL print_value(' r_dust_af13     = ', r_dust_af13)
    CALL print_value(' r_dust_af14     = ', r_dust_af14)
    CALL print_value(' r_dust_af15     = ', r_dust_af15)
    CALL print_value(' r_dust_af16     = ', r_dust_af16)
    CALL print_value(' r_dust_af17     = ', r_dust_af17)
    CALL message('', '', level=em_param)
    CALL print_value(' k_dust_smst     = ', k_dust_smst)
    CALL print_value(' k_dust_easo     = ', k_dust_easo)

    CALL message('',separator)

  END SUBROUTINE setdust

!====================================================================================================

  SUBROUTINE get_dust_namelist_defaults

    ! Description:
    ! ------------
    ! *get_dust_namelist_defaults* applies context-specific defaults for the dust namelist parameters
    !
    ! Authors:
    ! --------
    ! Sylvaine Ferrachat (ETH Zurich), January 2016

    USE mo_control, ONLY: nn, lnudge, lcouple

    !--- Configuration for dust scheme: preconfigured settings
    SELECT CASE (ndust)
       CASE (2)
           !Tiantao's version with modified scale factor for threshold velocity
           ndurough        = 0.000_dp
           r_dust_lai      = 1.E-10_dp
           r_dust_umin     = 21._dp
           r_dust_z0s      = 0.001_dp
           r_dust_scz0     = 1._dp
           r_dust_z0min    = 1.E-05_dp
           k_dust_smst     = 0
           k_dust_easo     = 0
           nduscale_reg(:) = 0.68_dp
       CASE (3)
           !Stier et al. (2005)
           ndurough     = 0.001_dp
           r_dust_lai   = 1.E-10_dp
           r_dust_umin  = 21._dp
           r_dust_z0s   = 0.001_dp
           r_dust_scz0  = 1._dp
           r_dust_z0min = 1.E-05_dp
           k_dust_smst  = 1
           k_dust_easo  = 1
    
           ! Choose the parameter for the threshold wind friction velocity based on
           ! the horizontal resolution: This parameterization was obtained by tuning
           ! nduscale_reg at the horizontal resolutions T21 and T42 to reproduce the
           ! total annual dust emissions at T63 with nduscale_reg(:)= 0.86 in nudged year
           ! 2000 runs.
           !
           ! THIS PARAMETERIZATION WILL LIKELY PRODUCE ERRONEOUS RESULTS FOR
           ! HORIZONTAL RESOLUTIONS > T63, FOR WHICH nduscale_reg MUST BE RE-TUNED!
           !
           ! jan.kazil@noaa.gov 2009-02-26 21:34:11 -07:00
           nduscale_reg(:)= -7.93650E-05_dp*dble(nn)**2.0_dp + 0.0095238_dp*dble(nn) + 0.575_dp
    
           IF(nn.gt.63) nduscale_reg(:)= 0.86_dp

       CASE (4, 5)
           ! Stier et al. (2005) + East Asia soil properties 
           !
           ! OR 
           !
           ! Stier et al. (2005) with Saharan dust source
           ! activation map + East Asia soil properties

           ndurough     = 0.001_dp
           r_dust_lai   = 1.E-1_dp
           r_dust_umin  = 21._dp
           r_dust_z0s   = 0.001_dp
           r_dust_scz0  = 1._dp
           r_dust_z0min = 1.E-05_dp
           r_dust_af14  = 1.E-06_dp
           r_dust_sf13  = 0.6
           k_dust_smst  = 1
           k_dust_easo  = 2

           SELECT CASE (nn)
               CASE (63)
                   IF (lnudge) THEN
                      nduscale_reg(:) = (/ &
                                        0.95_dp, & ! All the other locations than the following regions
                                        1.25_dp, & ! North America
                                        1.25_dp, & ! South America
                                        0.95_dp, & ! North Africa
                                        0.95_dp, & ! South Africa
                                        0.95_dp, & ! Middle East
                                        1.25_dp, & ! Asia
                                        0.95_dp  & ! Australia
                                       /)
                   ELSE
                      nduscale_reg(:) = (/ &
                                        1.05_dp, & ! All the other locations than the following regions
                                        1.45_dp, & ! North America
                                        1.45_dp, & ! South America
                                        1.05_dp, & ! North Africa
                                        1.05_dp, & ! South Africa
                                        1.05_dp, & ! Middle East
                                        1.45_dp, & ! Asia
                                        1.05_dp  & ! Australia
                                       /)
                   ENDIF
               CASE DEFAULT
                   nduscale_reg(:)= 0.86_dp
           END SELECT

           !>>SF obsolete code. Kept as example for potential resurrection
           !IF (lcouple) THEN
           !   nduscale_reg(:) = 0.92_dp  !checked with T63L47, by Lorenzo Tomassini (ZMAW) on 2013-01
           !ENDIF
           !
           !IF (lnudge) THEN
           !   nduscale_reg(:)= -7.93650E-05_dp*dble(nn)**2.0_dp + 0.0095238_dp*dble(nn) + 0.575_dp
           !   IF(nn.ge.63) nduscale_reg(:)= 0.86_dp
           !ENDIF    
           !<<SF

    END SELECT

  END SUBROUTINE get_dust_namelist_defaults

!====================================================================================================


  SUBROUTINE bgc_dust_initialize

    ! Description:
    ! ------------
    ! *bgc_dust_initialize* initializes several parameters and variables of the MPI BGC dust emission scheme.
    ! It also performs the wind independent calculations of threshold friction velocity and soil particle
    ! distribution. The effective soil fraction for dust deflation is determined, too.
    !
    ! Authors:
    ! --------
    ! M. Werner, MPI BGC, November 2002
    ! P. Stier, MPI MET,  June     2004: added scale factor for wind stress threshold
    ! M. Werner, MPI BGC, October  2004: changed emission cheme from 0.5x0.5 degree grid to 
    !                                    default Echam5-HAM grid (to save computational costs)
    !
    ! Interface:
    ! ----------
    ! *bgc_dust_initialize* is called from *call_read_forcing* in file *call_submodels.f90*.
    !
    ! This routine is based on the offline dust emission scheme by Ina Tegen, MPI BGC.
    ! 
    
    USE mo_math_constants,     ONLY: pi
    USE mo_physical_constants, ONLY: grav
    USE mo_control,       ONLY: nlon, ngl
    USE mo_decomposition, ONLY: lc => local_decomposition
    
    IMPLICIT NONE
    
    REAL(dp), PARAMETER :: gravi=grav*100._dp              ! gravity acceleration in cm/s2

    REAL(dp)            :: AAA, BB, CCC, DDD, EE, FF    ! tool variables
    
    REAL(dp)            :: dlastj(nclass)
    REAL(dp)            :: rdp

    INTEGER             :: i, j                         ! do indices
    INTEGER             :: kk                           ! do index
    INTEGER             :: nd, nsi, np
    INTEGER             :: ns                           ! do index

    REAL(dp)            :: rsize(nclass)                ! granulometric class size
    REAL(dp)            :: stotal                       ! total surface
    REAL(dp)            :: stotalV
    REAL(dp)            :: su                           ! surface 
    REAL(dp)            :: suV
    REAL(dp)            :: su_loc                       ! surface
    REAL(dp)            :: su_locV
    REAL(dp)            :: su_class(nclass)             ! surface occupied by each granulometric class
    REAL(dp)            :: su_classV(nclass)
    REAL(dp)            :: sum_srel(nats,nclass)        ! srel accumulated

    REAL(dp)            :: utest(nats)

    INTEGER             :: vgtp        

    REAL(dp)            :: xnV
    REAL(dp)            :: xk,xl,xm,xn                  ! tool variables

!---

    ! --- former mo_bgc_dust_data
    CALL set_dust_data

    ! --- initialisation now called from mo_ham_init
    CALL bgc_prestat                  ! read in dust emission statistics
    CALL bgc_read_annual_fields

    ! --- original code from here on
    IF (.NOT. ALLOCATED(flux_ann))  ALLOCATE (flux_ann(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(flux_a1))   ALLOCATE (flux_a1(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(flux_a2))   ALLOCATE (flux_a2(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(flux_a10))  ALLOCATE (flux_a10(lc%nproma, lc%ngpblks))
    IF (.NOT. ALLOCATED(idust))     ALLOCATE (idust(nlon, ngl))

    flux_ann(:,:) = 0._dp
    flux_a1 (:,:) = 0._dp
    flux_a2 (:,:) = 0._dp
    flux_a10(:,:) = 0._dp
    rsize     (:) = 0._dp
    srel    (:,:) = 0._dp
    srelV   (:,:) = 0._dp
    sum_srel(:,:) = 0._dp
    su_srelV(:,:) = 0._dp
    Uth       (:) = 0._dp
    utest     (:) = 0._dp

    i = 0

    rdp = Dmin

    ! calculation of the threshold friction velocity Uth, mat-95 formula (6)-(7)

    DO WHILE(rdp.lE.Dmax + 1.E-5_dp)

       i = i + 1

       rsize(i) = rdp

       BB = a_rnolds * (rdp ** x_rnolds) + b_rnolds  !formula (5) in mat-95

       AAA = SQRT(rop * gravi * rdp / roa)           !formula after (4) in mat-95

       CCC = SQRT(1._dp + d_thrsld /(rdp ** 2.5_dp)) !formula after (4) in mat-95

       IF (BB.LT.10._dp) THEN                        !formula (6) in mat-95
          DDD=SQRT(1.928_dp * (BB ** 0.092_dp) - 1._dp)
          Uth(i) = 0.129_dp * AAA * CCC / DDD
       ELSE                                          !formula (7) in mat-95
          EE = -0.0617_dp * (BB - 10._dp) 
          FF = 1._dp -0.0858_dp * EXP(EE)
!gf          Uth(i) = 0.12_dp * AAA * CCC * FF       
          Uth(i) = 0.129_dp * AAA * CCC * FF      
       ENDIF

       rdp = rdp * EXP(Dstep)   

    END DO      

    ! calculation of the soil particle distribution and related surfaces

    DO ns = 1,nats                                      ! loop over all soil types
       
       rdp = Dmin                                       
       kk = 0
       stotal  = 0._dp
       stotalV = 0._dp
       su_class (:) = 0._dp
       su_classV(:) = 0._dp

       DO WHILE (rdp.LE.Dmax+1.E-5_dp)                  ! surface calculations

          kk = kk + 1
          su  = 0._dp
          suV = 0._dp
          
          DO i = 1, nmode
          
             nd  = ((i - 1) *3 ) + 1
             nsi = nd + 1
             np  = nd + 2

             IF (solspe(ns,nd).EQ.0._dp) THEN            
                su_loc = 0._dp
                su_locV= 0._dp
             ELSE
                xk = solspe(ns,np)/(sqrt(2._dp* pi)*log(solspe(ns,nsi)))
                xl = ((log(rdp)-log(solspe(ns,nd)))**2)/(2._dp*(log(solspe(ns,nsi)))**2)
                xm = xk * exp(-xl)   !mass size distribution, formula (29) in mat-95
                xn = rop*(2._dp/3._dp)*(rdp/2._dp)  !surface, formula (30) in mat-95 
                xnV= 1._dp !volume
                su_loc  = (xm*Dstep/xn)       
                su_locV = (xm*Dstep/xnV)     
             ENDIF !

             su = su + su_loc
             suV = suV + su_locV

          END DO !Nmode

          su_class (kk) = su
          su_classV(kk) = suV
          stotal  = stotal  + su
          stotalV = stotalV + suV

          dlastj(kk) = rdp*5000._dp
          rdp = rdp * exp(Dstep)

       END DO !rdp
           
       DO j = 1,nclass !formula (32) in mat-95 

          IF (stotal.eq.0._dp) THEN
             srel(ns,j) = 0._dp
             srelV(ns,j) = 0._dp
          ELSE
             srel (ns,j) = su_class (j)/stotal
             srelV(ns,j) = su_classV(j)/stotalV
             utest(ns)   = utest(ns)+srelV(ns,j)
             su_srelV(ns,j) = utest(ns) !sum of nclasses
          ENDIF

       END DO !j=1,nclass
         
    END DO !ns (soil type)


    !currently not used in this code 

    DO j = 1,ngl                       ! loop over all latitudes
    DO i = 1,nlon                      ! loop over all longitudes
       idust(i,j) = 0                  ! check if vegetation type allows dust deflation
       vgtp = int(biome(i,j))          ! (idust=0: no deflation allowed; idust=1: deflation allowed)
       IF (vgtp.GT.0) idust(i,j) = active(vgtp)
    END DO !i
    END DO !j
    

    !>>Kai Zhang, 2009-02

    !the factor of F (horizontal flux) / G (vertical flux)

    solspe(13,nmode*3+1) = r_dust_af13
    solspe(14,nmode*3+1) = r_dust_af14
    solspe(15,nmode*3+1) = r_dust_af15
    solspe(16,nmode*3+1) = r_dust_af16
    solspe(17,nmode*3+1) = r_dust_af17

    !<<Kai Zhang, 2009-02

!!mgs, 2010-02!!    CALL bgc_dust_init_diag

  END SUBROUTINE bgc_dust_initialize


!====================================================================================================

  SUBROUTINE bgc_dust_init_diag(emi_stream)

  ! initialize dust emission diagnostics
  ! add field to emis stream (init_emi_stream in mo_emi_interface must have been called!)

    USE mo_linked_list,   ONLY: SURFACE
    USE mo_memory_base,   ONLY: t_stream, add_stream_element, SURFACE

    TYPE(t_stream), POINTER  :: emi_stream

    ! auxilliary variable ustar_acrit
!sschr: if not explicitly set (via lpost) ustar will not appear in output stream
    CALL add_stream_element (emi_stream, 'ustar_acrit', ustar_acrit, units='%',                   &
                             longname='time fraction with ustar > dust threshold', laccu=.TRUE.,  &
                             lpost=.TRUE., &
                             contnorest=.TRUE., leveltype=SURFACE, lrerun=.TRUE.)

  END SUBROUTINE bgc_dust_init_diag

!====================================================================================================

  SUBROUTINE bgc_dust_calc_emis(kproma, krow)

    ! Description:
    ! ------------
    ! *bgc_dust_calc_emis* calculates online global dust emission fields
    ! using the MPI BGC dust emissions scheme.
    !
    ! Authors:
    ! --------
    ! M. Werner, MPI BGC, October 2002
    ! P. Stier, MPI MET,  June    2004: added scale factor for wind stress threshold
    ! M. Werner, MPI BGC, October 2004: changed emission cheme from 0.5x0.5 degree grid to 
    !                                   default Echam5-HAM grid (to save computational costs)
    ! L. Kornblueh, MPI MET, July 2006: fix bug with respect to subscriptin in flux_6h
    !         
    !
    ! Interface:
    ! ----------
    ! The routine is called from *dust_emissions_bgc* in module *mo_bgc_dust_emis.f90*.
    !
    ! Method:
    ! -------
    ! The emissions are calculated every model time step using ECHAM 10m wind speeds.
    !
    ! This version of the emission scheme calculates dust emssissions for different bins, controlled by 
    ! dmin:   minimum particle diameter, 
    ! nbin:   number of binned size fractions, 
    ! ntrace: maximum number of bins, 
    ! dbmin(ntrace) & dbmax(ntrace): diameter limits of each bin
    !
    ! This routine is based on the offline dust emission scheme by Ina Tegen, MPI BGC.

    USE mo_physical_constants, ONLY: grav
    USE mo_memory_g3b,    ONLY: slm, glac, vlt, wl, ws, wsmx, sn, snc, orostd
    USE mo_physc2,        ONLY: cqsncr, cwlmax
    USE mo_time_control,  ONLY: delta_time
    USE mo_decomposition, ONLY: lc => local_decomposition
    USE mo_vphysc,        ONLY: vphysc
    USE mo_time_control,  ONLY: lstart,lresume !csld #433

    IMPLICIT NONE

! Parameter
    INTEGER, INTENT(in)          :: kproma, krow
    REAL(dp), PARAMETER          :: cd=1.00_dp*roa/(grav*100._dp) ! flux dimensioning parameter

    REAL(dp), PARAMETER          :: zepsec=1.E-12_dp    ! ECHAM parameter for snow cover calculation
    REAL(dp), PARAMETER          :: zsigfac=0.15_dp     ! ECHAM parameter for snow cover calculation
      
! Variables
    REAL(dp)           :: alpha

    REAL(dp)           :: dbmin(ntrace)              ! bin size limits
    REAL(dp)           :: dbmax(ntrace)              ! bin size limits
    REAL(dp)           :: dlast                      ! bin size limits
    REAL(dp)           :: rdp
    REAL(dp)           :: dpk(ntrace)
    INTEGER            :: dn(kproma)               ! number of dislocated particle classes
    INTEGER            :: dk(kproma,nclass)        ! dislocated particle classes
    
    REAL(dp)           :: du_snow(kproma)          ! ECHAM snow coverage interpolated to 0.5x0.5 grid resolution
    REAL(dp)           :: du_W1r(kproma)           ! ECHAM relative soil humidity interpolated to 0.5x0.5 grid resolution
    REAL(dp)           :: du_wind(kproma)          ! ECHAM wind field interpolated to 0.5x0.5 grid resolution

    INTEGER            :: dust_mask(kproma)        ! grid point mask for dust emissions (0: no emission)

!lk to reduce bank conflicts on SX-6 added 1 to kproma
    REAL(dp)           :: fluxtyp(kproma+1,nclass)
    REAL(dp)           :: fluxbin(kproma,ntrace)
    REAL(dp)           :: fdp1(nats), fdp2(nats)
    REAL(dp)           :: c_eff(kproma)             ! fraction efficiency
    
    REAL(dp)           :: fluxdiam1(kproma,nclass)  ! flux for soil type #1
    REAL(dp)           :: fluxdiam2(kproma,nclass)  ! flux for soil type #2
    REAL(dp)           :: fluxdiam3(kproma,nclass)  ! flux for soil type #3
    REAL(dp)           :: fluxdiam4(kproma,nclass)  ! flux for soil type #4
    REAL(dp)           :: fluxdiam6(kproma,nclass)  ! flux for soil type #6
    REAL(dp)           :: fluxdiam_pf(kproma,nclass)! flux for preferential sources soil type
    REAL(dp)           :: fluxdiam13(kproma,nclass) ! flux for soil type #13
    REAL(dp)           :: fluxdiam14(kproma,nclass) ! flux for soil type #14
    REAL(dp)           :: fluxdiam15(kproma,nclass) ! flux for soil type #15
    REAL(dp)           :: fluxdiam16(kproma,nclass) ! flux for soil type #16
    REAL(dp)           :: fluxdiam17(kproma,nclass) ! flux for soil type #17
    
    INTEGER            :: i, j                      ! loop index
    INTEGER            :: i_soil                    ! soil type
 
    REAL(dp)           :: AAA, BB, CCC, DDD, EE, FF !tool variables
    REAL(dp)           :: d1                        !distance between obstacles
    REAL(dp)           :: feff
    
    INTEGER            :: kk,kkk                    ! loop index
    INTEGER            :: kkmin
        
    INTEGER            :: n,nn                      ! loop index
    
    REAL(dp)           :: Ustar,Ustar_d(kproma)     ! threshold friction velocity for saltation
    REAL(dp)           :: uthp
    
    REAL(dp)           :: zw1r(kproma)             ! ECHAM relative soil moisture
    REAL(dp)           :: zcvs(kproma)             ! ECHAM snow cover fraction
    REAL(dp)           :: zcvsc(kproma)            ! ECHAM snow cover fraction at the canopy
    REAL(dp)           :: zwlmx(kproma)            ! auxiliary variable for ECHAM snow cover calculation
    REAL(dp)           :: zcvw(kproma)             ! auxiliary variable for ECHAM snow cover calculation
    REAL(dp)           :: zsn_mm, zsigh            ! auxiliary variable for ECHAM snow cover calculation

    !>>Kai Zhang, 2009-02
    REAL(dp)           :: mf13(kproma) 
    REAL(dp)           :: mf14(kproma) 
    REAL(dp)           :: mf15(kproma) 
    REAL(dp)           :: mf16(kproma) 
    REAL(dp)           :: mf17(kproma) 
    REAL(dp)           :: utsc(kproma) 
    REAL(dp)           :: tmpa, tmpb, tmpc         ! tmp vars 
    !<<Kai Zhang, 2009-02

!---

    !>>csld #433 computation of some quantities deriving from input files read with the BC scheme
    IF (lstart .OR. lresume) THEN
       CALL comp_nduscale_reg(kproma, krow)
    END IF
    !<<csld


    IF (.NOT. ALLOCATED(flux_6h)) ALLOCATE (flux_6h(lc%nproma, ntrace, lc%ngpblks))

    du_wind(:) = vphysc%velo10m(1:kproma,krow)

    !>>SF #458 (replacing where statements)
    du_wind(:) = MERGE( &
                       0._dp, & ! set missing values to zero
                       du_wind(:), &
                       (du_wind(:) < 0._dp))
    !<<SF #458 (replacing where statements)

    ! calculate ECHAM snow cover fraction *zcvs* (algorithm taken from routine *physc.f90*)
    ! and ECHAM relative soil moisture *zw1r* (ratio *wl* to *wlmax*)

    DO i=1,kproma
     
      IF ((slm(i,krow).GT.0.5_dp).AND.(glac(i,krow).LT.0.5_dp)) THEN    ! land surface, but not a glacier 

         zw1r (i)= MIN(ws(i,krow)/wsmx(i,krow),1._dp)
         zwlmx(i)= cwlmax*(1._dp+vlt(i,krow))
         zcvw (i)= MIN(wl(i,krow)/zwlmx(i),1._dp)
         zsn_mm=1000._dp*sn(i,krow)
         zsigh=SQRT(zsn_mm/(zsn_mm+zepsec+zsigfac*orostd(i,krow)))
         zcvs(i)=cqsncr*TANH(zsn_mm/10._dp)*zsigh
         zcvsc(i)=MIN(1._dp,snc(i,krow)/(zwlmx(i)-cwlmax+EPSILON(1._dp)))

         IF (zcvs(i).LT.EPSILON(1._dp) .AND. zcvsc(i).GE.EPSILON(1._dp)) THEN
            zcvs(i)=zcvsc(i)
         END IF

      ELSEIF((slm(i,krow).GT.0.5_dp).AND.(glac(i,krow).GT.0.5_dp)) THEN ! glacier on land surface
         zw1r(i)= 1._dp
         zcvs(i)= 1._dp
      ELSE                                                        ! ocean grid point
         zw1r(i)= -9.e9_dp
         zcvs(i)= -9.e9_dp      
      END IF
      
    END DO !i

    du_W1r (:) = zw1r(:)
    du_snow(:) = zcvs(:)

    !>>SF #458 (replacing where statements)
    du_W1r(:) = MERGE( &
                       0._dp, & ! set missing values to zero
                       du_W1r(:), &
                       (du_W1r(:) < 0._dp))

    du_snow(:) = MERGE( &
                       0._dp, & ! set missing values to zero
                       du_snow(:), &
                       (du_snow(:) < 0._dp))
    !<<SF #458 (replacing where statements)

    ! changed calculation for c_eff (T. Cheng)

      d1 = 0._dp
    feff = 0._dp
     AAA = 0._dp
      BB = 0._dp
     CCC = 0._dp
     DDD = 0._dp
      EE = 0._dp
      FF = 0._dp


    !get east-asia soil type fraction 

    mf13(1:kproma) = mat_s13(1:kproma,krow) !taklimakan
    mf14(1:kproma) = mat_s14(1:kproma,krow) !loess
    mf15(1:kproma) = mat_s15(1:kproma,krow) !gobi
    mf16(1:kproma) = mat_s16(1:kproma,krow) !other mixture soils
    mf17(1:kproma) = mat_s17(1:kproma,krow) !desert and sand land 

    !scale fator for Ut over East Asia desert (for sensitivity test) 

    utsc(:) = 1._dp

    !>>SF #458 (replacing where statements)
    utsc(:) = MERGE( &
                    r_dust_sf13, & !after r1060 (exp A0109)
                    utsc(:), &
                    (mf13(:) > 0._dp))

    utsc(:) = MERGE( &
                    r_dust_sf14, & !after r1060 (exp A0109)
                    utsc(:), &
                    (mf14(:) > 0._dp))

    utsc(:) = MERGE( &
                    r_dust_sf15, & !after r1060 (exp A0109)
                    utsc(:), &
                    (mf15(:) > 0._dp))

    utsc(:) = MERGE( &
                    r_dust_sf16, & !after r1060 (exp A0109)
                    utsc(:), &
                    (mf16(:) > 0._dp))

    utsc(:) = MERGE( &
                    r_dust_sf17, & !after r1060 (exp A0109)
                    utsc(:), &
                    (mf17(:) > 0._dp))
    !<<SF #458 (replacing where statements)

    !zero soil fraction for type 13-17 when k_dust_easo=1

    IF (k_dust_easo.eq.1) THEN 
       mf13(:) = 0._dp
       mf14(:) = 0._dp
       mf15(:) = 0._dp
       mf16(:) = 0._dp
       mf17(:) = 0._dp
    END IF


    !calculate the efficient friction velocity ratio feff, formula (17) in mat-95 
    !if use uniform z0 = 0.001cm, feff will be zero. 

    DO i=1,kproma
    
       IF (Z01(i,krow) .EQ. 0._dp .OR. Z02(i,krow) .EQ. 0._dp) THEN

          feff = 0._dp 

       ELSE

         AAA = log(Z01(i,krow)/r_dust_z0s)
          BB = log(aeff*(xeff/r_dust_z0s)**0.8_dp)
         CCC = 1._dp-AAA/BB

         IF (d1.eq.0._dp) THEN
            FF = 1._dp
         ELSE
            DDD = log(Z02(i,krow)/Z01(i,krow))
             EE = log(aeff * (d1/Z01(i,krow))**0.8_dp)
             FF = 1._dp- DDD/EE
         END IF

         feff = FF*CCC

         if (feff.lt.0._dp) feff=0._dp
         if (feff.gt.1._dp) feff=1._dp 

       END IF 

       c_eff(i) = feff
    
    END DO !DO i=1,kproma

    nn=1     
    rdp  =Dmin  
    dlast=Dmin
    dpk  (:)=0._dp
    dbmin(:)=0._dp
    dbmax(:)=0._dp

    DO kk=1,nclass                              ! assign fluxes to bins
       IF (mod(kk,nbin).eq.0) THEN
          dbmax(nn)=  rdp*10000._dp*0.5_dp      ! calculate bin minimum/maximum radius in um
          dbmin(nn)=dlast*10000._dp*0.5_dp     
          dpk(nn)=sqrt(dbmax(nn)*dbmin(nn))
          nn=nn+1
          dlast=rdp
       ENDIF
       rdp = rdp * exp(Dstep)
    ENDDO !kk      

    fluxbin(:,:)     =0._dp
    flux_6h(:,:,krow)=0._dp
       
    dust_mask(:)=0
     

    DO i = 1,kproma

    IF (c_eff(i).GT.0._dp) THEN    

       Ustar = (vk * du_wind(i) *100._dp)/(log(ZZ/Z02(i,krow)))  ! U*(cm/s): formular (15) in mat-95 

       ! check critical wind speed (wind stress threshold)

       IF ((Ustar .GT. 0._dp) .AND. (Ustar .GE. r_dust_umin*nduscale_2d(i,krow)/c_eff(i))) THEN  

          IF (k_fpar_eff(i,krow).GT.r_dust_lai) THEN   ! check if the grid cell is a potential dust source
             dust_mask(i)=1                            ! set grid point as a potential dust source point
             Ustar_d(i)=Ustar                          ! store wind speed of dust grid
          ENDIF ! k_fpar_eff.gt.0.

          ustar_acrit(i,krow)=ustar_acrit(i,krow)+delta_time

       ENDIF   ! Ustar

    ENDIF !IF (c_eff(i).GT.0._dp) THEN 

    ENDDO !DO i = 1,kproma
     

    fluxtyp(:,:)=0._dp
    fluxdiam1(:,:)=0._dp
    fluxdiam2(:,:)=0._dp
    fluxdiam3(:,:)=0._dp
    fluxdiam4(:,:)=0._dp
    fluxdiam6(:,:)=0._dp
    fluxdiam_pf(:,:)=0._dp
    fluxdiam13(:,:)=0._dp
    fluxdiam14(:,:)=0._dp
    fluxdiam15(:,:)=0._dp
    fluxdiam16(:,:)=0._dp
    fluxdiam17(:,:)=0._dp

    AAA = 0._dp
    BB  = 0._dp

    dn(:)=0


    DO kk=1,nclass


    DO i = 1,kproma 
   

    IF (c_eff(i).GT.0._dp) THEN 


    IF (dust_mask(i).EQ.1) THEN   

    !>>Kai Zhang, 2009-02

    !calculate the horizontal flux G(Dp), formula (28) in mat-95 

    IF (k_dust_smst.eq.0) THEN

       ! caculation of wet erosion threshold friction velocity (T. Cheng)
       ! replace ws with zmr 

       IF ((slm(i,krow).GT.0.5_dp).AND.(glac(i,krow).LT.0.5_dp)) THEN
          BB = MIN(ws(i,krow)/rop,1._dp) * 100._dp
       ELSE
          BB = 0._dp
       ENDIF
       
       DO j = 1, nats
       
          AAA = solspe(j,nspe)*100._dp
       
          IF ( BB .LE. AAA ) THEN
             uthp = Uth(kk)
          ELSE
             uthp = Uth(kk)*SQRT(1._dp+1.21_dp*(BB-AAA)**0.68_dp) 
          ENDIF
          uthp = uthp * nduscale_2d(i,krow)

           ! Marticorena:
           !
           !fdp1(j) = (1._dp+(uthp/(c_eff(i) * Ustar_d(i))))**2
           !fdp2(j) = (1._dp-(uthp/(c_eff(i) * Ustar_d(i))))
 
           !>>Kai Zhang, 2009-02
           fdp1(j) = (1._dp+(uthp*utsc(i)/(c_eff(i) * Ustar_d(i))))**2
           fdp2(j) = (1._dp-(uthp*utsc(i)/(c_eff(i) * Ustar_d(i))))
           !<<Kai Zhang, 2009-02

           ! Shao:
           ! 
           ! fdp1 = (1.-(Uthp/(c_eff(i) * Ustar_d(i)))**2)
           ! fdp2 = 1.
       
       ENDDO 


    ELSE 

       ! do not consider effect of soil moisture on threshold velocity 

       uthp = Uth(kk) * nduscale_2d(i,krow)
       
       DO j = 1, nats

         ! Marticorena:

           !>>Kai Zhang, 2009-02

           fdp1(j) = (1._dp+(uthp*utsc(i)/(c_eff(i) * Ustar_d(i))))**2
           fdp2(j) = (1._dp-(uthp*utsc(i)/(c_eff(i) * Ustar_d(i))))

           !<<Kai Zhang, 2009-02

         ! Shao:

         ! fdp1 = (1.-(Uthp/(c_eff(i) * Ustar_d(i)))**2)
         ! fdp2 = 1.

       ENDDO 


    ENDIF !IF (k_dust_smst.eq.0) THEN

    !<<Kai Zhang, 2009-02


    !calculate horizontal flux distribution as function of Dp
    !formular (33) in mat-95, formula (3) in tegen-2002

    tmpa = cd * Ustar_d(i)**3

    i_soil=1

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1) 
      fluxdiam1(i,kk) = srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha 
    ENDIF 

    i_soil=2

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1)
      fluxdiam2(i,kk) = srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha 
    ENDIF

    i_soil=3

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1)
      fluxdiam3(i,kk) = srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha 
    ENDIF

    i_soil=4

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1)
      fluxdiam4(i,kk) = srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha 
    ENDIF

    i_soil=6

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1)
      fluxdiam6(i,kk) = srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha 
    ENDIF

    i_soil=10

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1)
      fluxdiam_pf(i,kk)= srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha 
    ENDIF

    !add asian dust source emission (cheng)

    i_soil=13

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1)
      fluxdiam13(i,kk) = srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha
    ENDIF

    i_soil=14

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1)
      fluxdiam14(i,kk) = srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha
    ENDIF

    i_soil=15

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1)
      fluxdiam15(i,kk) = srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha
    ENDIF

    i_soil=16

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1)
      fluxdiam16(i,kk) = srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha
    ENDIF

    i_soil=17

    IF (fdp2(i_soil).gt.0._dp) THEN
      alpha= solspe(i_soil,nmode*3+1)
      fluxdiam17(i,kk) = srel(i_soil,kk) * fdp1(i_soil) * fdp2(i_soil) * tmpa *alpha
    ENDIF


    IF (kk.eq.1) THEN 
       
       tmpb = 1._dp-mat_psrc(i,krow)
       tmpc = (mat_s1(i,krow)-mat_s2(i,krow)-mat_s3(i,krow)-mat_s4(i,krow)-mat_s6(i,krow) & 
             - mf13(i)-mf14(i)-mf15(i)-mf16(i)-mf17(i))

       fluxtyp(i,kk) = fluxtyp(i,kk) &
         + fluxdiam1 (i,kk)*tmpb*tmpc &
         + fluxdiam2 (i,kk)*tmpb*mat_s2 (i,krow)                                      &
         + fluxdiam3 (i,kk)*tmpb*mat_s3 (i,krow)                                      &
         + fluxdiam4 (i,kk)*tmpb*mat_s4 (i,krow)                                      &
         + fluxdiam6 (i,kk)*tmpb*mat_s6 (i,krow)                                      &
         + fluxdiam13(i,kk)*tmpb*mf13(i)                                      &
         + fluxdiam14(i,kk)*tmpb*mf14(i)                                      &
         + fluxdiam15(i,kk)*tmpb*mf15(i)                                      &
         + fluxdiam16(i,kk)*tmpb*mf16(i)                                      &
         + fluxdiam17(i,kk)*tmpb*mf17(i)                                      &
         + fluxdiam_pf(i,kk)*mat_psrc(i,krow)

    ELSE

       dn(i)=dn(i)+1           ! increase number of dislocated particle classes
       dk(i,dn(i))=kk          ! store dislocated particle class

    ENDIF !IF (kk.eq.1) THEN 


    ENDIF !IF (dust_mask(i).EQ.1) THEN   

    ENDIF !IF (c_eff(i).GT.0._dp) THEN

    ENDDO !DO i = 1,kproma 

    ENDDO !DO kk=1,nclass

    !calculate horizontal flux distribution

    kkmin=1

    DO i = 1,kproma

    tmpa = 1._dp-mat_psrc(i,krow)

    tmpb = mat_s1(i,krow)-mat_s2(i,krow) &
          -mat_s3(i,krow)-mat_s4(i,krow)-mat_s6(i,krow) &
          -mf13(i)-mf14(i)-mf15(i)- mf16(i)-mf17(i)

    IF (dust_mask(i).EQ.1) THEN

       DO n=1,dn(i)                 ! loop over dislocated dust particle classes

       DO kkk=1,dk(i,n)             ! scaling with relative contribution of dust size  fraction


          i_soil=1
          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + tmpa*tmpb * fluxdiam1(i,dk(i,n))*tmpc

          i_soil=2
          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + tmpa*mat_s2(i,krow) * fluxdiam2(i,dk(i,n))*tmpc

          i_soil=3
          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + tmpa*mat_s3(i,krow) * fluxdiam3(i,dk(i,n))*tmpc

          i_soil=4
          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + tmpa*mat_s4(i,krow) * fluxdiam4(i,dk(i,n))*tmpc

          i_soil=6
          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + tmpa*mat_s6(i,krow) * fluxdiam6(i,dk(i,n))*tmpc

          IF(du_wind(i).gt.10._dp) THEN
             i_soil=11 ! flux from preferential source at high wind speeds
          ELSE
             i_soil=10 ! flux from preferential source at low wind speeds
          ENDIF

          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + mat_psrc(i,krow) * fluxdiam_pf(i,dk(i,n))*tmpc


          i_soil=13
          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + tmpa*mf13(i) * fluxdiam13(i,dk(i,n))*tmpc

          i_soil=14
          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + tmpa*mf14(i) * fluxdiam14(i,dk(i,n))*tmpc 

          i_soil=15
          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + tmpa*mf15(i) * fluxdiam15(i,dk(i,n))*tmpc

          i_soil=16
          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + tmpa*mf16(i) * fluxdiam16(i,dk(i,n))*tmpc

          i_soil=17
          tmpc = srelV(i_soil,kkk)/((su_srelV(i_soil,dk(i,n))-su_srelV(i_soil,kkmin)))
          fluxtyp(i,kkk) = fluxtyp(i,kkk) + tmpa*mf17(i) * fluxdiam17(i,dk(i,n))*tmpc


       ENDDO ! kkk
       ENDDO ! n

    ENDIF ! dust_mask.eq.1

    ENDDO ! i


    !calculate the vertical dust particle flux based on White (1979), formula (2) in tegen-2002

    DO nn=1,ntrace
       DO i = 1,kproma
          IF (dust_mask(i).EQ.1) THEN

!%%%%%%%%%%%%
!CDIR NOVETOR
!%%%%%%%%%%%%

          fluxbin(i,nn) = fluxbin(i,nn)+SUM(fluxtyp(i,(nn-1)*nbin+1:nn*nbin))

          ! mask out dust fluxes, where soil moisture threshold reached

          IF (du_W1r(i).gt.w0) fluxbin(i,nn)=0._dp             

          ! fluxbin: g/cm2/sec; flux_6h: g/m2/sec

          flux_6h(i,nn,krow)=fluxbin(i,nn)*10000._dp*(1._dp-du_snow(i))*k_fpar_eff(i,krow)

          ENDIF ! dust_mask.eq.1
       ENDDO ! i
    ENDDO ! nn

    DO nn=1,ntrace
       DO i=1,kproma
          tmpa = flux_6h(i,nn,krow)*delta_time
          flux_ann(i,krow)=flux_ann(i,krow)+tmpa                          ! flux_ann: g/m2/year
          IF (dpk(nn).le.10._dp) flux_a10(i,krow)=flux_a10(i,krow)+tmpa   ! flux_a10: g/m2/year
          IF (dpk(nn).le. 2._dp) flux_a2(i,krow) = flux_a2(i,krow)+tmpa   ! flux_a2:  g/m2/year
          IF (dpk(nn).le. 1._dp) flux_a1(i,krow) = flux_a1(i,krow)+tmpa   ! flux_a1:  g/m2/year
       ENDDO ! i                           
    ENDDO ! nn

  END SUBROUTINE bgc_dust_calc_emis


!====================================================================================================

  SUBROUTINE bgc_prestat

    ! Description:
    ! ------------
    ! *bgc_prestat* reads total global dust emission numbers from the text file "dust.diagnostics.dat"
    ! The purpose is to correctly continue monitoring total dust emissions after a (re)start of the model
    !
    ! Authors:
    ! --------
    ! M. Werner, MPI BGC, January 2003
    !
    ! Interface:
    ! ----------
    ! *bgc_prestat* is called from *call_submodels.f90*.
    !
    ! Method:
    ! -------
    ! Reads diagnostic dust sums from text file "dust.diagnostics.dat".
    ! If the text file does not exist, the dust sums are reset to zero (=initial start).

    USE mo_filename,          ONLY: find_next_free_unit

    IMPLICIT NONE

    INTEGER                           :: iunit   ! fortran unit for input file
    LOGICAL                           :: lex     ! flag to test if file exists
    INTEGER                           :: kyear, kmonth, kday, khour, kminute, ksecond

    INQUIRE(file='./dust.diagnostics.dat', exist=lex)  ! read dust diagnostic file
    IF (lex) THEN
     IF (p_parallel_io) THEN
      iunit = find_next_free_unit(10,99)
      OPEN(iunit,file='./dust.diagnostics.dat',status='old')      
20    READ(iunit,*,END=25) kyear, kmonth, kday, khour, kminute, ksecond, dustsum(1), dustsum(2), dustsum(3), dustsum(4) 
      GOTO 20
25    CLOSE(iunit)
      dustsum(:) = dustsum(:)*1.e12_dp              ! convert back from Mt/yr to g/yr
     ENDIF
    ELSE
     dustsum(:) = 0._dp                             ! if diagnostic file does not exist: initialize dust diagnostics
    ENDIF     

    !--- Broadcast field to all CPUs:

    CALL p_bcast(dustsum, p_io) 

  END SUBROUTINE bgc_prestat

  
!====================================================================================================


  SUBROUTINE bgc_dust_trastat(kyear, kmonth, kday, khour, kminute, ksecond)

    ! Description:
    ! ------------
    ! *bgc_trastat* writes total global dust emission numbers to the text file "dust.diagnostics.dat"
    ! The purpose is to correctly continue monitoring total dust emissions after a (re)start of the model
    !
    ! Authors:
    ! --------
    ! M. Werner, MPI BGC, January 2003
    ! M. Werner, MPI BGC, October 2004: changed emission cheme from 0.5x0.5 degree grid to 
    !                                    default Echam5-HAM grid (to save computational costs)
    !
    ! Interface:
    ! ----------
    ! *bgc_dust_trastat* is called from *call_submodels.f90*.
    !
    ! Method:
    ! -------
    ! Writes diagnostic dust sums to text file "dust.diagnostics.dat".
    ! If the text file does not exist, it is newly created (=cold start), 
    ! otherwise the output is appended (=restart).

    USE mo_math_constants,    ONLY: pi
    USE mo_control,           ONLY: nlon, ngl
    USE mo_io_units,          ONLY: nout
    USE mo_filename,          ONLY: find_next_free_unit
    USE mo_decomposition,     ONLY: gl_dc=> global_decomposition
    USE mo_transpose,         ONLY: gather_gp


    IMPLICIT NONE

    REAL(dp), POINTER :: zflux_ann(:,:)
    REAL(dp), POINTER :: zflux_a1(:,:)
    REAL(dp), POINTER :: zflux_a2(:,:)
    REAL(dp), POINTER :: zflux_a10(:,:)
    REAL(dp), POINTER :: zk_fpar_eff(:,:)

    INTEGER                           :: iunit   ! fortran unit for input file
    INTEGER, INTENT(IN)               :: kyear, kmonth, kday, khour, kminute, ksecond

    REAL(dp)           :: dustsumt, dustsum10, dustsum2, dustsum1  ! dummy values for diagnostics

    REAL(dp)           :: latitude    
    REAL(dp)           :: S
    REAL(dp)           :: surf(ngl)                  ! grid surface area
    REAL(dp)           :: tool
        
    INTEGER        :: x, y                       ! longitude / latitude index

    REAL(dp)             :: resol
    REAL(dp), PARAMETER  :: radius = 6371000._dp        ! radius of Earth


    resol = 180._dp/ngl

    dustsumt = dustsum(1)                        ! initialize dustsums with values from last (re)start
    dustsum10 = dustsum(2)                       ! initialize dustsums with values from last (re)start
    dustsum2 = dustsum(3)
    dustsum1 = dustsum(4)

    DO y = 1,ngl                                 ! calculate grid surface area for different latitudes
      latitude =(90._dp + resol * (0.5_dp - REAL(y,dp)) )
      S = radius * radius * pi / (180._dp/resol)
      tool = COS(latitude * pi/180._dp)
      surf(y) = S * abs(tool) * 2._dp*pi/(360._dp/resol)
    ENDDO

    IF (p_pe == p_io) THEN
       ALLOCATE (zflux_ann(nlon,ngl))
       ALLOCATE (zflux_a10(nlon,ngl))
       ALLOCATE (zflux_a2(nlon,ngl))
       ALLOCATE (zflux_a1(nlon,ngl))
       ALLOCATE (zk_fpar_eff(nlon,ngl))
    END IF

    CALL gather_gp (zflux_ann, flux_ann, gl_dc)
    CALL gather_gp (zflux_a10, flux_a10, gl_dc)
    CALL gather_gp (zflux_a2, flux_a2, gl_dc)
    CALL gather_gp (zflux_a1, flux_a1, gl_dc)
    CALL gather_gp (zk_fpar_eff, k_fpar_eff, gl_dc)

    IF (p_parallel_io) THEN                      ! write out diagnostics to standard output   
      DO y = 1,ngl
        DO x = 1,nlon
          IF (zk_fpar_eff(x,y).GT.1.E-10_dp) THEN        ! if the area is a potential dust source: sum up for diagnostic output
            dustsumt=dustsumt+zflux_ann(x,y)*surf(y)   ! total grams of dust per year
            dustsum10=dustsum10+zflux_a10(x,y)*surf(y) ! small dust particles < 10um
            dustsum2=dustsum2+zflux_a2(x,y)*surf(y)    ! small dust particles < 2um
            dustsum1=dustsum1+zflux_a1(x,y)*surf(y)    ! small dust particles < 1um
          ENDIF !! (k_fpar_eff.gt.1.e-10)
        END DO !! x
      END DO !! y

      DEALLOCATE (zflux_ann)
      DEALLOCATE (zflux_a10)
      DEALLOCATE (zflux_a2)
      DEALLOCATE (zflux_a1)
      DEALLOCATE (zk_fpar_eff)

     WRITE(nout,'(4(a,e10.4))')                                 &
                   ' DUST (Mt/yr): total: ',dustsumt*1.e-12_dp,    &
                   '     < 10um: ',dustsum10*1.e-12_dp,            &
                   '     < 2um: ',dustsum2*1.e-12_dp,              &
                   '     < 1um: ',dustsum1*1.e-12_dp
     iunit = find_next_free_unit(55,99)
     OPEN(iunit,file='dust.diagnostics.dat',status='unknown',position='append')
     WRITE(iunit,'(i0,i0,i0,i0,i0,i0,e25.15,e25.15,e25.15,e25.15)') kyear, kmonth, kday, khour, kminute, ksecond, &
                    dustsumt*1.e-12_dp, dustsum10*1.e-12_dp, dustsum2*1.e-12_dp, dustsum1*1.e-12_dp
     CLOSE(iunit)

    ENDIF

  END SUBROUTINE bgc_dust_trastat

  
!====================================================================================================


  SUBROUTINE bgc_dust_cleanup

    ! Description:
    ! ------------
    ! *bgc_dust_cleanup* deallocates memory needed for the dust scheme.
    !
    ! Authors:
    ! --------
    ! M. Werner, MPI BGC, November 2004
    !
    ! Interface:
    ! ----------
    ! *bgc_dust_cleanup* is called from *call_free_submodel_memory*
    
    IMPLICIT NONE

    IF (ALLOCATED(biome))       DEALLOCATE (biome)
    IF (ALLOCATED(mat_s1))      DEALLOCATE (mat_s1)
    IF (ALLOCATED(mat_s2))      DEALLOCATE (mat_s2)
    IF (ALLOCATED(mat_s3))      DEALLOCATE (mat_s3)
    IF (ALLOCATED(mat_s4))      DEALLOCATE (mat_s4)
    IF (ALLOCATED(mat_s6))      DEALLOCATE (mat_s6)
    IF (ALLOCATED(mat_psrc))    DEALLOCATE (mat_psrc)
    IF (ALLOCATED(mat_s13))     DEALLOCATE (mat_s13)
    IF (ALLOCATED(mat_s14))     DEALLOCATE (mat_s14)
    IF (ALLOCATED(mat_s15))     DEALLOCATE (mat_s15)
    IF (ALLOCATED(mat_s16))     DEALLOCATE (mat_s16)
    IF (ALLOCATED(mat_s17))     DEALLOCATE (mat_s17)
    IF (ALLOCATED(k_fpar_eff))  DEALLOCATE (k_fpar_eff)
    IF (ALLOCATED(Z01))         DEALLOCATE (Z01)
    IF (ALLOCATED(Z02))         DEALLOCATE (Z02) 
    IF (ALLOCATED(flux_ann))    DEALLOCATE (flux_ann)
    IF (ALLOCATED(flux_a1))     DEALLOCATE (flux_a1)
    IF (ALLOCATED(flux_a2))     DEALLOCATE (flux_a2)
    IF (ALLOCATED(flux_a10))    DEALLOCATE (flux_a10)
    IF (ALLOCATED(idust))       DEALLOCATE (idust)
    IF (ALLOCATED(flux_6h))     DEALLOCATE (flux_6h)
    IF (ALLOCATED(nduscale_2d)) DEALLOCATE (nduscale_2d) !csld #433

  END SUBROUTINE bgc_dust_cleanup

!====================================================================================================


END MODULE mo_ham_dust
