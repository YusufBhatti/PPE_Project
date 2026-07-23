!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_emi_biogenic.f90
!!
!! \brief
!! contains routines for the calculation of emissions of monoterpenes,
!! isoprene and other reactive volatile organic carbon (OVOC) species
!! from vegetation
!!
!! \author Martin Schultz (FZ Juelich)
!! \author Declan O'Donnell (FMI)
!!
!! \responsible_coder
!! Martin Schultz (FZ Juelich), m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# Adetutu Aghedo (MPI-Met) - original meganv1 code
!!   -# Declan O'Donnell (FMI) - original meganv2 code (2010-06)
!!   -# Martin Schultz (FZ Juelich) - rewrite to megan v2.10 using PFT specific EFs
!!                                  - also uses PAR radiation field from JSBACH now
!!
!! \limitations
!! A number of scientific issues remain which should be investigated for example in
!! a master or PhD study:
!! - MEGAN 2.10 uses a canopy model (EMPROC/canopy.f) where the temperature and light-
!!   dependent adjustment factors are computed for each PFT. In our model we neglect all
!!   canopy processes! Effectively we thus have MEGAN 2.04 implemented.
!! - emission factors and deposition: the 2012 paper mentions that MEGAN 2.10 emission factors
!!   are net primary emission factors which are corrected up to account for in-canopy deposition.
!!   It should be checked if the emission-deposition balance for "bidirectional" compounds
!!   (including methanol) is OK.
!! - vegetation map: it has been shown that MEGAN is very sensitive to the PFT distribution
!!   and gives different results when high-resolution, detailed plant species maps are used
!!   instead of broad-scale PFTs. This module should be extended to allow for input of more
!!   detailed maps, and it should be tested if the CLM4 PFTs could be replaced by the JSBACH
!!   PFTs (which are however variable).
!! - The emission adjustments for soil moisture and CO2 inhibition are not activated.
!! - equations 6, 9, 10 in Guenther et al., 2012, use 10-day averages in addition to daily
!!   averages. These have not been implemented yet.
!!
!! \details
!! The routines are using the MEGANv2.10 model by Alex Guenther, NCAR.
!!
!! \bibliographic_references
!! Guenther, A. et al., Geosci. Model Dev., 5, 1471-1492 (2012)
!! Guenther, A. et al., Atmos. Chem. Phys. 6, 3181-3210 (2006) + corrigendum from 2007
!! Guenther, A. et al., J. Geophys. Res. 100 D5, 8873-8892 (1995)
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

! ToDo:
!! -- diagnostic for "total gamma"
!! -- check for necessity of testing for loland(jl) !
!! -- testing !!
!! -- activate gamma_age routine
!!


MODULE mo_hammoz_emi_biogenic

  !---inherited types, data and functions
  USE mo_kind,                     ONLY: dp
  USE mo_tracdef,                  ONLY: ln, ll
  USE mo_linked_list,              ONLY: t_stream 
  USE mo_submodel_diag,            ONLY: vmem2d
  USE mo_decomposition,            ONLY: dc=>local_decomposition
  USE mo_external_field_processor, ONLY: EF_MODULE, EF_FILE
  
  IMPLICIT NONE

  !---public member functions
  PUBLIC :: start_biogenic_emissions, &   ! read namelist, define tracers and parameters
            init_emi_biogenic,        &   ! initialize emissions module (set boundary condition etc.)
            init_emi_biogenic_stream, &   ! initialize diagnostic output
            calc_biogenic_emissions       ! calculate vegetation emissions for current time step

  !---debugging flag
  LOGICAL, PUBLIC :: ldebug_bioemi = .FALSE.    ! additional diagnostic output

!>>SF
  !---public logical to keep track of whether biogenic emissions are dynamic (ie from megan scheme)
  !   or read from file
  !   this replaces the former megan switch, which was useless as the information it was carrying
  !   can be derived from the eftype value of the biogenic sector in the emi_spec file (IssueID #153)
  LOGICAL, PUBLIC :: lbioemi_dyn = .FALSE. !default: biogenic emissions are not dynamically computed 
                                           !         (no MEGAN scheme)
!<<SF

  INTEGER, PARAMETER, PRIVATE :: ncompounds = 32   ! total number of tracers in this module
                                                   ! note: original MEGAN has 150
  INTEGER, PARAMETER   :: NPFT = 16   ! number of PFTs from MEGAN 2.1 (actually from CLM4)
                                      ! PFTs are:
                                      !  - : bare
                                      !  1 : needleleaf evergreen temperate tree
                                      !  3 : needleleaf evergreen boreal tree
                                      !  2 : needleleaf deciduous boreal tree
                                      !  4 : broadleaf evergreen tropical tree
                                      !  5 : broadleaf evergreen temperate tree
                                      !  6 : broadleaf deciduous tropical tree
                                      !  7 : broadleaf deciduous temperate tree
                                      !  8 : broadleaf deciduous boreal tree
                                      !  9 : broadleaf evergreen temperate shrub
                                      ! 10 : broadleaf deciduous temperate shrub
                                      ! 11 : broadleaf deciduous boreal shrub
                                      ! 12 : arctic C3 grass
                                      ! 13 : cool C3 grass
                                      ! 14 : warm C4 grass
                                      ! 15 : crop1
                                      ! 16 : crop2 (presently always zero!)
                                      ! ATTENTION: index 2 and 3 are swapped - this is a "feature" of the MEGAN code!
                                      ! The PFT input file has bare soils as index 1 -> shift index by 1 when reading

  
  ! ### declaration of number of PFT exported from JSBACH to MEGAN following Colombe's work
  INTEGER, PARAMETER, PRIVATE :: npft_red_jsbach = 11    ! number of PFTs from JSBACH, CMIP5 setup, 11 PFTs-tiles
                                                         ! 1:Tropical evergreen
                                                         ! 2:Tropical deciduous
                                                         ! 3:Extra tropical evergreen (broadl+connif)
                                                         ! 4:Extra tropical deciduous (broadl+connif)
                                                         ! 5:Raingreen shrubs
                                                         ! 6:Deciduous shrubs
                                                         ! 7:C3 grasses
                                                         ! 8:C4 grasses
                                                         ! 9:C3 pasture
                                                         ! 10:C4 pasture
                                                         ! 11:C3+C4 crops

! ###

  !---derived data types
  !   compound properties
  TYPE, PRIVATE :: compound_properties
     CHARACTER(len=ll) :: longname              ! name of species
     CHARACTER(len=ln) :: shortname             ! chemical formula or other short name
     INTEGER           :: class                 ! MEGAN compound class (1..20)
     INTEGER           :: id_tr                 ! index of compound tracer in tracer list
     INTEGER           :: id_sp                 ! index of compound species in species list
     REAL(dp)          :: mweight               ! molecular weight
     INTEGER           :: cnumber               ! carbon number
     REAL(dp)          :: beta                  ! temperature adjustment factor (2012 paper, Table 4)
     REAL(dp)          :: ldf                   ! light-dependent fraction (2012 paper, Table 4)
     INTEGER           :: cafac                 ! relative emission rates for leaf age (2012 paper, Table 4)
     REAL(dp)          :: ef(NPFT)              ! PFT specific emission factor @ standard conditions
  END TYPE compound_properties

  !---linked list of surface fields
  TYPE, PRIVATE :: fldseries                        ! Type for 24-h records of surface field data
    REAL(dp), POINTER :: fld2d(:,:)                 ! Surface field
    CHARACTER(LEN=9)  :: fldname                    ! Name of the field in the rerun file
    TYPE(fldseries), POINTER :: nextfldseries         
  END TYPE fldseries
    
  !---module data
  REAL(dp), PARAMETER       :: zfac = 1.E-9_dp / 3600._dp   ! ug (m-2) h-1 --> kg (m-2) s-1

  !   definition of compounds
  TYPE (compound_properties), PRIVATE, ALLOCATABLE  :: compound(:)

  !   plant functional types map
  REAL(dp), ALLOCATABLE         :: pft(:,:,:)

  !   boundary conditions: leaf area index (LAI), emission factors (EF), and emission fluxes (EMI)
  INTEGER,  PRIVATE :: ibc_lai, ibc_ef(ncompounds), ibc_emi(ncompounds)
  !   boundary condition: CO2 atm. conc. (co2), relative soil water (sw)
  INTEGER, PRIVATE  :: ibc_co2
  INTEGER, PRIVATE  :: ibc_sw
  !### note: jsswpar (and possibly jsswpardif) from mo_surface_memory could also be implemented as boundary conditions


! ### declaration of array for JSBACH PFT fractions computation
  INTEGER, PRIVATE :: ibc_frpft(npft_red_jsbach), ibc_fr_brdle_ev, ibc_fr_brdle_de
  REAL(dp), ALLOCATABLE         :: vegfract_de(:,:,:,:)
  REAL(dp), ALLOCATABLE         :: vegfract_ev(:,:,:,:)
  REAL(dp), ALLOCATABLE         :: borealtr_limit(:,:,:)
  REAL(dp), ALLOCATABLE         :: borealgr_limit(:,:,:)


!>>gf #244
  INTEGER, PRIVATE :: nlai_biogenic_ef_type = EF_MODULE              !< Choice of lai external field type in the
                                                                     !< biogenic emissions module
                                                                     !<  = EF_FILE (2) from external input file
                                                                     !<  = EF_MODULE (3) online from jsbach
  INTEGER, PRIVATE :: nef_pft = EF_FILE                              !< choice of PFTs for calculation of emission factor
                                                                     !< = EF_FILE PFT fractions from MEGAN-CLM4
                                                                     !< = EF_MODULE PFT fractions from JSBACH
  CHARACTER(LEN=ln), SAVE :: emifact_files_species(ncompounds)       !< Choice of emission factors
                                                                     !< if the short name of the compound is given in
                                                                     !< the biogenic_emission namelist, then 
                                                                     !< emission factor is read from specific file,
                                                                     !< otherwise it is read from the list below
!<<gf

  !---biogenic emissions streams
  TYPE (t_stream), PUBLIC, POINTER :: bioemiav

  !   instantaneous and 24h-average temperature and surface SW radiation
  REAL(dp), POINTER :: avg_ppfd(:,:)         ! ppfd = photosynthetic photon flux density (=PAR)
  REAL(dp), POINTER :: avg_tm1(:,:)
  REAL(dp), POINTER :: avg_lai(:,:)
  REAL(dp), POINTER :: prt_lai(:,:)
  REAL(dp), POINTER :: LAI_daymean(:,:)
  REAL(dp), POINTER :: LAI_monthmean(:,:)
  REAL(dp), POINTER :: LAI_currentmonth(:,:)
  REAL(dp), POINTER :: LAI_sum(:,:)
  REAL(dp), POINTER :: LAI_sum_month(:,:)
 
  !   fields for the calculation of 24-h averages of temperature and surface SW radiation
  TYPE(fldseries), PRIVATE, ALLOCATABLE, TARGET  :: T24(:), ppfd24(:),lai24(:)
  TYPE(fldseries), PRIVATE, POINTER :: oldest_T24, newest_T24, oldest_ppfd24, newest_ppfd24
  TYPE(fldseries), PRIVATE, POINTER :: oldest_lai24, newest_lai24

  !   diagnostic pointers for efs and gamma functions
  TYPE(vmem2d) :: ef_diag(ncompounds)
  TYPE(vmem2d) :: gamma_diag(23)        ! gamma_lai, gamma_ppdf, gamma_temp_iso, gamma_temp_other,
                                        ! gamma_age, and combined gamma factors for:
                                        ! isoprene, a-pinene, methanol, acetone, ethene
  !   others
  INTEGER, PRIVATE :: nstepsperday                   ! Number of timesteps in ppfd24 and T24 arrays
  INTEGER, PRIVATE :: iavgstep                       ! Current step number for averaging
                                                     ! 1...nstepsperday
  !   for calculation of 24-h averages (see subroutine biogenic_averages)
  INTEGER,  PRIVATE, SAVE :: itimestep = -1


CONTAINS
  !-----------------------------------start_biogenic_emissions--------------------------------------------------
  SUBROUTINE start_biogenic_emissions

    !---inherited types, data and functions
    USE mo_mpi,                      ONLY: p_parallel_io, p_parallel, p_io, p_bcast
    USE mo_exception,                ONLY: finish, message, message_text, em_error, em_warn, em_info
    USE mo_submodel,                 ONLY: print_value, lbioemi_stdalone
    USE mo_tracdef,                  ONLY: trlist, CONSTANT, RESTART, GAS
    USE mo_tracer,                   ONLY: get_tracer, new_tracer
    USE mo_species,                  ONLY: new_species, query_species
    USE mo_namelist,                 ONLY: open_nml,close_nml, position_nml, POSITIONED, READ_ERROR
    USE mo_external_field_processor, ONLY: EF_MODULE, EF_FILE

    IMPLICIT NONE

    INCLUDE 'biogenic_emissionsctl.inc' !gf #244

    !---local parameters
    REAL(dp), PARAMETER :: nsphr  = 3600._dp     ! no of sec/hour
    REAL(dp), PARAMETER :: awC  = 12.011_dp      ! atomic weight of carbon

    !---local data
    INTEGER :: i, ierr, ispid, inml, iunit
    REAL(dp) :: rm
 
    !---executable procedure

! set default value of compound flag
    emifact_files_species(:) =''
    !>>csld #243 & #257 : Martin wants the reading of PFT averaged over grid box from the file 
    !                     ef_specifics.nc, when available
    emifact_files_species(1) = 'APIN'
    emifact_files_species(2) = 'BPIN'
    emifact_files_species(3) = 'CARENE3'
    emifact_files_species(4) = 'C5H8'
    emifact_files_species(5) = 'LIMON'
    emifact_files_species(6) = 'MYRC'
    emifact_files_species(7) = 'TBETAOCI'
    emifact_files_species(8) = 'SABIN'
    emifact_files_species(9) = 'MBO'
    emifact_files_species(10) = 'NO'
    !<<csld 

!>>gf #244
    !---1.1 Read biogenic_emissionsctl namelist
    IF (p_parallel_io) THEN
       inml = open_nml ('namelist.echam')
       iunit = position_nml ('biogenic_emissionsctl', inml, status=ierr)
       SELECT CASE (ierr)
          CASE (POSITIONED)
             READ(iunit,biogenic_emissionsctl)
          CASE (READ_ERROR)
             CALL finish ('start_biogenic_emissions','general read error in namelist.echam')
       END SELECT
       CALL close_nml(inml)  
    END IF

    IF (p_parallel) THEN
       CALL p_bcast (nlai_biogenic_ef_type,p_io)
       CALL p_bcast (nef_pft,p_io)
       CALL p_bcast (emifact_files_species,p_io)
       CALL p_bcast (ldebug_bioemi,p_io)
    END IF

    CALL message('start_biogenic_emissions','----------------------------------------------------------')
    CALL message('','--- Initialization of the biogenic emissions module')
    CALL message('', '---')
    CALL print_value('lbioemi_stdalone',lbioemi_stdalone)
    CALL message('', '---')
    CALL message('start_biogenic_emissions','----------------------------------------------------------')

    !---1.2 consistency check for nlai_biogenic_ef_type !SF + user information
    SELECT CASE(nlai_biogenic_ef_type)
        CASE(EF_FILE)
           WRITE(message_text,'(a)') 'Leaf area index read from file'
        CASE(EF_MODULE)
           WRITE(message_text,'(a)') 'Leaf area index from JSBACH'
        CASE default
         WRITE(message_text,'(a,i0,2a,2(i0,a))') 'nlai_biogenic_ef_type = ',nlai_biogenic_ef_type, &
                                ' --> this is not currently supported.', &
                                ' Only ',EF_FILE,' (from file) or ',EF_MODULE, &
                                ' (from another module) external field types are possible!'
         CALL finish('start_biogenic_emissions',message_text)
    END SELECT

    WRITE(message_text,'(a,a)') TRIM(message_text),' for usage in the biogenic emissions module.'
    CALL message('start_biogenic_emissions',message_text,level=em_info)
!<<gf
    SELECT CASE(nef_pft)
        CASE(EF_FILE)
           WRITE(message_text,'(a)') 'PFT fractions from MEGAN'
        CASE(EF_MODULE)
           WRITE(message_text,'(a)') 'PFT fractions from JSBACH'
        CASE default
         WRITE(message_text,'(a,i0,a)') 'nef_pft = ',nef_pft, &
                                ' --> this is not currently supported.'
         CALL finish('start_biogenic_emissions',message_text)
    END SELECT

    WRITE(message_text,'(a,a)') TRIM(message_text),' for usage in the biogenic emissions module.'
    CALL message('start_biogenic_emissions',message_text,level=em_info)

    ALLOCATE(compound(ncompounds))
    compound(:)%id_sp = 0
    compound(:)%id_tr = 0

    !---names of compound tracers
                                                    ! total emissions according to Guenther et al., 2012 (Tg/yr)
    compound( 1)%longname = 'isoprene'              ! 535.
    compound( 2)%longname = 'alpha pinene'          !  66.1
    compound( 3)%longname = 't beta ocimene'        !  19.4
    compound( 4)%longname = 'b pinene'              !  18.9
    compound( 5)%longname = 'limonene'              !  11.4
    compound( 6)%longname = 'sabinene'              !   9.0
    compound( 7)%longname = 'myrcene'               !   8.7
    compound( 8)%longname = '3-carene'              !   7.1
                                                    !  another 35 monoterpenes add another 217. Tg/yr
                                                    !  sesquiterpenes are ignored for now
    compound( 9)%longname = '2-methyl-3-buten-2-ol' !   2.2   (MBO)
    compound(10)%longname = 'methanol'              !  99.6
    compound(11)%longname = 'acetone'               !  43.7
    compound(12)%longname = 'ethanol'               !  20.7
    compound(13)%longname = 'acetaldehyde'          !  20.7
    compound(14)%longname = 'formaldehyde'          !   5.0
    compound(15)%longname = 'acetic acid'           !   3.7
    compound(16)%longname = 'formic acid'           !   3.7
    compound(17)%longname = 'ethene'                !  26.9
    compound(18)%longname = 'hydrogen cyanide'      !
    compound(19)%longname = 'toluene'               !
    compound(20)%longname = 'methyl bromide'        !
    compound(21)%longname = 'methyl chloride'       !
    compound(22)%longname = 'methyl iodide'         !
    compound(23)%longname = 'dimethyl sulfide'      !
    compound(24)%longname = 'methane'               !
    compound(25)%longname = 'ethane'                !
    compound(26)%longname = 'propane'               !
    compound(27)%longname = 'propene'               !  15.8
    compound(28)%longname = 'butene'                !   8.0
    compound(29)%longname = 'benzaldehyde'          !
                                                    ! plus others (in particular oxygenated VOC)
    compound(30)%longname = 'carbon monoxide'       !  81.6
    compound(31)%longname = 'nitricoxide'           !
    compound(32)%longname = 'beta caryophyllene'    ! 7.4 sesquiterpene

    compound( 1)%shortname = 'C5H8'               ! isoprene
    compound( 2)%shortname = 'APIN'               ! alpha pinene
    compound( 3)%shortname = 'TBETAOCI'           ! t beta ocimene    (add to a pinene!)
    compound( 4)%shortname = 'BPIN'               ! b pinene          (add to a pinene!)
    compound( 5)%shortname = 'LIMON'              ! limonene          (add to a pinene!)
    compound( 6)%shortname = 'SABIN'              ! sabinene          (add to a pinene!)
    compound( 7)%shortname = 'MYRC'               ! myrcene           (add to a pinene!)
    compound( 8)%shortname = 'CARENE3'            ! 3-carene          (add to a pinene!)
    compound( 9)%shortname = 'MBO'                ! MBO
    compound(10)%shortname = 'CH3OH'              ! methanol
    compound(11)%shortname = 'CH3COCH3'           ! acetone
    compound(12)%shortname = 'C2H5OH'             ! ethanol
    compound(13)%shortname = 'CH3CHO'             ! acetaldehyde
    compound(14)%shortname = 'CH2O'               ! formaldehyde
    compound(15)%shortname = 'CH3COOH'            ! acetic acid
    compound(16)%shortname = 'HCOOH'              ! formic acid
    compound(17)%shortname = 'C2H4'               ! ethene
    compound(18)%shortname = 'HCN'                ! hydrogen cyanide
    compound(19)%shortname = 'TOL'                ! toluene
    compound(20)%shortname = 'CH3BR'              ! methyl bromide
    compound(21)%shortname = 'CH3CL'              ! methyl chloride
    compound(22)%shortname = 'CH3I'               ! methyl iodide
    compound(23)%shortname = 'DMS'                ! dimethyl sulfide
    compound(24)%shortname = 'CH4'                ! methane
    compound(25)%shortname = 'C2H6'               ! ethane
    compound(26)%shortname = 'C3H8'               ! propane
    compound(27)%shortname = 'C3H6'               ! propene
    compound(28)%shortname = 'BIGENE'             ! butene (original name: C4H8)
    compound(29)%shortname = 'BZALD'              ! benzaldehyde
    compound(30)%shortname = 'CO'                 ! carbon monoxide
    compound(31)%shortname = 'NO'                 ! nitric oxide
    compound(32)%shortname = 'BCARY'             ! beta caryophyllene
 
    ! -- compound class
    !    from table MGN2MECH/INCLDIR/SPC_NOCONVER.EXT
    compound( 1)%class =  1                       ! isoprene
    compound( 2)%class =  8                       ! alpha pinene
    compound( 3)%class =  6                       ! t beta ocimene
    compound( 4)%class =  7                       ! b pinene
    compound( 5)%class =  4                       ! limonene
    compound( 6)%class =  3                       ! sabinene
    compound( 7)%class =  2                       ! myrcene
    compound( 8)%class =  5                       ! 3-carene
    compound( 9)%class = 13                       ! MBO
    compound(10)%class = 14                       ! methanol
    compound(11)%class = 15                       ! acetone
    compound(12)%class = 18                       ! ethanol
    compound(13)%class = 18                       ! acetaldehyde
    compound(14)%class = 18                       ! formaldehyde
    compound(15)%class = 18                       ! acetic acid
    compound(16)%class = 18                       ! formic acid
    compound(17)%class = 19                       ! ethene
    compound(18)%class = 19                       ! hydrogen cyanide
    compound(19)%class = 19                       ! toluene
    compound(20)%class = 20                       ! methyl bromide
    compound(21)%class = 20                       ! methyl chloride
    compound(22)%class = 20                       ! methyl iodide
    compound(23)%class = 20                       ! dimethyl sulfide
    compound(24)%class = 20                       ! methane
    compound(25)%class = 20                       ! ethane
    compound(26)%class = 20                       ! propane
    compound(27)%class = 20                       ! propene
    compound(28)%class = 20                       ! butene
    compound(29)%class = 20                       ! benzaldehyde
    compound(30)%class = 16                       ! carbon monoxide
    compound(31)%class = 17                       ! nitrogen monoxide
    compound(32)%class = 11                       ! beta caryophyllene
    !---molecular weights
    !    from table MGN2MECH/INCLDIR/SPC_NOCONVER.EXT
    compound( 1)%mweight =  68.12_dp              ! isoprene
    compound( 2)%mweight = 136.23_dp              ! alpha pinene
    compound( 3)%mweight = 136.23_dp              ! t beta ocimene
    compound( 4)%mweight = 136.23_dp              ! b pinene
    compound( 5)%mweight = 136.23_dp              ! limonene
    compound( 6)%mweight = 136.23_dp              ! sabinene
    compound( 7)%mweight = 136.23_dp              ! myrcene
    compound( 8)%mweight = 136.23_dp              ! 3-carene
    compound( 9)%mweight =  86.13_dp              ! MBO
    compound(10)%mweight =  32.04_dp              ! methanol
    compound(11)%mweight =  58.08_dp              ! acetone
    compound(12)%mweight =  46.07_dp              ! ethanol
    compound(13)%mweight =  44.05_dp              ! acetaldehyde
    compound(14)%mweight =  30.03_dp              ! formaldehyde
    compound(15)%mweight =  60.05_dp              ! acetic acid
    compound(16)%mweight =  46.03_dp              ! formic acid
    compound(17)%mweight =  28.05_dp              ! ethene
    compound(18)%mweight =  27.03_dp              ! hydrogen cyanide
    compound(19)%mweight =  92.14_dp              ! toluene
    compound(20)%mweight =  94.94_dp              ! methyl bromide
    compound(21)%mweight =  50.49_dp              ! methyl chloride
    compound(22)%mweight = 141.94_dp              ! methyl iodide
    compound(23)%mweight =  62.13_dp              ! dimethyl sulfide   (A_2met_s in ori MEGAN)
    compound(24)%mweight =  16.04_dp              ! methane
    compound(25)%mweight =  30.07_dp              ! ethane
    compound(26)%mweight =  44.10_dp              ! propane
    compound(27)%mweight =  42.08_dp              ! propene
    compound(28)%mweight =  56.11_dp              ! butene
    compound(29)%mweight = 106.12_dp              ! benzaldehyde
    compound(30)%mweight =  28.01_dp              ! carbon monoxide
    compound(31)%mweight =  30.01_dp              ! nitric oxide
    compound(32)%mweight =  204.35_dp             ! beta caryophyllene

    !---carbon number
    ! 0 is used for compounds where emissions should not be converted in diagnostics!
    compound( 1)%cnumber =  5                     ! isoprene
    compound( 2)%cnumber = 10                     ! alpha pinene
    compound( 3)%cnumber = 10                     ! t beta ocimene
    compound( 4)%cnumber = 10                     ! b pinene
    compound( 5)%cnumber = 10                     ! limonene
    compound( 6)%cnumber = 10                     ! sabinene
    compound( 7)%cnumber = 10                     ! myrcene
    compound( 8)%cnumber = 10                     ! 3-carene
    compound( 9)%cnumber =  4                     ! MBO
    compound(10)%cnumber =  1                     ! methanol
    compound(11)%cnumber =  3                     ! acetone
    compound(12)%cnumber =  2                     ! ethanol
    compound(13)%cnumber =  2                     ! acetaldehyde
    compound(14)%cnumber =  1                     ! formaldehyde
    compound(15)%cnumber =  2                     ! acetic acid
    compound(16)%cnumber =  1                     ! formic acid
    compound(17)%cnumber =  2                     ! ethene
    compound(18)%cnumber =  0                     ! hydrogen cyanide
    compound(19)%cnumber =  7                     ! toluene
    compound(20)%cnumber =  1                     ! methyl bromide
    compound(21)%cnumber =  1                     ! methyl chloride
    compound(22)%cnumber =  1                     ! methyl iodide
    compound(23)%cnumber =  3                     ! dimethyl sulfide
    compound(24)%cnumber =  1                     ! methane
    compound(25)%cnumber =  2                     ! ethane
    compound(26)%cnumber =  3                     ! propane
    compound(27)%cnumber =  3                     ! propene
    compound(28)%cnumber =  4                     ! butene
    compound(29)%cnumber =  6                     ! benzaldehyde
    compound(30)%cnumber =  0                     ! carbon monoxide
    compound(31)%cnumber =  0                     ! nitric oxide
    compound(32)%cnumber =  15                    ! beta caryophyllene

    !---classes for relative emission rates used in
    !   calculation of leaf age emission factor
    compound( 1)%cafac = 2                     ! isoprene
    compound( 2)%cafac = 3                     ! alpha pinene
    compound( 3)%cafac = 3                     ! t beta ocimene
    compound( 4)%cafac = 3                     ! b pinene
    compound( 5)%cafac = 3                     ! limonene
    compound( 6)%cafac = 3                     ! sabinene
    compound( 7)%cafac = 3                     ! myrcene
    compound( 8)%cafac = 3                     ! 3-carene
    compound( 9)%cafac = 5                     ! MBO
    compound(10)%cafac = 6                     ! methanol
    compound(11)%cafac = 1                     ! acetone
    compound(12)%cafac = 1                     ! ethanol
    compound(13)%cafac = 1                     ! acetaldehyde
    compound(14)%cafac = 1                     ! formaldehyde
    compound(15)%cafac = 1                     ! acetic acid
    compound(16)%cafac = 1                     ! formic acid
    compound(17)%cafac = 1                     ! ethene
    compound(18)%cafac = 1                     ! hydrogen cyanide
    compound(19)%cafac = 1                     ! toluene
    compound(20)%cafac = 1                     ! methyl bromide
    compound(21)%cafac = 1                     ! methyl chloride
    compound(22)%cafac = 1                     ! methyl iodide
    compound(23)%cafac = 1                     ! dimethyl sulfide
    compound(24)%cafac = 1                     ! methane
    compound(25)%cafac = 1                     ! ethane
    compound(26)%cafac = 1                     ! propane
    compound(27)%cafac = 1                     ! propene
    compound(28)%cafac = 1                     ! butene
    compound(29)%cafac = 1                     ! benzaldehyde
    compound(30)%cafac = 1                     ! carbon monoxide
    compound(31)%cafac = 1                     ! nitric oxide
    compound(32)%cafac = 4                     ! beta caryophyllene

    !---temperature adjustment factor
    compound( 1)%beta = 0.13_dp                ! isoprene
    compound( 2)%beta = 0.10_dp                ! alpha pinene
    compound( 3)%beta = 0.10_dp                ! t beta ocimene
    compound( 4)%beta = 0.10_dp                ! b pinene
    compound( 5)%beta = 0.10_dp                ! limonene
    compound( 6)%beta = 0.10_dp                ! sabinene
    compound( 7)%beta = 0.10_dp                ! myrcene
    compound( 8)%beta = 0.10_dp                ! 3-carene
    compound( 9)%beta = 0.13_dp                ! MBO
    compound(10)%beta = 0.08_dp                ! methanol
    compound(11)%beta = 0.10_dp                ! acetone
    compound(12)%beta = 0.13_dp                ! ethanol
    compound(13)%beta = 0.13_dp                ! acetaldehyde
    compound(14)%beta = 0.13_dp                ! formaldehyde
    compound(15)%beta = 0.13_dp                ! acetic acid
    compound(16)%beta = 0.13_dp                ! formic acid
    compound(17)%beta = 0.10_dp                ! ethene
    compound(18)%beta = 0.10_dp                ! hydrogen cyanide
    compound(19)%beta = 0.10_dp                ! toluene
    compound(20)%beta = 0.10_dp                ! methyl bromide
    compound(21)%beta = 0.10_dp                ! methyl chloride
    compound(22)%beta = 0.10_dp                ! methyl iodide
    compound(23)%beta = 0.10_dp                ! dimethyl sulfide
    compound(24)%beta = 0.10_dp                ! methane
    compound(25)%beta = 0.10_dp                ! ethane
    compound(26)%beta = 0.10_dp                ! propane
    compound(27)%beta = 0.10_dp                ! propene
    compound(28)%beta = 0.10_dp                ! butene
    compound(29)%beta = 0.10_dp                ! benzaldehyde
    compound(30)%beta = 0.08_dp                ! carbon monoxide
    compound(31)%beta = 0.00_dp                ! nitric oxide      ### CHECK !! ###
    compound(32)%beta = 0.17_dp                ! beta caryophyllene

    !---fraction of light dependence for correction factors
    compound( 1)%ldf = 1.00_dp                ! isoprene
    compound( 2)%ldf = 0.60_dp                ! alpha pinene
    compound( 3)%ldf = 0.80_dp                ! t beta ocimene
    compound( 4)%ldf = 0.20_dp                ! b pinene
    compound( 5)%ldf = 0.20_dp                ! limonene
    compound( 6)%ldf = 0.60_dp                ! sabinene
    compound( 7)%ldf = 0.60_dp                ! myrcene
    compound( 8)%ldf = 0.20_dp                ! 3-carene
    compound( 9)%ldf = 1.00_dp                ! MBO
    compound(10)%ldf = 0.80_dp                ! methanol
    compound(11)%ldf = 0.20_dp                ! acetone
    compound(12)%ldf = 0.80_dp                ! ethanol
    compound(13)%ldf = 0.80_dp                ! acetaldehyde
    compound(14)%ldf = 0.80_dp                ! formaldehyde
    compound(15)%ldf = 0.80_dp                ! acetic acid
    compound(16)%ldf = 0.80_dp                ! formic acid
    compound(17)%ldf = 0.80_dp                ! ethene
    compound(18)%ldf = 0.80_dp                ! hydrogen cyanide
    compound(19)%ldf = 0.80_dp                ! toluene
    compound(20)%ldf = 0.20_dp                ! methyl bromide
    compound(21)%ldf = 0.20_dp                ! methyl chloride
    compound(22)%ldf = 0.20_dp                ! methyl iodide
    compound(23)%ldf = 0.20_dp                ! dimethyl sulfide
    compound(24)%ldf = 0.20_dp                ! methane
    compound(25)%ldf = 0.20_dp                ! ethane
    compound(26)%ldf = 0.20_dp                ! propane
    compound(27)%ldf = 0.20_dp                ! propene
    compound(28)%ldf = 0.20_dp                ! butene
    compound(29)%ldf = 0.20_dp                ! benzaldehyde
    compound(30)%ldf = 1.00_dp                ! carbon monoxide
    compound(31)%ldf = 0.00_dp                ! nitric oxide      ### CHECK !! ###
    compound(32)%ldf = 0.50_dp                ! beta caryophyllene

    !--- MEGAN specific tracer emission factor ----
    ! -- PFT specific emission factors
    !    from table MGN2MECH/INCLDIR/EFS_PFT.EXT.womap
    !    use python script ef_hammoz.py to generate the following input data
       compound( 1)%ef(:) =  &                       ! isoprene (class: 1)
   (/     600.0000_dp,       1.0000_dp,    3000.0000_dp,    7000.0000_dp, &
        10000.0000_dp,    7000.0000_dp,   10000.0000_dp,   11000.0000_dp, &
         2000.0000_dp,    4000.0000_dp,    4000.0000_dp,    1600.0000_dp, &
          800.0000_dp,     200.0000_dp,      50.0000_dp,       1.0000_dp /)
       compound( 2)%ef(:) =  &                       ! pinene_a (class: 8)
   (/     500.0000_dp,     510.0000_dp,     500.0000_dp,     600.0000_dp, &
          400.0000_dp,     600.0000_dp,     400.0000_dp,     400.0000_dp, &
          200.0000_dp,     300.0000_dp,     200.0000_dp,       2.0000_dp, &
            2.0000_dp,       2.0000_dp,       2.0000_dp,       2.0000_dp /)
       compound( 3)%ef(:) =  &                       ! ocimene_t_b (class: 6)
   (/      70.0000_dp,      60.0000_dp,      70.0000_dp,     150.0000_dp, &
          120.0000_dp,     150.0000_dp,     120.0000_dp,     120.0000_dp, &
           90.0000_dp,     150.0000_dp,      90.0000_dp,       2.0000_dp, &
            2.0000_dp,       2.0000_dp,       2.0000_dp,       2.0000_dp /)
       compound( 4)%ef(:) =  &                       ! pinene_b (class: 7)
   (/     300.0000_dp,     200.0000_dp,     300.0000_dp,     120.0000_dp, &
          130.0000_dp,     120.0000_dp,     130.0000_dp,     130.0000_dp, &
          100.0000_dp,     150.0000_dp,     100.0000_dp,       1.5000_dp, &
            1.5000_dp,       1.5000_dp,       1.5000_dp,       1.5000_dp /)
       compound( 5)%ef(:) =  &                       ! limonene (class: 4)
   (/     100.0000_dp,     130.0000_dp,     100.0000_dp,      80.0000_dp, &
           80.0000_dp,      80.0000_dp,      80.0000_dp,      80.0000_dp, &
           60.0000_dp,     100.0000_dp,      60.0000_dp,       0.7000_dp, &
            0.7000_dp,       0.7000_dp,       0.7000_dp,       0.7000_dp /)
       compound( 6)%ef(:) =  &                       ! sabinene (class: 3)
   (/      70.0000_dp,      40.0000_dp,      70.0000_dp,      80.0000_dp, &
           50.0000_dp,      80.0000_dp,      50.0000_dp,      50.0000_dp, &
           50.0000_dp,      70.0000_dp,      50.0000_dp,       0.7000_dp, &
            0.7000_dp,       0.7000_dp,       0.7000_dp,       0.7000_dp /)
       compound( 7)%ef(:) =  &                       ! myrcene (class: 2)
   (/      70.0000_dp,      60.0000_dp,      70.0000_dp,      80.0000_dp, &
           30.0000_dp,      80.0000_dp,      30.0000_dp,      30.0000_dp, &
           30.0000_dp,      50.0000_dp,      30.0000_dp,       0.3000_dp, &
            0.3000_dp,       0.3000_dp,       0.3000_dp,       0.3000_dp /)
       compound( 8)%ef(:) =  &                       ! carene_3 (class: 5)
   (/     160.0000_dp,      80.0000_dp,     160.0000_dp,      40.0000_dp, &
           30.0000_dp,      40.0000_dp,      30.0000_dp,      30.0000_dp, &
           30.0000_dp,     100.0000_dp,      30.0000_dp,       0.3000_dp, &
            0.3000_dp,       0.3000_dp,       0.3000_dp,       0.3000_dp /)
       compound( 9)%ef(:) =  &                       ! MBO_3m2e1ol (class: 13)
   (/     700.0000_dp,       0.0100_dp,      60.0000_dp,       0.0100_dp, &
            0.0100_dp,       0.0100_dp,       0.0100_dp,       2.0000_dp, &
            0.0100_dp,       0.0100_dp,       0.0100_dp,       0.0100_dp, &
            0.0100_dp,       0.0100_dp,       0.0100_dp,       0.0100_dp /)
       compound(10)%ef(:) =  &                       ! methanol (class: 14)
   (/     900.0000_dp,     900.0000_dp,     900.0000_dp,     500.0000_dp, &
          900.0000_dp,     500.0000_dp,     900.0000_dp,     900.0000_dp, &
          900.0000_dp,     900.0000_dp,     900.0000_dp,     500.0000_dp, &
          500.0000_dp,     500.0000_dp,     900.0000_dp,     900.0000_dp /)
       compound(11)%ef(:) =  &                       ! acetone (class: 15)
   (/     240.0000_dp,     240.0000_dp,     240.0000_dp,     240.0000_dp, &
          240.0000_dp,     240.0000_dp,     240.0000_dp,     240.0000_dp, &
          240.0000_dp,     240.0000_dp,     240.0000_dp,      80.0000_dp, &
           80.0000_dp,      80.0000_dp,      80.0000_dp,      80.0000_dp /)
       compound(12)%ef(:) =  &                       ! ethanol (class: 18)
   (/     200.0000_dp,     200.0000_dp,     200.0000_dp,     200.0000_dp, &
          200.0000_dp,     200.0000_dp,     200.0000_dp,     200.0000_dp, &
          200.0000_dp,     200.0000_dp,     200.0000_dp,      20.0000_dp, &
           20.0000_dp,      20.0000_dp,      20.0000_dp,      20.0000_dp /)
       compound(13)%ef(:) =  &                       ! acetaldehyde (class: 18)
   (/     200.0000_dp,     200.0000_dp,     200.0000_dp,     200.0000_dp, &
          200.0000_dp,     200.0000_dp,     200.0000_dp,     200.0000_dp, &
          200.0000_dp,     200.0000_dp,     200.0000_dp,      20.0000_dp, &
           20.0000_dp,      20.0000_dp,      20.0000_dp,      20.0000_dp /)
       compound(14)%ef(:) =  &                       ! formaldehyde (class: 18)
   (/      40.0000_dp,      40.0000_dp,      40.0000_dp,      40.0000_dp, &
           40.0000_dp,      40.0000_dp,      40.0000_dp,      40.0000_dp, &
           40.0000_dp,      40.0000_dp,      40.0000_dp,      16.0000_dp, &
           16.0000_dp,      16.0000_dp,      16.0000_dp,      16.0000_dp /)
       compound(15)%ef(:) =  &                       ! acetic_acid (class: 18)
   (/      30.0000_dp,      30.0000_dp,      30.0000_dp,      30.0000_dp, &
           30.0000_dp,      30.0000_dp,      30.0000_dp,      30.0000_dp, &
           30.0000_dp,      30.0000_dp,      30.0000_dp,      12.0000_dp, &
           12.0000_dp,      12.0000_dp,      12.0000_dp,      12.0000_dp /)
       compound(16)%ef(:) =  &                       ! formic_acid (class: 18)
   (/      30.0000_dp,      30.0000_dp,      30.0000_dp,      30.0000_dp, &
           30.0000_dp,      30.0000_dp,      30.0000_dp,      30.0000_dp, &
           30.0000_dp,      30.0000_dp,      30.0000_dp,      12.0000_dp, &
           12.0000_dp,      12.0000_dp,      12.0000_dp,      12.0000_dp /)
       compound(17)%ef(:) =  &                       ! ethene (class: 19)
   (/     174.0000_dp,     174.0000_dp,     174.0000_dp,     174.0000_dp, &
          174.0000_dp,     174.0000_dp,     174.0000_dp,     174.0000_dp, &
          174.0000_dp,     174.0000_dp,     174.0000_dp,     174.0000_dp, &
          174.0000_dp,     174.0000_dp,     174.0000_dp,     174.0000_dp /)
       compound(18)%ef(:) =  &                       ! hydrogen_cyanide (class: 19)
   (/       4.5000_dp,       4.5000_dp,       4.5000_dp,       4.5000_dp, &
            4.5000_dp,       4.5000_dp,       4.5000_dp,       4.5000_dp, &
            4.5000_dp,       4.5000_dp,       4.5000_dp,       4.5000_dp, &
            4.5000_dp,       4.5000_dp,       4.5000_dp,       4.5000_dp /)
       compound(19)%ef(:) =  &                       ! toluene (class: 19)
   (/       9.0000_dp,       9.0000_dp,       9.0000_dp,       9.0000_dp, &
            9.0000_dp,       9.0000_dp,       9.0000_dp,       9.0000_dp, &
            9.0000_dp,       9.0000_dp,       9.0000_dp,       9.0000_dp, &
            9.0000_dp,       9.0000_dp,       9.0000_dp,       9.0000_dp /)
       compound(20)%ef(:) =  &                       ! met_bromide (class: 20)
   (/       0.2800_dp,       0.2800_dp,       0.2800_dp,       0.2800_dp, &
            0.2800_dp,       0.2800_dp,       0.2800_dp,       0.2800_dp, &
            0.2800_dp,       0.2800_dp,       0.2800_dp,       0.2800_dp, &
            0.2800_dp,       0.2800_dp,       0.2800_dp,       0.2800_dp /)
       compound(21)%ef(:) =  &                       ! met_chloride (class: 20)
   (/       1.4000_dp,       1.4000_dp,       1.4000_dp,       1.4000_dp, &
            1.4000_dp,       1.4000_dp,       1.4000_dp,       1.4000_dp, &
            1.4000_dp,       1.4000_dp,       1.4000_dp,       1.4000_dp, &
            1.4000_dp,       1.4000_dp,       1.4000_dp,       1.4000_dp /)
       compound(22)%ef(:) =  &                       ! met_iodide (class: 20)
   (/       0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp, &
            0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp, &
            0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp, &
            0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp /)
       compound(23)%ef(:) =  &                       ! carbonyl_s (class: 20)
   (/       0.4200_dp,       0.4200_dp,       0.4200_dp,       0.4200_dp, &
            0.4200_dp,       0.4200_dp,       0.4200_dp,       0.4200_dp, &
            0.4200_dp,       0.4200_dp,       0.4200_dp,       0.4200_dp, &
            0.4200_dp,       0.4200_dp,       0.4200_dp,       0.4200_dp /)
       compound(24)%ef(:) =  &                       ! methane (class: 20)
   (/       0.7000_dp,       0.7000_dp,       0.7000_dp,       0.7000_dp, &
            0.7000_dp,       0.7000_dp,       0.7000_dp,       0.7000_dp, &
            0.7000_dp,       0.7000_dp,       0.7000_dp,       0.7000_dp, &
            0.7000_dp,       0.7000_dp,       0.7000_dp,       0.7000_dp /)
       compound(25)%ef(:) =  &                       ! ethane (class: 20)
   (/       1.4000_dp,       1.4000_dp,       1.4000_dp,       1.4000_dp, &
            1.4000_dp,       1.4000_dp,       1.4000_dp,       1.4000_dp, &
            1.4000_dp,       1.4000_dp,       1.4000_dp,       1.4000_dp, &
            1.4000_dp,       1.4000_dp,       1.4000_dp,       1.4000_dp /)
       compound(26)%ef(:) =  &                       ! propane (class: 20)
   (/       0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp, &
            0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp, &
            0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp, &
            0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp /)
       compound(27)%ef(:) =  &                       ! propene (class: 20)
   (/      67.2000_dp,      67.2000_dp,      67.2000_dp,      67.2000_dp, &
           67.2000_dp,      67.2000_dp,      67.2000_dp,      67.2000_dp, &
           67.2000_dp,      67.2000_dp,      67.2000_dp,      67.2000_dp, &
           67.2000_dp,      67.2000_dp,      67.2000_dp,      67.2000_dp /)
       compound(28)%ef(:) =  &                       ! butene (class: 20)
   (/      33.6000_dp,      33.6000_dp,      33.6000_dp,      33.6000_dp, &
           33.6000_dp,      33.6000_dp,      33.6000_dp,      33.6000_dp, &
           33.6000_dp,      33.6000_dp,      33.6000_dp,      33.6000_dp, &
           33.6000_dp,      33.6000_dp,      33.6000_dp,      33.6000_dp /)
       compound(29)%ef(:) =  &                       ! benzaldehyde (class: 20)
   (/       0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp, &
            0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp, &
            0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp, &
            0.1400_dp,       0.1400_dp,       0.1400_dp,       0.1400_dp /)
       compound(30)%ef(:) =  &                       ! carbon_monoxide (class: 16)
   (/     600.0000_dp,     600.0000_dp,     600.0000_dp,     600.0000_dp, &
          600.0000_dp,     600.0000_dp,     600.0000_dp,     600.0000_dp, &
          600.0000_dp,     600.0000_dp,     600.0000_dp,     600.0000_dp, &
          600.0000_dp,     600.0000_dp,     600.0000_dp,     600.0000_dp /)
       compound(31)%ef(:) =  &                       ! nitric_OXD (class: 17)
   (/       2.0000_dp,       2.0000_dp,       2.0000_dp,       2.0000_dp, &
            2.0000_dp,       2.0000_dp,       2.0000_dp,       2.0000_dp, &
            2.0000_dp,       2.0000_dp,       2.0000_dp,      27.0000_dp, &
           27.0000_dp,      27.0000_dp,      40.0000_dp,      68.0000_dp /)
       compound(32)%ef(:) =  &                       ! beta_caryophyllene (class: 11)
   (/      80.0000_dp,      80.0000_dp,      80.0000_dp,      60.0000_dp, &
           40.0000_dp,      60.0000_dp,      40.0000_dp,      40.0000_dp, &
           50.0000_dp,      50.0000_dp,      50.0000_dp,       1.0000_dp, &
            1.0000_dp,       1.0000_dp,       2.0000_dp,       4.0000_dp /)


! ### old code - unit conversion factors were species-dependant ---
!!     !---conversion factor = 1.E-9*molecular weight/3600
!!     compound(1)%conv_ft = 1._dp
!!     compound(2)%conv_ft = 1._dp
!!     DO i=nocompoundstart,ncompounds
!!        compound(i)%conv_ft = 1.E-9_dp*(compound(i)%mweight/(awC*REAL(compound(i)%cnumber,dp)))/nsphr
!!     END DO
! ### end old code

       !---3. Set tracer and species ids for MEGAN compounds

    IF (lbioemi_stdalone) THEN

       !---define species and tracers in case of standalone megan
       DO i=1,ncompounds
          CALL query_species(shortname=TRIM(compound(i)%shortname), ierr=ierr)
          IF (ierr == 0) THEN
             ! hier auch mit _a, _b, _c, ...
             CALL message ('start_biogenic_emissions', 'species '//TRIM(compound(i)%shortname)//   &
                           ' already exists', level=em_warn )
             compound(i)%shortname = "MG_"//TRIM( compound(i)%shortname )
             CALL message ('', 'species renamed to '//TRIM(compound(i)%shortname), level=em_info )
          END IF
          ! subname bioemi not needed (because lham=lmoz=false for lbioemi_stdalone)!
          CALL new_species(GAS, compound(i)%longname, compound(i)%shortname, 'mole mole-1', compound(i)%mweight,      &
                           lburden=.TRUE., lemis=.TRUE., idx=compound(i)%id_sp)

          CALL get_tracer (TRIM(compound(i)%shortname)//'_bioemi',ierr=ierr)
          IF (ierr==0) THEN
             ! this happens for lumped species
             ! new spec with _a, _b, _c
          END IF
          CALL new_tracer (TRIM(compound(i)%shortname)//'_bioemi',           &
               'bioemi',                           &
               spid=compound(i)%id_sp,                &
               nphase=GAS,                        &
               units='mole mole-1',                       &
               ninit=CONSTANT+RESTART,            &
               moleweight=compound(i)%mweight, &
               idx=compound(i)%id_tr           )
       END DO

    ELSE
       !---get species and tracer ids in case of megan as HAMMOZ sub module
       !   species loop
       DO i=1,ncompounds
          CALL query_species(shortname=TRIM(compound(i)%shortname), ierr=ierr, index=ispid)
          IF (ierr == 0) THEN
             WRITE(message_text, '("Species ",a,"(",a,") found at index ",i4)') TRIM(compound(i)%shortname), &
                   TRIM(compound(i)%longname), ispid
             CALL message('start_biogenic_emissions', message_text, level=em_info)
             compound(i)%id_sp = ispid
          END IF
       END DO
       !   tracer loop
       DO i=1,ncompounds
          CALL get_tracer (TRIM(compound(i)%shortname), idx=compound(i)%id_tr, ierr=ierr)
          IF (ierr == 0) THEN
             CALL message('start_biogenic_emissions', 'Using tracer '//TRIM(compound(i)%shortname)//    &
                          ' from submodel '//TRIM(trlist%ti(compound(i)%id_tr)%modulename), level=em_info)
             ! overwrite moleweight of megan with moleweight of tracer already defined
             rm=trlist%ti(compound(i)%id_tr)%moleweight
             IF ( ABS(rm-compound(i)%mweight)/compound(i)%mweight*100._dp .LE. 2._dp ) THEN
                compound(i)%mweight = rm
             ELSE
                CALL message ('start_biogenic_emissions','moleweight of existing '//TRIM(compound(i)%shortname)//          &
                              ' tracer deviates by more than 2% from moleweight defined in MEGAN', level=em_error)
             END IF
          END IF
       END DO
    END IF    !--end IF standalone

    !---check if any tracers were found/defined and abort if not
    IF (ANY(compound(:)%id_tr == 0)) THEN
        CALL message('start_biogenic_emissions', 'No tracer found or defined for the following MEGAN compounds:', &
                     level=em_info)
        DO i=1,ncompounds
           IF (compound(i)%id_tr == 0) THEN
               WRITE(message_text, '("     ",a,"(",a,")")') TRIM(compound(i)%shortname), TRIM(compound(i)%longname)
               CALL message('', message_text, level=em_info)
           END IF
       END DO
    END IF
    IF (ALL(compound(:)%id_tr == 0)) THEN
       CALL message('start_biogenic_emissions', 'No tracers found! Check submodel tracer names or '  &
                                  //'set lbioemi_stdalone=.false.', level=em_error)
    END IF

  END SUBROUTINE start_biogenic_emissions
  !-----------------------------------end start_biogenic_emissions-----------------------------------

  !-----------------------------------init_emi_biogenic----------------------------------------------
  !   set-up boundary conditions:
  !   - default: read plant functional types and assign to ibc_ef
  !   - some compounds may overwrite these ibc_ef with compound-specific emission factor maps
  !   compound-specific maps are available (from acd.ucar.edu/~guenther/MEGAN/MEGANv2.10_beta/EFglobal/)
  !   and have been interpolated (T31 or T63) for:
  !   - alpha pinene
  !   - beta pinene
  !   - 3-carene
  !   - isoprene
  !   - limonene
  !   - MBO
  !   - myrcene
  !   - nox
  !   - ocimene
  !   - sabinene

  SUBROUTINE init_emi_biogenic

    !---inherited types, data and functions
    USE mo_exception,                ONLY: finish, message, message_text, em_warn, em_info
    USE mo_control,                  ONLY: nlon, ngl
    USE mo_mpi,                      ONLY: p_parallel_io, p_io, p_pe
    USE mo_read_netcdf77,            ONLY: read_var_nf77_3d, read_var_nf77_4d
    USE mo_decomposition,            ONLY: lc => local_decomposition, global_decomposition
    USE mo_transpose,                ONLY: scatter_gp
    USE mo_boundary_condition,       ONLY: bc_find, bc_query, bc_define, bc_nml, BC_REPLACE, &
                                           BC_BOTTOM
    USE mo_external_field_processor, ONLY: EF_FILE, EF_MODULE, EF_LONLAT, EF_IGNOREYEAR, EF_CONSTANT, &
                                           EF_NOINTER
    USE mo_string_utls,             ONLY: st1_in_st2_proof

    !---local data
    INTEGER, PARAMETER            :: npftfile  = 17   ! 16 vegetation types plus bare soil
    CHARACTER(LEN=32), PARAMETER  :: pft_file  ='megan_clm_pft_map.nc            '
    CHARACTER(LEN=32), PARAMETER  :: efac_file ='megan_emission_factors.nc       '
    CHARACTER(LEN=64)             :: bc_name
    CHARACTER(LEN=32), PARAMETER  :: broadlde_file  ='megan_fraction_broadl_de.nc            '
    CHARACTER(LEN=32), PARAMETER  :: broadlev_file  ='megan_fraction_broadl_ev.nc            '
    CHARACTER(LEN=32), PARAMETER  :: borealtr_file  ='megan_boreal_trees_limit.nc            '
    CHARACTER(LEN=32), PARAMETER  :: borealgr_file  ='megan_boreal_grasses_limit.nc          '

    REAL(dp), ALLOCATABLE, TARGET :: zin(:,:,:)
    REAL(dp), POINTER             :: gl_data(:,:)
    REAL(dp), ALLOCATABLE, TARGET :: zinde(:,:,:,:)
    REAL(dp), POINTER             :: gl_datade(:,:)
    REAL(dp), ALLOCATABLE, TARGET :: zinev(:,:,:,:)
    REAL(dp), POINTER             :: gl_dataev(:,:)
    REAL(dp), ALLOCATABLE, TARGET :: zinbotr(:,:,:)
    REAL(dp), POINTER             :: gl_databotr(:,:)
    REAL(dp), ALLOCATABLE, TARGET :: zinbogr(:,:,:)
    REAL(dp), POINTER             :: gl_databogr(:,:)

    !>>csld #404 & #366 : renamed bc_ef1 and bc_ef2 into bc_ef_mod (for "ef from module")
    !                     and bc_ef_fil (for "ef from file") 
    TYPE(bc_nml)     :: bc_ef_mod, bc_ef_fil  
    !csld<<
    TYPE(bc_nml)     :: bc_frpft
    TYPE(bc_nml)     :: bc_fr_brdle
    TYPE(bc_nml)     :: bc_co2
    TYPE(bc_nml)     :: bc_sw
    INTEGER          :: i, ii, ires, ierr, ieftype, j
    LOGICAL          :: lex

    CHARACTER(len=2) :: tilnum
    
    LOGICAL :: emifact_flag(ncompounds)
    CHARACTER(LEN=ln) :: emifact_name(ncompounds)

    !---executable procedure

    IF (nef_pft == EF_FILE) THEN

       !---read PFT distribution from file
       !csld this is the default (pft distribution of clm4)

       IF (.NOT.ALLOCATED( pft )) ALLOCATE (pft(lc%nproma, lc%ngpblks, npftfile-1))   ! exclude bare soil
       IF (p_parallel_io) THEN
          ALLOCATE (zin(lc%nlon,lc%nlat,npftfile))
          INQUIRE (file=TRIM(pft_file), exist=lex)
          IF (lex) THEN
             CALL read_var_nf77_3d (pft_file, "lon", "lat", "pft", "PCT_PFT", zin, ierr)
             !---values of pft: 1..100 (percent) ==> calculate as percentage
             zin=zin/100.0
          ELSE
             CALL finish('init_emi_biogenic','missing input file '//TRIM(pft_file))
          END IF
       END IF
       NULLIFY (gl_data)

       ! start loop at index 2, because index 1 belongs to "bare soil"
       
       DO i=2,npftfile
          IF (p_pe == p_io) gl_data => zin(:,:,i)
          ii = i
          !   need to swap index 2 and 3
          IF (i == 2) ii = 3
          IF (i == 3) ii = 2
          CALL scatter_gp (gl_data, pft(:,:,ii-1), global_decomposition)
       END DO
       IF (p_parallel_io) THEN
          DEALLOCATE (zin)
       ENDIF
 
    ELSEIF (nef_pft == EF_MODULE) THEN    

           !csld the files broadlde_file, broadlev_file, borealtr_file and borealgr_file are necessary
           !     to later compute the PFT's. 
 
       !---read BL fractions distribution from file
       IF (.NOT.ALLOCATED( vegfract_de )) ALLOCATE (vegfract_de(lc%nproma, lc%ngpblks,1,1))  
       IF (.NOT.ALLOCATED( vegfract_ev )) ALLOCATE (vegfract_ev(lc%nproma, lc%ngpblks,1,1))  
       IF (p_parallel_io) THEN
          ALLOCATE (zinde(lc%nlon,lc%nlat,1,1))
          ALLOCATE (zinev(lc%nlon,lc%nlat,1,1))
          INQUIRE (file=TRIM(broadlde_file), exist=lex)
          IF (lex) THEN
             CALL read_var_nf77_4d (broadlde_file, "lon", "lat", "vegtype", "time", "vegfract", zinde, ierr)
          ELSE
             CALL finish('init_emi_biogenic','missing input file '//TRIM(broadlde_file))
          ENDIF
          INQUIRE (file=TRIM(broadlev_file), exist=lex)
          IF (lex) THEN
             CALL read_var_nf77_4d (broadlev_file, "lon", "lat", "vegtype", "time", "vegfract", zinev, ierr)
          ELSE
             CALL finish('init_emi_biogenic','missing input file '//TRIM(broadlev_file))
          ENDIF
       END IF
       NULLIFY (gl_datade)
       NULLIFY (gl_dataev)
       IF (p_pe == p_io) gl_datade => zinde(:,:,1,1)
       CALL scatter_gp (gl_datade, vegfract_de(:,:,1,1), global_decomposition)
       IF (p_pe == p_io) gl_dataev => zinev(:,:,1,1)
       CALL scatter_gp (gl_dataev, vegfract_ev(:,:,1,1), global_decomposition)
       IF (p_parallel_io) THEN
          DEALLOCATE (zinde)
          DEALLOCATE (zinev)
       ENDIF
       
       !    !---read boreal/temperate climatic threshold from file
       IF (.NOT.ALLOCATED( borealtr_limit )) ALLOCATE (borealtr_limit(lc%nproma, lc%ngpblks, 1))  
       IF (.NOT.ALLOCATED( borealgr_limit )) ALLOCATE (borealgr_limit(lc%nproma, lc%ngpblks, 1))  
       IF (p_parallel_io) THEN
          ALLOCATE (zinbotr(lc%nlon,lc%nlat,1))
          ALLOCATE (zinbogr(lc%nlon,lc%nlat,1))
          INQUIRE (file=TRIM(borealtr_file), exist=lex)
          IF (lex) THEN
             CALL read_var_nf77_3d (borealtr_file, "lon", "lat", "time", "borealtr_limit", zinbotr, ierr)
          ELSE
             CALL finish('init_emi_biogenic','missing input file '//TRIM(borealtr_file))
          ENDIF
          INQUIRE (file=TRIM(borealgr_file), exist=lex)
          IF (lex) THEN
             CALL read_var_nf77_3d (borealgr_file, "lon", "lat", "time", "borealgr_limit", zinbogr, ierr)
          ELSE
             CALL finish('init_emi_biogenic','missing input file '//TRIM(borealgr_file))
          ENDIF
       END IF
       NULLIFY (gl_databotr)
       NULLIFY (gl_databogr)
       IF (p_pe == p_io) gl_databotr => zinbotr(:,:,1)
       CALL scatter_gp (gl_databotr, borealtr_limit(:,:,1), global_decomposition)
       IF (p_pe == p_io) gl_databogr => zinbogr(:,:,1)
       CALL scatter_gp (gl_databogr, borealgr_limit(:,:,1), global_decomposition)
       IF (p_parallel_io) THEN
          DEALLOCATE (zinbotr)
          DEALLOCATE (zinbogr)
       ENDIF
    ENDIF 

    !---Set-up templates for EF boundary conditions
    !   EF means emission factor here and denotes compound-specific emission factor maps
    !   The default is to calculate emission factor maps based on PFT-specific emission factors (compound(i)%ef)
    bc_ef_mod%ef_type      = EF_MODULE                  ! calculated online at first time step
    bc_ef_mod%bc_mode      = BC_REPLACE
    bc_ef_mod%bc_domain    = BC_BOTTOM

    bc_ef_fil%ef_type      = EF_FILE
    bc_ef_fil%ef_template  = efac_file
    bc_ef_fil%ef_varname   = 'ef_<COMPOUND%LONGNAME>'   
    bc_ef_fil%ef_geometry  = EF_LONLAT
    bc_ef_fil%ef_timedef   = EF_CONSTANT                ! ### may need to be revised for long-term climate change runs
    bc_ef_fil%ef_timeindex = 1                          ! ### may need to be revised for long-term climate change runs
    bc_ef_fil%bc_mode      = BC_REPLACE
    bc_ef_fil%bc_domain    = BC_BOTTOM

    !---Define boundary conditions for vegetation emissions
    ibc_lai = -1
    ibc_ef(:) = -1
    ibc_emi(:) = -1
    ibc_frpft(1:npft_red_jsbach) = -1
    ibc_fr_brdle_ev = -1
    ibc_fr_brdle_de = -1
    ibc_co2 = -1
    ibc_sw = -1

    !-- locate boundary condition for emission mass flux (ibc_emi) of biogenic tracers
    !   only do this for compounds that are associated with tracer ids
    !   bc_define was called from init_emissions in emi_interface
    emifact_flag(:) = .FALSE.

    DO i=1,ncompounds
       emifact_name(i) = compound(i)%shortname
       emifact_flag(i) = st1_in_st2_proof(compound(i)%shortname, emifact_files_species)
    END DO

    DO i=1,ncompounds
       bc_name = 'BIOGENIC emissions of '//compound(i)%shortname
       IF (compound(i)%id_tr > 0) THEN
          CALL bc_find(bc_name, ires, ierr=ierr)
          IF (ierr == 0) THEN
             !---report use of boundary condition from emi_interface
             WRITE(message_text,'(a,i0)') 'Located boundary condition "'//TRIM(bc_name)//   &
                   '" for compound '//TRIM(compound(i)%longname)//' at index ', ires
             CALL message('init_emi_biogenic', message_text, level=em_info)
             ibc_emi(i) = ires
             !---define boundary condition for emission factor map
             IF (emifact_flag(i)) THEN 
                bc_ef_fil%ef_varname = 'ef_'//TRIM(compound(i)%shortname)
                ibc_ef(i) = bc_define(TRIM(compound(i)%longname)//' emission factors', bc_ef_fil, 2, .TRUE.)
             ELSE 
                ibc_ef(i) = bc_define(TRIM(compound(i)%longname)//' emission factors', bc_ef_mod, 2, .TRUE.)
             ENDIF

          ELSE     ! (ierr /= 0)
             WRITE(message_text, '(a)') 'Failed to find boundary condition "'//TRIM(bc_name)//'"'
             CALL message('init_emi_biogenic', message_text, level=em_warn)
          END IF
       END IF
    END DO

! ### modified from Colombe's work
! DEFINE PFT fractions from JSBACH
    bc_frpft%bc_domain = BC_BOTTOM
    bc_frpft%bc_mode = BC_REPLACE
    bc_frpft%ef_type = EF_MODULE
    bc_frpft%ef_actual_unit = '- (% of gridox)'

    DO i=1, npft_red_jsbach
        WRITE(tilnum,'(i2)') i
        ibc_frpft(i) = bc_define('Tile number '//TRIM(tilnum), bc_frpft, 2, .TRUE.)
    END DO

    bc_co2%ef_type      = EF_MODULE                  ! calculated online at first time step
    bc_co2%bc_mode      = BC_REPLACE
    bc_co2%bc_domain    = BC_BOTTOM
    bc_co2%ef_geometry  = EF_LONLAT
    bc_co2%ef_timedef       = EF_IGNOREYEAR

    ibc_co2 = bc_define('CO2 concentration ', bc_co2, 2, .TRUE.)

    bc_sw%ef_type      = EF_MODULE                  ! calculated online at first time step
    bc_sw%bc_mode      = BC_REPLACE
    bc_sw%bc_domain    = BC_BOTTOM
    bc_sw%ef_geometry  = EF_LONLAT
    bc_sw%ef_timedef   = EF_IGNOREYEAR
         
    ibc_sw = bc_define('Relative soil moisture', bc_sw, 2, .TRUE.)
   

  END SUBROUTINE init_emi_biogenic
  !-----------------------------------end init_emi_biogenic-----------------------------------------

  !-----------------------------------init_emi_biogenic_stream--------------------------------------
  SUBROUTINE init_emi_biogenic_stream

    ! MEGAN requires 24-hour averages of temperature and surface SW radiation.
    ! This averaging is done every radiation timestep.
    ! It is implemented using a circular linked list: each element in the list records
    ! the temperature / SW at one averaging time. 
    ! An index to the oldest and to the newest list element is maintained. 
    ! At each averaging step, the new 24-hour average is computed, the new temperature
    ! and SW fields are saved by overwriting the oldest field and the indexes to
    ! oldest / newest are updated, see subroutine biogenic_averages.
    ! This is a computationally efficient way of computing the averages, but it 
    ! makes handling of rerun and of the first day of an initial run a bit tricky. 

    !---inherited types, data and functions
    USE mo_memory_base,   ONLY: new_stream, default_stream_setting, &
                                add_stream_element
    USE mo_linked_list,   ONLY: NETCDF, SURFACE
    USE mo_time_control,  ONLY: delta_time, lstart, lresume, get_time_step, trigrad
    USE mo_time_event,    ONLY: TIME_INC_HOURS
    USE mo_exception,     ONLY: finish, message, message_text, em_warn, em_info
    USE mo_submodel,      ONLY: print_value, lbioemi_stdalone
    USE mo_boundary_condition,       ONLY: bc_find, bc_query, bc_define, bc_nml, BC_REPLACE, &
                                           BC_BOTTOM
    USE mo_external_field_processor, ONLY: EF_FILE, EF_LONLAT, EF_IGNOREYEAR, EF_CONSTANT, &
                                           EF_NOINTER, EF_MODULE, EF_VALUE

    IMPLICIT NONE

    !---subroutine interface ---
    !    -

    !---local data ---
    INTEGER                       :: i, j, k, istep, ierr, ix, ndx, len, ief_type
    CHARACTER(LEN=3)              :: ichar
    CHARACTER(len=ll+3)           :: helpname
    TYPE(bc_nml)                  :: bc_lai, bc_frpft

    !---executable procedure ---

    !---create fields used for 24-h averaging of temperature and surface SW radiation

    IF (trigrad%unit == TIME_INC_HOURS) THEN
       nstepsperday = 24/trigrad%counter
       IF (MOD(24,trigrad%counter) /= 0) nstepsperday = nstepsperday + 1
    ELSE
       CALL finish('init_emi_biogenic_stream', &
                   'Sorry, the biogenic emissions module only supports TRIGRAD intervals in hours')
    END IF

    ALLOCATE(T24(nstepsperday))
    ALLOCATE(ppfd24(nstepsperday))
    ALLOCATE(lai24(nstepsperday))

    CALL message ('init_emi_biogenic_stream','no. of steps per day:')
    CALL print_value('nstepsperday',nstepsperday)

    CALL new_stream(bioemiav, 'bioemi', filetype=NETCDF, lpost=ldebug_bioemi)  

    CALL default_stream_setting (bioemiav,lrerun  = .TRUE., &
                                         laccu     = .FALSE.,   &
                                         leveltype = SURFACE      )

    CALL add_stream_element(bioemiav, 'AVGTEMP', avg_tm1, leveltype=SURFACE, lrerun=.TRUE.,  &
         longname = '24h average land temperature', units = 'K')

    CALL add_stream_element(bioemiav, 'AVGPPFD', avg_ppfd,  leveltype=SURFACE, lrerun=.TRUE., &
         longname = '24h average photosynthetically active radiation', units = 'W m-2')

    CALL add_stream_element(bioemiav, 'LAIJSBACH', prt_lai,  leveltype=SURFACE, lrerun=.TRUE., &
         longname = 'LAi from JSBACH', units = 'm2 m-2')

    CALL add_stream_element(bioemiav, 'AVGLAI', avg_lai,  leveltype=SURFACE, lrerun=.TRUE., &
         longname = '24h average leaf area index', units = 'm2 m-2')

    CALL add_stream_element(bioemiav, 'LAI_daymean', LAI_daymean, leveltype=SURFACE, lrerun=.TRUE., &
         longname = 'daily average leaf area index', units = 'm2 m-2')

    CALL add_stream_element(bioemiav, 'LAI_monthmean', LAI_monthmean, leveltype=SURFACE, lrerun=.TRUE., &
         longname = 'monthly average leaf area index', units = 'm2 m-2')

    CALL add_stream_element(bioemiav, 'LAI_currentmonth', LAI_currentmonth, leveltype=SURFACE, lrerun=.TRUE., &
         longname = 'current month average leaf area index', units = 'm2 m-2')

    CALL add_stream_element(bioemiav, 'LAI_sum', LAI_sum,  leveltype=SURFACE, lrerun=.TRUE., &
         longname = 'daily sum leaf area index', units = 'm2 m-2')

    CALL add_stream_element(bioemiav, 'LAI_sum_month', LAI_sum_month, leveltype=SURFACE, lrerun=.TRUE., &
         longname = 'monthly sum leaf area index', units = 'm2 m-2')

    DO i=1,nstepsperday

       WRITE(ichar,'(i0)') i                           ! translate integer to character

       T24(i)%fldname = 'T24_'//TRIM(ADJUSTL(ichar))   ! field names for the rerun output stream
       ppfd24(i)%fldname = 'PPFD24_'//TRIM(ADJUSTL(ichar))
       lai24(i)%fldname = 'LAI24_'//TRIM(ADJUSTL(ichar))

       ! Add the fields to the output stream
       CALL message('init_emi_biogenic_stream','adding stream elements '//TRIM(ADJUSTL(T24(i)%fldname))// &
                                                                   ', '//TRIM(ADJUSTL(ppfd24(i)%fldname))// &
                                                                   ','//TRIM(ADJUSTL(lai24(i)%fldname)))

       CALL add_stream_element(bioemiav, T24(i)%fldname, T24(i)%fld2d, leveltype=SURFACE, units='K')
       CALL add_stream_element(bioemiav, ppfd24(i)%fldname, ppfd24(i)%fld2d, leveltype=SURFACE, units='Wm-2')
       CALL add_stream_element(bioemiav, lai24(i)%fldname, lai24(i)%fld2d, leveltype=SURFACE, units='m2 m-2')

    END DO

    !---add other diagnostic elements in debug mode
    IF (ldebug_bioemi) THEN
       DO i=1,ncompounds
          IF (compound(i)%id_tr > 0) THEN
             helpname='ef_'//TRIM(compound(i)%longname)
             len = LEN_TRIM(helpname)
             ndx=INDEX(helpname,' ')
             DO WHILE ((ndx > 0) .AND. (ndx < len))
               helpname(ndx:ndx)='_'
               ndx=INDEX(helpname,' ')
             ENDDO
             CALL add_stream_element(bioemiav, helpname,    &
                                     ef_diag(i)%ptr, leveltype=SURFACE, units='ug m-2 h-1')
          ENDIF
       END DO
       CALL add_stream_element(bioemiav, 'gamma_lai', gamma_diag(1)%ptr, leveltype=SURFACE, units='unitless')
       CALL add_stream_element(bioemiav, 'gamma_ppfd', gamma_diag(2)%ptr, leveltype=SURFACE, units='unitless')
       CALL add_stream_element(bioemiav, 'gamma_temp_iso', gamma_diag(3)%ptr, leveltype=SURFACE, units='unitless')
       CALL add_stream_element(bioemiav, 'gamma_temp_other', gamma_diag(4)%ptr, leveltype=SURFACE, units='unitless')
       CALL add_stream_element(bioemiav, 'gamma_age', gamma_diag(5)%ptr, leveltype=SURFACE, units='unitless')
       CALL add_stream_element(bioemiav, 'gamma_isoprene', gamma_diag(6)%ptr, leveltype=SURFACE, units='unitless')
       CALL add_stream_element(bioemiav, 'gamma_apinene', gamma_diag(7)%ptr, leveltype=SURFACE, units='unitless')
       CALL add_stream_element(bioemiav, 'flux_isop', gamma_diag(8)%ptr, leveltype=SURFACE, units='unitless')
!      CALL add_stream_element(bioemiav, 'gamma_methanol', gamma_diag(8)%ptr, leveltype=SURFACE, units='unitless')  
       IF (nef_pft == EF_MODULE) THEN  !csld #404 & #366 : zfrpft_jsbach is not assigned if nef_pft != EF_MODULE
          CALL add_stream_element(bioemiav, 'frac_PFT_1', gamma_diag(9)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_2', gamma_diag(10)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_3', gamma_diag(11)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_4', gamma_diag(12)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_5', gamma_diag(13)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_6', gamma_diag(14)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_7', gamma_diag(15)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_8', gamma_diag(16)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_9', gamma_diag(17)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_10', gamma_diag(18)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_11', gamma_diag(19)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_12', gamma_diag(20)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_13', gamma_diag(21)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_14', gamma_diag(22)%ptr, leveltype=SURFACE, units='unitless')
          CALL add_stream_element(bioemiav, 'frac_PFT_15', gamma_diag(23)%ptr, leveltype=SURFACE, units='unitless')
       ENDIF  !csld
       !CALL add_stream_element(bioemiav, 'tm1_sum', gamma_diag(23)%ptr, leveltype=SURFACE, units='unitless')
       !CALL add_stream_element(bioemiav, 'tm1_test', gamma_diag(24)%ptr, leveltype=SURFACE, units='unitless')
!       CALL add_stream_element(bioemiav, 'tcmax', gamma_diag(20)%ptr, leveltype=SURFACE, units='degrees')
!      CALL add_stream_element(bioemiav, 'gamma_ethene', gamma_diag(10)%ptr, leveltype=SURFACE, units='unitless')
    END IF

    !---initialise the linked lists at start and rerun
    IF (lstart .OR. lresume) THEN
       avg_tm1(:,:) = 0_dp
       avg_ppfd(:,:) = 0._dp
       avg_lai(:,:) = 0._dp
       LAI_daymean(:,:) = 0._dp
       LAI_monthmean(:,:) = 0._dp
       LAI_currentmonth(:,:) = 0._dp
       LAI_sum(:,:) = 0._dp
       LAI_sum_month(:,:) = 0._dp

       oldest_T24 => T24(1)
       oldest_ppfd24 => ppfd24(1)
       oldest_lai24 => lai24(1)
       newest_T24 => T24(nstepsperday)
       newest_ppfd24 => ppfd24(nstepsperday)
       newest_lai24 => lai24(nstepsperday)

       iavgstep = 0
    END IF

    DO i = 1, nstepsperday
       IF (i == nstepsperday) THEN
          T24(i)%nextfldseries => T24(1)               ! Circular linked lists
          ppfd24(i)%nextfldseries => ppfd24(1)
          lai24(i)%nextfldseries => lai24(1)
       ELSE
          T24(i)%nextfldseries => T24(i+1)
          ppfd24(i)%nextfldseries => ppfd24(i+1)
          lai24(i)%nextfldseries => lai24(i+1)
       END IF
    END DO

    IF (lresume) THEN

       istep = get_time_step()                   ! Beware, zero numbering applies

       k = (trigrad%counter*3600)/delta_time     ! no of timesteps per rad time step

       j = 1 + (istep / k)                       ! radiation calc. step no

       iavgstep = j

       ix = MOD(j, nstepsperday)
       IF (ix == 0) ix = nstepsperday

       CALL message('init_emi_biogenic_stream','settings')
       CALL print_value('resuming at time step',istep)
       CALL print_value('no. of steps per day',nstepsperday)
       CALL print_value('current averaging step',iavgstep)
       CALL print_value('index to latest step',ix)

       newest_T24 => T24(ix)
       newest_ppfd24 => ppfd24(ix)

       IF (iavgstep .GE. nstepsperday) THEN            ! more than 1 day elapsed before rerun
          oldest_T24 => newest_T24%nextfldseries
          oldest_ppfd24 => newest_ppfd24%nextfldseries
       ELSE                                            ! rerun from 1st day
          oldest_T24 => T24(1)
          oldest_ppfd24 => ppfd24(1)
       END IF

    END IF

    !---find / define LAI boundary condition
    ! Note: cannot be done in init_emi_biogenic, because "leaf area index" is defined in
    ! drydep_lg_init, which is called from subm_init_memory
    CALL bc_find('leaf area index', ibc_lai, ierr=ierr)
    IF (ierr /= 0) THEN
       bc_lai%bc_domain = BC_BOTTOM
       bc_lai%bc_mode = BC_REPLACE
       bc_lai%ef_type = nlai_biogenic_ef_type
       bc_lai%ef_template = 'surface_properties.nc'  !gf #226
       bc_lai%ef_varname = 'lai'
       bc_lai%ef_geometry = EF_LONLAT
       bc_lai%ef_timedef = EF_IGNOREYEAR
       bc_lai%ef_factor = 1._dp
       bc_lai%ef_interpolate = EF_NOINTER   ! none
       bc_lai%ef_actual_unit = 'm2 m-2'
       ibc_lai = bc_define('leaf area index', bc_lai, 2, .TRUE.)
    ELSE
       WRITE(message_text, '(a,i4)') 'Located boundary condition "leaf area index" at index ',ibc_lai
       CALL message('init_emi_biogenic_stream', message_text, level=em_info)
       CALL bc_query(ibc_lai,ef_type=ief_type)
       IF (ief_type /= nlai_biogenic_ef_type) THEN
         WRITE(message_text, '(a,i4)') 'boundary condition "leaf area index" overwritten by previous definition: ',ief_type
         CALL message('init_emi_biogenic_stream', message_text, level=em_warn)
       ENDIF
    END IF

  END SUBROUTINE init_emi_biogenic_stream
  !-----------------------------------end init_emi_biogenic_stream---------------------------------

  !-----------------------------------calculate adjustment factors: gamma -------------------------
  ! (from original MEGAN 2.10 source code)
  !     Scientific algorithm
  !
  !           Emission = [EF][GAMMA][RHO]
  !         where [EF]    = emission factor (ug/m2h)
  !               [GAMMA] = emission activity factor (non-dimension)
  !               [RHO]   = production and loss within plant canopies
  !                         (non-dimensino)
  !               Assumption: [RHO] = 1 (11/27/06) (See PDT_LOT_CP.EXT)
  !
  !           GAMMA  = [GAMMA_CE][GAMMA_age][GAMMA_SM]
  !         where [GAMMA_CE]  = canopy correction factor
  !               [GAMMA_age] = leaf age correction factor
  !               [GAMMA_SM]  = soil moisture correction factor
  !               Assumption: [GAMMA_SM]  = 1 (11/27/06)
  !
  !           GAMMA_CE = [GAMMA_LAI][GAMMA_P][GAMMA_T]
  !         where [GAMMA_LAI] = leaf area index factor
  !               [GAMMA_P]   = PPFD emission activity factor
  !               [GAMMA_T]   = temperature response factor
  !
  !           Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
  !      Derivation:
  !           Emission = [EF][GAMMA_etc](1-LDF) + [EF][GAMMA_etc][LDF][GAMMA_P]
  !           Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] }
  !           Emission = [EF][GAMMA_ect]{ (1-LDF) + [LDF][GAMMA_P] }
  !         where LDF = light dependent function (non-dimension)
  !
  !   For ISOPRENE
  !               Assumption: LDF = 1 for isoprene            (11/27/06)
  !
  !      Final Equation
  !           Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
  !
  !   For NON-ISOPRENE
  !      Final Equation
  !           Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age][GAMMA_SM]*
  !                      { (1-LDF) + [LDF][GAMMA_P] }
  !
  !------------------------------------------------------------------------------------------------

  ! ----- leaf area index adjustment
  ! (unchanged from v2.04)
  SUBROUTINE gamma_lai(kproma, kbdim, plai, pgamma)
  INTEGER, INTENT(in)   :: kproma, kbdim
  REAL(dp), INTENT(in)  :: plai(kbdim)
  REAL(dp), INTENT(out) :: pgamma(kbdim)

  ! ### assumes that LAI is zero where not valid ! ###
  pgamma(1:kproma) = (0.49_dp*plai(1:kproma)) / SQRT(1.0_dp+0.2_dp*(plai(1:kproma)**2))
  RETURN

  END SUBROUTINE gamma_lai

  ! (new in v2.10: extra routine for "bidirectional VOC" (class 18))
  !-----------------------------------------------------------------------
  !From Alex Guenther 2010-01-26
  !If lai < 2 Then
  !gammaLAIbidir= 0.5 * lai
  !ElseIf lai <= 6 Then
  !gammaLAIbidir= 1 - 0.0625 * (lai - 2)
  !Else
  !gammaLAIbidir= 0.75
  !End If
  !
  !    Xuemei Wang-2010-01-28
  !
  !-----------------------------------------------------------------------
!!     SUBROUTINE GAMMA_LAIbidir(NCOLS, NROWS,LAI,GAM_LAIbidir)
!!
!!        IMPLICIT NONE
!!
!!     INTEGER,INTENT(IN) :: NCOLS, NROWS
!!     INTEGER :: I,J
!!     REAL,DIMENSION(NCOLS, NROWS),INTENT(IN) ::  LAI
!!     REAL,DIMENSION(NCOLS, NROWS),INTENT(OUT) :: GAM_LAIbidir
!!      DO I = 1,NCOLS
!!      DO J = 1, NROWS
!!       IF(LAI(I,J) < 2) THEN
!!      GAM_LAIbidir =  0.5 * LAI
!!      ELSEIF (LAI(I,J) .LE. 6 .AND. LAI(I,J) .GE. 2) THEN
!!      GAM_LAIbidir = 1 - 0.0625 * (LAI(I,J) - 2)
!!      ELSE
!!      GAM_LAIbidir = 0.75
!!      ENDIF
!!
!!      ENDDO
!!      ENDDO
!!
!!      RETURN
!!    END  SUBROUTINE GAMMA_LAIbidir


  ! ----- photosynthetically active radiation adjustment
  ! (unchanged from v2.04)
  ! Modified by A. Henrot (december 2013-january 2014)

  SUBROUTINE gamma_ppfd(kproma, kbdim, pzenith, ppfd, tppfd, avg_ppfd, pgamma)

!             GAMMA_P = 0.0         a<=0, a>=180, sin(a) <= 0.0
!
!             GAMMA_P = sin(a)[ 2.46*(1+0.0005(Pdaily-400))*PHI - 0.9*PHI^2 ]
!                                   0<a<180, sin(a) > 0.0
!           where PHI    = above canopy PPFD transmission (non-dimension)
!                 Pdaily = daily average above canopy PPFD (umol/m2s)
!                 a      = solar angle (degree)
!                 PPFD   = photosynthetic flux density
!
!                 Note: AAA = 2.46*BBB*PHI - 0.9*PHI^2
!                       BBB = (1+0.0005(Pdaily-400))
!                       GAMMA_P = sin(a)*AAA
!
!                       Pac
!             PHI = -----------
!                   sin(a)*Ptoa
!           where Pac  = above canopy PPFD (umol/m2s)
!                 Ptoa = PPFD at the top of atmosphere (umol/m2s)
!
!             Pac =  SRAD * 4.766 mmmol/m2-s * 0.5
!
!             Ptoa = 3000 + 99*cos[2*3.14-( DOY-10)/365 )]
!           where DOY = day of year

  INTEGER, INTENT(in)   :: kproma, kbdim
  REAL(dp), INTENT(in)  :: pzenith(kbdim)
  REAL(dp), INTENT(in)  :: ppfd(kbdim), tppfd(kbdim), avg_ppfd(kbdim)
  REAL(dp), INTENT(out) :: pgamma(kbdim)

  REAL(dp)              :: zsinbeta(kbdim),PHI(kbdim)
  REAL(dp)              :: AAA(kbdim), BBB(kbdim)
  INTEGER               :: ji

  REAL(dp), PARAMETER :: conv_ppfd = 4.766_dp          ! factor to convert avg_ppfd in W/m2 to umol/m2/s
 
 ! ### assumes zenith angle in degrees
  pgamma(:) = 0._dp
  DO ji = 1, kproma
    IF (pzenith(ji) < 0.0_dp) THEN   !pzenith = SIN(solar angle)
      pgamma(ji) = 0.0_dp
    ELSEIF (pzenith(ji) > 0.0_dp) THEN

!  ### PPFD and PPFD at the top of the atmosphere are adjusted for solar angle in ECHAM 
    IF(tppfd(ji) <= 0.0_dp) THEN 
       PHI(ji) = 0.0_dp
    ELSEIF(tppfd(ji) > 0.0_dp) THEN
       PHI(ji) = ppfd(ji)/tppfd(ji)
    ENDIF

      BBB(ji) = 1._dp + 0.0005_dp*( avg_ppfd(ji)*conv_ppfd-400._dp )

      AAA(ji) = ( 2.46_dp * BBB(ji) * PHI(ji) ) - ( 0.9_dp * PHI(ji)**2 )

      pgamma(ji) = pzenith(ji) * AAA(ji)
!   ELSE
!     MESG = 'Error: Solar angle is invalid'
!     CALL M3EXIT(FUNCNAME,JDATE,JTIME,MESG,2)
    ENDIF

    IF (ASIN(pzenith(ji))*180._dp/3.14159_dp < 1.0_dp .AND. pgamma(ji) > 0.1_dp) THEN
      pgamma(ji) = 0.0_dp
    ENDIF
  END DO

  RETURN

!### OLD code
!!        DO jl = 1,kproma
!!           IF (loland(jl) .AND. vphysc%sw_flux_toa(jl,krow) > zeps) THEN
!!
!!              !---2.2 Sunlight dependency (eqn 11-12 in MEGAN paper)
!!              !       We already have ToA flux adjusted for solar angle, so (12) in MEGAN
!!              !       paper reduces to:
!!              phi = vphysc%sw_flux_surf(jl,krow) / vphysc%sw_flux_toa(jl,krow)
!!
!!              !--- sin(a) of the solar angle according to MEGAN is the cosine of the zenith angle
!!              !    amu0_x
!!              !    SW radiation in eqn. 11b in MEGAN paper is in umol m-2 s-1, we have surface fluxes
!!              !    in Wm-2. Convert using Wm2toppfd and cvrad. Factor of 5E-4 is from eqn 11b.
!!              !    See also Corrigendum to MEGAN paper in ACP 7, 4327, 2007.
!!
!!              zfac_p = zWm2_ppfd*cvrad
!!              gamma_p(jl) = amu0_x(jl,krow) * (2.46_dp*(1._dp + 5.E-4_dp * &
!!                            (zfac_p * avg_sw_flux_surf(jl,krow)-400._dp) * phi) - 0.9_dp*phi**2)
!!              [...]
!###

  END SUBROUTINE gamma_ppfd


  ! ----- temperature adjustment for isoprene
  ! (unchanged from v2.04)
  SUBROUTINE gamma_temp_iso(kproma, kbdim, ptemp, ptemp24, pgamma)
!                          Eopt*CT2*exp(CT1*x)
!             GAMMA_T =  ------------------------
!                        [CT2-CT1*(1-exp(CT2*x))]
!           where x      = [ (1/Topt)-(1/Thr) ] / 0.00831
!                 Eopt   = 1.75*exp(0.08(Tdaily-297)
!                 CT1    = 80
!                 CT2    = 200
!                 Thr    = hourly average air temperature (K)
!                 Tdaily = daily average air temperature (K)
!                 Topt   = 313 + 0.6(Tdaily-297)
!
!                 Note: AAA = Eopt*CT2*exp(CT1*x)
!                       BBB = [CT2-CT1*(1-exp(CT2*x))]
!                       GAMMA_T = AAA/BBB

  USE mo_physical_constants,          ONLY: argas

  INTEGER, INTENT(in)   :: kproma, kbdim
  REAL(dp), INTENT(in)  :: ptemp(kbdim), ptemp24(kbdim)  ! current and 24-h average air temperature
                                                         ! at canopy level (original: "hourly" instead of current)
  REAL(dp), INTENT(out) :: pgamma(kbdim)

  REAL(dp), PARAMETER   :: CT1 = 80._dp, CT2 = 200._dp
  REAL(dp), PARAMETER   :: invR = 1.E3_dp/argas          ! (1/R) expressed in (mol K / kJ)
  REAL(dp)              :: Eopt(kbdim), Topt(kbdim), X(kbdim)
  REAL(dp)              :: AAA(kbdim), BBB(kbdim)

  Eopt(1:kproma) = 1.75_dp * exp(0.08_dp*(ptemp24(1:kproma)-297.0_dp))
  Topt(1:kproma) = 313.0_dp + ( 0.6_dp*(ptemp24(1:kproma)-297.0_dp) )
!!  X(1:kproma) = ( (1.0_dp/Topt(1:kproma))-(1.0_dp/ptemp(1:kproma)) ) / 0.00831_dp
  X(1:kproma) = ( (1.0_dp/Topt(1:kproma))-(1.0_dp/ptemp(1:kproma)) ) * invR

  AAA(1:kproma) = Eopt(1:kproma)*CT2*exp(CT1*X(1:kproma))
  BBB(1:kproma) = (  CT2-CT1*( 1.0_dp-exp(CT2*X(1:kproma)) )  )

  pgamma(1:kproma) = AAA(1:kproma)/BBB(1:kproma)
  RETURN

  END SUBROUTINE gamma_temp_iso


  ! ----- temperature adjustment for non isoprene species
  ! (unchanged from v2.04)
  SUBROUTINE gamma_temp_other(kproma, kbdim, ptemp, pbeta, pgamma)
!             GAMMA_T =  exp[BETA*(T-Ts)]
!           where BETA   = temperature dependent parameter
!                 Ts     = standard temperature (normally 303K, 30C)

  INTEGER, INTENT(in)   :: kproma, kbdim
  REAL(dp), INTENT(in)  :: ptemp(kbdim)        ! current canopy level air temperature
  REAL(dp), INTENT(in)  :: pbeta               ! species dependent parameter
  REAL(dp), INTENT(out) :: pgamma(kbdim)

  REAL(dp), PARAMETER   :: Ts = 303._dp

  pgamma(1:kproma) = exp(pbeta * (ptemp(1:kproma)-Ts) )
  RETURN

  END SUBROUTINE gamma_temp_other

  !----- adjustment for leaf age
  !### need to verify !
  ! (unchanged from v2.04)
  SUBROUTINE gamma_age(kproma, kbdim, kclass, plaic, plaip, ptempd, pgamma)
!
!             GAMMA_age = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
!           where Fnew = new foliage fraction
!                 Fgro = growing foliage fraction
!                 Fmat = mature foliage fraction
!                 Fold = old foliage fraction
!                 Anew = relative emission activity for new foliage
!                 Agro = relative emission activity for growing foliage
!                 Amat = relative emission activity for mature foliage
!                 Aold = relative emission activity for old foliage
!
!
!             For foliage fraction
!             Case 1) LAIc = LAIp
!             Fnew = 0.0  , Fgro = 0.1  , Fmat = 0.8  , Fold = 0.1
!
!             Case 2) LAIp > LAIc
!             Fnew = 0.0  , Fgro = 0.0
!             Fmat = 1-Fold
!             Fold = (LAIp-LAIc)/LAIp
!
!             Case 3) LAIp < LAIc
!             Fnew = 1-(LAIp/LAIc)                       t <= ti
!                  = (ti/t) * ( 1-(LAIp/LAIc) )          t >  ti
!
!             Fmat = LAIp/LAIc                           t <= tm
!                  = (LAIp/LAIc) +
!                      ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )  t >  tm
!
!             Fgro = 1 - Fnew - Fmat
!             Fold = 0.0
!
!           where
!             ti = 5 + (0.7*(300-Tt))                   Tt <= 303
!                = 2.9                                  Tt >  303
!             tm = 2.3*ti
!
!             t  = length of the time step (days)
!             ti = number of days between budbreak and the induction of
!                  emission
!             tm = number of days between budbreak and the initiation of
!                  peak emissions rates
!             Tt = average temperature (K) near top of the canopy during
!                  current time period (daily ave temp for this case)
!             LAIc = current month's LAI
!             LAIp = previous month's LAI
!
!             Relative emission activity
!             Class           Anew       Agro       Amat       Aold
!             Constant         1.0        1.0        1.0        1.0
!             Monoterpenes     2.0        1.8        0.95       1.0
!             Sesquiterpenes   0.4        0.6        1.075      1.0
!             Methanol         3.0        2.6        0.85       1.0
!             Isoprene         0.05       0.6        1.125      1.0
!
!             Species:
!             Constant:  'ACTO','ACTA','FORM','CH4','NO','CO'
!             Monoterpenes: 'MYRC','SABI','LIMO','3CAR','OCIM','BPIN','APIN','OMTP'
!             Sesquiterpenes: 'FARN','BCAR','OSQT'
!             Methanol: 'MEOH'
!             Isoprene: 'C5H8'
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in)   :: kproma, kbdim
  INTEGER, INTENT(in)   :: kclass              ! compound class
  REAL(dp), INTENT(in)  :: plaic(kbdim), plaip(kbdim)  ! LAI of current and previous month
  REAL(dp), INTENT(in)  :: ptempd(kbdim)       ! daily average above canopy temperature
  REAL(dp), INTENT(out) :: pgamma(kbdim)


  INTEGER, PARAMETER    :: NCLASS = 6          ! constant,isoprene,monoterpenes,sesquiterpenes,MBO,methanol
                                               ! updated from Guenther et al. 2012
  INTEGER, PARAMETER    :: NAGE = 4            ! new, grown, mature, old

! Fortran 90/95 reference: For 2D matrices, values are listed column wise.

  REAL(dp), PARAMETER   :: afac(NCLASS, NAGE) = &
                           RESHAPE(SOURCE = (/  1.0_dp,  0.05_dp,  2.0_dp,   0.4_dp,   0.05_dp,  3.5_dp,     &
                                                1.0_dp,  0.6_dp,   1.8_dp,   0.6_dp,   0.6_dp,   3.0_dp,     &
                                                1.0_dp,  1.0_dp,   1.0_dp,   1.0_dp,   1.0_dp,   1.0_dp,     &
                                                1.0_dp,  0.9_dp,   1.05_dp,  0.95_dp,  0.9_dp,   1.2_dp /), &
                                  SHAPE = (/NCLASS, NAGE/))
  REAL(dp), PARAMETER   :: t = 30._dp  ! original time step in MEGAN (set to TSTLEN = 30.)
                                       ! ### this is unclear! See formulae for fnew and fmat below

  REAL(dp)              :: Fnew, Fgro, Fmat, Fold     ! fraction of new, grown, mature and aged leafs

  REAL(dp)              :: ti          ! number of days between budbreak
                                       ! and the induction of emission
  REAL(dp)              :: tm          ! number of days between budbreak
                                       ! and the initiation of peak
                                       ! emissions rates

  INTEGER               :: ji          ! loop index


  ! ### this code may break if either plaic or plaip is == 0  !!! ###
  DO ji = 1, kproma
!...  Calculate foliage fraction
    IF (plaip(ji) == plaic(ji)) THEN     ! mature state
       Fnew = 0.0_dp
       Fgro = 0.1_dp
       Fmat = 0.8_dp
       Fold = 0.1_dp

    ELSEIF (plaip(ji) > plaic(ji)) THEN  ! senescing state
       Fnew = 0.0_dp
       Fgro = 0.0_dp
       Fold = ( plaip(ji)-plaic(ji) ) / plaip(ji)
       Fmat = 1.0_dp-Fold

    ELSEIF (plaip(ji) < plaic(ji)) THEN  ! growing season
!      Calculate ti and tm
       IF (ptempd(ji) <= 303.0_dp) THEN
          ti = 5.0_dp + 0.7_dp*(300._dp-ptempd(ji))
       ELSEIF (ptempd(ji) > 303.0_dp) THEN
          ti = 2.9_dp
       ENDIF
       tm = 2.3_dp*ti

!      Calculate Fnew and Fmat, then Fgro and Fold
! ### questionable algorithm! This seems to implicitly assume that time step for lai and "t" is the same (= 1 month)!
! ### also: if t <= tm (by definition tm>ti) then there will never be grown leaves, because Fnew+Fmat = 1
!      Fnew
       IF (t <= ti) THEN
          Fnew = 1.0_dp - (plaip(ji)/plaic(ji))
       ELSEIF (t > ti) THEN
          Fnew = (ti/t) * ( 1.0_dp-(plaip(ji)/plaic(ji)) )
       ENDIF

!      Fmat
       IF (t <= tm) THEN
          Fmat = plaip(ji)/plaic(ji)
       ELSEIF (t > tm) THEN
          Fmat = (plaip(ji)/plaic(ji)) + ( (t-tm)/t ) * ( 1-(plaip(ji)/plaic(ji)) )
       ENDIF

       Fgro = 1.0 - Fnew - Fmat
       Fold = 0.0

    ENDIF

!...  Calculate GAMMA_A
    pgamma(ji) = Fnew*afac(kclass,1) + Fgro*afac(kclass,2) + Fmat*afac(kclass,3) + Fold*afac(kclass,4)

  END DO

  RETURN
  END SUBROUTINE gamma_age

!-------------------- 
  SUBROUTINE gamma_sm(kproma,kbdim,psoilm,pgamma)

!! SUBROUTINE GAMMA_SM(NCOLS, NROWS, SM, WP, GAM_SM)
!!
  ! --- adjustment for soil moisture
  !  gamma_sm = 1 in standard version of MEGAN v2.10
!! From Guenther et al., 2012, MEGANv2.1
!! soil moisture dependence algorithm for isoprene emission:
!! 
!! GAM_SM = 1              IF SM > SM1
!! GAM_SM = (SM-WP)/DSM1   IF WP < SM < SM1
!! GAM_SM = 0              IF SM < WP
!!
!! where SM = soil moisture (m3/m3)
!!       WP = wilting point (m3/m3)
!!       DSM1 = 0.04 empirical parameter (from Pegoraro et al., 2004; modified in MEGANv2.1)
!!       SM1 = WP + DSM1
!----------------------------------------------------------------------------------------------
  INTEGER,  INTENT(in)  :: kproma, kbdim
  REAL(dp), INTENT(in)  :: psoilm(kbdim)
  REAL(dp), INTENT(out) :: pgamma(kbdim)
 
  REAL(dp) :: wiltp(kbdim)
  REAL(dp) :: soilm1(kbdim)
  REAL(dp), PARAMETER :: delta1 = 0.06_dp
  INTEGER  :: ji

  wiltp(:) = 0.35_dp
  soilm1(:) = wiltp(:) + delta1
 
  DO ji = 1, kproma
    
     IF(psoilm(ji) >= soilm1(ji)) THEN 
   
        pgamma(ji) = 1.0_dp
   
     ELSEIF(psoilm(ji) < soilm1(ji) .AND. psoilm(ji) >= wiltp(ji)) THEN
 
          pgamma(ji) = (psoilm(ji)-wiltp(ji))/delta1 
 
     ELSEIF(psoilm(ji) < wiltp(ji)) THEN
 
          pgamma(ji) = 0.0_dp
     
     ENDIF
   
  END DO
 
  RETURN
 
  END SUBROUTINE gamma_sm


  ! --- adjustment for CO2 inhibition
  !           GAMMA_CO2 =     1.0   (non-dimension)
  !           When CO2 =400ppm
  !
  !   SUBROUTINE GAM_CO2 returns the GAMMA_CO2 values
  !  Xuemei Wang-2009-06-22
!-----------------------------------------------------------------------
!!    SUBROUTINE GAMMA_CO2( NCOLS, NROWS, CO2, GAM_CO2 )
!!
!!    IMPLICIT NONE
!!
!!    INTEGER :: NCOLS, NROWS
!!    REAL,DIMENSION(NCOLS,NROWS) :: GAM_CO2,CO2,Ci
!!    REAL, PARAMETER :: ISmax = 1.344, h=1.4614
!!    REAL, PARAMETER :: Cstar =585
!!
!!
!!    Ci = 0.7* CO2
!!    WHERE (CO2.eq.400.0)
!!       GAM_CO2 = 1.0
!!    ELSEWHERE
!!       GAM_CO2 = ISmax- ((ISmax*Ci**h) /(Cstar**h+Ci**h))
!!    ENDWHERE
!!
!!    RETURN
!!    END SUBROUTINE GAMMA_CO2
!!----------------------------------------------------------------------

  SUBROUTINE gamma_co2(kproma,kbdim,pco2,pgamma)
   
   INTEGER,  INTENT(in)  :: kproma, kbdim
   REAL(dp), INTENT(in)  :: pco2(kbdim)

   REAL(dp), INTENT(out) :: pgamma(kbdim)

   REAL(dp) :: Ci(kbdim)
   REAL(dp) :: CO2_ppm(kbdim)
   REAL(dp), PARAMETER :: Ismax = 1.344_dp
   REAL(dp), PARAMETER :: h = 1.4614_dp
   REAL(dp), PARAMETER :: Cstar = 585.0_dp
   REAL(dp), PARAMETER :: molmassCO2_kg = 44.011E-3_dp  ! mass of 1 mol of CO2 in kg
   REAL(dp), PARAMETER :: molmassAir_kg = 28.97E-3_dp   ! mass of 1 mol of CO2 in kg
   INTEGER :: ji

   CO2_ppm(:) = pco2(:)*(molmassAir_kg/molmassCO2_kg)*1E6_dp     ! conversion of mass mixing ratio to volume 
   DO ji= 1, kproma                                              ! mixing ratio and then ppmv
  
   Ci(1:kproma) = 0.7*CO2_ppm(1:kproma)
   
   IF (pco2(ji) == 400.0_dp) THEN
      pgamma(ji) = 1.0_dp
   ELSE
      pgamma(ji) = Ismax - ((Ismax*Ci(ji)**h)/(Cstar**h + Ci(ji)**h))
   ENDIF

   END DO
  RETURN

  END SUBROUTINE gamma_co2


  !-----------------------------------calc_biogenic_emissions------------------------------------------------
  SUBROUTINE calc_biogenic_emissions(kproma, kbdim, krow, loland, ptm1)
    !
    ! Calculates the isoprene, monoterpene and other biogenic VOC emission fluxes
    ! for the given time step.

    USE mo_exception,           ONLY: finish, message_text
    USE mo_geoloc,              ONLY: amu0_x
    USE mo_time_control,        ONLY: lfirst_day, l_trigrad, lstart,lresume
    USE mo_vphysc,              ONLY: vphysc
    USE mo_hammoz_vegetation,   ONLY: cvrad
    USE mo_surface_memory,      ONLY: jsswpar      ! photosynthetically active radiation (PAR = PPFD)
    USE mo_external_field_processor, ONLY: EF_MODULE, EF_FILE
    USE mo_boundary_condition,  ONLY: bc_query, bc_apply, bc_set, bc_find
    
    IMPLICIT NONE

    INTEGER,  INTENT(in) :: kproma             ! geographic block number of locations
    INTEGER,  INTENT(in) :: kbdim              ! geographic block maximum number of locations
    INTEGER,  INTENT(in) :: krow               ! geographic block number
    LOGICAL,  INTENT(in) :: loland(kbdim)      ! land mask
    REAL(dp), INTENT(in) :: ptm1(kbdim)        ! temperature at t-dt

    !--- local declarations
    !    Parameters
    REAL(dp), PARAMETER :: zWm2_ppfd = 4.766_dp          ! factor to convert W to umol s-1
    REAL(dp), PARAMETER :: zrho_iso = 0.96_dp            ! MEGAN model parameter
    REAL(dp), PARAMETER :: cce = 1.00_dp                 ! MEGAN canopy model factor ### canopy model not yet implemented! ###

    !--- local variables ---

    REAL(dp) :: zlai(kbdim)                     ! leaf area index
    REAL(dp) :: zppfd(kbdim)                    ! photosynthetically active radiation for gamma_ppfd
    REAL(dp) :: ztppfd(kbdim)                   ! top SW radiation for gamma_ppfd
    REAL(dp) :: zamu(kbdim)                     ! cos(zenith angle) = sin(zenith angle) in MEGAN geometry
    REAL(dp) :: zgamma_t(kbdim),             &  ! temperature dependency
                zgamma_ppfd(kbdim),          &  ! light dependency
                zgamma_lai(kbdim),           &  ! leaf area index dependency
                zgamma_age(kbdim),           &  ! leaf area index dependency
                zgamma_co2(kbdim),           &  ! CO2 inhibition factor for isoprene
                zgamma_sm(kbdim),            &  ! soil moisture dependency
                zgamma(kbdim) ,              &  ! combined iadjustment factor
                zfrpft_jsbach(kbdim,NPFT),   &
                ztest(kbdim)

               
    REAL(dp) :: zco2conc(kbdim)               
    REAL(dp) :: zpco2_4emi(kbdim)

    REAL(dp) :: zsw(kbdim)               
    REAL(dp) :: zrelsw(kbdim)

    REAL(dp) :: zef(kbdim)                      ! compound emission factor
    REAL(dp) :: zflux(kbdim)                    ! compound emission flux in kg m-2 s-1

    REAL(dp) :: zeps, qnsteps, zfac_p

    INTEGER :: i, j,k, jl, itrc, ieftype, nm, ierr

   
    !--- executable procedure

    ! ### vegetation fraction from JSBACH
    IF (nef_pft == EF_MODULE) CALL calc_emissions_factors(kproma,kbdim,krow,zfrpft_jsbach)

    !--- at first time step calculate appropriate emission factors
    IF (lstart .OR. lresume) THEN
       DO i=1,ncompounds
          IF (ibc_ef(i) > 0) THEN
             CALL bc_query(ibc_ef(i), ef_type=ieftype)
             ! calculate emission factors based on PFT-specific values
             ! and set as boundary condition
             ! (nothing to be done if those are replaced by compound-specific file input)
             IF (ieftype == EF_MODULE) THEN
                zef(1:kproma) = 0._dp

                IF (nef_pft == EF_FILE) THEN !PFT fractions from MEGAN-CLM4
                   
                   DO j=1,NPFT
                      zef(1:kproma) = zef(1:kproma) + compound(i)%ef(j)*pft(1:kproma,krow,j)   ! units: ug m-2 h-1
                   END DO
                   
                ELSEIF (nef_pft == EF_MODULE) THEN
                   
                   DO j=1,NPFT          !PFT fractions from JSBACH
                      zef(1:kproma) = zef(1:kproma) + compound(i)%ef(j)*zfrpft_jsbach(1:kproma,j) ! units: ug m-2 h-1
                   END DO
                   
                END IF
                CALL bc_set(ibc_ef(i), kproma, krow, zef(1:kproma))
             END IF
          END IF
       END DO
    END IF


    !---initialisations
    zlai(:) = 0._dp
    zgamma_lai(:) = 0._dp
    zgamma_ppfd(:) = 0._dp
    zgamma_t(:) = 0._dp
    zgamma_age(:) = 0._dp
    zgamma_co2(:) = 0._dp
    zgamma_sm(:) = 0._dp
    zef(:) = 0._dp
    zflux(:) = 0._dp
    zppfd(:) = 0._dp
    ztppfd(:) =0._dp
    zamu(:) = 0._dp
    ztest(:) = 0._dp
    zco2conc(:) = 0._dp
    zsw(:) = 0._dp

!!    zppfd(1:kproma) = jsswpar(1:kproma, krow)
!! calculate normalized PAR in gamma_ppfd subroutine
    zppfd(1:kproma)  = vphysc%sw_flux_surf(1:kproma,krow)
    ztppfd(1:kproma) = vphysc%sw_flux_toa(1:kproma,krow) 
    zamu(1:kproma) = amu0_x(1:kproma, krow)



    !---Get LAI values
    CALL bc_apply(ibc_lai, kproma, krow, zlai)
    prt_lai(1:kproma,krow) = zlai(1:kproma)

!>>AH    
    !---biogenic averages...calculated each radiation timestep
    IF (l_trigrad) CALL biogenic_averages(kbdim, kproma, krow, ptm1, zlai)
!<<AH

    !---Get CO2 concentration
    IF (lstart .OR. lresume) THEN
       CALL bc_find ('CO2 concentration', ibc_co2, ierr=ierr)
       IF (ierr /= 0 .OR. ibc_co2 <=0) CALL finish ('calc bioemi', &
          'cannot locate boundary condition for CO2 concentration')
    ENDIF
    CALL bc_apply(ibc_co2,kbdim, krow, zpco2_4emi(:))
 
    zco2conc(:) = zpco2_4emi(:)
    
    !---Get relative soil moisture
    
    IF (lstart .OR. lresume) THEN
       CALL bc_find ('Relative soil moisture', ibc_sw, ierr=ierr)
       IF (ierr /= 0 .OR. ibc_sw <=0) CALL finish ('calc bioemi', &
          'cannot locate boundary condition for relative soil moisture')
    ENDIF
    CALL bc_apply(ibc_sw,kbdim, krow, zrelsw(:))
    zsw(:) = zrelsw(:)
    
    
    !--- calculate vegetation emissions
    !    Isoprene:
    !         Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
    !
    !    non Isoprene:
    !         Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age][GAMMA_SM]*
    !                    { (1-LDF) + [LDF][GAMMA_P] }
    !--- first compute gamma factors independent of compound
    CALL gamma_lai( kproma, kbdim, zlai, zgamma_lai)
    CALL gamma_ppfd(kproma, kbdim, zamu, zppfd, ztppfd, avg_ppfd(:,krow), zgamma_ppfd) !SF #483 bugfix

    !-- compound loop: compute gamma and fluxes for each MEGAN compound that shall be emitted
    !   in HAMMOZ
    DO i=1,ncompounds
      IF (ibc_emi(i) > 0) CALL bc_set(ibc_emi(i), kproma, krow, zflux(1:kproma))    ! set to zero (see above)
    END DO

    DO i=1,ncompounds
       IF (ibc_emi(i) > 0) THEN
          !### safety measure: should be possible to remove in production code!
          IF (ibc_ef(i) <= 0) THEN
             WRITE(message_text, '(a,i4)') 'ibc_emi > 0, but ibc_ef <= 0 ! Compound index ', i
             CALL finish('calc_biogenic_emissions', message_text)
          END IF
          ! retrieve emission factor for standard conditions
          CALL bc_apply(ibc_ef(i), kproma, krow, zef)
          CALL bc_apply(ibc_emi(i), kproma, krow, zflux)
          CALL gamma_age(kproma, kbdim,compound(i)%cafac, LAI_daymean(:,krow), LAI_monthmean(:,krow), avg_tm1(:,krow), zgamma_age)
          IF (i == 1) THEN
             ! isoprene special case
             CALL gamma_temp_iso(kproma, kbdim, ptm1, avg_tm1(:,krow), zgamma_t) !SF #483 bugfix
             CALL gamma_co2(kproma,kbdim,zco2conc,zgamma_co2)
             CALL gamma_sm(kproma,kbdim,zsw,zgamma_sm)
             IF (lfirst_day) THEN
                zflux(1:kproma) = 0._dp     ! need to wait for 24h temperature average
             ELSE
                zflux(1:kproma) = zflux(1:kproma) + zfac * zef(1:kproma) * zgamma_lai(1:kproma) * zgamma_ppfd(1:kproma)  &
                                  * zgamma_t(1:kproma) * zgamma_age(1:kproma) !* zgamma_co2(1:kproma) !* zgamma_sm(1:kproma)
             END IF
!### OLD code
!!              zisop_flux(jl) = gamma_p(jl)*gamma_t(jl)*gamma_lai(jl)  &
!!                               * zef_isop(jl) * zrho_iso * zisop2oc       !>>dod removed dead code <<dod
!!            ## note: zisop2oc not needed --> zfac in bc_set !
!!            ## WHAT ABOUT zrho_iso ???
!###
          ELSE
             CALL gamma_temp_other(kproma, kbdim, ptm1, compound(i)%beta, zgamma_t)
!>>AH
!             zflux(1:kproma) = zflux(1:kproma) + zfac * zef(1:kproma) * zgamma_lai(1:kproma) * zgamma_ppfd(1:kproma)  &
!                               * zgamma_t(1:kproma) * zgamma_age(1:kproma)
             zflux(1:kproma) = zflux(1:kproma) + zfac * zef(1:kproma) * zgamma_lai(1:kproma) &
                               *((1-compound(i)%ldf)+compound(i)%ldf * zgamma_ppfd(1:kproma))  &
                               * zgamma_t(1:kproma)   &
                               * zgamma_age(1:kproma)
!<<AH

!### OLD code
!!           mea2(jl) = zlai(jl) / 5._dp
!!           hea2(jl) = EXP(0.09*(ptm1(jl) - 303.15))
!!           zmterp_flux(jl) =  hea2(jl) * mea2(jl) * zef_mterp(jl) * zmterp2oc     !>>dod removed dead code <<dod
!###
          END IF
          ! ### WARNING!! This should ADD data to the current BC values, not overwrite them!!
          ! ### This is important for monoterpene emissions which are added from individual compounds !!!
          CALL bc_set(ibc_emi(i), kproma, krow, zflux(1:kproma))    ! convert to kg m-2 s-1
          ! store diagnostics
          IF (ldebug_bioemi) THEN
             ef_diag(i)%ptr(1:kproma,krow) = zef(1:kproma)
             IF (i == 1) THEN
                gamma_diag(1)%ptr(1:kproma,krow) = zgamma_lai(1:kproma)
                gamma_diag(2)%ptr(1:kproma,krow) = zgamma_ppfd(1:kproma)
                gamma_diag(3)%ptr(1:kproma,krow) = zgamma_t(1:kproma)
                gamma_diag(5)%ptr(1:kproma,krow) = zgamma_age(1:kproma)
                gamma_diag(6)%ptr(1:kproma,krow) = zgamma_lai(1:kproma) * zgamma_ppfd(1:kproma) * zgamma_t(1:kproma) 
                gamma_diag(8)%ptr(1:kproma,krow) = zflux(1:kproma) 
                IF (nef_pft == EF_MODULE) THEN 
                   gamma_diag(9)%ptr(1:kproma,krow)  = zfrpft_jsbach(1:kproma,1) 
                   gamma_diag(10)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,2) 
                   gamma_diag(11)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,3) 
                   gamma_diag(12)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,4) 
                   gamma_diag(13)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,5) 
                   gamma_diag(14)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,6) 
                   gamma_diag(15)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,7) 
                   gamma_diag(16)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,8) 
                   gamma_diag(17)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,9) 
                   gamma_diag(18)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,10) 
                   gamma_diag(19)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,11) 
                   gamma_diag(20)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,12) 
                   gamma_diag(21)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,13) 
                   gamma_diag(22)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,14) 
                   gamma_diag(23)%ptr(1:kproma,krow) = zfrpft_jsbach(1:kproma,15)
                END IF
             ELSE
                gamma_diag(4)%ptr(1:kproma,krow) = zgamma_t(1:kproma)
                IF (i == 2) &
                   gamma_diag(7)%ptr(1:kproma,krow) = zgamma_lai(1:kproma) * zgamma_ppfd(1:kproma) * zgamma_t(1:kproma)
             END IF
          END IF
       END IF
    END DO

  END SUBROUTINE calc_biogenic_emissions
!---------------------------------end calc_biogenic_emissions------------------------------------------

!-------------------------------------------biogenic_averages------------------------------------------

  SUBROUTINE biogenic_averages(kbdim,kproma,krow,ptm1,plai)

    ! Computes 24h moving averages for net surface irradiance and lowest model level temperature

    ! MEGAN requires 24-hour averages of temperature and surface SW radiation.
    ! This averaging is done every radiation timestep.
    ! It is implemented using a circular linked list: each element in the list records
    ! the temperature / SW at one averaging time.
    ! An index to the oldest and to the newest list element is maintained.
    ! At each averaging step, the new 24-hour average is computed, the new temperature
    ! and SW fields are saved by overwriting the oldest field and the indexes to
    ! oldest / newest are updated
    ! This is a computationally efficient way of computing the averages, but it
    ! makes handling of the first day of an initial run a bit tricky.

  USE mo_time_control,   ONLY: get_time_step, &
                               get_date_components, current_date, previous_date, lstart
  USE mo_vphysc,         ONLY: vphysc
  USE mo_debugs
  USE mo_exception,           ONLY: finish, message_text
  USE mo_jsbach,         ONLY: new_day, new_month, new_year

  IMPLICIT NONE

! compulsory arguments

    INTEGER, INTENT(in) :: kproma                     ! geographic block number of locations
    INTEGER, INTENT(in) :: kbdim                      ! geographic block maximum number of locations

! optional arguments

    INTEGER,  INTENT(in),    OPTIONAL :: krow                      ! geographic block number
    REAL(dp), INTENT(in),    OPTIONAL :: ptm1(kbdim)               ! temperature at t-dt
    REAL(dp), INTENT(in),    OPTIONAL :: plai(kbdim)               ! LAI at t-dt

  !---local data
  !   Parameters
  !   -

  !---local variables
  REAL(dp) :: qnsteps
  LOGICAL  :: lnewtimestep
  INTEGER  :: istep
  INTEGER  :: day_of_month
  INTEGER  :: day_of_current_month
  INTEGER  :: month_lai
  LOGICAL  :: first_day_lai

  !---executable procedure

  qnsteps = 1._dp / REAL(nstepsperday,dp)

  istep = get_time_step()                 ! Model timestep number

  IF (istep /= itimestep) THEN            ! cycle the linked lists only at new time step
     lnewtimestep = .TRUE.                ! (this subroutine is called for every row and
     itimestep = istep                    !  radiation timestep)
     iavgstep = iavgstep + 1              ! averaging step number (=no. of radiation calculation
  ELSE                                    ! steps)
     lnewtimestep = .FALSE.
  END IF

  IF (iavgstep > nstepsperday) THEN        ! if TRUE then linked list complete ...tried using
                                           ! lfirst_day (mo_time_control) but that is not reset
                                           ! until the second time step of the second day....
     !---compute new averages
     avg_ppfd(1:kproma,krow) = avg_ppfd(1:kproma,krow) + qnsteps *&
                                       (vphysc%sw_flux_surf(1:kproma,krow) - oldest_ppfd24%fld2d(1:kproma,krow))

     avg_tm1(1:kproma,krow) = avg_tm1(1:kproma,krow) + qnsteps *&
                              (ptm1(1:kproma) - oldest_T24%fld2d(1:kproma,krow))

     avg_lai(1:kproma,krow) = avg_lai(1:kproma,krow) + qnsteps *&
                              (plai(1:kproma) - oldest_lai24%fld2d(1:kproma,krow))

     !--- cycle the circular linked lists
     IF (lnewtimestep) THEN
        oldest_T24 => oldest_T24%nextfldseries
        newest_T24 => newest_T24%nextfldseries
        oldest_ppfd24 => oldest_ppfd24%nextfldseries
        newest_ppfd24 => newest_ppfd24%nextfldseries
        oldest_lai24 => oldest_lai24%nextfldseries
        newest_lai24 => newest_lai24%nextfldseries
     END IF

     newest_ppfd24%fld2d(1:kproma,krow) = vphysc%sw_flux_surf(1:kproma,krow)

     newest_T24%fld2d(1:kproma,krow) = ptm1(1:kproma)

     newest_lai24%fld2d(1:kproma,krow) = plai(1:kproma)

  ELSE              ! first day

     ppfd24(iavgstep)%fld2d(1:kproma,krow) = vphysc%sw_flux_surf(1:kproma,krow)

     T24(iavgstep)%fld2d(1:kproma,krow) = ptm1(1:kproma)

     lai24(iavgstep)%fld2d(1:kproma,krow) = plai(1:kproma)

     avg_ppfd(1:kproma,krow) = avg_ppfd(1:kproma,krow) + &
                                       qnsteps*vphysc%sw_flux_surf(1:kproma,krow)

     avg_tm1(1:kproma,krow) = avg_tm1(1:kproma,krow) + qnsteps*ptm1(1:kproma)

     avg_lai(1:kproma,krow) = avg_lai(1:kproma,krow) + qnsteps*plai(1:kproma)

     IF (lnewtimestep) THEN
        newest_T24 => T24(iavgstep)
        newest_ppfd24 => ppfd24(iavgstep)
        newest_lai24 => lai24(iavgstep)
     END IF

  END IF

  !LAI daily and monthly mean calculation, needed for gamma_age
  !based on mo_climbuff.f90

  ! get first day
  first_day_lai = (istep < nstepsperday)

  IF( .NOT. new_day .OR. lstart) THEN

    LAI_sum(1:kproma,krow) = LAI_sum(1:kproma,krow) + prt_lai(1:kproma,krow)

    IF(first_day_lai) THEN !Put a value for LAI first day
      LAI_daymean(1:kproma,krow) = prt_lai(1:kproma,krow)
      LAI_sum_month(1:kproma,krow) = LAI_daymean(1:kproma,krow)
      LAI_monthmean(1:kproma,krow) = LAI_daymean(1:kproma,krow)
      LAI_currentmonth(1:kproma,krow) = LAI_daymean(1:kproma,krow)
    ENDIF

  ELSE

    LAI_daymean(1:kproma,krow) = LAI_sum(1:kproma,krow)*qnsteps

    LAI_sum(1:kproma,krow) = prt_lai(1:kproma,krow)

    !monthly mean LAI
    LAI_sum_month(1:kproma,krow) = LAI_sum_month(1:kproma,krow) + LAI_daymean(1:kproma,krow)

    CALL get_date_components(current_date, month=month_lai)
    CALL get_date_components(current_date, DAY=day_of_current_month)

    LAI_currentmonth(1:kproma,krow) = LAI_sum_month(1:kproma,krow)/day_of_current_month

    IF (.NOT. new_year .AND. month_lai < 2) THEN ! Put a value for first month
      LAI_monthmean(1:kproma,krow) = LAI_currentmonth(1:kproma,krow)
    ENDIF

    ! add LAI_prevmonth=LAI_monthmean before if new_month
    ! add if first month then LAI_monthmean=LAI_prevmonth=LAI_sum/number of day of
    ! current month???
    IF (new_month) THEN

      CALL get_date_components(previous_date, DAY=day_of_month)
      CALL get_date_components(current_date, DAY=day_of_current_month)
      LAI_monthmean(1:kproma,krow) = LAI_sum_month(1:kproma,krow)/day_of_month

      IF (day_of_current_month < 2) THEN
        LAI_currentmonth(1:kproma,krow) = LAI_daymean(1:kproma,krow)
      ELSE
        LAI_currentmonth(1:kproma,krow) = LAI_sum_month(1:kproma,krow)/day_of_current_month
      ENDIF

      LAI_sum_month(1:kproma,krow) = 0._dp

    ENDIF !new month

  END IF !new_day
  END SUBROUTINE biogenic_averages

!  
!-------------------------------------end biogenic_averages--------------------------------------------
! New subroutine for interactive computation of emission factor from JSBACH PFTs from Colombe's work
    
  SUBROUTINE calc_emissions_factors(kproma,kbdim,krow,pfrpft_jsbach)

    ! Called by calc_biogenic_emissions 
    ! Mapping of JSBACH PFTs into MEGANv2.1/CLM4 PFTs
    ! Gives the corresponding fractions of PFTs for NPFT


    USE mo_boundary_condition,     ONLY: bc_find, bc_apply,bc_find  
    USE mo_exception,              ONLY: finish
    USE mo_time_control,           ONLY: lstart, lresume
    
    INTEGER,  INTENT(in)  :: kproma                     ! geographic block number of locations
    INTEGER,  INTENT(in)  :: kbdim                      ! geographic block maximum number of locations 
    INTEGER,  INTENT(in)  :: krow                       ! geographic block number
    REAL(dp), INTENT(out) :: pfrpft_jsbach(kbdim,NPFT)  ! Fraction PFT JSBACH
    
    ! local variables

    REAL(dp)         :: zfrpft_jsbach(kbdim,npft_red_jsbach)
    REAL(dp)         :: zfr_brdle_ev(kbdim),zfr_brdle_de(kbdim)
    REAL(dp)         :: zborealtr_limit(kbdim),zborealgr_limit(kbdim)
    REAL(dp)         :: zfrpft_jsbach_sep(kbdim,NPFT)   !mapping of JSBACH PFTs to MEGAN-CLM4 PFTs
    REAL(dp)         :: ztcmax(kbdim)   
    INTEGER          :: i,j,ij,ierr
    CHARACTER(len=2) :: tilnum
   

    !initialisation
    pfrpft_jsbach(:,:) = 0.0_dp
    zfrpft_jsbach_sep(:,:) = 0.0_dp

    ! BL/NL fractions from potential vegetation maps
    zfr_brdle_de(1:kproma)=vegfract_de(1:kproma,krow,1,1)
    zfr_brdle_ev(1:kproma)=vegfract_ev(1:kproma,krow,1,1)
    ! BOREAL/TEMPERATE separation from maps
    zborealtr_limit(1:kproma)=borealtr_limit(1:kproma,krow,1)
    zborealgr_limit(1:kproma)=borealgr_limit(1:kproma,krow,1)

    ! get the fraction of each vegetation class (as defined by JSBACH) with the help of boundary conditions
    DO i = 1, npft_red_jsbach  
       IF (lstart .OR. lresume) THEN
          write(tilnum,'(i2)') i
          CALL bc_find('Tile number '//TRIM(tilnum), ibc_frpft(i), ierr=ierr)
          IF (ierr /= 0 .OR. ibc_frpft(i)<= 0) CALL finish('init_emi_biogenic_stream',         &
             'Cannot locate boundary condition for fraction coverage land cover type'//char(i)) 
       ENDIF
       CALL bc_apply(ibc_frpft(i), kproma, krow,  zfrpft_jsbach(:,i))
    END DO

   ! Mapping of JSBACH PFTs into MEGAN/CLM4 PFTs to calculate BVOC emissions from JSBACH PFT fractions and emission factors  
   ! using JSBACH 11 PFTs(tiles) scheme and MEGAN/CLM4 16 PFTs
 
   ! Tropical PFTs (1 and 2 in JSBACH to 4 and 6 in MEGAN/CLM4)
   zfrpft_jsbach_sep(:,4) = zfrpft_jsbach(:,1) !Tropical BL EV
   zfrpft_jsbach_sep(:,6) = zfrpft_jsbach(:,2) !Tropical BL DE

   ! Extra-tropical PFTs ( 3 and 4 in JSBACH to 1,2,3,5,7,8 in MEGAN/CLM4)
   ! BOREAL/TEMPERATE separation: based on the maximum limit (-2 deg. Celcius) 
   ! of the minimum coldest month temperature from a 20-years run (1990-2010) with climate_buffer
   ! BROADLEAF/NEEDLELEAF separation: based on potential vegetation maps ( based on the work of J. Pongratz et al., 2008)

   ! Extra-tropcial EV
   DO ij = 1, kproma 
     IF (zborealtr_limit(ij) < 0.5_dp) THEN
       zfrpft_jsbach_sep(ij,1) = zfrpft_jsbach(ij,3)*(1._dp-zfr_brdle_ev(ij))  !Temperate NL EV
       zfrpft_jsbach_sep(ij,5) = zfrpft_jsbach(ij,3)*zfr_brdle_ev(ij)          !Temperate BL EV
     ELSE
       zfrpft_jsbach_sep(ij,3) = zfrpft_jsbach(ij,3)*(1._dp-zfr_brdle_ev(ij))  !Boreal NL EV ! Attention PFT 2 and 3 swapped in MEGAN data
                                                                               !Normally should be PFT 2 in MEGAN/CLM4 classification
     END IF
   END DO

   ! Extra-tropical DE
   DO ij = 1, kproma 
     IF (zborealtr_limit(ij) < 0.5_dp) THEN
       zfrpft_jsbach_sep(ij,7) = zfrpft_jsbach(ij,4)*zfr_brdle_de(ij)         !Temperate BL DE
     ELSE 
       zfrpft_jsbach_sep(ij,2) = zfrpft_jsbach(ij,4)*(1._dp-zfr_brdle_de(ij)) !Boreal NL DE ! Attention PFT2 and 3 swapped in MEGAN data
                                                                              !Normally should be PFT 3 in MEGAN/CLM4 classification
       zfrpft_jsbach_sep(ij,8) = zfrpft_jsbach(ij,4)*zfr_brdle_de(ij)         !Boreal BL DE
     END IF
   END DO

   ! Shrubs (treated like trees, same temperature threshold for BOREAL/TEMPERATE separation)
   DO ij = 1, kproma 
     zfrpft_jsbach_sep(:,9) = 0.0_dp                  !Temperate EV shrubs
     zfrpft_jsbach_sep(ij,10) = zfrpft_jsbach(ij,5)   !Temperate DE shrubs
   !ELSE
     zfrpft_jsbach_sep(ij,11) = zfrpft_jsbach(ij,6)   !Boreal DE shrubs
   !END IF
   END DO

   ! Grasses 
   ! Arctic/Cool C3 separation: based on the maximum limit (-17 deg. Celcius, from CLM4 Tech. Note, 2010, Table 16.1, p.232) 
   ! of the minimum coldest month temperature from a 20-years run (1990-2010) with climate_buffer
   ! Clustering of grasses and pasture from JSBACH PFTs (C3 grasses + C3 pastures, C4 grasses + C4 pastures)
   DO ij = 1, kproma 
     IF (zborealgr_limit(ij) < 0.5_dp) THEN
       zfrpft_jsbach_sep(ij,13) = zfrpft_jsbach(ij,7)+zfrpft_jsbach(ij,9) ! C3 cool grasses
     ELSE
       zfrpft_jsbach_sep(ij,12) = zfrpft_jsbach(ij,7)+zfrpft_jsbach(ij,9) ! C3 arctic grasses
     END IF
   END DO
   zfrpft_jsbach_sep(:,14) = zfrpft_jsbach(:,8)+zfrpft_jsbach(:,10)   ! C4 grasses

   ! Crops
   zfrpft_jsbach_sep(:,15) = zfrpft_jsbach(:,11)           !crops
   zfrpft_jsbach_sep(:,16) = 0._dp                         !currently not used

   ! Corrections of error values???    
   DO ij=1,kproma
     DO i =1,NPFT
       pfrpft_jsbach(ij,i)=zfrpft_jsbach_sep(ij,i)
       IF(pfrpft_jsbach(ij,i)<=0._dp) pfrpft_jsbach(ij,i)=0._dp
       IF(pfrpft_jsbach(ij,i)>100._dp) pfrpft_jsbach(ij,i)=0._dp
     END DO
   END DO

  END SUBROUTINE calc_emissions_factors

END MODULE mo_hammoz_emi_biogenic
