!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham.f90
!!
!! \brief
!! mo_ham contains physical switches and other control parameters for 
!! ECHAM aerosol models (particularly HAM).
!!
!! \author P. Stier (MPI-Met)
!!
!! \responsible_coder
!! Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# P. Stier (MPI-Met) - original version - (2002-12-xx) 
!!   -# K. Zhang (MPI-Met) - changes for ECHAM6 - (2009-08-11)
!!   -# M. Schultz (FZ Juelich) - cleanup (2009-09-23)
!!   -# H. Kokkola (FMI) - definition of aerocomp (former mo_ham mode) (2011-12-12)
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

MODULE mo_ham

  ! *mo_ham* contains physical switches and parameters 
  !           for the ECHAM-HAM aerosol model.
  !
  ! Author:
  ! -------
  ! Philip Stier, MPI-MET                    12/2002
  !

  USE mo_kind,                ONLY: dp
  USE mo_species,             ONLY: t_species, nmaxspec
  USE mo_external_field_processor, ONLY: EF_FILE, EF_MODULE !SF #244

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setham !TB initialise ham
  PUBLIC :: new_aerocomp ! create aerosol component
  PUBLIC :: print_aerocomp_info
  PUBLIC :: nham_subm ! choice for aerosol microphysics submodel
  PUBLIC :: naeroclass, naerosol ! arrays of possible number of classes (resp. soluble classes)
  PUBLIC :: nclass ! number of classes in current model
  PUBLIC :: nsol ! number of soluble classes in current model
  PUBLIC :: naerocomp, naerorad, nraddiag, nrad, nradmix, nseasalt, nmaxclass
!>>NAJS: optical properties of dry aerosol
  PUBLIC :: nraddry
!<<NAJS
!>>NAJS: lidar backscatter
  PUBLIC :: nradbeta
!<<NAJS
  PUBLIC :: HAM_BULK, HAM_M7, HAM_SALSA
  PUBLIC :: npist, laerocom_diag, ndrydep, nwetdep, lomassfix
  PUBLIC :: ndust !SF #479
  PUBLIC :: burden_keytype
  PUBLIC :: mw_s, mw_so2, mw_so4, mw_dms, mw_oc
  PUBLIC :: lhetfreeze
  PUBLIC :: nlai_drydep_ef_type !gf #244
  PUBLIC :: lscond ! Condensation of H2SO4
  PUBLIC :: lscoag ! Coagulation
  PUBLIC :: lgcr ! galactic cosmic ray ionization
  PUBLIC :: nsolact ! solar activity parameter
  PUBLIC :: lmass_diag ! Mass diagnostics switch
  PUBLIC :: nccndiag ! (C)CN diagnostics at fixed supersaturations
  PUBLIC :: ibc_dust, ibc_seasalt, idsec_dust, idsec_biogenic !alaak moved from mo_ham_m7/salsa_emissions
  PUBLIC :: sigma_fine, sigma_coarse !SF #320

  INTEGER, PARAMETER     :: nmaxclass = 20    ! maximum number of aerosol modes or bins
  ! flags for aerosol microphysics scheme
  INTEGER, PARAMETER     :: HAM_BULK  = 1
  INTEGER, PARAMETER     :: HAM_M7    = 2
  INTEGER, PARAMETER     :: HAM_SALSA = 3

  !--- derived types
  TYPE, PUBLIC :: t_aeroclass
     INTEGER                  :: iclass          ! Aerosol mode or bin number
     TYPE(t_species), POINTER :: species        ! Aerosol species 
     INTEGER                  :: spid           ! Index in species list
     INTEGER                  :: aero_idx       ! Index to aerosol species list 
     INTEGER                  :: tracer_type    ! Tracer type (none/diagnostic/prognostic)
     INTEGER                  :: idt            ! Tracer identity
  END TYPE t_aeroclass

  TYPE, PUBLIC :: t_sizeclass 
     CHARACTER (LEN=32) :: classname         ! Long mode name, e.g. "nucleation soluble"
     CHARACTER (LEN=4)  :: shortname         ! Short mode name, e.g. "NS"
     INTEGER            :: self              ! =mode index, for quick comparisons, etc
     LOGICAL            :: lsoluble          ! Mode soluble (T/F)
     LOGICAL            :: lactivation       ! Mode can activate to cloud droplets !dod #377
     LOGICAL            :: lsed              ! Sedimentation for this mode (T/F)
     LOGICAL            :: lsoainclass       ! Secondary organics occur in this mode (T/F)
     INTEGER            :: idt_no            ! Tracer identity for aerosol number
  END TYPE t_sizeclass                       ! sigma, sigmaln and the conversion factors are
                                             ! kept separate to avoid too many impacts on 
                                             ! other modules and subroutines
  INTEGER, SAVE          :: nclass=7           ! number of aerosol modes or size bins (e.g. 7 for M7)
  INTEGER, SAVE          :: nsol=4           ! number of aerosol modes or size bins (e.g. 4 for M7)
  INTEGER                :: naerocomp          ! number of aerocomps defined (see mo_ham_init)

  TYPE(t_sizeclass), TARGET, PUBLIC :: sizeclass(nmaxclass)
  !--- aerocomp: linear list of class*species

  TYPE(t_aeroclass), PUBLIC, ALLOCATABLE :: aerocomp(:)
  TYPE(t_aeroclass), PUBLIC, ALLOCATABLE :: aerowater(:)

  INTEGER, PUBLIC :: subm_ngasspec  = 0
  INTEGER, PUBLIC :: subm_naerospec = 0
  INTEGER, PUBLIC :: subm_gasspec(nmaxspec)      ! gas species indices for microphysical processes
  INTEGER, PUBLIC :: subm_aerospec(nmaxspec)     ! aero species indices for microphysical processes
  !>>SF for convenience: special mapping for all but water aero species
  INTEGER, PUBLIC :: subm_naerospec_nowat = 0
  INTEGER, PUBLIC :: subm_aerospec_nowat(nmaxspec) ! all but water aero species indices for microphysical processes
  !<<SF
  INTEGER, PUBLIC :: subm_aero_idx(nmaxspec)     ! mapping from speclist to subm_aerospec
  INTEGER, PUBLIC :: subm_gasunitconv(nmaxspec)  ! unit conversion flag for subm_gasspec
  INTEGER, PUBLIC :: subm_aerounitconv(nmaxspec) ! unit conversion flag for subm_aerospec

  PUBLIC :: immr2ug, immr2molec, ivmr2molec  ! unit conversion for some species

  INTEGER, PARAMETER     :: immr2ug    = 1        ! Mass mixing ratio to ug m-3
  INTEGER, PARAMETER     :: immr2molec = 2        ! Mass mixing ratio to molecules cm-3
  INTEGER, PARAMETER     :: ivmr2molec = 3        ! Volume mixing ratio to molecules cm-3

  !--- 1) Switches:

  !--- 1.0) Logical:

  INTEGER :: nham_subm = 2             ! Switch for aerosol microphysics scheme:
                                       !
                                       ! nham_subm = 1  Bulk scheme
                                       !           = 2  Modal scheme (M7) (default)
                                       !           = 3  Sectional scheme (SALSA)

  !--- 1.1) Physical:

  !--- Define control variables and pre-set with default values: 

  INTEGER :: nseasalt    = 2           ! Sea Salt emission scheme:
                                       ! 
                                       !    nseasalt = 1  Monahan (1986)
                                       !             = 2  Schulz et al. (2002)
                                       !             = 3  Reserved (Martensson)
                                       !             = 4  Monahan (1986), bin scheme
                                       !             = 5  Guelle (2001)
                                       !             = 6  Gong (2003)
                                       !             = 7  Long (2011) (SST dep)
                                       !             = 8  Gong + SST dependence 
                                       !
  INTEGER :: npist       = 3           ! DMS emission scheme:
                                       !
                                       !    npist = 1 Liss & Merlivat (1986) 
                                       !          = 2 Wanninkhof (1992)
                                       !          = 3 Nightingale (2000)
                                       !
  INTEGER :: naerorad    = 1           ! HAM aerosols are radiatively active
                                       !    
                                       !    naerorad = 0 HAM aerosol radiation deactivated
                                       !                 (requires iaero/=1)
                                       !             = 1 HAM aerosol radiation prognostic
                                       !             = 2 HAM aerosol radiation diagnostic only
                                       !
  INTEGER :: ndrydep     = 2           ! dry deposition scheme (default = interactive (Ganzeveld))
                                       !
  INTEGER :: nwetdep     = 1           ! default wetdep scheme (?????)
                                       ! this flag is a simple, user-friendly interface
                                       ! which drives a more detailed scavenging setup
                                       ! see details in mo_wetdep_interface and
                                       ! related routines.
                                       ! nwetdep = 0 wetdep (scavenging) off
                                       !           1 standard mode-wise prescribed scavenging parameters
                                       !           2 standard in-cloud scav + aerosol size-dep below-cloud scav 
                                       !           3 size-dep in-cloud and below-cloud scav
                                       !             WARNING: !!!size-dep IC not yet implemented!!!
  INTEGER :: ndust = 5                 ! choice of parameter set for the BGC dust scheme !SF #479
                                       ! 2: Cheng et al. (2008)
                                       ! 3: Stier et al. (2005) default
                                       ! 4: Stier et al. (2005) + East Asia soil properties
                                       ! 5: Stier et al. (2005) + MSG-based Saharan dust sources
                                       !    (Schepanski et al., GRL 2007; RSE 2012) + East Asia soil properties
                                       !    (!BH #382)
                                       !SF ToDo: shift of ndust values / cleanup needed: 
                                       !         currently ndust=[0,1] not implemented
                                       !
  LOGICAL :: laerocom_diag = .FALSE.   ! Extended diagnostics
                                       !
  LOGICAL :: lhetfreeze  = .FALSE.     ! switch to set heterogeneous freezing below 235K (cirrus scheme)
                  
!>>gf #244
  INTEGER :: nlai_drydep_ef_type = EF_MODULE  ! Choice of lai external field type in the drydep scheme
                                              !   = EF_FILE (2) from external input file
                                              !   = EF_MODULE (3) online from jsbach
!<<gf
  LOGICAL :: lscond     = .TRUE.    ! Condensation of H2SO4
  LOGICAL :: lscoag     = .TRUE.    ! Coagulation
  

  LOGICAL :: lgcr       = .TRUE.    ! Calculate ionization due to galactic cosmic rays
  
  REAL(dp):: nsolact    = -99.99_dp ! Solar activity parameter [-1,1]; if outside of
                                    ! this range (as per default), then the model will
                                    ! determine the solar activity based on the model
                                    ! calendar date; otherwise, it will use the user
                                    ! set solar activity parameter throughout the run.
                                    ! -1 is solar minimum, 1 solar maximum.
  
  LOGICAL :: lmass_diag = .FALSE.   ! Mass balance check in m7_interface

  INTEGER :: nccndiag = 0           ! (C)CN diagnostics at fixed supersaturations
  


!++mgs: changed dimension and initialisation of nrad and nradmix to allocatable to avoid fixed dependency on nmod
!  unfortunately, ALLOCATABLE is not possible because of namelist
  INTEGER :: nrad(nmaxclass)                      !
!!  INTEGER :: nrad(nmod)    = (/ 0, &   ! Radiation calculation (for each mode)
!!                                1, &   !
!!                                1, &   !    nrad = 0 NO    radiation calculation
!!                                1, &   !         = 1 SW    radiation calculation
!!                                1, &   !         = 2 LW    radiation calculation
!!                                1, &   !         = 3 SW+LW radiation calculation
!!                                1  /)  !

  INTEGER :: nradmix(nmaxclass)
!!  INTEGER :: nradmix(nmod) = (/ 1, &   ! Mixing scheme for refractive indices
!!                                1, &   ! (for each mode)
!!                                1, &   !
!!                                1, &   !    nradmix = 1 volume weighted mixing
!!                                1, &   !            = 2 Maxwell-Garnet mixing
!!                                1, &   !            = 3 Bruggeman mixing
!!                                1  /)  !
                                       !
  INTEGER :: nraddiag    = 1           ! Extended radiation diagnostics
                                       !
                                       !    nraddiag = 0 off
                                       !             = 1 2D diagnostics
                                       !             = 2 2D+3D diagnostics
                                       !
!>>NAJS: optical properties of dry aerosol
  INTEGER :: nraddry     = 0           ! 0: no dry diagnostics
                                       ! 1: only dry AOT
                                       ! 2: dry AOT and AAOT
                                       ! 3: dry AOt, AAOT and extinction
!<<NAJS
!>>NAJS: lidar backscatter
  INTEGER :: nradbeta    = 0           ! 0: no LIDAR diagnostics
                                       ! 1: basic LIDAR diagnostics: alfa_aer (aerosol extinction)
                                       ! 2: (1) + aerosol & gas backscatter & lidar ratio
                                       ! 3: (2) + attenuated backscatter from space & Earth
!<<NAJS
  LOGICAL :: lomassfix   = .TRUE.      ! Mass fixer in convective scheme

  INTEGER, PUBLIC :: nsoa = 0          ! Choice for the secondary organics scheme 
                                       !! 0: no SOA scheme
                                       !! 1: SOA scheme from O'Donnell et all, ACP 2011
                                       !! 2: SOA scheme with VBS approach from Farina et al, JGR 2010
                                       !!    (curr. SALSA only)
  INTEGER, PUBLIC :: nsoaspec          ! number of SOA species 

  INTEGER, PUBLIC :: nsoalumping = 0   ! SOA lumping scheme
                                       ! 0: no lumping (DEFAULT: WILL BE CHANGED)
                                       ! 1: lump anthropogenic non-volatile SOA
                                       ! 2: lump anthropogenic non-volatile SOA and map onto OC
                                       ! 3: lump anthropogenic non-volatile SOA and map onto OC
                                       !    and lump all anthropogenic precursors

  !--- 1.1) output control...>>dod deleted <<dod

  ! Output of optional diagnostic fields
  ! ------------------------------------
  ! key_types (see mo_submodel_diag.f90)
  !         0 = no output (in which case no memory is allocated(?))
  !         1 = by tracer
  !         2 = by species
  !         3 = by mode
  ! Not all switches can be set for all key_types!
  !
  ! Gas phase species are always output by tracer
  ! Aerosol numbers are always output per mode, with the obvious exception of 'no output'.

  INTEGER         :: burden_keytype          ! options: OFF, BYTRACER, BYSPECIES
  
  !>>dod deleted iwritetrac <<dod

  !--- 1.2) boundary conditions 
  INTEGER, PUBLIC :: ibc_oh, ibc_o3, & ! boundary condition indices for tracer fields
                     ibc_h2o2, ibc_no2, ibc_no3


  !--- 2) Parameters:

  INTEGER, PARAMETER :: iaeroham = 1      ! iaero switch for submodel aerosol (mo_radiation_parameters)

  INTEGER, PARAMETER :: ntype=6
  CHARACTER(LEN=3), PARAMETER :: ctype(ntype)=(/'SO4','BC ','OC ','SS ','DU ','WAT'/)
  !<<dod

  !--- Molecular weight of SO2 and SO4

  INTEGER             :: ibc_dust = 0, ibc_seasalt = 0 !alaak moved from mo_ham_m7/salsa_emissions
  INTEGER             :: idsec_dust, idsec_biogenic

  INTEGER, PARAMETER  :: naeroclass(3) = (/1, 7, 17/)  ! number of size classes (bulk, M7, SALSA)

  INTEGER, PARAMETER  :: naerosol(3) = (/1, 4, 10/)    ! number of soluble compounds alaak: nsol-> nsol(3)
 
  REAL(dp), PARAMETER :: mw_s   = 32.0655_dp, & ! molecular weight S
                         mw_so2 = 64.0643_dp, & ! molecular weight SO2
                         mw_so4 = 96.0631_dp, & ! molecular weight SO4
                         mw_dms = 62.1345_dp, & ! molecular weight DMS
                         mw_oc  = 180._dp       ! molecular weight OC 
  !--- Ratio of sulfate mass to mass of sulfur:
  !    Depends on the assumed degree of neutralization.
  !    >>dod deleted. 

  !>>SF #320: these are std's for fine and coarse modes. This is placed here and not in M7-only code,
  !           because it's needed outside of M7-specific code (see radiation, sea salt emisions...)
  REAL(dp), PARAMETER :: sigma_fine   = 1.59_dp ! std for fine modes
  REAL(dp), PARAMETER :: sigma_coarse = 2._dp   ! std for coarse modes
  !<<SF #320

CONTAINS 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! setham modifies pre-set switches of the hamctl namelist for the 
!! configuration of the ECHAM/HAM aerosol model
!!  
!! @author see above
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! init_ham
!!
!! @par Externals:
!! <ol>
!! <li>None
!! </ol>
!!
!! @par Notes
!! 
!!
!! @par Responsible coder
!! kai.zhang@zmaw.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  SUBROUTINE setham

    ! *setham* modifies pre-set switches of the hamctl
    !           namelist for the configuration of the 
    !           ECHAM/HAM aerosol model
    ! 
    ! Authors:
    ! --------
    ! Philip Stier, MPI-MET                        12/2002
    ! Martin Schultz, FZ Juelich, Oct 2009, added nmod parameter
    !

    USE mo_mpi,                 ONLY: p_parallel, p_parallel_io, p_bcast, p_io
    USE mo_namelist,            ONLY: open_nml, position_nml, POSITIONED, MISSING
    USE mo_radiation_parameters,ONLY: iaero !SF obsolete to remove , l_srtm, l_lrtm
    USE mo_exception,           ONLY: finish, message, message_text, em_info, em_error, em_param,  &
                                      em_info, em_warn
    USE mo_util_string,         ONLY: separator
    USE mo_submodel,            ONLY: print_value
    USE mo_submodel_diag,       ONLY: OFF, ON, BYTRACER, BYSPECIES
    USE mo_param_switches,      ONLY: nic_cirrus, ncd_activ
    USE mo_control,             ONLY: nn, nlev
    USE mo_boundary_condition,  ONLY: bc_nml, bc_define, p_bcast_bc,  &
                                      BC_REPLACE, BC_EVERYWHERE
    ! other available options (bc_domain=BC_EVERYWHERE): BC_ALTITUDE, BC_PRESSURE, BC_LEVEL, BC_TOP, BC_BOTTOM
    !                         (bc_mode=BC_REPLACE): BC_SPECIAL, BC_RELAX, BC_ADD
    USE mo_external_field_processor, ONLY: EF_FILE, &
                                           EF_3D, EF_IGNOREYEAR, EF_NOINTER
    ! other available options: EF_CONSTANT, EF_TIMERESOLVED, EF_LATLEV, EF_LONLAT, EF_UNDEFINED
    USE mo_submodel,            ONLY: lham, lhammoz
    USE mo_physical_constants,  ONLY: amd
    USE mo_util_string,         ONLY: separator

    IMPLICIT NONE

    TYPE(bc_nml)     ::  bc_oh, bc_o3, bc_h2o2, bc_no2, &   ! tracer boundary conditions/input specs
                         bc_no3                             ! SOA oxidant

    INCLUDE 'hamctl.inc'

!    INTEGER       :: nclass         ! number of modes/bins

    !--- Local variables

    CHARACTER(len=24)         :: csubmname    ! name of aerosol sub model
    INTEGER                   :: jj, ierr, inml, iunit


    !--- 0) Set defaults
    ibc_oh    = 0
    ibc_o3    = 0
    ibc_h2o2  = 0
    ibc_no2   = 0
    ibc_no3   = 0

    burden_keytype    = BYSPECIES
    
    !--- Set defaults for tracer boundary conditions
    ! If lhammoz is active, default is to use tracer fields from gas-phase chemistry module
    ! else the fields are read from file
    ! define entries for first tracer, then copy
    ! explanation: %bc_... = how to apply
    !              %ef_... = how to get values
    bc_oh%bc_domain = BC_EVERYWHERE
    bc_oh%bc_mode = BC_REPLACE
    bc_oh%ef_type = EF_FILE
    bc_oh%ef_template = 'ham_oxidants_monthly.nc' 
    bc_oh%ef_varname = 'OH_VMR_avrg'
    bc_oh%ef_geometry = EF_3D
    bc_oh%ef_timedef = EF_IGNOREYEAR
    !! additional options (just for illustration) - see mo_boundary_condition
    !! bc_oh%ef_timeoffset = 0._dp
    !! bc_oh%ef_timeindex = 1
    !! bc_oh%ef_value = 0._dp
    bc_oh%ef_factor = 17._dp/amd        ! convert from VMR to MMR
    bc_oh%ef_interpolate = EF_NOINTER   ! none
    bc_oh%ef_actual_unit = 'VMR'

    bc_o3   = bc_oh
    bc_o3%ef_varname   = 'O3_VMR_avrg'
    bc_o3%ef_factor = 48._dp/amd        ! convert from VMR to MMR
    bc_h2o2 = bc_oh
    bc_h2o2%ef_varname = 'H2O2_VMR_avrg'
    bc_h2o2%ef_factor = 34._dp/amd      ! convert from VMR to MMR
    bc_no2  = bc_oh
    bc_no2%ef_varname  = 'NO2_VMR_avrg'
    bc_no2%ef_factor = 46._dp/amd       ! convert from VMR to MMR
    bc_no3  = bc_oh
    bc_no3%ef_varname  = 'NO3_VMR_avrg'
    bc_no3%ef_factor = 62._dp/amd       ! convert from VMR to MMR

    nrad(1)      = 0
    nrad(2:nmaxclass) = 3
    nradmix(:)   = 1

    subm_gasspec(:)  = 0
    subm_aerospec(:) = 0
    subm_aero_idx(:) = 0
 
    !--- 1) Read namelist:

    CALL message('',separator)
    CALL message('setham', 'Reading namelist hamctl...', level=em_info)

    IF (p_parallel_io) THEN
       inml = open_nml('namelist.echam')
       iunit = position_nml ('HAMCTL', inml, status=ierr)
       SELECT CASE (ierr)
       CASE (POSITIONED)
          READ (iunit, hamctl)
       CASE (MISSING)
          WRITE(message_text,'(a,i0,a)') 'Namelist HAMCTL not found. Will use default values'
          CALL message('setham', message_text, level=em_warn)
       CASE DEFAULT ! LENGTH_ERROR or READ_ERROR
          WRITE(message_text,'(a,i0)') 'Namelist HAMCTL not correctly read! ierr = ', ierr
          CALL finish('setham', message_text)
       END SELECT
    ENDIF

    !--- 2) Broadcast over processors:
    IF (p_parallel) THEN
       CALL p_bcast (nham_subm,     p_io)
       CALL p_bcast (nseasalt,      p_io)
       CALL p_bcast (npist,         p_io)
       CALL p_bcast (naerorad,      p_io)       
       CALL p_bcast (ndrydep,       p_io)
       CALL p_bcast (nwetdep,       p_io)
       CALL p_bcast (ndust,         p_io) !SF #479
       CALL p_bcast (laerocom_diag, p_io)
       CALL p_bcast (nrad,          p_io)
       CALL p_bcast (nradmix,       p_io)
       CALL p_bcast (nraddiag,      p_io)
!>>NAJS: optical properties of dry aerosol
       CALL p_bcast (nraddry,       p_io)
!<<NAJS
!>>NAJS: lidar backscatter
       CALL p_bcast (nradbeta,      p_io)
!<<NAJS
       CALL p_bcast (lomassfix,     p_io)
       CALL p_bcast (nsoa,          p_io)                     
       CALL p_bcast (nsoalumping,   p_io) 
       CALL p_bcast (lhetfreeze,    p_io)
       CALL p_bcast (burden_keytype,p_io)
       CALL p_bcast_bc (bc_oh,      p_io)
       CALL p_bcast_bc (bc_o3,      p_io)
       CALL p_bcast_bc (bc_h2o2,    p_io)
       CALL p_bcast_bc (bc_no2,     p_io)
       CALL p_bcast_bc (bc_no3,     p_io)
       CALL p_bcast(nlai_drydep_ef_type, p_io) !gf #244
       CALL p_bcast (lscond,     p_io)
       CALL p_bcast (lscoag,     p_io)
       CALL p_bcast (lgcr,       p_io)
       CALL p_bcast (nsolact,    p_io)
       CALL p_bcast (lmass_diag, p_io)
       CALL p_bcast (nccndiag,   p_io)
  
    END IF

    !---------------------------------------------------------------------------------------------------
    !--- 3) Consistency and dependency checks:

!### normally this test for lham is not necessary, but perhaps needed when we support additional modules
    IF (lham) THEN

      !--- consistency check for microphysics scheme
      IF (nham_subm < 1 .OR. nham_subm > 3) THEN
         WRITE(message_text,'(a,i0)') 'Illegal option for aerosol microphysics scheme nham_subm = ', &
              nham_subm
         CALL message('setham',message_text, level=em_error)
      END IF

      IF (nham_subm == 1) THEN
         WRITE(message_text,'(a,i0,a)') 'Bulk microphysics nham_subm = ',nham_subm,'has not been implemented'
         CALL message('setham',message_text, level=em_error)
      END IF

      IF (nham_subm == HAM_SALSA) THEN
         IF (nwetdep > 1) THEN
            WRITE(message_text,'(a,i0)') &
                 'SALSA microphysics does not include wet deposition option nwetdep = ', nwetdep
            CALL message('setham',message_text, level=em_error)
         END IF

         IF (nccndiag > 0) THEN
            WRITE(message_text,'(a,i0)') &
                 'SALSA microphysics does not support detailed CCN diagnostics nccndiag = ', nccndiag
            CALL message('setham',message_text, level=em_error)
         END IF

      ENDIF

!>>SF #299
      !--- Check authorized values for burden_keytype
      SELECT CASE(burden_keytype)
         CASE(OFF,BYTRACER,BYSPECIES)
             CONTINUE
         CASE default
             WRITE(message_text,'(a,i0,a,i0,a,i0,a,i0)') 'invalid value for burden_keytype: ',&
                   burden_keytype,'. Authorized values: ',OFF,' ',BYTRACER,' ',BYSPECIES
             CALL finish('setham',message_text)
      END SELECT
!<<SF

!>>SF #479
      !--- Check authorized values for ndust
      IF (ndust < 2 .OR. ndust > 5) THEN
         WRITE(message_text,'(a,i0)') 'Illegal option for the parameter set of the BGC dust scheme  ndust = ', &
              ndust
         CALL message('setham',message_text, level=em_error)
      ENDIF
!<<SF

!>>gf #244
      !--- consistency check for lai in the drydep scheme !SF + user information
      SELECT CASE(nlai_drydep_ef_type)
          CASE(EF_FILE)
             WRITE(message_text,'(a)') 'Leaf area index read from file'
          CASE(EF_MODULE)
             WRITE(message_text,'(a)') 'Leaf area index from JSBACH'
          CASE default
           WRITE(message_text,'(a,i0,2a,2(i0,a))') 'nlai_drydep_ef_type = ',nlai_drydep_ef_type, &
                                  ' --> this is not currently supported.', &
                                  ' Only ',EF_FILE,' (from file) or ',EF_MODULE, &
                                  ' (from another module) external field types are possible!'
           CALL finish('setham',message_text)
      END SELECT
      WRITE(message_text,'(a,a)') TRIM(message_text),' for usage by the drydep scheme.'
      CALL message('setham',message_text,level=em_info)
!<<gf
      nclass=naeroclass(nham_subm)
      nsol = naerosol(nham_subm)
      IF (nclass > nmaxclass) CALL finish('setham', 'Maximum number of aerosol modes/bins exceeded!')


      !--- consistency checks for radiation
      IF(naerorad>0) THEN
        IF (ANY(nrad(:)>3) .OR. ANY(nrad(:)<0) ) THEN
          CALL message('setham','nrad>3 or nrad<0 not supported', level=em_warn)
          CALL message('', 'Will reset all nrad and nradmix values to default values!', level=em_info)
          nrad(1)         = 0
          nrad(2:nclass)    = 1
          nradmix(1:nclass) = 1
        END IF
        IF (iaero/=iaeroham .AND. ANY(nrad(:)>0) .AND. naerorad/=2) THEN
          CALL message('setham','inconsistent setting of iaero, nrad and naerorad', level=em_warn)
          CALL message('', 'Will run ECHAM with selected iaero and switch HAM radiation to diagnostic mode', level=em_info)
          naerorad = 2
        END IF
        IF (iaero==iaeroham .AND. .NOT. ANY(nrad(:)>0)) THEN
          CALL message('setham','HAM radiation requested but all nrad set to zero.', level=em_warn)
          CALL message('', 'Will reset all nrad and nradmix values to default values!', level=em_info)
          nrad(1)         = 0
          nrad(2:nclass)    = 1
          nradmix(1:nclass) = 1
        ENDIF

      END IF

      !--- consistency checks for SOA scheme choices
      IF ((nham_subm .NE. HAM_SALSA) .AND. nsoa == 2) THEN

         CALL message('setham','nsoa = 2 currently only supported with SALSA (nham_subm = 3)', &
                      level=em_error)

      ENDIF

!>>SF override the following check with -DWITH_LHET in your Makefile
!     note: lhet-dependent code has not be fully tested yet
!           (implemented here to avoid losing existing code in earlier versions)
#ifndef WITH_LHET
      IF (lhetfreeze) CALL message('setham','lhetfreeze=.TRUE. is not currently supported', level=em_error)
#endif
!<<SF

    ELSE      
       ! reset to default values for lham=.FALSE.
       laerocom_diag = .FALSE.
       lomassfix     = .FALSE.
       nrad          = 0
       nraddiag      = 0
       lhetfreeze    = .FALSE.
       nham_subm     = 0
    END IF

    !---------------------------------------------------------------------------------------------------
    !--- define the boundary conditions based on the default values and namelist input
    ibc_oh   = bc_define('OH mass mixing ratio', bc_oh, 3, .TRUE.)
    ibc_o3   = bc_define('O3 mass mixing ratio', bc_o3, 3, .TRUE.)
    ibc_h2o2 = bc_define('H2O2 mass mixing ratio', bc_h2o2, 3, .TRUE.)
    ibc_no2  = bc_define('NO2 mass mixing ratio', bc_no2, 3, .TRUE.)

!gf see #146
    ibc_no3  = bc_define('NO3 mass mixing ratio', bc_no3, 3, .TRUE.)
!gf

    !--- Display parameter settings 
  
    csubmname = 'UNKNOWN' 
    IF (lham) csubmname = 'HAM'
    IF (lhammoz) csubmname = 'HAMMOZ' 

    CALL message('','')
    CALL message('', separator)
    CALL message('setham','Parameter settings for the ECHAM-'//TRIM(csubmname)//' aerosol model', &
                 level=em_info)

    CALL message('','---')
    CALL print_value('Seasalt emissions (nseasalt)              = ', nseasalt)
    SELECT CASE(nseasalt)
!!    CASE (0)
!!      CALL message('aero_initialize','No seasalt emissions (nseasalt=0). Are you sure?',level=em_warn)
      CASE (1)
        CALL message('','                            -> Monahan (1986)',level=em_param)
      CASE (2)
        CALL message('','                            -> Schulz et al. (2002)',level=em_param)
      !>>dod (redmine #44)
      CASE (4)
        CALL message('','                            -> Monahan (1986), bin scheme',level=em_param)
      CASE (5)
        CALL message('','                            -> Guelle et al. (2001)',level=em_param)
      CASE (6)
        CALL message('','                            -> Gong et al. (2003)',level=em_param)
      CASE (7)
        CALL message('','                            -> Long et al.(2011)',level=em_param)
      CASE (8)
        CALL message('','                            -> Gong et al.(2003)+T-dep.',level=em_param)
      !<<dod
      CASE DEFAULT
        WRITE (message_text,'(a,i0,a)') '   nseasalt = ',nseasalt,' not supported.'
        CALL message('setham',message_text, level=em_error)
    END SELECT

    CALL message('','---')
    CALL print_value('Air-sea exchange parameterisation for DMS emissions (npist) = ', npist)
    SELECT CASE(npist)
      CASE (1)
        CALL message('','                            -> Liss & Merlivat (1986)',level=em_param)
      CASE (2)
        CALL message('','                            -> Wanninkhof (1992)',level=em_param)
      CASE (3)
        CALL message('','                            -> Nightingale (2000)',level=em_param)
      CASE DEFAULT
        WRITE (message_text,'(a,i0,a)') '   npist = ',npist,' not supported.'
        CALL message('setham',message_text, level=em_error)
    END SELECT

    CALL message('','---')
    CALL print_value('             Aerosol diagnostics (laerocom_diag) = ', laerocom_diag)
!! ### is there a check that any(nrad > 0) if naerorad=true ???
    CALL print_value('Aerosol feedback with radiation (naerorad) = ', naerorad)
    IF (naerorad>0) THEN
      WRITE (message_text,'(A,25i3)') '                     nrad = ', nrad  !### hardcoded 25
      CALL message ('', message_text, level=em_param)
      CALL message ('', '        mode/bin   SW   LW  nradmix', level=em_param)
      DO jj=1, nclass
        IF (nrad(jj)==0) THEN
           WRITE (message_text,'(i0,t12,a1,t17,a1,t21,i0)') jj, 'F', 'F', nradmix(jj)
        ELSE IF (nrad(jj)==1) THEN
           WRITE (message_text,'(i0,t12,a1,t17,a1,t21,i0)') jj, 'T', 'F', nradmix(jj)
        ELSE IF (nrad(jj)==2) THEN
           WRITE (message_text,'(i0,t12,a1,t17,a1,t21,i0)') jj, 'F', 'T', nradmix(jj)
        ELSE IF (nrad(jj)==3) THEN
           WRITE (message_text,'(i0,t12,a1,t17,a1,t21,i0)') jj, 'T', 'T', nradmix(jj)
        END IF
        CALL message ('', message_text, level=em_param)
        IF (nrad(jj) > 3 .OR. nradmix(jj) > 3) THEN
          WRITE (message_text,'(a,i0,a,i0,a,i0)') 'Invalid nrad or nradmix: jbin=',jj,' nrad(jbin)=', &
                 nrad(jj),' nradmix(jbin)=',nradmix(jj)
          CALL message('setham',message_text, level=em_error)
        END IF
      END DO
      CALL message('','Refractive index mixing rule (nradmix):', level=em_param)
      CALL message('','     0 = ---, 1 = volume weighted, 2 = Maxwell-Garnett, 3 = Bruggeman', &
                   level=em_param)

      CALL message('','---')

      SELECT CASE(nraddiag)
        CASE (0)
          CALL message ('','No extended radiation diagnostics (nraddiag=0).',level=em_param)
        CASE (1)
          CALL message ('','Extended radiation diagnostics 2D (nraddiag=1).',level=em_param)
        CASE (2)
          CALL message ('','Extended radiation diagnostics 2D and 3D (nraddiag=2).',level=em_param)
        CASE DEFAULT
          WRITE (message_text,'(a,i0)') 'Invalid value for nraddiag: ', nraddiag
          CALL message ('setham', message_text, level=em_error)
      END SELECT

!>>NAJS: optical properties of dry aerosol
      SELECT CASE(nraddry)
        CASE (0)
          CALL message ('','No dry radiation diagnostics (nraddry=0).',level=em_param)
        CASE (1)
          CALL message ('','Dry radiation diagnostics: AOT  (nraddry=1).',level=em_param)
        CASE (2)
          CALL message ('','Dry radiation diagnostics: AOT & AAOT (nraddry=2).',level=em_param)
        CASE (3)
          CALL message ('','Dry radiation diagnostics: AOT, AAOT & extinction (nraddry=3).',level=em_param)
        CASE DEFAULT
          WRITE (message_text,'(a,i0)') 'Invalid value for nraddry: ', nraddry
          CALL message ('setham', message_text, level=em_error)
      END SELECT
!<<NAJS
!>>NAJS: lidar backscatter
      SELECT CASE(nradbeta)
        CASE (0)
          CALL message ('','No LIDAR diagnostics (nradbeta=0).',level=em_param)
        CASE (1)
          CALL message ('','LIDAR diagnostics: alfa_aer  (nradbeta=1).',level=em_param)
        CASE (2)
          CALL message ('','LIDAR diagnostics: alfa_aer, alfa_gas, beta_tot & lidar_ratio (nradbeta=2).',level=em_param)
        CASE (3)
          CALL message ('','Full LIDAR diagnostics: incl. att. backscatter from Earth and space (nradbeta=3).',level=em_param)
        CASE DEFAULT
          WRITE (message_text,'(a,i0)') 'Invalid value for nradbeta: ', nradbeta
          CALL message ('setham', message_text, level=em_error)
      END SELECT
!<<NAJS

    ELSE       ! naerorad==0
      CALL message('setham','Aerosol radiation interactions deactivated!', level=em_warn)
    END IF

    CALL message('','---')
    CALL print_value (' Mass fixer in convection : ', lomassfix)

    CALL message('','---')
    CALL message('', separator)

    !>>dod ndrydep = 1 is not tested. Override this check (at your own risk) with -DDRYDEP1 in your Makefile
#ifndef DRYDEP1
    IF (ndrydep == 1) THEN
       CALL message('', 'ndrydep = 1 in namelist hamctl: this value is not tested!')
       CALL message('', 'to proceed with ndrydep=1 set flag -DDRYDEP1 in the fortran 90 flags in your Makefile and re-compile')
       CALL message('', 'then run the model again at your own risk')
       CALL message('', 'else select another value for ndrydep')
       CALL finish('setham', 'unsupported namelist value ndrydep=1')
    END IF
#endif    
    !<<dod

    CALL print_value('Condensation of H2SO4 (lscond)', lscond)
    CALL print_value('Coagulation (lscoag)', lscoag)
    CALL print_value('M7 mass conservation check (lmass_diag)', lmass_diag)
    CALL print_value('(C)CN diagnostics (nccndiag)', nccndiag)
    SELECT CASE(nccndiag)
      CASE (0)
        CALL message('','-> OFF',level=em_param)
      CASE (1)
        CALL message('','-> 2D CCN diagnostics',level=em_param)
      CASE (2)
        CALL message('','-> 3D CCN diagnostics',level=em_param)
      CASE (3)
        CALL message('','-> 2D CCN + CN diagnostics',level=em_param)
      CASE (4)
        CALL message('','-> 3D CCN + CN diagnostics',level=em_param)
      CASE (5)
        CALL message('','-> 2D CCN + CN diagnostics + burdens',level=em_param)
      CASE (6)
        CALL message('','-> 3D CCN + CN diagnostics + burdens',level=em_param)
      CASE (7)
        CALL message('','-> 2D CCN (Surface) + CN Burden',level=em_param)
      CASE DEFAULT
        WRITE (message_text,*) 'nccndiag = ',nccndiag,' not supported.'
        CALL message('setham',message_text, level=em_error)
    END SELECT

  END SUBROUTINE setham


  INTEGER FUNCTION new_aerocomp(iclass, ispec, itrtype)

    USE mo_species,      ONLY: speclist
    USE mo_tracdef,      ONLY: itrprog

    IMPLICIT NONE

    !---function interface
    INTEGER, INTENT(IN) :: iclass
    INTEGER, INTENT(IN) :: ispec
    INTEGER, INTENT(IN), OPTIONAL :: itrtype

    !---local variables
    INTEGER :: i

    naerocomp = naerocomp + 1
          i = naerocomp

    aerocomp(i)%iclass   =  iclass
    aerocomp(i)%spid    =  ispec
    aerocomp(i)%species => speclist(ispec)

    IF (PRESENT(itrtype)) THEN
       aerocomp(i)%tracer_type = itrtype
    ELSE
       aerocomp(i)%tracer_type = itrprog
    END IF

    new_aerocomp = i

  END FUNCTION new_aerocomp

  SUBROUTINE print_aerocomp_info

    USE mo_tracdef,   ONLY: trlist
    USE mo_exception, ONLY: message_text, message, em_param, em_info
    USE mo_util_string, ONLY: separator

    INTEGER :: jt

    CALL message('',separator)
    WRITE(message_text,'(a)') 'Aerosol component information'
    CALL message('',message_text,level=em_info)
    CALL message('','',level=em_param)
    WRITE(message_text,'(a)') ': aerocomp id  iclass    spid     idt    tracer name'
    CALL message('',message_text,level=em_param)
    DO jt=1,naerocomp
       WRITE(message_text,'(a,5x,4(i7,1x),3x,a)') ':',jt, aerocomp(jt)%iclass, &
                                                  aerocomp(jt)%spid, aerocomp(jt)%idt, &
                                                  trlist%ti(aerocomp(jt)%idt)%fullname
       CALL message('',message_text,level=em_param)
    ENDDO
    CALL message('','',level=em_param)
    CALL message('',separator)

  END SUBROUTINE print_aerocomp_info

END MODULE mo_ham
