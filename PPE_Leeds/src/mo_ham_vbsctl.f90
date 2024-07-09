!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_vbsctl.f90
!!
!! \brief
!! mo_ham_vbsctl contains parameters, switches and initialization routines 
!! for the volatility basis set (VBS)
!!
!! \responsible_coder
!! Thomas Kühn (thk), thomas.h.kuhn@uef.fi
!!
!! \revision_history
!!   -# T. Kühn (UEF) - initial setup
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
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any 
!! kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ham_vbsctl

  !--- inherited types, data and functions
  USE mo_kind, ONLY: &
       dp       

  IMPLICIT NONE
  PRIVATE

  ! switch for basic setup
  INTEGER, PUBLIC :: nvbs_setup = 0 

  ! 3 VBS class setup with OC as non-volatile
  INTEGER, PUBLIC, PARAMETER :: vbs_setup_3class_OC = 0 

  ! 3 VBS class setup with separate non-volatile class
  INTEGER, PUBLIC, PARAMETER :: vbs_setup_3class = 1 

  ! switch for aqueous phase soa scheme (default = .true.)
  LOGICAL, PUBLIC :: laqsoa = .TRUE.

  ! Collection of parameters for the volatility base precursers
  TYPE, PUBLIC :: t_voc_prec
     !cross-referencing
     INTEGER :: spid        ! species index
     INTEGER :: idt         ! gas phase tracer index (set in mo_ham_init)
     INTEGER :: id_gasspec  ! index for HAM zgas array

     !physical properties

     !oxidation rates and their temperature dependence (Arrhenius)
     !OH
     REAL(dp) :: k_0_OH     ! pre-factor [mol/(m3*s)]
     REAL(dp) :: Eact_p_OH  ! reduced activation energy [K]: Eact_p = Eact/R
     !O3
     REAL(dp) :: k_0_O3     ! pre-factor [mol/(m3*s)]
     REAL(dp) :: Eact_p_O3  ! reduced activation energy [K]: Eact_p = Eact/R
     !NO3
     REAL(dp) :: k_0_NO3    ! pre-factor [mol/(m3*s)]
     REAL(dp) :: Eact_p_NO3 ! reduced activation energy [K]: Eact_p = Eact/R

     ! define what fraction goes where 
     ! (vector length depends on size of basis set)
     REAL(dp), ALLOCATABLE :: stoich_coeff(:) !Stoichiometric Coefficients

  END TYPE t_voc_prec

  ! Collection of the volatility base parameters
  TYPE, PUBLIC :: t_vbs_group
     ! cross-referencing
     INTEGER :: spid        ! species index
     INTEGER :: idt         ! gas phase tracer index (set in mo_ham_init)
     INTEGER, ALLOCATABLE :: &
          idx(:)            ! list of aerosol tracer indices (set in mo_ham_init)
     INTEGER :: id_vols     ! index in vols structure (for SALSA)
     INTEGER :: id_gasspec  ! index for HAM zgas array (set in mo_ham_subm_species)
     INTEGER :: id_aerospec ! index for HAM zaerml

     LOGICAL :: lcreateaero ! whether aerosol tracers should be created in ham_define_bins
     INTEGER :: spid_aero   ! needed if lcreateaero == .FALSE.

     ! effective saturation concentration -- Farina et al (2010)
     REAL(dp) :: C0         ! saturation concentration at reference temperature [kg/m3]
     REAL(dp) :: T0         ! reference temperature [K]
     REAL(dp) :: Hvap_eff   ! effective heat of vaporization Hvap_eff=Hvap/R  [K]

     ! physical properties
     REAL(dp) :: mv         ! molecular volume [m3]
  END TYPE t_vbs_group


  ! Collection of the aqsoa parameters
  TYPE, PUBLIC :: t_aq_soa
     ! cross-referencing
     INTEGER :: spid        ! species index
     INTEGER :: idt         ! gas phase tracer index (set in mo_ham_init)
     INTEGER, ALLOCATABLE :: &
          idx(:)            ! list of aerosol tracer indices (set in mo_ham_init)
     INTEGER :: id_aqsoa    ! index in aqsoa structure (for SALSA)
     INTEGER :: id_gasspec  ! index for HAM zgas array (set in mo_ham_subm_species)
     INTEGER :: id_aerospec ! index for HAM zaerml

     LOGICAL :: lcreateaero ! whether aerosol tracers should be created in ham_define_bins
     INTEGER :: spid_aero   ! needed if lcreateaero == .FALSE.

     !oxidation rates and their temperature dependence (Arrhenius)
     !OH
     REAL(dp) :: k_0_OH     ! pre-factor [mol/(m3*s)]
     REAL(dp) :: Eact_p_OH  ! reduced activation energy [K]: Eact_p = Eact/R
     !O3
     REAL(dp) :: k_0_O3     ! pre-factor [mol/(m3*s)]
     REAL(dp) :: Eact_p_O3  ! reduced activation energy [K]: Eact_p = Eact/R
     !NO3
     REAL(dp) :: k_0_NO3    ! pre-factor [mol/(m3*s)]
     REAL(dp) :: Eact_p_NO3 ! reduced activation energy [K]: Eact_p = Eact/R

     !Photodissociation 
     REAL(dp) :: photodis   ! pre-factor [rate s-1]    

     !Effective Henry coefficients  
     REAL(dp) :: Eff_henry_aerosol   ! pre-factor [rate s-1]   
     REAL(dp) :: Eff_henry_water   ! pre-factor [rate s-1] 

     ! physical properties
     REAL(dp) :: mv         ! molecular volume
  END TYPE t_aq_soa


  ! the following four variable will be set in
  ! *vbs_species* in *mo_ham_vbs*
  TYPE(t_voc_prec), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:) :: vbs_voc_prec
  TYPE(t_vbs_group), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:) :: vbs_set
  TYPE(t_aq_soa), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:) :: aqsoa_set

  INTEGER, PUBLIC :: vbs_nvocs = 0    ! number of VBS precursors
  INTEGER, PUBLIC :: vbs_ngroup = 0   ! number of VBS groups
  INTEGER, PUBLIC :: nclass_vbs = 0   ! number of modes/bins containing VBS soa (will
                                      ! be counted in *setham_vbs* in *mo_ham_vbs*)
  INTEGER, PUBLIC :: aqsoa_ngroup = 0 ! number of wet SOA groups
  INTEGER, PUBLIC :: nclass_aqsoa = 0 ! number of modes/bins containing wet soa (will
                                      ! be counted in *setham_vbs* in *mo_ham_vbs*)

  PUBLIC :: setham_vbs


CONTAINS

  SUBROUTINE setham_vbs
    
    ! -----------------------------------------------------------------------
    !
    ! SUBROUTINE setham_vbs
    !
    ! modifies pre-set switches of the ham_vbsctl
    ! namelist for the configuration of the volatility 
    ! basis set 'addon' for ham
    ! 
    ! Authors:
    ! --------
    ! Thomas Kühn, UEF    6/2015 --
    !
    ! setham_vbs is called from start_ham in mo_ham_init)
    ! 
    ! -----------------------------------------------------------------------
   
    USE mo_mpi,         ONLY: &
         p_parallel,    &
         p_parallel_io, &
         p_bcast,       &
         p_io

    USE mo_namelist,    ONLY: &
         open_nml,     &
         position_nml, &
         POSITIONED

    USE mo_exception,   ONLY: &
         message,      &
         message_text, &
         em_info,      &
         finish
    
    USE mo_util_string, ONLY: &
         separator

    IMPLICIT NONE
    
    INCLUDE 'ham_vbsctl.inc'
    
    ! Local variables:
    
    INTEGER :: ierr, inml, iunit

    ! -----------------------------------------------------------------------
    ! executable procedure
    ! -----------------------------------------------------------------------


    ! Read the namelist with the switches:
    CALL message('',separator)
    CALL message('setham_vbs', 'Reading namelist ham_vbsctl...', level=em_info)
    IF (p_parallel_io) THEN

      inml = open_nml('namelist.echam') 
      iunit = position_nml ('HAM_VBSCTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
         READ (iunit, ham_vbsctl)
      CASE DEFAULT
         CALL message('setham_vbs', 'Namelist HAM_VBSCTL not found -- using default values', level=em_info)
      END SELECT
    ENDIF
    
    ! Broadcast the switches over the processors:
    CALL p_bcast (nvbs_setup,     p_io)
    CALL p_bcast (laqsoa, p_io)

    ! write the values of the switches:
    CALL sethamvbs_log(nvbs_setup, laqsoa)

  END SUBROUTINE setham_vbs



  SUBROUTINE sethamvbs_log(nvbs_setup, laqsoa)

    USE mo_exception,    ONLY: &
         finish,       &
         message,      &
         message_text, &
         em_param

    USE mo_submodel,     ONLY: &
         print_value

    USE mo_util_string,  ONLY: &
         separator 
    
    IMPLICIT NONE

    INTEGER :: nvbs_setup
    LOGICAL :: laqsoa

    ! -----------------------------------------------------------------------
    ! executable procedure
    ! -----------------------------------------------------------------------

    CALL message('', separator) 
    CALL message(&
         'sethamvbs_log',&
         'Initialization of settings for volatility base',&
         level=em_param&
         ) 

    IF (laqsoa) THEN
       call message(&
            '',&
            'aqueous phase SOA on', &
            level=em_param&
            )
    END IF


    SELECT CASE(nvbs_setup)
    CASE (vbs_setup_3class_OC)
       CALL message(&
            '',&
            'two VBS classes plus OC selected',  &
            level=em_param&
            )
    CASE (vbs_setup_3class)
       CALL message(&
            '',&
            'three VBS classes selected',  &
            level=em_param&
            )
    CASE DEFAULT
       WRITE (message_text,'(a,i0)') &
            'ERROR: No scheme is implemented for nvbs_setup =',&
            nvbs_setup

       CALL message(&
            '',&
            message_text&
            )

       CALL finish(&
            'sethamvbs_log',&
            'run terminated'&
            )
    END SELECT

    CALL message('', separator)

  END SUBROUTINE sethamvbs_log



END MODULE mo_ham_vbsctl
