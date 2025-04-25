!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! This module contains all subroutines necessary to handle a CCN climatology
!!
!! @author 
!! <ol>
!! <li>S. Ferrachat (ETHZ)
!! </ol>
!!
!! @par Revision History
!! <ol>
!! <li>S. Ferrachat   (ETHZ) -  original code structure - (2010-03-xx) 
!!                                  
!! </ol>
!!
!! @par This module is used by
!! and to_be_added
!! 
!! @par Responsible coder
!! sylvaine.ferrachat@env.ethz.ch
!!
!! @par Copyright
!! 2010 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ECHAM is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!! violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!! copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!! an according license agreement with MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ccnclim

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream,  HYBRID, SURFACE
  USE mo_decomposition, ONLY: gc => global_decomposition, &
                              dc => local_decomposition
  USE mo_boundary_condition, ONLY: bc_nml, bc_define, p_bcast_bc,  &
                                   BC_REPLACE, BC_EVERYWHERE, BC_ALTITUDE, BC_LEVEL, BC_BOTTOM
  USE mo_external_field_processor, ONLY: EF_FILE, EF_LONLAT, EF_3D, EF_IGNOREYEAR, EF_NOINTER, &
                                         EF_VALUE, EF_UNDEFINED

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_ccnclim_submodel
  PUBLIC :: ccnclim_define_tracer
  PUBLIC :: set_ccn_3d
  PUBLIC :: ccnclim_avail_activ_lin_leaitch
  PUBLIC :: ccnclim_IN_setup
  PUBLIC :: ccnclim_free_memory

  INTEGER, PUBLIC :: idt_cdnc_ccnclim, idt_icnc_ccnclim !tracer indice

  !scalar quantities related to IN, in case the IN input is chosen 0D:
  REAL(dp), PUBLIC :: fracdusol_0d   = 0.010_dp  !fraction of dust particules in all soluble modes
  REAL(dp), PUBLIC :: fracduai_0d    = 0.094_dp  !fraction of dust particules in the insoluble accumulation mode
  REAL(dp), PUBLIC :: fracduci_0d    = 0.229_dp  !fraction of dust particules in the insoluble coarse mode
  REAL(dp), PUBLIC :: fracbcsol_0d   = 0.042_dp  !fraction of BC particules in all soluble modes
  REAL(dp), PUBLIC :: fracbcinsol_0d = 0._dp     !fraction of BC particules in all insoluble modes
  REAL(dp), PUBLIC :: rwetki_0d      = 1._dp       !wet radius for the insoluble aitken mode (m) dummy value
  REAL(dp), PUBLIC :: rwetai_0d      = 2.59e-08_dp !wet radius for the insoluble accumulation mode (m)
  REAL(dp), PUBLIC :: rwetci_0d      = 7.97e-08_dp !wet radius for the insoluble coarse mode (m)
  REAL(dp), PUBLIC :: rwetas_0d      = 1.09e-07_dp !wet radius for the soluble accumulation mode (m)

  !scalar quantities related to IN, in case the IN input is chosen 0D and with 0
  !as content (ninput_type_CCN = : IN_00D = -1
  REAL(dp), PUBLIC :: fracdusol_00d   = 0._dp !fraction of dust particules in all soluble modes
  REAL(dp), PUBLIC :: fracduai_00d    = 0._dp !fraction of dust particules in the insoluble accumulation mode
  REAL(dp), PUBLIC :: fracduci_00d    = 0._dp !fraction of dust particules in the insoluble coarse mode
  REAL(dp), PUBLIC :: fracbcsol_00d   = 0._dp !fraction of BC particules in all soluble modes
  REAL(dp), PUBLIC :: fracbcinsol_00d = 0._dp !fraction of BC particules in all insoluble modes
  REAL(dp), PUBLIC :: rwetki_00d      = 0._dp !wet radius for the insoluble aitken mode (m) dummy value
  REAL(dp), PUBLIC :: rwetai_00d      = 0._dp !wet radius for the insoluble accumulation mode (m)
  REAL(dp), PUBLIC :: rwetci_00d      = 0._dp !wet radius for the insoluble coarse mode (m)
  REAL(dp), PUBLIC :: rwetas_00d      = 0._dp !wet radius for the soluble accumulation mode (m)

  LOGICAL :: l_variable_ss = .true. !to take into account variable supersat (-> variable updraft velocity)
                                    ! when false, diagnostics for stratiform and convective cl are used
                                    ! therefore in this later case nss is used as a convenient proxy
                                    ! (though a bit dirty!) to distinguish between strat and conv (nss=2)

  INTEGER, PARAMETER, PUBLIC :: CCN_2D = 2
  INTEGER, PARAMETER, PUBLIC :: CCN_3D = 3
  INTEGER, PARAMETER, PUBLIC :: CCN_00D = -1
  INTEGER, PARAMETER, PUBLIC :: CCN_0CCN3D = 301
  INTEGER, PARAMETER, PUBLIC :: CCN_0SD3D = 302
  INTEGER, PUBLIC :: ninput_type_CCN = CCN_2D   !parameter to describe the CCN input type
                                                ! ninput_type_CCN = 2 (2D CCN input @1km and 8km)
                                                ! ninput_type_CCN = 3 (3D CCN input) 
                                                ! ninput_type_CCN = -1 (0 as CCN input) 
                                                ! ninput_type_CCN = 301 (3D CCN input, but CCN that enter CD activation artificially set to 0) 
                                                ! ninput_type_CCN = 302 (3D CCN input, but CCN that enter SD freezing artificially set to 0) 
  INTEGER :: nss = 4 !number of different supersaturation levels 
  REAL(dp), PARAMETER :: zupdraft_limits(1:3) = (/ 1.e-02_dp, 1.e-01_dp, 1._dp /) ! m.s-1

  INTEGER :: ninput_CCN !number of input field per physical CCN input variable 
                        ! (ninput_CCN=1 if 3D, ninput_CCN=2 if 2D for 'surf' and 'up')

  INTEGER, PARAMETER, PUBLIC :: IN_0D = 0
  INTEGER, PARAMETER, PUBLIC :: IN_3D = 3
  INTEGER, PARAMETER, PUBLIC :: IN_00D = -1
  INTEGER, PUBLIC :: ninput_type_IN = IN_0D !parameter to describe the IN input type
                                            ! ninput_type_IN = 0 (0D IN input, ie scalar)
                                            ! ninput_type_IN = 3 (3D IN input)
                                            ! ninput_type_IN = -1 (0 values in 0D IN input)

  INTEGER, ALLOCATABLE, PUBLIC :: ibc_ccn_an(:,:), ibc_ccn_pi(:,:), & !boundary condition indices for
                                  ibc_ccn_du(:), ibc_ccn_co(:)        !the CCN input 
  TYPE(bc_nml), ALLOCATABLE :: bc_ccn_an(:,:), bc_ccn_pi(:,:), bc_ccn_du(:), bc_ccn_co(:)

  INTEGER, PUBLIC :: ibc_fracdusol, ibc_fracduai, ibc_fracduci, &
                     ibc_fracbcsol, ibc_fracbcinsol, &
                     ibc_rwetki, ibc_rwetai, ibc_rwetci, ibc_rwetas

  TYPE(bc_nml) :: bc_fracdusol, bc_fracduai, bc_fracduci, &
                  bc_fracbcsol, bc_fracbcinsol, &
                  bc_rwetki, bc_rwetai, bc_rwetci, bc_rwetas 

  TYPE (t_stream), PUBLIC, POINTER :: ccnclim

  REAL(dp), ALLOCATABLE :: zccn_du(:,:)        !input CCN for dust
  REAL(dp), ALLOCATABLE :: zccn_co(:,:,:)      !input CCN for coarse mode (natural)
  REAL(dp), ALLOCATABLE :: zccn_an(:,:,:,:)    !input CCN for anthropogenic aerosols (fine)
  REAL(dp), ALLOCATABLE :: zccn_pi(:,:,:,:)    !input CCN for pre-industrial aerosols (fine)
  REAL(dp), ALLOCATABLE :: zccn_tot_2d(:,:,:)  !utility variable

  REAL(dp), ALLOCATABLE :: zccn_tot_3d(:,:,:,:)    !3d tot ccn (3d field + supersat dimension) 

  REAL(dp), PUBLIC, POINTER :: ccn_tot_strat(:,:,:) !final tot ccn field for stratiform case 
                                                    ! (supersat case chosen upon local strat updraft)
  REAL(dp), PUBLIC, POINTER :: ccn_tot_conv(:,:,:)  !final tot ccn field for convective case
                                                    ! (supersat case chosen upon local conv updraft)
  REAL(dp), PUBLIC, POINTER :: ccn_dust(:,:,:)      !final dust ccn field


  !SF: useless??
  !REAL(dp), PUBLIC, POINTER :: fracdusol(:,:,:)   !fraction of dust particules in all soluble modes
  !REAL(dp), PUBLIC, POINTER :: fracduai(:,:,:)    !fraction of dust particules in the insoluble accumulation mode
  !REAL(dp), PUBLIC, POINTER :: fracduci(:,:,:)    !fraction of dust particules in the insoluble coarse mode
  !REAL(dp), PUBLIC, POINTER :: fracbcsol(:,:,:)   !fraction of BC particules in all soluble modes
  !REAL(dp), PUBLIC, POINTER :: fracbcinsol(:,:,:) !fraction of BC particules in all insoluble modes

  !REAL(dp), PUBLIC, POINTER :: rwetki(:,:,:) !wet radius for the insoluble aitken mode
  !REAL(dp), PUBLIC, POINTER :: rwetai(:,:,:) !wet radius for the insoluble accumulation mode
  !REAL(dp), PUBLIC, POINTER :: rwetci(:,:,:) !wet radius for the insoluble coarse mode
  !REAL(dp), PUBLIC, POINTER :: rwetas(:,:,:) !wet radius for the soluble accumulation mode

CONTAINS

  !---------------------------------------------------------------------------
  !>
  !! @brief Initializes the CCN climatology submodel
  !!
  !! @remarks This subroutine reads into the ccnclim namelist
  !! and constructs the CCN climatology stream
  !! It also sets all the boundary conditions
  !! 

  SUBROUTINE init_ccnclim_submodel

    USE mo_memory_base,        ONLY: new_stream, add_stream_element,     &
                                     default_stream_setting, add_stream_reference, AUTO
    USE mo_submodel,           ONLY: print_value
    USE mo_mpi,                ONLY: p_io, p_parallel_io, p_parallel, p_bcast
    USE mo_namelist,           ONLY: open_nml, position_nml, POSITIONED, &
                                     LENGTH_ERROR, READ_ERROR
    USE mo_filename,           ONLY: out_filetype
    USE mo_exception,          ONLY: finish, message, message_text, em_info, em_error
    USE mo_activ,              ONLY: nfrzmod
    USE mo_param_switches,     ONLY: ncd_activ

    INTEGER :: ierr, inml, iunit, kinput, kss, kef_geom_ccn, kef_geom_in, kef_type_in, kef_dims_in, &
               kef_type_ccn, kef_dims_ccn
    CHARACTER(LEN=512) :: kef_template_in, kef_template_ccn
    CHARACTER(LEN=:), ALLOCATABLE :: c_ss(:), c_input_suffix(:)

    INCLUDE 'ccnclimctl.inc'

    !--- Namelist reading:
    IF (p_parallel_io) THEN
       inml = open_nml('namelist.echam')
       iunit = position_nml ('CCNCLIMCTL', inml, status=ierr)
       SELECT CASE (ierr)
       CASE (POSITIONED)
          READ (iunit, ccnclimctl)
      CASE (LENGTH_ERROR)
        CALL finish ('init_ccnclim_submodel', &
             'length error in namelist ccnclimctl')
      CASE (READ_ERROR)
        CALL finish ('init_ccnclim_submodel', &
             'read error in namelist.echam')
       END SELECT
    END IF

    IF (p_parallel) THEN
       CALL p_bcast (fracdusol_0d, p_io)
       CALL p_bcast (fracduai_0d, p_io)
       CALL p_bcast (fracduci_0d, p_io)
       CALL p_bcast (fracbcsol_0d, p_io)
       CALL p_bcast (fracbcinsol_0d, p_io)
       CALL p_bcast (rwetki_0d, p_io)
       CALL p_bcast (rwetai_0d, p_io)
       CALL p_bcast (rwetci_0d, p_io)
       CALL p_bcast (rwetas_0d, p_io)
       CALL p_bcast (l_variable_ss, p_io)
       CALL p_bcast (ninput_type_CCN, p_io)
       CALL p_bcast (ninput_type_IN, p_io)
    END IF

    !--- security checks
    IF (ncd_activ == 2) THEN
        ncd_activ = 1
        CALL message('init_ccn_submodel','ncd_activ reset to 1 because 2 is not supported with ccn climatology', &
                   level=em_info)
    ENDIF 

    SELECT CASE(ninput_type_ccn)
       CASE(CCN_2D,CCN_3D,CCN_00D,CCN_0CCN3D,CCN_0SD3D)
           CONTINUE
       CASE DEFAULT
           WRITE(message_text,'(a,i0,a)') 'Unsupported value for ninput_type_ccn (currently: ', &
                                        ninput_type_ccn, ').'
           CALL message('init_ccn_submodel',message_text,level=em_error)
    END SELECT

    SELECT CASE(ninput_type_in)
       CASE(IN_0D,IN_3D,IN_00D)
           CONTINUE
       CASE DEFAULT
           WRITE(message_text,'(a,i0,a)') 'Unsupported value for ninput_type_in (currently: ', &
                                        ninput_type_in, ').'
           CALL message('init_ccn_submodel',message_text,level=em_error)
    END SELECT

    !--- Adjust number of supersat levels
    IF (.NOT. l_variable_ss) THEN
       nss=2
    ENDIF

    ALLOCATE(CHARACTER(LEN=32) :: c_ss(nss))
    DO kss=1,nss
       IF (l_variable_ss) THEN
          WRITE(c_ss(kss),'(i0)') kss
       ELSE
          SELECT CASE (kss)
             CASE(1)
                WRITE(c_ss(kss),'(a)') 'strat'
             CASE(2)
                WRITE(c_ss(kss),'(a)') 'conv'
          END SELECT
       ENDIF
    ENDDO

    !--- Adjust the number of input fields per physical CCN variable, and set corrollary parameters
    SELECT CASE (ninput_type_CCN)
       CASE(CCN_2D)
           kef_type_ccn     = EF_FILE
           kef_template_ccn = 'ccn%Y0.nc'
           kef_geom_ccn = EF_LONLAT
           ninput_CCN = 2
           ALLOCATE(CHARACTER(len=32) :: c_input_suffix(ninput_CCN))
           c_input_suffix(1) = '_surf'
           c_input_suffix(2) = '_up'
       CASE(CCN_3D,CCN_0CCN3D,CCN_0SD3D)
           kef_type_ccn     = EF_FILE
           kef_template_ccn = 'ccn%Y0.nc'
           kef_geom_ccn = EF_3D
           ninput_CCN = 1
           ALLOCATE(CHARACTER(len=32) :: c_input_suffix(ninput_CCN))
          c_input_suffix(1) = ''
       CASE(CCN_00D)
           kef_type_ccn     = EF_VALUE
           kef_template_ccn = ''
           kef_geom_ccn     = EF_UNDEFINED
           kef_dims_ccn     = 0
           ninput_CCN       = 1
           ALLOCATE(CHARACTER(len=32) :: c_input_suffix(ninput_CCN))
          c_input_suffix(1) = ''
    END SELECT

    !--- Set parameters related to IN input
    SELECT CASE (ninput_type_IN)
       CASE(IN_0D,IN_00D)
           kef_type_in     = EF_VALUE
           kef_template_in = ''
           kef_geom_in     = EF_UNDEFINED
           kef_dims_in        = 0
       CASE(IN_3D)
           kef_type_in     = EF_FILE
           kef_template_in = 'in%Y0.nc' 
           kef_geom_in     = EF_3D
           kef_dims_in     = 3
    END SELECT

    !--- Boundary condition allocation:
    ALLOCATE(ibc_ccn_an(nss,ninput_CCN), ibc_ccn_pi(nss,ninput_CCN), &
             ibc_ccn_du(ninput_CCN), ibc_ccn_co(ninput_CCN) )
    ALLOCATE(bc_ccn_an(nss,ninput_CCN),  bc_ccn_pi(nss,ninput_CCN),  &
             bc_ccn_du(ninput_CCN),  bc_ccn_co(ninput_CCN)  )

    !--- Set the boundary conditions:
    ibc_ccn_an(1:nss,1:ninput_CCN) = 0
    ibc_ccn_pi(1:nss,1:ninput_CCN) = 0
    ibc_ccn_du(1:ninput_CCN) = 0
    ibc_ccn_co(1:ninput_CCN) = 0

    bc_ccn_an(1:nss,1:ninput_CCN)%ef_type        = kef_type_ccn
    bc_ccn_an(1:nss,1:ninput_CCN)%ef_template    = kef_template_ccn
    bc_ccn_an(1:nss,1:ninput_CCN)%ef_timedef     = EF_IGNOREYEAR
    bc_ccn_an(1:nss,1:ninput_CCN)%ef_interpolate = EF_NOINTER
    bc_ccn_an(1:nss,1:ninput_CCN)%ef_geometry    = kef_geom_ccn  
    bc_ccn_an(1:nss,1:ninput_CCN)%ef_actual_unit = 'number m-3'  

    bc_ccn_pi(1:nss,1:ninput_CCN)%ef_type        = kef_type_ccn
    bc_ccn_pi(1:nss,1:ninput_CCN)%ef_template    = kef_template_ccn
    bc_ccn_pi(1:nss,1:ninput_CCN)%ef_timedef     = EF_IGNOREYEAR
    bc_ccn_pi(1:nss,1:ninput_CCN)%ef_interpolate = EF_NOINTER
    bc_ccn_pi(1:nss,1:ninput_CCN)%ef_geometry    = kef_geom_ccn  
    bc_ccn_pi(1:nss,1:ninput_CCN)%ef_actual_unit = 'number m-3'  

    bc_ccn_du(1:ninput_CCN)%ef_type        = kef_type_ccn
    bc_ccn_du(1:ninput_CCN)%ef_template    = kef_template_ccn
    bc_ccn_du(1:ninput_CCN)%ef_timedef     = EF_IGNOREYEAR
    bc_ccn_du(1:ninput_CCN)%ef_interpolate = EF_NOINTER
    bc_ccn_du(1:ninput_CCN)%ef_geometry    = kef_geom_ccn  
    bc_ccn_du(1:ninput_CCN)%ef_actual_unit = 'number m-3'  

    bc_ccn_co(1:ninput_CCN)%ef_type        = kef_type_ccn
    bc_ccn_co(1:ninput_CCN)%ef_template    = kef_template_ccn
    bc_ccn_co(1:ninput_CCN)%ef_timedef     = EF_IGNOREYEAR
    bc_ccn_co(1:ninput_CCN)%ef_interpolate = EF_NOINTER
    bc_ccn_co(1:ninput_CCN)%ef_geometry    = kef_geom_ccn  
    bc_ccn_co(1:ninput_CCN)%ef_actual_unit = 'number m-3'  

    DO kinput=1,ninput_CCN
       DO kss=1,nss

          SELECT CASE(ninput_type_CCN)
          CASE(CCN_00D)
              bc_ccn_an(kss,kinput)%ef_value   = 0._dp
              bc_ccn_pi(kss,kinput)%ef_value   = 0._dp
           
          CASE(CCN_3D,CCN_0CCN3D,CCN_0SD3D)
             bc_ccn_an(kss,kinput)%ef_varname = 'CCN_an_'//TRIM(ADJUSTL(c_ss(kss)))// &
                                                '_mo'//TRIM(ADJUSTL(c_input_suffix(kinput))) 
             bc_ccn_pi(kss,kinput)%ef_varname = 'CCN_pi_'//TRIM(ADJUSTL(c_ss(kss)))// &
                                                '_mo'//TRIM(ADJUSTL(c_input_suffix(kinput))) 
          END SELECT

          ibc_ccn_an(kss,kinput) = bc_define('Fine anthropogenic CCN input',bc_ccn_an(kss,kinput), &
                                             ninput_type_CCN, .TRUE.) 
          ibc_ccn_pi(kss,kinput) = bc_define('Fine pre-industrial CCN input',bc_ccn_pi(kss,kinput), &
                                             ninput_type_CCN, .TRUE.) 
       ENDDO
          SELECT CASE(ninput_type_CCN)
          CASE(CCN_00D)
              bc_ccn_du(kinput)%ef_value   = 0._dp
              bc_ccn_co(kinput)%ef_value   = 0._dp
           
          CASE(CCN_3D,CCN_0CCN3D,CCN_0SD3D)
               bc_ccn_du(kinput)%ef_varname = 'CCN_dust'// &
                                              '_mo'//TRIM(ADJUSTL(c_input_suffix(kinput))) 
               bc_ccn_co(kinput)%ef_varname = 'CCN_coarse'// &
                                              '_mo'//TRIM(ADJUSTL(c_input_suffix(kinput))) 
          END SELECT

       ibc_ccn_du(kinput) = bc_define('Dust CCN input',bc_ccn_du(kinput), &
                                      ninput_type_CCN, .TRUE.) 
       ibc_ccn_co(kinput) = bc_define('Coarse (natural) CCN input',bc_ccn_co(kinput), &
                                      ninput_type_CCN, .TRUE.) 
    ENDDO

    ibc_fracdusol   = 0
    ibc_fracduai    = 0
    ibc_fracduci    = 0
    ibc_fracbcsol   = 0
    ibc_fracbcinsol = 0
    ibc_rwetki = 0
    ibc_rwetai = 0
    ibc_rwetci = 0
    ibc_rwetas = 0

    bc_fracdusol%ef_type        = kef_type_in
    bc_fracdusol%ef_template    = kef_template_in
    bc_fracdusol%ef_timedef     = EF_IGNOREYEAR
    bc_fracdusol%ef_interpolate = EF_NOINTER
    bc_fracdusol%ef_geometry    = kef_geom_in
    bc_fracdusol%ef_actual_unit = '1'
  
    bc_fracduai%ef_type        = kef_type_in
    bc_fracduai%ef_template    = kef_template_in
    bc_fracduai%ef_timedef     = EF_IGNOREYEAR
    bc_fracduai%ef_interpolate = EF_NOINTER
    bc_fracduai%ef_geometry    = kef_geom_in
    bc_fracduai%ef_actual_unit = '1'
  
    bc_fracduci%ef_type        = kef_type_in
    bc_fracduci%ef_template    = kef_template_in
    bc_fracduci%ef_timedef     = EF_IGNOREYEAR
    bc_fracduci%ef_interpolate = EF_NOINTER
    bc_fracduci%ef_geometry    = kef_geom_in
    bc_fracduci%ef_actual_unit = '1'
  
    bc_fracbcsol%ef_type        = kef_type_in
    bc_fracbcsol%ef_template    = kef_template_in
    bc_fracbcsol%ef_timedef     = EF_IGNOREYEAR
    bc_fracbcsol%ef_interpolate = EF_NOINTER
    bc_fracbcsol%ef_geometry    = kef_geom_in
    bc_fracbcsol%ef_actual_unit = '1'
  
    bc_fracbcinsol%ef_type        = kef_type_in
    bc_fracbcinsol%ef_template    = kef_template_in
    bc_fracbcinsol%ef_timedef     = EF_IGNOREYEAR
    bc_fracbcinsol%ef_interpolate = EF_NOINTER
    bc_fracbcinsol%ef_geometry    = kef_geom_in
    bc_fracbcinsol%ef_actual_unit = '1'
  
    bc_rwetki%ef_type        = kef_type_in
    bc_rwetki%ef_template    = kef_template_in
    bc_rwetki%ef_timedef     = EF_IGNOREYEAR
    bc_rwetki%ef_interpolate = EF_NOINTER
    bc_rwetki%ef_geometry    = kef_geom_in
    bc_rwetki%ef_actual_unit = 'm'
  
    bc_rwetai%ef_type        = kef_type_in
    bc_rwetai%ef_template    = kef_template_in
    bc_rwetai%ef_timedef     = EF_IGNOREYEAR
    bc_rwetai%ef_interpolate = EF_NOINTER
    bc_rwetai%ef_geometry    = kef_geom_in
    bc_rwetai%ef_actual_unit = 'm'
  
    bc_rwetci%ef_type        = kef_type_in
    bc_rwetci%ef_template    = kef_template_in
    bc_rwetci%ef_timedef     = EF_IGNOREYEAR
    bc_rwetci%ef_interpolate = EF_NOINTER
    bc_rwetci%ef_geometry    = kef_geom_in
    bc_rwetci%ef_actual_unit = 'm'
  
    bc_rwetas%ef_type        = kef_type_in
    bc_rwetas%ef_template    = kef_template_in
    bc_rwetas%ef_timedef     = EF_IGNOREYEAR
    bc_rwetas%ef_interpolate = EF_NOINTER
    bc_rwetas%ef_geometry    = kef_geom_in
    bc_rwetas%ef_actual_unit = 'm'

    SELECT CASE (ninput_type_IN)
       CASE(IN_0D)
           bc_fracdusol%ef_value   = fracdusol_0d
           bc_fracduai%ef_value    = fracduai_0d
           bc_fracduci%ef_value    = fracduci_0d
           bc_fracbcsol%ef_value   = fracbcsol_0d
           bc_fracbcinsol%ef_value = fracbcinsol_0d
    
           bc_rwetki%ef_value = rwetki_0d
           bc_rwetai%ef_value = rwetai_0d
           bc_rwetci%ef_value = rwetci_0d
           bc_rwetas%ef_value = rwetas_0d
       CASE(IN_00D)
           bc_fracdusol%ef_value   = fracdusol_00d
           bc_fracduai%ef_value    = fracduai_00d
           bc_fracduci%ef_value    = fracduci_00d
           bc_fracbcsol%ef_value   = fracbcsol_00d
           bc_fracbcinsol%ef_value = fracbcinsol_00d
    
           bc_rwetki%ef_value = rwetki_00d
           bc_rwetai%ef_value = rwetai_00d
           bc_rwetci%ef_value = rwetci_00d
           bc_rwetas%ef_value = rwetas_00d
       CASE(IN_3D)
           bc_fracdusol%ef_varname    = 'FRACDUSOL'
           bc_fracduai%ef_varname     = 'FRACDUAI'
           bc_fracduci%ef_varname     = 'FRACDUCI'
           bc_fracbcsol%ef_varname    = 'FRACBCSOL'
           bc_fracbcinsol%ef_varname  = 'FRACBCINSOL'

           bc_rwetki%ef_varname  = 'RWETKI'
           bc_rwetai%ef_varname  = 'RWETAI'
           bc_rwetci%ef_varname  = 'RWETCI'
           bc_rwetas%ef_varname  = 'RWETAS'
    END SELECT

    ibc_fracdusol = bc_define('Fraction of dust particules in all soluble modes', &
                              bc_fracdusol, kef_dims_in, .TRUE.)
    ibc_fracduai = bc_define('Fraction of dust particules in the insoluble accumulation mode', &
                             bc_fracduai, kef_dims_in, .TRUE.)
    ibc_fracduci = bc_define('Fraction of dust particules in the insoluble coarse mode', &
                             bc_fracduci, kef_dims_in, .TRUE.)
    ibc_fracbcsol = bc_define('Fraction of BC particules in all soluble modes', &
                              bc_fracbcsol, kef_dims_in, .TRUE.)
    ibc_fracbcinsol = bc_define('Fraction of BC particules in all insoluble modes', &
                                bc_fracbcinsol, kef_dims_in, .TRUE.)
    ibc_rwetki = bc_define('Wet radius for the insoluble aitken mode', &
                           bc_rwetki, kef_dims_in, .TRUE.)
    ibc_rwetai = bc_define('Wet radius for the insoluble accumulation mode', &
                           bc_rwetai, kef_dims_in, .TRUE.)
    ibc_rwetci = bc_define('Wet radius for the insoluble coarse mode', &
                           bc_rwetci, kef_dims_in, .TRUE.)
    ibc_rwetas = bc_define('Wet radius for the soluble accumulation mode', &
                           bc_rwetas, kef_dims_in, .TRUE.)
  
    !--- Print info:
    IF (p_parallel_io) THEN
       CALL print_value('init_ccnclim_submodel, variable supersat. taken into account', &
                       l_variable_ss)
       CALL print_value('init_ccnclim_submodel, CCN input type', &
                       ninput_type_CCN)
       CALL print_value('init_ccnclim_submodel, IN input type', &
                       ninput_type_IN)
    ENDIF

    !--- set the number of freezing modes:
    nfrzmod = 1

    !--- creates the ccnclim stream
    CALL new_stream (ccnclim, 'ccnclim', filetype = out_filetype, &
         post_suf = '_ccnclim', rest_suf = '_ccnclim')

    !--- add standard fields for post-processing:

    CALL add_stream_reference (ccnclim, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (ccnclim, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (ccnclim, 'aps'     ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (ccnclim, 'gboxarea','geoloc',lpost=.TRUE.)

    !--- add stream elements:

    CALL default_stream_setting (ccnclim, units     = 'number m-3',&
                                          lrerun    = .FALSE. ,     &
                                          laccu     = .FALSE. ,    &
                                          lpost     = .TRUE. ,     &
                                          leveltype = SURFACE,     &
                                          code      = AUTO         )

    !----- 3D fields:
    CALL default_stream_setting (ccnclim, leveltype = HYBRID)

    CALL add_stream_element (ccnclim, 'CCN_TOT_STRAT', ccn_tot_strat,    &
                            longname='total CCN for stratiform clouds')

    CALL add_stream_element (ccnclim, 'CCN_TOT_CONV', ccn_tot_conv,      &
                            longname='total CCN for convective clouds')   

    CALL add_stream_element (ccnclim, 'CCN_DUST', ccn_dust,    &
                            longname='dust CCN (fraction of coarse natural)')

    !--- allocate memory for module variables:

    SELECT CASE(ninput_type_CCN)
           CASE(CCN_2D)
               ALLOCATE(zccn_du(dc%nproma,ninput_CCN))
               ALLOCATE(zccn_co(dc%nproma,1,ninput_CCN))
               ALLOCATE(zccn_an(dc%nproma,1,nss,ninput_CCN))
               ALLOCATE(zccn_pi(dc%nproma,1,nss,ninput_CCN))
               ALLOCATE(zccn_tot_2d(dc%nproma,nss,ninput_CCN))
           CASE(CCN_3D,CCN_00D,CCN_0CCN3D,CCN_0SD3D)
               ALLOCATE(zccn_co(dc%nproma,dc%nlev,ninput_CCN))
               ALLOCATE(zccn_an(dc%nproma,dc%nlev,nss,ninput_CCN))
               ALLOCATE(zccn_pi(dc%nproma,dc%nlev,nss,ninput_CCN))
    END SELECT
    ALLOCATE(zccn_tot_3d(dc%nproma, dc%nlev, nss, dc% ngpblks))

  END SUBROUTINE init_ccnclim_submodel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! ccnclim_define_tracer: create ECHAM tracers for CDNC and ICNC
!!
!! @author see module info
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info
!!
!! @par This subroutine is called by
!! init_submodels
!!
!! @par Externals:
!! <ol>
!! <li>none
!! </ol>
!!
!! @par Responsible coder
!! sylvaine.ferrachat@env.ethz.ch
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ccnclim_define_tracer

    USE mo_exception,         ONLY: message, message_text, em_debug, em_error
    USE mo_tracer,            ONLY: new_tracer, new_diag_burden
    USE mo_tracdef,           ONLY: OFF, ON, OK

    IMPLICIT NONE

    INTEGER :: iburdenid, ierr

    !---executable procedure

    CALL message('ccnclim_define_tracer', 'Starting ...', level=em_debug)

    iburdenid = new_diag_burden('CDNC', itype=1, lclobber=.false.)

    CALL new_tracer('CDNC', 'CCNCLIM', ierr=ierr,      &
                     idx=idt_cdnc_ccnclim,         &
                     units='1 kg-1',               &
                     table=131,                    &
                     code=131,                     &
                     nwrite=1,                     &
                     burdenid=iburdenid,           &
                     nconv=OFF,                    &
              !!     nconvmassfix=OFF,             &
                     nvdiff=ON,                    &
                     nint=ON,                      &
                     longname='cloud droplet number concentration')

    IF (ierr /= OK) THEN
      WRITE(message_text,*) 'new_tracer CDNC returned error code ',ierr
      CALL message('ham_define_tracer', message_text, level=em_error)
    END IF

    iburdenid = new_diag_burden('ICNC', itype=1, lclobber=.false.)

    CALL new_tracer('ICNC', 'CCNCLIM', ierr=ierr,      &
                     idx=idt_icnc_ccnclim,         &
                     units='1 kg-1',               &
                     table=131,                    &
                     code=132,                     &
                     nwrite=1,                     &
                     burdenid=iburdenid,           &
                     nconv=OFF,                    &
              !!     nconvmassfix=OFF,             &
                     nvdiff=ON,                    &
                     nint=ON,                      &
                     longname='ice crystal number concentration')

    IF (ierr /= OK) THEN
      WRITE(message_text,*) 'new_tracer ICNC returned error code ',ierr
      CALL message('ham_define_tracer', message_text, level=em_error)
    END IF

  END SUBROUTINE ccnclim_define_tracer

  !---------------------------------------------------------------------------
  !>
  !! @brief Computes 3D CCN values
  !! 
  !! @remarks This routine puts together the necessary CCN boundary
  !! conditions, either 2D or 3D, and creates 3D CCN fields out of that.
  !!

  SUBROUTINE set_ccn_3d(kproma, kbdim, klev, krow, papm1)

    USE mo_boundary_condition, ONLY: bc_apply

    INTEGER,  INTENT(in) :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(in) :: papm1(kbdim,klev)
    
    !--- Local vars:
    INTEGER  :: kss, kinput

    SELECT CASE(ninput_type_CCN)
       CASE(CCN_2D) !2D input --> needs to reconstruct the vertical profiles based on surf and up inputs 

           DO kinput=1,ninput_CCN

              CALL bc_apply(ibc_ccn_du(kinput), kproma, krow, zccn_du(1:kproma,kinput))
              CALL bc_apply(ibc_ccn_co(kinput), kproma, krow, zccn_co(1:kproma,1,kinput))

              DO kss=1,nss
                 CALL bc_apply(ibc_ccn_an(kss,kinput), kproma, krow, zccn_an(1:kproma,1,kss,kinput))
                 CALL bc_apply(ibc_ccn_pi(kss,kinput), kproma, krow, zccn_pi(1:kproma,1,kss,kinput))

                 zccn_tot_2d(1:kproma,kss,kinput) = zccn_an(1:kproma,1,kss,kinput)  &
                                                  + zccn_pi(1:kproma,1,kss,kinput)  &
                                                  + zccn_co(1:kproma,1,kinput)
              ENDDO
           ENDDO

           DO kss=1,nss
              CALL ccn_profile(kproma, kbdim, klev, papm1,                 &
                               zccn_tot_2d(1:kproma,kss,1), zccn_tot_2d(1:kproma,kss,2), &
                               zccn_tot_3d(1:kproma,:,kss,krow))
           ENDDO

           CALL ccn_profile(kproma, kbdim, klev, papm1,               &
                            zccn_du(1:kproma,1), zccn_du(1:kproma,2), &
                            ccn_dust(1:kproma,:,krow))

       CASE(CCN_3D,CCN_00D,CCN_0CCN3D,CCN_0SD3D) !3D input --> CCN input ready to use

           CALL bc_apply(ibc_ccn_du(1), kproma, krow, ccn_dust(1:kproma,:,krow))
           CALL bc_apply(ibc_ccn_co(1), kproma, krow, zccn_co(1:kproma,:,1))

           DO kss=1,nss
              CALL bc_apply(ibc_ccn_an(kss,1), kproma, krow, zccn_an(1:kproma,:,kss,1))
              CALL bc_apply(ibc_ccn_pi(kss,1), kproma, krow, zccn_pi(1:kproma,:,kss,1))

              zccn_tot_3d(1:kproma,:,kss,krow) = zccn_an(1:kproma,:,kss,1) &
                                               + zccn_pi(1:kproma,:,kss,1) &
                                               + zccn_co(1:kproma,:,1)
           ENDDO
    END SELECT

  END SUBROUTINE set_ccn_3d

  !---------------------------------------------------------------------------
  !>
  !! @brief Utility routine to compute the vertical profiles
  !! 
  !! @remarks The calculation reproduces what is implemented for
  !! acdnc setup in case of fully prescribed CDNC

  SUBROUTINE ccn_profile(kproma, kbdim, klev, papm1, pccnsurf, pccnup, pccn)

    INTEGER,  INTENT(in)  :: kproma, kbdim, klev
    REAL(dp), INTENT(in)  :: papm1(kbdim,klev)
    REAL(dp), INTENT(in)  :: pccnsurf(kbdim), pccnup(kbdim)
    REAL(dp), INTENT(out) :: pccn(kbdim,klev)

    !--- Local variables:

    INTEGER  :: jexp, & !exponent
                jk

    REAL(dp) :: zpthresh,          & !pressure threshold (boundary layer)
                zprat(kbdim,klev), &
                ztmp1(kbdim,klev), &
                ztmp2(kbdim,klev)

    LOGICAL  :: lobl(kbdim,klev)

    jexp = 2
    zpthresh = 80000._dp

    zprat(1:kproma,:) = ( MIN(8._dp, zpthresh/papm1(1:kproma,:)) )**jexp

    lobl(1:kproma,:)  = (papm1(1:kproma,:) < zpthresh)

    DO jk=1,klev
       ztmp1(1:kproma,jk) = pccnup(1:kproma) &
                          + (pccnsurf(1:kproma)-pccnup(1:kproma))*(EXP(1._dp-zprat(1:kproma,jk)))

       ztmp2(1:kproma,jk) = pccnsurf(1:kproma)
    ENDDO

    pccn(1:kproma,:) = MERGE(ztmp1(1:kproma,:), ztmp2(1:kproma,:), lobl(1:kproma,:))

  END SUBROUTINE ccn_profile

  !---------------------------------------------------------------------------
  !>
  !! @brief Computes available particules for activation
  !! 
  !! @remarks Preparatory routine for Lin&Leaitch activation scheme. 
  !! Available CCN's is computed depending on pre-defined values at nss different supersaturation levels
  !! WARNING: 
  !! nss=1 means that supersaturation levels *are not taken into account*
  !! If nss/=1, then a given supersaturation is chosen according to the local updraft velocity

  SUBROUTINE ccnclim_avail_activ_lin_leaitch(kproma, kbdim, klev, krow, lstrat, pw)

    USE mo_activ, ONLY: na
    USE mo_conv,  ONLY: na_cv

    INTEGER, INTENT(IN)  :: kproma, kbdim, klev, krow
    LOGICAL, INTENT(IN)  :: lstrat !switch for stratiform or convective case
    REAL(dp), INTENT(IN) :: pw(kbdim,klev)

    !--- local variables:
    INTEGER  :: iss
    LOGICAL  :: ll(kbdim,klev)
    REAL(dp) :: zupdraft_matrix(kbdim,klev,nss) ! zupdraft_matrix(1:kproma,:,iss) = 1 where 
                                                ! the updraft velocity satisfies the criteria for
                                                ! the iss-th' supersaturation level
    REAL(dp) :: zccn(kbdim,klev)

    !--- Compute the updraft matrix

    IF (l_variable_ss) THEN ! variable supersaturation levels taken into account

       ll(1:kproma,:) = (pw(1:kproma,:) <= zupdraft_limits(1)) !lower limit
       zupdraft_matrix(1:kproma,:,1) = MERGE(1, 0, ll(1:kproma,:))
    
       ll(1:kproma,:) = (pw(1:kproma,:) > zupdraft_limits(nss-1)) !upper limit
       zupdraft_matrix(1:kproma,:,nss) = MERGE(1, 0, ll(1:kproma,:))
    
       DO iss=2,nss-1
          ll(1:kproma,:) = (pw(1:kproma,:) <= zupdraft_limits(iss)) .AND. &
                           (pw(1:kproma,:) >  zupdraft_limits(iss-1))         !intermediate ranges
          zupdraft_matrix(1:kproma,:,iss) = MERGE(1, 0, ll(1:kproma,:))
       ENDDO

       !--- Fill the CCN final field with CCN at the relevant supersat
       !    (matrix product of updraft_matrix times CCN's at different supersat vector)
       zccn(1:kproma,:) = 0._dp
   
       DO iss=1,nss
          zccn(1:kproma,:) = zccn(1:kproma,:) + zupdraft_matrix(1:kproma,:,iss)*zccn_tot_3d(1:kproma,:,iss,krow)
       ENDDO

       IF (lstrat) THEN
          !-- stratiform case
          ccn_tot_strat(1:kproma,:,krow) = zccn(1:kproma,:)
          na(1:kproma,:,krow)            = zccn(1:kproma,:) ![m-3]
       ELSE
          !-- convective case
          ccn_tot_conv(1:kproma,:,krow) = zccn(1:kproma,:)
          na_cv(1:kproma,:,krow)        = zccn(1:kproma,:) ![m-3] 
       ENDIF

    ELSE !no variable supersaturation levels (strat / conv instead)

       IF (lstrat) THEN
          !-- stratiform case
          ccn_tot_strat(1:kproma,:,krow) = zccn_tot_3d(1:kproma,:,1,krow)
          SELECT CASE(ninput_type_CCN)
          CASE(CCN_0CCN3D)
                  na(1:kproma,:,krow)    = 0._dp ![m-3]
          CASE DEFAULT 
                  na(1:kproma,:,krow)    = zccn_tot_3d(1:kproma,:,1,krow) ![m-3]
          END SELECT
       ELSE
          !-- convective case
          ccn_tot_conv(1:kproma,:,krow) = zccn_tot_3d(1:kproma,:,2,krow)
          SELECT CASE(ninput_type_CCN)
          CASE(CCN_0CCN3D)
                  na_cv(1:kproma,:,krow) = 0._dp ![m-3] 
          CASE DEFAULT 
                  na_cv(1:kproma,:,krow) = zccn_tot_3d(1:kproma,:,2,krow) ![m-3] 
          END SELECT
       ENDIF

    ENDIF

  END SUBROUTINE ccnclim_avail_activ_lin_leaitch

  !---------------------------------------------------------------------------
  !>
  !! @brief Routine to set up necessary vars for mixed-phase and cirrus freezing calculations
  !! 
  !! @remarks Here, some quantities are surrogates for their m7 modes counterparts
  !! and for this reason they keep names refering to m7 modes.
  !! There is no dependency to HAM and m7.
  !! These quantities were diagnosed in a preliminary full HAM run
  !!

  SUBROUTINE ccnclim_IN_setup(kproma, kbdim, klev, krow,                                  &
                              prho, prho_rcp,                                             &
                              prwetki, prwetai, prwetci,                                  &
                              pfracdusol, pfracduai, pfracduci, pfracbcsol, pfracbcinsol, &
                              pascs, papnx, paprx, papsigx, ld_het                        )

    USE mo_activ, ONLY: nfrzmod
    USE mo_boundary_condition, ONLY: bc_apply

    INTEGER,  INTENT(in)  :: kproma, kbdim, klev, krow
    REAL(dp), INTENT(in)  :: prho(kbdim,klev)     ! air density
    REAL(dp), INTENT(in)  :: prho_rcp(kbdim,klev) ! reciprocal air density
    !-- for mixed-phase freezing:
    REAL(dp), INTENT(out) :: prwetki(kbdim,klev)  ! wet radius, aitken insoluble mode
    REAL(dp), INTENT(out) :: prwetai(kbdim,klev)  ! wet radius, accumulation insoluble mode
    REAL(dp), INTENT(out) :: prwetci(kbdim,klev)  ! wet radius, coarse insoluble mode
    REAL(dp), INTENT(out) :: pfracdusol(kbdim,klev)   ! total fraction of dust particules in all soluble modes
    REAL(dp), INTENT(out) :: pfracduai(kbdim,klev)    ! fraction of dust particules in the accum. soluble mode 
    REAL(dp), INTENT(out) :: pfracduci(kbdim,klev)    ! fraction of dust particules in the coarse soluble mode 
    REAL(dp), INTENT(out) :: pfracbcsol(kbdim,klev)   ! total fraction of BC particules in all soluble modes 
    REAL(dp), INTENT(out) :: pfracbcinsol(kbdim,klev) ! total fraction of BC particules in all insoluble modes 
    !-- for cirrus freezing:
    REAL(dp), INTENT(out) :: pascs(kbdim,klev)           ! soluble aerosol number conc. 
    REAL(dp), INTENT(out) :: papnx(kbdim,klev,nfrzmod)   ! aerosol number available for freezing [1/cm3] 
    REAL(dp), INTENT(out) :: paprx(kbdim,klev,nfrzmod)   ! radius of aerosols avail. for freezing  [cm] 
    REAL(dp), INTENT(out) :: papsigx(kbdim,klev,nfrzmod) ! std. dev. of aerosols available for freezing
    LOGICAL, INTENT(out)  :: ld_het  !switch to set heterogeneous freezing below 235K (cirrus scheme)

    REAL(dp) :: zrwetas(kbdim,klev) ! wet radius, soluble accumulation mode

    !--- Mixed-phase freezing setup:

     !-- Wet radii:
     CALL bc_apply(ibc_rwetki, kproma, krow, prwetki(1:kproma,:))
     CALL bc_apply(ibc_rwetai, kproma, krow, prwetai(1:kproma,:))
     CALL bc_apply(ibc_rwetci, kproma, krow, prwetci(1:kproma,:))
     CALL bc_apply(ibc_rwetas, kproma, krow, zrwetas(1:kproma,:))

     !-- Various dust and BC fractions:
!SFold     pfracdusol(1:kproma,:)   = ccn_dust(1:kproma,:,krow) / MAX(ccn_tot_strat(1:kproma,:,krow), EPSILON(1.0_dp))

     CALL bc_apply(ibc_fracdusol, kproma, krow, pfracdusol(1:kproma,:))
     CALL bc_apply(ibc_fracduai, kproma, krow, pfracduai(1:kproma,:))
     CALL bc_apply(ibc_fracduci, kproma, krow, pfracduci(1:kproma,:))
     CALL bc_apply(ibc_fracbcsol, kproma, krow, pfracbcsol(1:kproma,:))
     CALL bc_apply(ibc_fracbcinsol, kproma, krow, pfracbcinsol(1:kproma,:))

    !--- Cirrus freezing setup:

     !-- Soluble aerosol number concentration:
     pascs(1:kproma,:) = ccn_tot_strat(1:kproma,:,krow)*prho_rcp(1:kproma,:) ![kg-1]
     SELECT CASE(ninput_type_CCN)
     CASE(CCN_0SD3D)
             pascs(1:kproma,:) = 0._dp
             !UP: I'm conciously removing the MAX here, and seeing whether that
             !gives sensible results at all, but one may want to include the MAX
             !again later
     CASE DEFAULT
             pascs(1:kproma,:) = MAX(pascs(1:kproma,:), 10.E6_dp)
     END SELECT

     !-- Number, radius and std. dev. of aerosols available for cirrus freezing
     !   note: the number will be later reduced by the actual icnc in the cloud routine

     ld_het = .false.

     papnx(1:kproma,:,1) = prho(1:kproma,:) * pascs(1:kproma,:) ![1/m3]

     paprx(1:kproma,:,1) = 100._dp * zrwetas(1:kproma,:) ![cm]
     SELECT CASE(ninput_type_CCN)
     CASE(CCN_0SD3D)
             paprx(1:kproma,:,1) = 0._dp
     CASE DEFAULT
             paprx(1:kproma,:,1) = MAX(paprx(1:kproma,:,1), 0.05E-4_dp)
     END SELECT

     papsigx(1:kproma,:,1) = 1._dp

  END SUBROUTINE ccnclim_IN_setup

  !---------------------------------------------------------------------------
  !>
  !! @brief Routine to cleanup memory related to the CCN clim submodel
  !! 

  SUBROUTINE ccnclim_free_memory

    DEALLOCATE(zccn_co)
    DEALLOCATE(zccn_an)
    DEALLOCATE(zccn_pi)
    DEALLOCATE(zccn_tot_3d)

    IF (ninput_type_CCN == CCN_2D) THEN
       DEALLOCATE(zccn_du)
       DEALLOCATE(zccn_tot_2d)
    ENDIF     
  END SUBROUTINE ccnclim_free_memory

END MODULE mo_ccnclim
