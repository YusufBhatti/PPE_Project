!! \filename
!! mo_ham_aeroclim.f90
!!
!! \brief
!! 
!! \author Sylvaine Ferrachat (ETHZ)
!!
!! \responsible_coder
!! Sylvaine Ferrachat, sylvaine.ferrachat@env.ethz.ch
!!
!! \revision_history
!!   -# S. Ferrachat (ETHZ) - 
!!
!! \limitations
!! None
!!
!! \details
!!
!! \bibliographic_references
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

MODULE mo_ham_aeroclim

  USE mo_kind,               ONLY: wp

  USE mo_boundary_condition, ONLY: bc_nml, bc_define, p_bcast_bc,  &
                                   BC_REPLACE, BC_EVERYWHERE, BC_ALTITUDE, BC_LEVEL, BC_BOTTOM
  USE mo_external_field_processor, ONLY: EF_FILE, EF_LONLAT, EF_3D, EF_IGNOREYEAR, EF_NOINTER, &
                                         EF_VALUE, EF_UNDEFINED

  !dbg USE mo_ham,                ONLY: naerocomp
  USE mo_tracdef,            ONLY: trlist, ntrac, t_trinfo
  USE mo_boundary_condition, ONLY: bc_nml, bc_define, p_bcast_bc,  &
                                   BC_REPLACE, BC_EVERYWHERE, BC_ALTITUDE, BC_LEVEL, BC_BOTTOM, &
                                   bc_apply
  USE mo_external_field_processor, ONLY: EF_FILE, EF_LONLAT, EF_3D, EF_IGNOREYEAR, EF_NOINTER, &
                                         EF_VALUE, EF_UNDEFINED


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_aeroclim_submodel, set_prescribed_tracer_list, set_prescribed_tracer_bc, & 
            reset_tracer, prescribe_tracer

  !dbg INTEGER                   :: nmax_presc_trc = naerocomp
  INTEGER, PARAMETER        :: nmax_presc_trc = 100
  INTEGER                   :: n_presc_trc
  INTEGER                   :: map_presc_trc(nmax_presc_trc)
  INTEGER, ALLOCATABLE      :: ibc_presc_trc(:)
  TYPE(bc_nml), ALLOCATABLE :: bc_presc_trc(:)
  CHARACTER(LEN=32)         :: presc_trc_nml(nmax_presc_trc)

CONTAINS

  !---------------------------------------------------------------------------
  !>
  !! @brief Initializes the HAM climatology submodel
  !!
  !! @remarks This subroutine reads into the hamclim namelist
  !! and constructs the hamclim stream
  !! It also sets all the relevant boundary conditions
  !! 

  SUBROUTINE init_aeroclim_submodel

    USE mo_memory_base,        ONLY: new_stream, add_stream_element,     &
                                     default_stream_setting, add_stream_reference, AUTO
    USE mo_submodel,           ONLY: print_value
    USE mo_mpi,                ONLY: p_io, p_parallel_io, p_parallel, p_bcast
    USE mo_namelist,           ONLY: open_nml, position_nml, POSITIONED, &
                                     LENGTH_ERROR, READ_ERROR
    USE mo_filename,           ONLY: out_filetype
    USE mo_exception,          ONLY: finish, message, message_text, em_info, em_error

    INTEGER :: ierr, inml, iunit

    INCLUDE 'ham_aeroclimctl.inc'

    !-- Set defaults
    presc_trc_nml(:) = ''
    presc_trc_nml(1) = 'all'

    !-- Namelist reading:
    IF (p_parallel_io) THEN
       inml = open_nml('namelist.echam')
       iunit = position_nml ('HAM_AEROCLIMCTL', inml, status=ierr)
       SELECT CASE (ierr)
       CASE (POSITIONED)
          READ (iunit, ham_aeroclimctl)
      CASE (LENGTH_ERROR)
        CALL finish ('init_ham_aeroclim_submodel', &
             'length error in namelist hamclimctl')
      CASE (READ_ERROR)
        CALL finish ('init_ham_aeroclim_submodel', &
             'read error in namelist.echam')
      END SELECT

    END IF

    IF (p_parallel) THEN
       CALL p_bcast (presc_trc_nml, p_io)
    ENDIF

    !ToDo Implement security checks

  END SUBROUTINE init_aeroclim_submodel

  SUBROUTINE set_prescribed_tracer_list
     
     INTEGER :: jtrac
     TYPE(t_trinfo) ,POINTER :: ti

     map_presc_trc(:) = -1
     n_presc_trc = 0

     !ToDo Extend possibilities of tracer selection from namelist
     ! So far: only 'all' is possible
     !DBG buggy IF (ANY(TRIM(presc_trc_nml)) == 'all') THEN
     !>>SF tmp
     IF (TRIM(presc_trc_nml(1)) == 'all') THEN
     !<<SF tmp
         DO jtrac=1,ntrac
             ti => trlist% ti (jtrac)
             IF ( TRIM(ti %modulename) == 'HAM' .AND.                     &
                 (TRIM(ti %basename) /= 'CDNC' .AND. TRIM(ti %basename) /= 'ICNC') &
                ) THEN
                n_presc_trc = n_presc_trc + 1
                map_presc_trc(n_presc_trc) = jtrac
             ENDIF
         ENDDO
     ENDIF
  END SUBROUTINE set_prescribed_tracer_list

  SUBROUTINE set_prescribed_tracer_bc

     INTEGER :: jp, jtrac
     TYPE(t_trinfo) ,POINTER :: ti

     ALLOCATE(ibc_presc_trc(n_presc_trc))
     ALLOCATE(bc_presc_trc(n_presc_trc))

     DO jp=1,n_presc_trc
         jtrac = map_presc_trc(jp)
         ti => trlist% ti (jtrac)

         bc_presc_trc(jp) %ef_type        = EF_FILE
         bc_presc_trc(jp) %ef_template    = 'presc_tracer_%Y4%M2.nc'
         bc_presc_trc(jp) %ef_timedef     = EF_IGNOREYEAR
         bc_presc_trc(jp) %ef_interpolate = EF_NOINTER
         bc_presc_trc(jp) %ef_geometry    = EF_3D
         bc_presc_trc(jp) %ef_actual_unit = ti %units
         bc_presc_trc(jp) %ef_varname     = ti %fullname

         ibc_presc_trc(jp) = bc_define(ti %longname, bc_presc_trc(jp), 3, .TRUE.)
     ENDDO

  END SUBROUTINE set_prescribed_tracer_bc

  SUBROUTINE reset_tracer

     USE mo_tracdef,   ONLY: ON, OFF

     INTEGER :: jp, jtrac
     TYPE(t_trinfo) ,POINTER :: ti

     DO jp=1,n_presc_trc
         jtrac = map_presc_trc(jp)
         ti => trlist% ti (jtrac)

         ti %ninit   = -1
         ti %ntran   = OFF
         ti %nconv   = OFF
         ti %nvdiff  = OFF
         ti %nint    = OFF
         ti %ndrydep = OFF
         ti %nwetdep = OFF
         ti %nsedi   = OFF
     ENDDO

  END SUBROUTINE reset_tracer

  SUBROUTINE prescribe_tracer(kproma, kbdim, klev, krow, ktrac, pxtm1, pxtte)
      INTEGER,  INTENT(in)    :: kproma                 ! geographic block number of locations
      INTEGER,  INTENT(in)    :: kbdim                  ! geographic block maximum number of locations 
      INTEGER,  INTENT(in)    :: klev                   ! number of levels
      INTEGER,  INTENT(in)    :: ktrac                  ! number of tracers
      INTEGER,  INTENT(in)    :: krow                   ! geographic block number
    
      REAL(wp), INTENT(inout) :: pxtm1 (kbdim,klev,ktrac) ! tracer mass/number mixing ratio (t-dt)
      REAL(wp), INTENT(inout) :: pxtte (kbdim,klev,ktrac) ! tracer mass/number mixing ratio tendency

      INTEGER :: jp, jtrac
      TYPE(t_trinfo) ,POINTER :: ti

      DO jp=1,n_presc_trc
          jtrac = map_presc_trc(jp)
          ti => trlist% ti (jtrac)

          CALL bc_apply(ibc_presc_trc(jp), kproma, krow, pxtm1(1:kproma,:,jtrac))

          pxtte(1:kproma,:,jtrac) = 0._wp
      ENDDO

  END SUBROUTINE prescribe_tracer

END MODULE mo_ham_aeroclim
