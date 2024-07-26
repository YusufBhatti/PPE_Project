!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_memory_cfdiag
!
! module to write out 3-D radiative and convective mass fluxes 
!
! CFMIP
!
! uses separate output files, so that it can more 
! conveniently be switched on and off during a run as requested by CFMIP
! 
!------------------------------------------------------------------------------
  !
  ! Modules used
  !

  USE mo_kind,         ONLY : wp
  USE mo_linked_list,  ONLY : t_stream, HYBRID_H
  USE mo_memory_base,  ONLY : new_stream, delete_stream, add_stream_element, &
                              default_stream_setting
  USE mo_filename,     ONLY : ftyp => out_filetype
  USE mo_time_control, ONLY : trigrad
  USE mo_exception,    ONLY : finish
  USE mo_time_event,   ONLY : io_time_event

  USE mo_mpi,          ONLY : p_parallel, p_bcast, p_io
  USE mo_namelist,     ONLY : open_nml, position_nml, POSITIONED

  IMPLICIT NONE
!------------------------------------------------------------------------------

  PRIVATE

  PUBLIC :: construct_cfdiag ! construct the cfdiag table
  PUBLIC :: destruct_cfdiag  ! destruct  the cfdiag table
  PUBLIC :: setup_cfdiag     ! read namelist

  PUBLIC :: cfdiag           ! the cfdiag table

  PUBLIC :: calc_avps, calc_cfdiag

  LOGICAL, PUBLIC :: locfdiag = .FALSE.

  ! declaration of predefined fields within this module 

  REAL(wp), PUBLIC :: odrad
  REAL(wp), PUBLIC :: omrad

   INTEGER, PUBLIC :: ncfdiacall = 1
!$OMP THREADPRIVATE(ncfdiacall)

! radiative fluxes
  REAL(wp), POINTER, PUBLIC :: rlu(:,:,:)
  REAL(wp), POINTER, PUBLIC :: rsu(:,:,:)
  REAL(wp), POINTER, PUBLIC :: rld(:,:,:)
  REAL(wp), POINTER, PUBLIC :: rsd(:,:,:)
  REAL(wp), POINTER, PUBLIC :: rlucs(:,:,:)
  REAL(wp), POINTER, PUBLIC :: rsucs(:,:,:)
  REAL(wp), POINTER, PUBLIC :: rldcs(:,:,:)
  REAL(wp), POINTER, PUBLIC :: rsdcs(:,:,:)

  REAL(wp), POINTER, PUBLIC :: irlu(:,:,:)
  REAL(wp), POINTER, PUBLIC :: irsu(:,:,:)
  REAL(wp), POINTER, PUBLIC :: irld(:,:,:)
  REAL(wp), POINTER, PUBLIC :: irsd(:,:,:)
  REAL(wp), POINTER, PUBLIC :: irlucs(:,:,:)
  REAL(wp), POINTER, PUBLIC :: irsucs(:,:,:)
  REAL(wp), POINTER, PUBLIC :: irldcs(:,:,:)
  REAL(wp), POINTER, PUBLIC :: irsdcs(:,:,:)

  REAL(wp), POINTER, PUBLIC :: srsu(:,:,:)
  REAL(wp), POINTER, PUBLIC :: srsd(:,:,:)
  REAL(wp), POINTER, PUBLIC :: srsucs(:,:,:)
  REAL(wp), POINTER, PUBLIC :: srsdcs(:,:,:)

! convective mass fluxes
  REAL(wp), POINTER, PUBLIC :: mcu(:,:,:)
  REAL(wp), POINTER, PUBLIC :: mcd(:,:,:)
  REAL(wp), POINTER, PUBLIC :: mc(:,:,:)
  REAL(wp), POINTER, PUBLIC :: smc(:,:,:)
  REAL(wp), POINTER, PUBLIC :: dmc(:,:,:)

  REAL(wp), POINTER, PUBLIC :: imc(:,:,:)

  REAL(wp), POINTER, PUBLIC :: avps(:,:)
  
! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: cfdiag

  INCLUDE 'cfdiag_ctl.inc'
  
CONTAINS
!------------------------------------------------------------------------------
  SUBROUTINE construct_cfdiag

    ! fraction of day between time steps (needed for accumulated variables)
     IF ( trigrad%unit .EQ. 'hours' ) THEN
       odrad =  3600._wp * REAL(trigrad%counter, wp) / 86400._wp
     ELSE
       CALL finish('mo_memory_cfdiag ','expecting different trigrad unit')
     END IF


    ! construct the cfdiag table
    !
    ! all information specific to this table is set in this subroutine


    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

 
    ! do not change output frequency
    CALL new_stream (cfdiag, 'cfdiag', ftyp, rest_suf= '_echam', post_suf='_cfdiag', &
                        interval=io_time_event(1,'months','first',0) )


    CALL default_stream_setting (cfdiag, lrerun=.TRUE., lpost=.TRUE., &
                                   contnorest=.TRUE., &
                                  leveltype = HYBRID_H, reset = 1.e-50_wp )


! averaged

    CALL add_stream_element (cfdiag,'avps',  avps       ,code=134,          &
          longname='surface pressure'  , laccu=.TRUE.  ,units='Pa'      )
    CALL add_stream_element (cfdiag,'rlu',   rlu,  code=1 , laccu=.TRUE.,    &
          longname='upwelling_longwave_flux_in_air'  ,units='W/m2')
    CALL add_stream_element (cfdiag,'rsu',   rsu,  code=2 , laccu=.TRUE.,  &
          longname='upwelling_shortwave_flux_in_air'   ,units='W/m2')
    CALL add_stream_element (cfdiag,'rld',   rld,  code=3 , laccu=.TRUE.,  &
          longname='downwelling_longwave_flux_in_air'  ,units='W/m2')
    CALL add_stream_element (cfdiag,'rsd',   rsd,  code=4 , laccu=.TRUE.,  &
       longname='downwelling_shortwave_flux_in_air' ,units='W/m2')
    CALL add_stream_element (cfdiag,'rlucs', rlucs,   code=5 , laccu=.TRUE., &
          longname='upwelling_longwave_flux_in_air_assuming_clear_sky' , &
          units='W/m2')
    CALL add_stream_element (cfdiag,'rsucs', rsucs,  code=6 , laccu=.TRUE., &
         longname='upwelling_shortwave_flux_in_air_assuming_clear_sky' , &
         units='W/m2')
    CALL add_stream_element (cfdiag,'rldcs', rldcs, code=7 ,  laccu=.TRUE., &
          longname='downwelling_longwave_flux_in_air_assuming_clear_sky', &
          units='W/m2')
    CALL add_stream_element (cfdiag,'rsdcs',  rsdcs, code=8 , laccu=.TRUE., &
          longname='downwelling_shortwave_flux_in_air_assuming_clear_sky', &
          units='W/m2')
    CALL add_stream_element (cfdiag,'mc',  mc, code=9 , laccu=.FALSE., &
          longname='atmosphere_net_upward_convective_mass_flux', &
          units='kg/s/m2')
    CALL add_stream_element (cfdiag,'mcu',  mcu, code=10 , laccu=.FALSE., &
          longname='atmosphere_updraft_convective_mass_flux', &
          units='kg/s/m2')
    CALL add_stream_element (cfdiag,'mcd',  mcd, code=11 , laccu=.FALSE., &
          longname='atmosphere_downdraft_convective_mass_flux', &
          units='kg/s/m2')
    CALL add_stream_element (cfdiag,'dmc',  dmc, code=12 , laccu=.FALSE., &
          longname='atmosphere_net_upward_deep_convective_mass_flux', &
          units='kg/s/m2')
    CALL add_stream_element (cfdiag,'smc',  smc, code=13 , laccu=.FALSE., &
          longname='atmosphere_net_upward_shallow_convective_mass_flux', &
          units='kg/s/m2')


! save 

    CALL add_stream_element (cfdiag,'srsu',   srsu,    units='W/m2', &
                               lpost=.FALSE., lrerun=.TRUE. )
    CALL add_stream_element (cfdiag,'srsd',   srsd,    units='W/m2', &
                               lpost=.FALSE., lrerun=.TRUE. )
    CALL add_stream_element (cfdiag,'srsucs', srsucs,  units='W/m2', &
                               lpost=.FALSE., lrerun=.TRUE. )
    CALL add_stream_element (cfdiag,'srsdcs', srsdcs,  units='W/m2', &
                               lpost=.FALSE., lrerun=.TRUE. )

! instantaneous

    CALL add_stream_element (cfdiag,'irlu',   irlu,    units='W/m2', &
                               lpost=.FALSE., lrerun=.TRUE. )
    CALL add_stream_element (cfdiag,'irsu',   irsu,    units='W/m2', &
                               lpost=.FALSE., lrerun=.FALSE. )
    CALL add_stream_element (cfdiag,'irld',   irld,    units='W/m2', &
                               lpost=.FALSE., lrerun=.TRUE. )
    CALL add_stream_element (cfdiag,'irsd',   irsd,    units='W/m2', &
                               lpost=.FALSE., lrerun=.FALSE. )
    CALL add_stream_element (cfdiag,'irlucs', irlucs,  units='W/m2', &
                               lpost=.FALSE., lrerun=.TRUE. )
    CALL add_stream_element (cfdiag,'irsucs', irsucs,  units='W/m2', &
                               lpost=.FALSE., lrerun=.FALSE. )
    CALL add_stream_element (cfdiag,'irldcs', irldcs,  units='W/m2', &
                               lpost=.FALSE., lrerun=.TRUE. )
    CALL add_stream_element (cfdiag,'irsdcs', irsdcs,  units='W/m2', &
                               lpost=.FALSE., lrerun=.FALSE. )
    CALL add_stream_element (cfdiag,'imc', imc,  units='kg/s/m2', &
                               lpost=.FALSE., lrerun=.FALSE. )


  END SUBROUTINE construct_cfdiag
!------------------------------------------------------------------------------
  SUBROUTINE destruct_cfdiag

    CALL delete_stream (cfdiag)
  
  END SUBROUTINE destruct_cfdiag
!------------------------------------------------------------------------------
  SUBROUTINE setup_cfdiag

    INTEGER :: ierr, inml, iunit

      !
      ! --- Read NAMELIST
      ! 
      inml = open_nml('namelist.echam')
      iunit = position_nml ('CFDIAGCTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
        READ (iunit, cfdiagctl)
      END SELECT

     IF (p_parallel) THEN
        CALL p_bcast (locfdiag, p_io)
     END IF 

  END  SUBROUTINE setup_cfdiag
!------------------------------------------------------------------------------

  SUBROUTINE calc_avps

      USE mo_time_control,      ONLY: delta_time
      USE mo_decomposition, ONLY: ldc=>local_decomposition
      USE mo_memory_g3b,    ONLY: aps

     INTEGER :: nproma, ngpblks
     INTEGER :: jrow, jl

     ngpblks = ldc % ngpblks

     DO jrow = 1, ngpblks        
       IF ( jrow == ngpblks ) THEN
         nproma = ldc% npromz
       ELSE
         nproma = ldc% nproma
       END IF

       DO jl=1,nproma      
           avps(jl,jrow) = avps(jl,jrow) + aps(jl,jrow) * delta_time
       END DO
     END DO

  END SUBROUTINE calc_avps

!------------------------------------------------------------------------------

  SUBROUTINE calc_cfdiag( kproma, kbdim, pi0, krow )

    USE mo_time_control,   ONLY : delta_time

    INTEGER, INTENT(IN ) :: kproma, kbdim, krow 
    REAL(wp), INTENT(IN )  :: pi0(kbdim)

    INTEGER :: jl

     DO jl=1,kproma      
       irsu(jl,:,krow)   = pi0(jl) * srsu(jl,:,krow) 
       irsd(jl,:,krow)   = pi0(jl) * srsd(jl,:,krow) 
       irsucs(jl,:,krow) = pi0(jl) * srsucs(jl,:,krow) 
       irsdcs(jl,:,krow) = pi0(jl) * srsdcs(jl,:,krow) 

       rlu(jl,:,krow)   =  rlu(jl,:,krow) + irlu(jl,:,krow) * delta_time
       rsu(jl,:,krow)   =  rsu(jl,:,krow) + irsu(jl,:,krow) * delta_time
       rld(jl,:,krow)   =  rld(jl,:,krow) + irld(jl,:,krow) * delta_time
       rsd(jl,:,krow)   =  rsd(jl,:,krow) + irsd(jl,:,krow) * delta_time
       rlucs(jl,:,krow) =  rlucs(jl,:,krow) + irlucs(jl,:,krow) * delta_time
       rsucs(jl,:,krow) =  rsucs(jl,:,krow) + irsucs(jl,:,krow) * delta_time
       rldcs(jl,:,krow) =  rldcs(jl,:,krow) + irldcs(jl,:,krow) * delta_time
       rsdcs(jl,:,krow) =  rsdcs(jl,:,krow) + irsdcs(jl,:,krow) * delta_time

     END DO

  END SUBROUTINE  calc_cfdiag


END MODULE mo_memory_cfdiag
