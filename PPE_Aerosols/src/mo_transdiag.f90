!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_transdiag

!   Diagnose atmospheric energy transport. 
!
!   This module was created to diagnose horizontal energy 
!   transport online while ECHAM is running. The four 
!   components of the energy transport is calculated
!   in both the zonal and meridional directions. In 
!   addition hydrological transport of water vapor
!   cloud liquid water and cloud ice water is diagnosed.
!
!   The module consists of two subroutines. The first 
!   creates a separate output stream which should be
!   called from 'init_subm_memory' in the submodel
!   interface. The second does the diagnostic calculations
!   of energy transport components, and is called from
!   'vdiff_subm'. Results are output in 2D horizontal 
!   fields. Accumulation is done.
!
!   To activate the module one needs to set the switch
!   'ltransdiag' to .TRUE. in runctl.
!
!   Thorsten Mauritsen, 24/6-2010: Created
!   Thorsten Mauritsen, 3/11-2010: Added cloud water- and ice transports, 
!                                  fetch constants from mo_physical_constants.f90
!

USE mo_linked_list,    ONLY: t_stream 
USE mo_kind,           ONLY: wp

IMPLICIT NONE

PUBLIC :: construct_transdiag_stream

!-------------------------

! Output transdiag stream:

TYPE (t_stream),    PUBLIC, POINTER :: transdiag_stream

  REAL(wp),           PUBLIC, POINTER :: ptr_ju_pot(:,:)         ! Zonal transport of potential energy
  REAL(wp),           PUBLIC, POINTER :: ptr_ju_int(:,:)         ! Zonal transport of internal energy
  REAL(wp),           PUBLIC, POINTER :: ptr_ju_lat(:,:)         ! Zonal transport of latent energy
  REAL(wp),           PUBLIC, POINTER :: ptr_ju_kin(:,:)         ! Zonal transport of kinetic energy
  REAL(wp),           PUBLIC, POINTER :: ptr_ju_q(:,:)           ! Zonal transport of water vapor
  REAL(wp),           PUBLIC, POINTER :: ptr_ju_xl(:,:)          ! Zonal transport of cloud liquid
  REAL(wp),           PUBLIC, POINTER :: ptr_ju_xi(:,:)          ! Zonal transport of cloud ice

  REAL(wp),           PUBLIC, POINTER :: ptr_jv_pot(:,:)         ! Meridional transport of potential energy
  REAL(wp),           PUBLIC, POINTER :: ptr_jv_int(:,:)         ! Meridional transport of internal energy
  REAL(wp),           PUBLIC, POINTER :: ptr_jv_lat(:,:)         ! Meridional transport of latent energy
  REAL(wp),           PUBLIC, POINTER :: ptr_jv_kin(:,:)         ! Meridional transport of kinetic energy
  REAL(wp),           PUBLIC, POINTER :: ptr_jv_q(:,:)           ! Meridional transport of water vapor
  REAL(wp),           PUBLIC, POINTER :: ptr_jv_xl(:,:)          ! Meridional transport of cloud liquid
  REAL(wp),           PUBLIC, POINTER :: ptr_jv_xi(:,:)          ! Meridional transport of cloud ice

CONTAINS

!-----------------------------------------------
!-----------------------------------------------

  SUBROUTINE construct_transdiag_stream

    ! Here we create an output stream for transdiag.

    USE mo_exception,      ONLY: message
    USE mo_memory_base,    ONLY: new_stream, add_stream_element
    USE mo_linked_list,    ONLY: NETCDF

    IMPLICIT NONE

    call message('','Create: construct_transdiag_stream')
    CALL new_stream(transdiag_stream, 'trdiag',filetype=NETCDF, &
                    lpost=.TRUE.,lpout=.TRUE.,lrerun=.FALSE.,linit=.FALSE.)

    CALL add_stream_element(transdiag_stream, 'ju_pot', ptr_ju_pot, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., longname='Zonal potential energy transport', units='W/m')
    CALL add_stream_element(transdiag_stream, 'ju_int', ptr_ju_int, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., longname='Zonal internal energy transport', units='W/m')
    CALL add_stream_element(transdiag_stream, 'ju_lat', ptr_ju_lat, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., longname='Zonal latent energy transport', units='W/m')
    CALL add_stream_element(transdiag_stream, 'ju_kin', ptr_ju_kin, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., longname='Zonal kinetic energy transport', units='W/m')
    CALL add_stream_element(transdiag_stream, 'ju_q', ptr_ju_q, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., longname='Zonal water vapor transport', units='kg/m/s')
    CALL add_stream_element(transdiag_stream, 'ju_xl', ptr_ju_xl, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., longname='Zonal cloud liquid water transport', units='kg/m/s')
    CALL add_stream_element(transdiag_stream, 'ju_xi', ptr_ju_xi, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., longname='Zonal cloud ice water transport', units='kg/m/s')

    CALL add_stream_element(transdiag_stream, 'jv_pot', ptr_jv_pot, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., &
                            longname='Meridional potential energy transport', &
                            units='W/m')
    CALL add_stream_element(transdiag_stream, 'jv_int', ptr_jv_int, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., &
                            longname='Meridional internal energy transport', &
                            units='W/m')
    CALL add_stream_element(transdiag_stream, 'jv_lat', ptr_jv_lat, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., &
                            longname='Meridional latent energy transport', &
                            units='W/m')
    CALL add_stream_element(transdiag_stream, 'jv_kin', ptr_jv_kin, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., &
                            longname='Meridional kinetic energy transport', &
                            units='W/m')
    CALL add_stream_element(transdiag_stream, 'jv_q', ptr_jv_q, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., &
                            longname='Meridional water vapor transport', &
                            units='kg/m/s')
    CALL add_stream_element(transdiag_stream, 'jv_xl', ptr_jv_xl, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., &
                            longname='Meridional cloud liquid water transport', &
                            units='kg/m/s')
    CALL add_stream_element(transdiag_stream, 'jv_xi', ptr_jv_xi, lpost=.TRUE., lrerun=.FALSE., &
                            laccu=.TRUE., &
                            longname='Meridional cloud ice water transport', &
                            units='kg/m/s')

  END SUBROUTINE construct_transdiag_stream

!-----------------------------------------------
!-----------------------------------------------

  SUBROUTINE transdiag(kproma,kbdim,klev,klevp1,krow,pum1,pvm1,ptm1,pqm1,pxlm1,pxim1,paphm1,pgeom1)

  ! Here energy transport components are calculated.
  !
  ! This subroutine is called from mo_submodel_interface.f90 from the 
  ! subroutine 'vdiff_subm'. This is necessary in order to have winds
  ! passed on. 

  USE mo_time_control,   ONLY: delta_time
  USE mo_physical_constants, ONLY: rgrav, alv, cpd    ! inverse gravity
                                                      ! latent heat for vaporisation
                                                      ! inverse specific heat of dry air,

  IMPLICIT NONE

  INTEGER,  INTENT(in) :: kproma                      ! geographic block number of locations
  INTEGER,  INTENT(in) :: kbdim                       ! geographic block maximum number of locations
  INTEGER,  INTENT(in) :: klev                        ! number of levels
  INTEGER,  INTENT(in) :: klevp1                      ! number of levels + 1
  INTEGER,  INTENT(in) :: krow                        ! geographic block number

  REAL(wp), INTENT(in) :: pum1     (kbdim,klev)       ! u-wind (t-dt)
  REAL(wp), INTENT(in) :: pvm1     (kbdim,klev)       ! v-wind (t-dt)
  REAL(wp), INTENT(in) :: ptm1     (kbdim,klev)       ! temperature (t-dt)
  REAL(wp), INTENT(in) :: pqm1     (kbdim,klev)       ! specific humidity (t-dt)
  REAL(wp), INTENT(in) :: pxlm1    (kbdim,klev)       ! cloud liquid water (t-dt)
  REAL(wp), INTENT(in) :: pxim1    (kbdim,klev)       ! cloud ice water (t-dt)
  REAL(wp), INTENT(in) :: paphm1   (kbdim,klevp1)     ! air pressure at layer interface (t-dt)
  REAL(wp), INTENT(in) :: pgeom1   (kbdim,klev)       ! geopotential (t-dt)

  REAL(wp)             :: dpr      (klev)             ! pressure thicknesses
  INTEGER              :: k,l


  !---------------------------------------------------

  do k=1,kproma

    ! Calculate layer pressure thicknesses for gridpoint:

    do l=1,klev
      dpr(l) = paphm1(k,l+1) - paphm1(k,l)
    end do

    ! Perform vertical integration and time accumulation of transport components:

    ptr_ju_pot(k,krow) = ptr_ju_pot(k,krow) + &
                         rgrav*sum(pum1(k,:)*pgeom1(k,:)*dpr)*delta_time
    ptr_ju_int(k,krow) = ptr_ju_int(k,krow) + &
                         rgrav*cpd*sum(pum1(k,:)*ptm1(k,:)*dpr)*delta_time
    ptr_ju_lat(k,krow) = ptr_ju_lat(k,krow) + &
                         rgrav*alv*sum(pum1(k,:)*pqm1(k,:)*dpr)*delta_time
    ptr_ju_kin(k,krow) = ptr_ju_kin(k,krow) + &
                         rgrav*0.5_wp*sum(pum1(k,:)* &
                         (pum1(k,:)*pum1(k,:)+pvm1(k,:)*pvm1(k,:))*dpr)*delta_time
    ptr_ju_q(k,krow)   = ptr_ju_q(k,krow) + &
                         rgrav*sum(pum1(k,:)*pqm1(k,:)*dpr)*delta_time
    ptr_ju_xl(k,krow)  = ptr_ju_xl(k,krow) + &
                         rgrav*sum(pum1(k,:)*pxlm1(k,:)*dpr)*delta_time
    ptr_ju_xi(k,krow)  = ptr_ju_xi(k,krow) + &
                         rgrav*sum(pum1(k,:)*pxim1(k,:)*dpr)*delta_time


    ptr_jv_pot(k,krow) = ptr_jv_pot(k,krow) + &
                         rgrav*sum(pvm1(k,:)*pgeom1(k,:)*dpr)*delta_time
    ptr_jv_int(k,krow) = ptr_jv_int(k,krow) + &
                         rgrav*cpd*sum(pvm1(k,:)*ptm1(k,:)*dpr)*delta_time
    ptr_jv_lat(k,krow) = ptr_jv_lat(k,krow) + &
                         rgrav*alv*sum(pvm1(k,:)*pqm1(k,:)*dpr)*delta_time
    ptr_jv_kin(k,krow) = ptr_jv_kin(k,krow) + &
                         rgrav*0.5_wp*sum(pvm1(k,:)* & 
                         (pum1(k,:)*pum1(k,:)+pvm1(k,:)*pvm1(k,:))*dpr)*delta_time
    ptr_jv_q(k,krow)   = ptr_jv_q(k,krow) + &
                         rgrav*sum(pvm1(k,:)*pqm1(k,:)*dpr)*delta_time
    ptr_jv_xl(k,krow)  = ptr_jv_xl(k,krow) + &
                         rgrav*sum(pvm1(k,:)*pxlm1(k,:)*dpr)*delta_time
    ptr_jv_xi(k,krow)  = ptr_jv_xi(k,krow) + &
                         rgrav*sum(pvm1(k,:)*pxim1(k,:)*dpr)*delta_time


  end do

  END SUBROUTINE transdiag

!-----------------------------------------------
!-----------------------------------------------

END MODULE mo_transdiag
