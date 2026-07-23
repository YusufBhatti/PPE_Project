!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_emi_fire.f90
!!
!! \brief
!! Module for handling of fire emissions
!!
!! \author M. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (2010-02)
!!   -# G. Frontoso (C2SM-ETHZ) - Fire emissions injection accounts for PBL height (2012-07)
!!                                (M. Val Martin, 2010), see #161 on redmine
!!   -# A. Veira (MPI-M) - Modification of vertical emission distribution:
!!                         Constant mass mixing ratio throughout PBL instead of simple
!!                         pro-rata distribution to all PBL model layers (see #310)
!! \limitations
!! None
!!
!! \details
!! - vertical distribution of fire emissions
!! - planned extension: derive emission factors based on MCE and read only
!!   CO and CO2 data from file.
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
MODULE mo_hammoz_emi_fire

  USE mo_kind,       ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: distribute_emi_fire

  CONTAINS

!gf  SUBROUTINE distribute_emi_fire(kproma, kbdim, klev, zbc2d, zbc3d)
  SUBROUTINE distribute_emi_fire(kproma,kbdim,klev,krow,ihpbl,zbc2d,zbc3d)

  USE mo_time_control,     ONLY: lstart
  USE mo_vphysc,           ONLY: vphysc
  USE mo_control,          ONLY: vct, nvclev

  INTEGER, INTENT(in)     :: kproma
  INTEGER, INTENT(in)     :: kbdim
  INTEGER, INTENT(in)     :: klev
  INTEGER, INTENT(in)     :: krow
  INTEGER, INTENT(in)     :: ihpbl(kbdim)

  REAL(dp), INTENT(in)    :: zbc2d(kbdim)
  REAL(dp), INTENT(inout) :: zbc3d(kbdim, klev)

! Local variables

  INTEGER   :: jl,jk,ilower,iupper
  INTEGER   :: lev_blw_pbl
  REAL(dp)  :: zfac_k(klev),zbc2d_abv1(kbdim),zbc2d_abv2(kbdim),zbc2d_blw(kbdim)
  REAL(dp)  :: zdepth, p_pblh
  REAL(dp)  :: as(nvclev), bs(nvclev), ph(nvclev), zf(klev)

! This is to by-pass the fact that the variable grheightm1 is not available at the first time step. 
! It has been taken from mo_hammoz_emi_volcano.f90

  as(:) = vct(1:nvclev)
  bs(:) = vct(nvclev+1:2*nvclev)
  ph(:) = as(:)+bs(:)*1.e5_dp    ! nominal pressure values
  ph(1) = 1.e-19_dp

!!$! convert to "height above surface"
  DO jk=1,klev
     zf(jk) = 8000._dp * (LOG(ph(klev+1)/1.e5_dp) - LOG(ph(jk)/1.e5_dp))
  END DO

  ilower = klev       ! lower level

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Amount of fire emissions below and above top of PBL

  zbc2d_abv1(1:kproma)=0.17_dp*zbc2d(1:kproma)
  zbc2d_abv2(1:kproma)=0.08_dp*zbc2d(1:kproma)
  zbc2d_blw(1:kproma)=0.75_dp*zbc2d(1:kproma)

! Vertical profile

  DO jl=1,kproma
     iupper=ihpbl(jl) ! upper level - top of the PBL

! Compute PBL depth in m for the 1st timestep
     zdepth=zf(iupper)
      IF (.NOT. lstart) THEN
         zdepth = 0._dp
         DO jk=ilower,iupper,-1
            zdepth=zdepth+vphysc%grheightm1(jl,jk,krow)
         ENDDO
! Height at the middle of the level
        zdepth=zdepth-vphysc%grheightm1(jl,iupper,krow)/2.
      ENDIF

! Check if the depth of the PBL is higher than 4km
! Emit only within the PBL

      IF (zdepth .GT. 4000._dp) THEN
         zbc2d_abv1(jl)=0._dp
         zbc2d_abv2(jl)=0._dp
         zbc2d_blw(jl)=zbc2d(jl)
      ENDIF

     lev_blw_pbl=ilower-ihpbl(jl)+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Distribute emission in PBL with constant mass mixing ratio (!AV #310)
     zfac_k(:) = 0._dp
     p_pblh = ph(nvclev)-ph(iupper)  ! pressure difference surface - top of the pbl level
     DO jk=ilower, iupper, -1
         zfac_k(jk) = (ph(jk+1)-ph(jk))/p_pblh   
         zbc3d(jl,jk)= zfac_k(jk)*zbc2d_blw(jl)
     ENDDO

! Distribute emissions in the FT
     zbc3d(jl,iupper-1)     = zbc2d_abv1(jl)
     zbc3d(jl,iupper-2)     = zbc2d_abv2(jl)
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Old way of fire injections

!  zbc3d(1:kproma,klev)   = 0.33 * zbc2d(1:kproma)
!  zbc3d(1:kproma,klev-1) = 0.26 * zbc2d(1:kproma)
!  zbc3d(1:kproma,klev-2) = 0.19 * zbc2d(1:kproma)
!  zbc3d(1:kproma,klev-3) = 0.15 * zbc2d(1:kproma)
!  zbc3d(1:kproma,klev-4) = 0.07 * zbc2d(1:kproma)

  ! note: zbc2d will be reset to zero in calling routine.

  END SUBROUTINE distribute_emi_fire

END MODULE mo_hammoz_emi_fire
