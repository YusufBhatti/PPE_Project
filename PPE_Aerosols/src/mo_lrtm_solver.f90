!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to compute longwave radiative transfer (no scattering)
!!
!! @remarks
!!   This module contains routines that compute longwave radiative transfer. 
!!   The algorithms follow those used by RRTMG but have been simplified and vectorized
!!   accross columns. 
!!
!! @author Robert Pincus, U. Colorado, visiting MPI-M, Hamburg (2011-07)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Major segments of this code combines and rewrites (for the ICON standard) 
!!   code previously provided by AER and copyrighed by them.  The authors of the
!!   original AER code are: Eli J. Mlawer, Jennifer S. Delamere, Michael J. 
!!   Iacono and Shepard A. Clough with acknowledgments to Steven J. Taubman, 
!!   Karen Cady-Pereira, Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert 
!!   Bergstrom. The rewrites were designed to better interface with the structure
!!   of the ICON family of models and elements of the ICON programming standard.
!!
!! @par Copyright
!!   The AER copyright
!!   on the original code is as follows: Copyright 2002-2009, Atmospheric and
!!   Environmental Research, Inc. (AER). This software may be used, copied, or
!!   redistributed as long as it is not sold and this copyright notice is
!!   reproduced on each copy made.  This model is provided as is without any
!!   express or implied warranties. (http://www.rtweb.aer.com/)               
!! 
!
MODULE mo_lrtm_solver
  USE mo_kind,             ONLY : wp
  USE mo_math_constants,   ONLY : pi
  USE mo_rrtm_params,      ONLY : nbndlw
  USE mo_rad_fastmath,     ONLY : transmit, tautrans 

  IMPLICIT NONE

  REAL(wp), PARAMETER :: fluxfac = 2.0e+04_wp * pi
  
CONTAINS
  ! -------------------------------------------------------------------------------
   SUBROUTINE lrtm_solver(kproma, kbdim,   klev      , &
     & tau       , layPlnk, levPlnk, weights, secdiff, &
     & surfPlanck, surfEmis   ,                        & 
     & fluxUp    , fluxDn)
    !
    ! Compute IR (no scattering) radiative transfer for a set of columns 
    !   Based on AER code RRTMG_LW_RTNMC, including approximations used there
    ! Layers are ordered from botton to top (i.e. tau(1) is tau in lowest layer)
    ! Computes all-sky RT given a total optical thickness in each layer
    !
    INTEGER,  INTENT(in) :: &
      & kproma,             & !< Number of columns 
      & kbdim,              & !< Maximum number of columns as declared in calling (sub)program
      & klev                  !< number of layers (one fewer than levels) 
      
    REAL(wp), INTENT(in) :: & !< dimension (kbdim, klev) 
      & tau(kbdim,klev),           & !< Longwave optical thickness
      & layPlnk(kbdim,klev),       & !< Planck function at layer centers
      & weights(kbdim,klev)          !< Fraction of total Planck function for this g-point 
      
    REAL(wp), INTENT(in) :: &
      & levPlnk(kbdim, 0:klev)        !< Planck function at layer edges, level i is the top of layer i
      
    REAL(wp), INTENT(in) :: & !< dimension (kbdim)
      & surfPlanck(kbdim),      & !< Planck function at surface
      & surfEmis(kbdim),        & !< Surface emissivity
      & secdiff(kbdim)            !< secant of integration angle - depends on band, column water vapor
      
    REAL(wp), INTENT(out) :: & !< dimension (kbdim, 0:klev)
      & fluxUp(kbdim, 0:klev),       & !< Fluxes at the interfaces 
      & fluxDn(kbdim, 0:klev) 
                                         
    ! -----------
    INTEGER  :: &
      & jk     !< Loop index for layers
      
    REAL(wp) ::               &
      & trans(kbdim,klev),   & !< Layer transmissivity
      & tfn(kbdim),          & !< TFN_TBL 
		!< Tau transition function; i.e. the transition of the Planck
		!< function from that for the mean layer temperature to that for
		!< the layer boundary temperature as a function of optical depth.
		!< The "linear in tau" method is used to make the table.
      & dPlnkup(kbdim,klev), & !< Upward derivative of Planck function 
      & dPlnkdn(kbdim,klev), & !< Downward derivative of Planck function 
      & bbdn(kbdim,klev),    & !< Interpolated downward emission 
      & bbup(kbdim,klev),    & !< Interpolated upward emission 
      & odepth(kbdim,klev)   !< Effective IR optical depth of layer

    REAL(wp) ::                &
      & rad_dn(kbdim,0:klev), & !< Radiance down at propagation angle
      & rad_up(kbdim,0:klev)    !< Radiance down at propagation angle

    ! This secant and weight corresponds to the standard diffusivity 
    ! angle.  The angle is redefined for some bands.
    REAL(wp), PARAMETER :: wtdiff = 0.5_wp, rec_6 = 0.166667_wp
    ! -----------
    
    !    
    ! 1.0 Initial preparations 
    ! Weight optical depth by 1/cos(diffusivity angle), which depends on band 
    ! This will be used to compute layer transmittance
    ! -----
    
!IBM* ASSERT(NODEPS)
    DO jk = 1, klev
      odepth(1:kproma,jk) = max(0._wp, secdiff(1:kproma) * tau(1:kproma,jk))
    END DO 
    !
    ! 2.0 Radiative transfer
    !
    ! -----
    !
    ! Plank function derivatives and total emission for linear-in-tau approximation
    !
!IBM* ASSERT(NODEPS)
    DO jk = 1, klev
      tfn(1:kproma) = tautrans(odepth(:,jk), kproma)
      dPlnkup(1:kproma,jk) = levPlnk(1:kproma,jk  ) - layPlnk(1:kproma,jk)
      dPlnkdn(1:kproma,jk) = levPlnk(1:kproma,jk-1) - layPlnk(1:kproma,jk)

      bbup(1:kproma,jk) = weights(1:kproma,jk) * (layPlnk(1:kproma,jk) + dPlnkup(1:kproma,jk) * tfn(1:kproma))
      bbdn(1:kproma,jk) = weights(1:kproma,jk) * (layPlnk(1:kproma,jk) + dPlnkdn(1:kproma,jk) * tfn(1:kproma))
    END DO         
    ! -----
    ! 2.1 Downward radiative transfer
    !
    ! Level 0 is closest to the ground
    ! 
    rad_dn(:, klev) = 0. ! Upper boundary condition - no downwelling IR

    DO jk = klev, 1, -1
      trans(1:kproma,jk) = transmit(odepth(:,jk), kproma)
      ! RHS is a rearrangment of rad_dn(:,jk) * (1._wp - trans(:,jk)) + trans(:,jk) * bbdn(:)
      rad_dn(1:kproma,jk-1) = rad_dn(1:kproma,jk) + (bbdn(1:kproma,jk) - rad_dn(1:kproma,jk)) * trans(1:kproma,jk)
    END DO 
    !
    ! 2.2 Surface contribution, including reflection
    !
    rad_up(1:kproma, 0) = weights(1:kproma, 1) * surfEmis(1:kproma) * surfPlanck(1:kproma) &
                 + (1._wp - surfEmis(1:kproma)) * rad_dn(1:kproma, 0)
    
    !
    ! 2.3 Upward radiative transfer
    !
    DO jk = 1, klev
      rad_up(1:kproma,jk) = rad_up(1:kproma,jk-1) * (1._wp - trans(1:kproma,jk)) + trans(1:kproma,jk) * bbup(1:kproma,jk)
    END DO 
    
    !
    ! 3.0 Covert intensities at diffusivity angles to fluxes
    !
    ! -----
    fluxUp(1:kproma, 0:klev) = rad_up(1:kproma,:) * wtdiff * fluxfac
    fluxDn(1:kproma, 0:klev) = rad_dn(1:kproma,:) * wtdiff * fluxfac
    
  END SUBROUTINE lrtm_solver
  ! -------------------------------------------------------------------------------
  ELEMENTAL FUNCTION find_secdiff(iband, pwvcm) 
    INTEGER,  INTENT(in) :: &
      & iband !< RRTMG LW band number
    REAL(wp), INTENT(in) :: &
      & pwvcm !< Precipitable water vapor (cm) 
    REAL(wp)             :: &
      & find_secdiff

    ! Compute diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
    ! and 1.80) as a function of total column water vapor.  The function
    ! has been defined to minimize flux and cooling rate errors in these bands
    ! over a wide range of precipitable water values.
    REAL(wp), DIMENSION(nbndlw), PARAMETER :: &
      a0 = (/ 1.66_wp,  1.55_wp,  1.58_wp,  1.66_wp, 1.54_wp, 1.454_wp,  1.89_wp,  1.33_wp,     &
              1.668_wp, 1.66_wp,  1.66_wp,  1.66_wp, 1.66_wp,  1.66_wp,  1.66_wp,  1.66_wp /),  &
      a1 = (/ 0.00_wp,  0.25_wp,  0.22_wp,  0.00_wp,  0.13_wp, 0.446_wp, -0.10_wp,  0.40_wp,    &
             -0.006_wp, 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /), &
      a2 = (/ 0.00_wp, -12.0_wp, -11.7_wp,  0.00_wp, -0.72_wp,-0.243_wp,  0.19_wp,-0.062_wp,    &
              0.414_wp, 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp, 0.00_wp /)

    if (iband == 1 .or. iband == 4 .or. iband >= 10) then
      find_secdiff = 1.66_wp
    else
      find_secdiff = MAX(MIN(a0(iband) + a1(iband) * exp(a2(iband)*pwvcm), 1.8_wp), 1.5_wp) 
    endif

  END FUNCTION find_secdiff
  ! -------------------------------------------------------------------------------
END MODULE mo_lrtm_solver
