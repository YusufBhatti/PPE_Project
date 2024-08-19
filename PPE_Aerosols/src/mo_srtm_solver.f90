#ifdef __xlC__
@PROCESS HOT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to do shortwave radiative transfer calculation
!!
!! @remarks
!!   This module contains routines to compute shortwave radiative transfer 
!!     given optical properties. Two-stream methods provide
!!     layer transmittance and reflectance; adding computing flux profiles
!! 
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2010-08)
!!         Robert Pincus, U. Colorado, visiting MPI-M, Hamburg (2011-07)
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
!!   Many of the comments are in the original. 
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
MODULE mo_srtm_solver

  USE mo_kind,        ONLY : wp
  USE mo_rad_fastmath,ONLY : expon, inv_expon, transmit

  IMPLICIT NONE
  
  INTERFACE delta_scale
    MODULE PROCEDURE delta_scale_1d, delta_scale_2d
  END INTERFACE delta_scale

  INTERFACE two_stream
    MODULE PROCEDURE two_stream_ec
  END INTERFACE two_stream

  INTERFACE srtm_solver
    MODULE PROCEDURE srtm_solver_opt, srtm_solver_tr
  END INTERFACE srtm_solver

  REAL(WP), PARAMETER :: w0Min      = 0.9999995_wp ! Min single scattering albedo for conservative approx.
  REAL(WP), PARAMETER :: thinLimit  = 1.E-4_wp     ! Max optical path (tau/mu0) for thin limit approx. 
  REAL(WP), PARAMETER :: thickLimit = 100._wp, &   ! Min optical thickness (tau) for semi-infinite limit approx. 
                         thickw0Limit = 0.1        ! w0 limit for total absorption in thick case
  PRIVATE 
  PUBLIC :: srtm_solver, delta_scale, two_stream

CONTAINS
! ---------------------------------------------------------------------------
!>
!! @brief Solver using optical properties of each layer as input
!!
  SUBROUTINE srtm_solver_opt(kproma, kbdim, klev, &
                          &  palbd, palbp, prmu0, &
                          &  ptau,  pasy,  pomg,  &
                          &  zfd, zfu, ztdbt_sfc)
    INTEGER,  INTENT(in) :: kproma, kbdim, klev
    REAL(wp), INTENT(in) :: palbd(kbdim)      !< surface albedo (diffuse)
    REAL(wp), INTENT(in) :: palbp(kbdim)      !< surface albedo (direct)
    REAL(wp), INTENT(in) :: prmu0(kbdim)      !< cosine of solar zenith angle
    !
    !   Dimensions: (klev)
    !
    REAL(wp), INTENT(in) :: ptau(kbdim,klev) !<  optical depth
    REAL(wp), INTENT(in) :: pasy(kbdim,klev) !<  asymmetry parameter
    REAL(wp), INTENT(in) :: pomg(kbdim,klev) !<  single scattering albedo
    !
    !   Dimensions: (klev+1)
    !
    REAL(wp), INTENT(out) :: zfd(kbdim,klev+1), zfu(kbdim,klev+1) !< Flux down, up
    REAL(wp), INTENT(out) :: ztdbt_sfc(kbdim)

    INTEGER  :: jk
    REAL(wp) :: wf(kbdim)                                        ! Delta-scaling variable
    REAL(wp) :: g(kbdim,klev), w0(kbdim,klev), tau(kbdim,klev) ! Optical parameters, delta-scaled
    REAL(wp) :: Rdir(kbdim,klev), Rdif(kbdim,klev) ! Direct and diffuse reflectance
    REAL(wp) :: Tdir(kbdim,klev), Tdif(kbdim,klev) ! Direct and diffuse transmittance

    ! Following are comments from the AER code from which this is derived
    ! ---------------------------------------------------------------------------
    !
    ! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
    !          using the two-stream method of H. Barker. 
    !
    ! Method:
    !    Adapted from two-stream model of H. Barker;
    !    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
    !        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
    !
    ! Modifications:
    !
    ! Original: H. Barker
    ! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
    ! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
    ! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
    ! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
    ! Revision: Code modified so that delta scaling is not done in cloudy profiles    
    !           if routine cldprop is used; delta scaling can be applied by swithcing
    !           code below if cldprop is not used to get cloud properties. 
    !           AER, Jan 2005
    ! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
    ! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
    !           Aug 2007 
    !
    ! ------------------------------------------------------------------

    ! Delta-scale input properties
    !   (Doing this in line instead of calling the function below saves copying input arrays)
    DO jk = 1, klev
	  wf(1:kproma) = pomg(1:kproma,jk)*pasy(1:kproma,jk)*pasy(1:kproma,jk) ! w0*f with f = g**2
      tau(1:kproma,jk) = (1.0_wp - wf(1:kproma)) * ptau(1:kproma, jk)
      w0(1:kproma,jk) = (pomg(1:kproma,jk) - wf(1:kproma)) / (1.0_wp - wf(1:kproma))
      ! g = (g - f) / (1 - f) ; f = g**w
      g (1:kproma,jk) =  pasy(1:kproma,jk) / ( 1.0_wp + pasy(1:kproma,jk))
    END DO 

    CALL two_stream(kproma         ,kbdim    ,klev     ,           &
                 &  prmu0(:)       ,tau      ,w0       ,g        , &
                 &  Rdir(:,:)      ,Rdif(:,:),Tdir(:,:),Tdif(:,:)  )

    !
    ! Boundary conditions and adding to get flux layer-by-layer
    !
    CALL srtm_solver_tr(kproma, kbdim, klev,      &
                        palbp, palbd, prmu0, tau, &
                        Rdir, Rdif, Tdir, Tdif,   &
                        zfd, zfu, ztdbt_sfc)
                        
  END SUBROUTINE srtm_solver_opt
! ---------------------------------------------------------------------------
!>
!! @brief Solver using optical properties and pre-computed direct/diffuse
!!   reflectanc/transmittance to compute fluxes layer-by-layer 
!!
  SUBROUTINE srtm_solver_tr(kproma, kbdim, klev,    &
                      &  palbd, palbp, prmu0, ptau, &
                      &  Rdir, Rdif, Tdir, Tdif,    &
                      &  zfd, zfu, ztdbt_sfc)
    INTEGER,  INTENT(in) :: kproma, kbdim, klev
    REAL(wp), INTENT(in) :: palbd(kbdim)      !< surface albedo (diffuse)
    REAL(wp), INTENT(in) :: palbp(kbdim)      !< surface albedo (direct)
    REAL(wp), INTENT(in) :: prmu0(kbdim)      !< cosine of solar zenith angle
    !
    !   Dimensions: (klev)
    !
    REAL(wp), INTENT(in) :: ptau(kbdim,klev)            !<  optical depth
    REAL(wp), INTENT(in) :: Rdir(kbdim,klev), Rdif(kbdim,klev) !<  Diffuse and direct reflectance
    REAL(wp), INTENT(in) :: Tdir(kbdim,klev), Tdif(kbdim,klev) !<  Diffuse and direct transmittance
    !
    !   Dimensions: (klev+1)
    !
    REAL(wp), INTENT(out) :: zfd(kbdim,klev+1), zfu(kbdim,klev+1) !< Flux down, up
    REAL(wp), INTENT(out) :: ztdbt_sfc(kbdim)

    INTEGER  :: jk, jkr
    REAL(wp) :: zref(kbdim,klev+1), zrefd(kbdim,klev+1) ! Direct and diffuse reflectance, reordered and with BCs
    REAL(wp) :: ztra(kbdim,klev+1), ztrad(kbdim,klev+1) ! Direct and diffuse transmittance
    REAL(wp) :: zdbt(kbdim,klev+1), ztdbt(kbdim,klev+1) ! Direct beam transmittance per layer and total to each layer
    ! ---------------------------------------------------------------------------

    ! Note: two-stream calculations proceed from top to bottom; 
    !   RRTMG_SW quantities are given bottom to top and are reversed here
    zref (1:kproma,1:klev) = Rdir(1:kproma,klev:1:-1)
    zrefd(1:kproma,1:klev) = Rdif(1:kproma,klev:1:-1)
    ztra (1:kproma,1:klev) = Tdir(1:kproma,klev:1:-1)
    ztrad(1:kproma,1:klev) = Tdif(1:kproma,klev:1:-1)


    ! Surface values
    zref (1:kproma,klev+1) = palbp(1:kproma)
    zrefd(1:kproma,klev+1) = palbd(1:kproma)
    ztra (1:kproma,klev+1) = 0.0_wp
    ztrad(1:kproma,klev+1) = 0.0_wp
    
    ! Direct beam transmission -- need to reverse vertical ordering of optical depth array
!IBM* ASSERT(NODEPS)
	DO jk = 1, klev
      jkr = klev+1-jk                       
      zdbt(1:kproma,jk) = inv_expon(ptau(1:kproma,jkr)/prmu0(1:kproma), kproma)
    END DO 
    zdbt(1:kproma,klev+1) = 0.0_wp
    
    ! Accumulated direct beam transmission 
    ztdbt(1:kproma,1) = 1._wp
    DO jk=1,klev
      ztdbt(1:kproma,jk+1) = zdbt(1:kproma,jk)*ztdbt(1:kproma,jk)
    ENDDO
    ztdbt_sfc(1:kproma)  = ztdbt(1:kproma,klev+1)

    ! Vertical quadrature for cloudy fluxes
    CALL adding (kproma, kbdim, klev, &
                &  zref(:,:), zrefd(:,:), ztra(:,:), ztrad(:,:), &      
                &  zdbt(:,:), ztdbt(:,:), palbp( :), palbd(  :), &
                &  zfd (:,:), zfu  (:,:)) 
  END SUBROUTINE srtm_solver_tr
    ! ---------------------------------------------------------------------------

!>
!! @brief "vertical quadrature" i.e. flux profiles based on layer-by-layer direct
!!    and diffuse reflectance and transmittance
!!
!! @ remarks Equations are developed in doi:10.1002/qj.49712555316. 

  SUBROUTINE adding(kproma, kbdim, klev, &
       pref, prefd, ptra, ptrad, &
       pdbt, ptdbt, palbp, palbd,&
       pfd, pfu)
    ! --------------------------------------------------------------------------

    ! Purpose: This routine performs the vertical quadrature integration
    !
    ! Modifications.
    ! 
    ! Original: H. Barker
    ! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
    ! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
    !
    !-----------------------------------------------------------------------
    INTEGER, INTENT (in) :: kproma, kbdim, klev  ! number of columns, max. number of col., model layers

    ! 
    ! All vectors are dimensioned kproma, klev+1; last value indicates surface
    ! 
    REAL(wp), INTENT(in) :: pref(kbdim,klev+1), prefd(kbdim,klev+1)  !< direct amd diffuse beam reflectivity
    REAL(wp), INTENT(in) :: ptra(kbdim,klev+1), ptrad(kbdim,klev+1)  !< direct and diffuse beam transmissivity
    REAL(wp), INTENT(in) :: pdbt(kbdim,klev+1), ptdbt(kbdim,klev+1)  !< direct beam transmission, by layer and accumulated
    REAL(wp), INTENT(in) :: palbp(kbdim),  palbd(kbdim)    !< Surface albedo for direct and diffuse incidence

    REAL(wp), INTENT(out) :: pfd(kbdim,klev+1)              ! downwelling flux (W/m2)
    REAL(wp), INTENT(out) :: pfu(kbdim,klev+1)              ! upwelling   flux (W/m2)

    INTEGER  :: ikp, ikx, jk
    REAL(wp) :: zreflect(kbdim)
    REAL(wp) :: ztdn(kbdim,klev+1)  
    REAL(wp) :: zrup(kbdim,klev+1), zrupd(kbdim,klev+1)
    REAL(wp) ::                     zrdnd(kbdim,klev+1)
    !-----------------------------------------------------------------------------

    ! Link lowest layer with surface
    zrup (1:kproma,klev+1) = palbp(1:kproma)
    zrupd(1:kproma,klev+1) = palbd(1:kproma)

    zreflect(1:kproma) = 1._wp / (1._wp - prefd(1:kproma,klev+1) * prefd(1:kproma,klev))
    zrup (1:kproma,klev) = pref(1:kproma,klev)                                                    &
                       & + (ptrad(1:kproma,klev) * ((ptra(1:kproma,klev) - pdbt(1:kproma,klev)) * & 
                           prefd(1:kproma,klev+1)                                                 & 
                       & + pdbt (1:kproma,klev) * pref(1:kproma,klev+1))) * zreflect(1:kproma)
    zrupd(1:kproma,klev) = prefd(1:kproma,klev) &
                       & + ptrad(1:kproma,klev) * ptrad(1:kproma,klev) * &
                       &   prefd(1:kproma,klev+1) * zreflect(1:kproma)

    ! Pass from bottom to top 

    DO jk = 1,klev-1
      ikp = klev+1-jk                       
      ikx = ikp-1
      zreflect(1:kproma) = 1._wp / (1._wp - zrupd(1:kproma,ikp) * prefd(1:kproma,ikx))
      zrup(1:kproma,ikx) = pref(1:kproma,ikx)                                                &
                       & + (ptrad(1:kproma,ikx) *                                            &
                            ((ptra(1:kproma,ikx) - pdbt(1:kproma,ikx)) * zrupd(1:kproma,ikp) & 
                       &   +  pdbt(1:kproma,ikx) * zrup(1:kproma,ikp))) * zreflect(1:kproma)
      zrupd(1:kproma,ikx) = prefd(1:kproma,ikx) &
                        & + ptrad(1:kproma,ikx) * ptrad(1:kproma,ikx) * zrupd(1:kproma,ikp) * &
                            zreflect(1:kproma)
    ENDDO

    ! Upper boundary conditions

    ztdn (1:kproma,1) = 1._wp
    zrdnd(1:kproma,1) = 0._wp
    ztdn (1:kproma,2) = ptra(1:kproma,1)
    zrdnd(1:kproma,2) = prefd(1:kproma,1)

    ! Pass from top to bottom

    DO jk = 2,klev
      ikp = jk+1
      zreflect(1:kproma) = 1._wp / (1._wp - prefd(1:kproma,jk) * zrdnd(1:kproma,jk))
      ztdn(1:kproma,ikp) = ptdbt(1:kproma,jk) * ptra(1:kproma,jk)                            &
                       & + ( ptrad(1:kproma,jk) * ((ztdn(1:kproma,jk) - ptdbt(1:kproma,jk))  &
                       &   + ptdbt(1:kproma,jk) * pref(1:kproma,jk) * zrdnd(1:kproma,jk))) * &
                    &   zreflect(1:kproma)
      zrdnd(1:kproma,ikp) = prefd(1:kproma,jk) &
                          + ptrad(1:kproma,jk) * ptrad(1:kproma,jk) * zrdnd(1:kproma,jk) * &
                            zreflect(1:kproma)
    ENDDO

    ! Up and down-welling fluxes at levels

!IBM* ASSERT(NODEPS)
    DO jk = 1,klev+1
      zreflect(1:kproma) = 1._wp / (1._wp - zrdnd(1:kproma,jk) * zrupd(1:kproma,jk))
      pfu(1:kproma,jk) = (ptdbt(1:kproma,jk) * zrup(1:kproma,jk) &
                     & + (ztdn(1:kproma,jk) - ptdbt(1:kproma,jk)) * zrupd(1:kproma,jk)) * & 
                     &    zreflect(1:kproma)
      pfd(1:kproma,jk) = ptdbt(1:kproma,jk) &
                     & + (ztdn(1:kproma,jk) - ptdbt(1:kproma,jk) &
                    &    + ptdbt(1:kproma,jk) * zrup(1:kproma,jk) * zrdnd(1:kproma,jk)) * &
                         zreflect(1:kproma)
    END DO
  END SUBROUTINE adding

  ! --------------------------------------------------------------------
!>
!! @brief Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
!!    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.  
!!
!! @ remarks Equations are developed in Meador and Weaver, 1980, 
!!    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
!!
  SUBROUTINE two_stream_rp(kproma, kbdim, klev, mu0, tau, w0, g, Rdir, Rdif, Tdir, Tdif, update) 
    INTEGER,  INTENT( IN) :: kproma, kbdim, klev            !< number of columns, max. number of col., number of levels
    REAL(WP), INTENT( IN) :: mu0(kbdim),               & !< Cosine of solar zenith angle
                             tau(kbdim,klev),             & !< Optical depth
                             w0(kbdim,klev),              & !< Single scattering albedo
                             g(kbdim,klev)                  !< Asymmetery parameter
    REAL(WP), INTENT(OUT) :: Rdir(kbdim,klev), Tdir(kbdim,klev), & !< Reflectance and transmittance for direct beam
                             Rdif(kbdim,klev), Tdif(kbdim,klev)    !< Reflectance and transmittance for diffuse illumination
    LOGICAL,  INTENT( IN), OPTIONAL :: update(kbdim,klev)   !< Update this cell? 
    
    !
    ! Take a 2D problem, collect three solution classes (thin/conservative/general), and solve in 
    !   1D batches of size kproma
    ! 
    ! Masks for transparent, thick, reflective, and general solutions
    LOGICAL  :: isThin(kbdim,klev), isThck(kbdim,klev), isRefl(kbdim,klev), doThis(kbdim,klev) 
    INTEGER  :: nUse
    REAL(WP) :: Rdir_l(kbdim*klev), Tdir_l(kbdim*klev), Rdif_l(kbdim*klev), Tdif_l(kbdim*klev), & 
                 mu0_l(kbdim*klev),  tau_l(kbdim*klev),   w0_l(kbdim*klev),    g_l(kbdim*klev)

!    ! Eddington approximation (Joseph et al., 1976; doi: 10.1175/1520-0469(1976)033<2452:TDEAFR>2.0.CO;2)
!    gamma1(:)= (7._wp - w0(1:kproma)*(4._wp + 3._wp*g(1:kproma))) / 4._wp 
!    gamma2(:)=-(1._wp - w0(1:kproma)*(4._wp - 3._wp*g(1:kproma))) / 4._wp
!    gamma3(:)= (2._wp - 3._wp * g(1:kproma) * mu0(1:kproma) ) / 4._wp

!    ! Discrete-ordinates quadrature (Liou, 1973; doi: 10.1175/1520-0469(1973)030<1303:ANEOCD>2.0.CO;2)
!    gamma1(:)= SQRT(3._wp)/2._wp * (2._wp - w0(1:kproma)*(1._wp + g(1:kproma))) 
!    gamma2(:)= SQRT(3._wp)/2._wp * w0(1:kproma) * (1._wp - g(1:kproma)) 
!    gamma3(:)= (1._wp - SQRT(3._wp) * g(1:kproma) * mu0(1:kproma) )/2._wp

    ! Zdunkowski "PIFM"  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
!    gamma1(:)= (8._wp - w0(1:kproma) * (5._wp + 3._wp*g(1:kproma))) * .25_wp
!    gamma2(:)=  3._wp *(w0(1:kproma) * (1._wp -       g(1:kproma))) * .25_wp
!    gamma3(:)= (2._wp - 3._wp * mu0(1:kproma) *       g(1:kproma) ) * .25_wp
!    gamma4(:)=  1._wp - gamma3(:) 

    !
    ! Is the sun above the horizon, and do we want to update this cell 
    !
    doThis(1:kproma,1:klev) = SPREAD(mu0(1:kproma) > 0._wp, DIM = 2, NCOPIES = klev)
    IF (PRESENT(update)) THEN 
      doThis(1:kproma,1:klev) = doThis(1:kproma,1:klev) .and. update(1:kproma,1:klev)      
    ELSE
      WHERE(.not. doThis) 
        Rdir(1:kproma,1:klev) = 1._wp
        Rdif(1:kproma,1:klev) = 1._wp
        Tdir(1:kproma,1:klev) = 0._wp
        Tdif(1:kproma,1:klev) = 0._wp
      END WHERE
    END IF 
    !
    ! Thin medium 
    !
    isThin(1:kproma,1:klev) = tau(1:kproma,1:klev)/SPREAD(mu0(1:kproma), DIM = 2, NCOPIES = klev) <= thinLimit .and. &
      doThis(1:kproma,1:klev)
    nUse = COUNT(isThin(1:kproma,1:klev))
    IF(nUse > 0) THEN 
      mu0_l(:nUse) = PACK(SPREAD(mu0(1:kproma), DIM = 2, NCOPIES = klev), isThin)
      tau_l(:nUse) = PACK(tau(1:kproma,1:klev), isThin)
       w0_l(:nUse) = PACK( w0(1:kproma,1:klev), isThin)
        g_l(:nUse) = PACK(  g(1:kproma,1:klev), isThin)

      CALL two_stream_thin(nUse, mu0_l, tau_l, w0_l, g_l, Rdir_l, Rdif_l, Tdir_l, Tdif_l)

      !
      ! UNPACK into 2D arrays
      ! 
      Rdir(1:kproma,1:klev) = UNPACK(Rdir_l(:nUse), MASK = isThin(1:kproma,1:klev), FIELD = Rdir(1:kproma,1:klev))
      Rdif(1:kproma,1:klev) = UNPACK(Rdif_l(:nUse), MASK = isThin(1:kproma,1:klev), FIELD = Rdif(1:kproma,1:klev))
      Tdir(1:kproma,1:klev) = UNPACK(Tdir_l(:nUse), MASK = isThin(1:kproma,1:klev), FIELD = Tdir(1:kproma,1:klev))
      Tdif(1:kproma,1:klev) = UNPACK(Tdif_l(:nUse), MASK = isThin(1:kproma,1:klev), FIELD = Tdif(1:kproma,1:klev))
    END IF  
    WHERE(isThin(1:kproma,1:klev)) doThis = .false. 

    !
    ! Thick medium 
    !
    isThck(:,:) = .false. ! tau(:,:) > ThickLimit .and. doThis(:,:)
    nUse = COUNT(isThck(1:kproma,1:klev))
    IF(nUse > 0) THEN 
      mu0_l(:nUse) = PACK(SPREAD(mu0(1:kproma), DIM = 2, NCOPIES = klev), isThck)
      tau_l(:nUse) = PACK(tau(1:kproma,1:klev), isThck)
       w0_l(:nUse) = PACK( w0(1:kproma,1:klev), isThck)
        g_l(:nUse) = PACK(  g(1:kproma,1:klev), isThck)

      CALL two_stream_thick(nUse, mu0_l, tau_l, w0_l, g_l, Rdir_l, Rdif_l, Tdir_l, Tdif_l)

      !
      ! UNPACK into 2D arrays
      ! 
      Rdir(1:kproma,1:klev) = UNPACK(Rdir_l(:nUse), MASK = isThck(1:kproma,1:klev), FIELD = Rdir(1:kproma,1:klev))
      Rdif(1:kproma,1:klev) = UNPACK(Rdif_l(:nUse), MASK = isThck(1:kproma,1:klev), FIELD = Rdif(1:kproma,1:klev))
      Tdir(1:kproma,1:klev) = UNPACK(Tdir_l(:nUse), MASK = isThck(1:kproma,1:klev), FIELD = Tdir(1:kproma,1:klev))
      Tdif(1:kproma,1:klev) = UNPACK(Tdif_l(:nUse), MASK = isThck(1:kproma,1:klev), FIELD = Tdif(1:kproma,1:klev))
    END IF  
    WHERE(isThck(1:kproma,1:klev)) doThis = .false. 

    !
    ! Conservative scattering solution
    !
    isRefl(1:kproma,1:klev) = doThis(1:kproma,1:klev) .and. w0(1:kproma,1:klev) > w0Min ! 1._wp - 100._wp * EPSILON(1._wp) 
    nUse = COUNT(isRefl)
    IF(nUse > 0) THEN 
      mu0_l(:nUse) = PACK(SPREAD(mu0(1:kproma), DIM = 2, NCOPIES = klev), isRefl)
      tau_l(:nUse) = PACK(tau(1:kproma,1:klev), isRefl)
       w0_l(:nUse) = PACK( w0(1:kproma,1:klev), isRefl)
        g_l(:nUse) = PACK(  g(1:kproma,1:klev), isRefl)

      CALL two_stream_conserv(nUse, mu0_l, tau_l, w0_l, g_l, Rdir_l, Rdif_l, Tdir_l, Tdif_l)

      Rdir(1:kproma,1:klev) = UNPACK(Rdir_l(:nUse), MASK = isRefl(1:kproma,1:klev), FIELD = Rdir(1:kproma,1:klev))
      Rdif(1:kproma,1:klev) = UNPACK(Rdif_l(:nUse), MASK = isRefl(1:kproma,1:klev), FIELD = Rdif(1:kproma,1:klev))
      Tdir(1:kproma,1:klev) = UNPACK(Tdir_l(:nUse), MASK = isRefl(1:kproma,1:klev), FIELD = Tdir(1:kproma,1:klev))
      Tdif(1:kproma,1:klev) = UNPACK(Tdif_l(:nUse), MASK = isRefl(1:kproma,1:klev), FIELD = Tdif(1:kproma,1:klev))
    END IF  
    WHERE(isRefl(1:kproma,1:klev)) doThis = .false. 

    !
    ! General solution
    !
    nUse = COUNT(doThis(1:kproma,1:klev))
    IF(nUse > 0) THEN 
      mu0_l(:nUse) = PACK(SPREAD(mu0(1:kproma), DIM = 2, NCOPIES = klev), doThis)
      tau_l(:nUse) = PACK(tau(1:kproma,1:klev), doThis)
       w0_l(:nUse) = PACK( w0(1:kproma,1:klev), doThis)
        g_l(:nUse) = PACK(  g(1:kproma,1:klev), doThis)
        
      CALL two_stream_general(nUse, mu0_l, tau_l, w0_l, g_l, Rdir_l, Rdif_l, Tdir_l, Tdif_l)
      
      Rdir(1:kproma,1:klev) = UNPACK(Rdir_l(:nUse), MASK = doThis(1:kproma,1:klev), FIELD = Rdir(1:kproma,1:klev))
      Rdif(1:kproma,1:klev) = UNPACK(Rdif_l(:nUse), MASK = doThis(1:kproma,1:klev), FIELD = Rdif(1:kproma,1:klev))
      Tdir(1:kproma,1:klev) = UNPACK(Tdir_l(:nUse), MASK = doThis(1:kproma,1:klev), FIELD = Tdir(1:kproma,1:klev))
      Tdif(1:kproma,1:klev) = UNPACK(Tdif_l(:nUse), MASK = doThis(1:kproma,1:klev), FIELD = Tdif(1:kproma,1:klev))
    END IF  
  END SUBROUTINE two_stream_rp
   ! -----------------
  SUBROUTINE two_stream_thin(kproma, mu0, tau, w0, g, Rdir, Rdif, Tdir, Tdif) 
    INTEGER, INTENT(in)  :: kproma
    REAL(wp),INTENT(in)  :: mu0(:), tau(:), w0(:), g(:)
    REAL(wp),INTENT(out) :: Rdir(:), Rdif(:),Tdir(:), Tdif(:)
  
    INTEGER  :: i
    REAL(wp) :: gamma1(kproma), gamma2(kproma), gamma3(kproma)
    REAL(wp) :: tau_over_mu(kproma)
    
!IBM* ASSERT(NODEPS)
    DO i = 1, kproma  ! Perhaps later we'll add the capability to loop over blocks of a set size
	  ! Zdunkowski "PIFM"  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
	  gamma1(i)= (8._wp - w0(i) * (5._wp + 3._wp * g(i))) * .25_wp
	  gamma2(i)=  3._wp *(w0(i) * (1._wp -         g(i))) * .25_wp
	  gamma3(i)= (2._wp - 3._wp * mu0(i) *         g(i) ) * .25_wp

	  tau_over_mu(i) = tau(i)/mu0(i)
	  
	  Rdir(i) =  w0(i)  * tau_over_mu(i) * gamma3(i)                 ! Eq. 19
	  Tdir(i) = 1._wp - (1._wp - w0(i)) * tau_over_mu(i) - Rdir(i)   ! Eq. 20
	  Rdif(i) = gamma2(i) * tau(i)                                   ! Eq. 27
	  Tdif(i) = 1._wp - Rdif(i) - (gamma1(i) - gamma2(i)) * tau(i)

    END DO 
  END SUBROUTINE two_stream_thin

   ! -----------------
  SUBROUTINE two_stream_thick(kproma, mu0, tau, w0, g, Rdir, Rdif, Tdir, Tdif) 
    !
    ! Semi-infinite limit
    !
    INTEGER, INTENT(in)  :: kproma
    REAL(wp),INTENT(in)  :: mu0(:), tau(:), w0(:), g(:)
    REAL(wp),INTENT(out) :: Rdir(:), Rdif(:),Tdir(:), Tdif(:)
  
    INTEGER  :: i
    REAL(wp) :: gamma1(kproma), gamma2(kproma), gamma3(kproma), gamma4(kproma)
    REAL(wp) :: k(kproma), alpha2(kproma), k_plus_gamma1_inv(kproma)
    
!IBM* ASSERT(NODEPS)
    DO i = 1, kproma  ! Perhaps later we'll add the capability to loop over blocks of a set size
      IF(w0(i) > thickw0Limit) THEN 
        ! Zdunkowski "PIFM"  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
        gamma1(i)= (8._wp - w0(i) * (5._wp + 3._wp * g(i))) * .25_wp
        gamma2(i)=  3._wp *(w0(i) * (1._wp -         g(i))) * .25_wp
        gamma3(i)= (2._wp - 3._wp * mu0(i) *         g(i) ) * .25_wp
        gamma4(i)=  1._wp - gamma3(i)

        k     (i) = SQRT((gamma1(i) - gamma2(i)) * (gamma1(i) + gamma2(i))) ! Eq 18;  k = SQRT(gamma1**2 - gamma2**2)
        alpha2(i) = gamma1(i) * gamma3(i) + gamma2(i) * gamma4(i)           ! Eq. 17
        k_plus_gamma1_inv(i) = 1._wp/(k(i) + gamma1(i)) 
      
        Rdir(i) = w0(i)*(alpha2(i) + k(i)*gamma3(i)) * k_plus_gamma1_inv(i) / (1._wp - k(i)*mu0(i))
                                                   ! Eq. 22
        Rdif(i) = gamma2(i) * k_plus_gamma1_inv(i) ! Eq. 25, large tau
	  ELSE
        Rdir(i) = 0._wp
        Rdif(i) = 0._wp
	  END IF 
    END DO 
    Tdir(:kproma) = 0._wp ! Eq. 15, limit of large tau
    Tdif(:kproma) = 0._wp ! Eq. 26, large tau
  END SUBROUTINE two_stream_thick
     ! -----------------
  SUBROUTINE two_stream_conserv(kproma, mu0, tau, w0, g, Rdir, Rdif, Tdir, Tdif) 
    INTEGER, INTENT(in)  :: kproma
    REAL(wp),INTENT(in)  :: mu0(:), tau(:), w0(:), g(:)
    REAL(wp),INTENT(out) :: Rdir(:), Rdif(:),Tdir(:), Tdif(:)
  
    INTEGER  :: i
    REAL(wp) :: gamma1(kproma),  gamma3(kproma) 
    REAL(wp) :: tau_over_mu(kproma), gamma1_tau(kproma), transDir(kproma)
    
!IBM* ASSERT(NODEPS)
    DO i = 1, kproma  ! Perhaps later we'll add the capability to loop over blocks of a set size
	  ! Zdunkowski "PIFM"  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
	  gamma1(i)= (8._wp - w0(i) * (5._wp + 3._wp * g(i))) * .25_wp
	  gamma3(i)= (2._wp - 3._wp * mu0(i) *         g(i) ) * .25_wp

	  tau_over_mu(i) = tau(i)/mu0(i)
	  gamma1_tau(i) = gamma1(i) * tau(i)
    END DO 
    
    transDir(:) = transmit(tau_over_mu, kproma) ! 1. - exp(-tau/mu) 
    
!IBM* ASSERT(NODEPS)
    DO i = 1, kproma
      ! Eq 24
      Rdir(i) = (gamma1_tau(i) + (gamma3(i) - gamma1(i)*mu0(i))*transdir(i))/(1.0_wp + gamma1_tau(i))
      Tdir(i) = 1._wp - Rdir(i)
      ! Eq 29
      Rdif(i) = gamma1_tau(i) / (1.0_wp + gamma1_tau(i))
      Tdif(i) = 1._wp - Rdif(i)
      ! Original RRTMG code also contains a check (when using lookup tables) for transDir == 1. 
      !  and sets Tdi[rf] = 1, Rdi[rf] = 0 where that's true
    END DO 
  END SUBROUTINE two_stream_conserv

   ! -----------------
  SUBROUTINE two_stream_general(kproma, mu0, tau, w0, g, Rdir, Rdif, Tdir, Tdif) 
    INTEGER, INTENT(in)  :: kproma
    REAL(wp),INTENT(in)  :: mu0(:), tau(:), w0(:), g(:)
    REAL(wp),INTENT(out) :: Rdir(:), Rdif(:),Tdir(:), Tdif(:)
    
    INTEGER  :: i
    ! Variables used in Meador and Weaver
    REAL(wp) :: gamma1(kproma), gamma2(kproma), gamma3(kproma), gamma4(kproma)
    REAL(wp) :: alpha1(kproma), alpha2(kproma), k(kproma)
    ! Ancillary variables
    REAL(wp) :: tau_over_mu(kproma)
    REAL(WP) :: one_minus_kmu_sqr(kproma), k_tau(kproma), k_mu(kproma), k_gamma3(kproma), k_gamma4(kproma)
    REAL(wp) :: exp_ktau(kproma), exp_2ktau(kproma), exp_minus_tau_over_mu0(kproma) 
    REAL(wp) :: RT_denom(kproma), RTdif_denom(kproma)
    REAL(wp) :: T_coeff_1(kproma), T_coeff_2(kproma), T_coeff_3(kproma)   
    
!IBM* ASSERT(NODEPS)
    DO i = 1, kproma  ! Perhaps later we'll add the capability to loop over blocks of a set size
      ! Zdunkowski "PIFM"  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
      gamma1(i)= (8._wp - w0(i) * (5._wp + 3._wp * g(i))) * .25_wp
      gamma2(i)=  3._wp *(w0(i) * (1._wp -         g(i))) * .25_wp
      gamma3(i)= (2._wp - 3._wp * mu0(i) *         g(i) ) * .25_wp
      gamma4(i)=  1._wp - gamma3(i)

      alpha1(i) = gamma1(i) * gamma4(i) + gamma2(i) * gamma3(i)           ! Eq. 16
      alpha2(i) = gamma1(i) * gamma3(i) + gamma2(i) * gamma4(i)           ! Eq. 17

      k     (i) = SQRT(MAX((gamma1(i) - gamma2(i)) * (gamma1(i) + gamma2(i)), 1.E-12_wp)) ! Eq 18;  k = SQRT(gamma1**2 - gamma2**2)
      k_tau (i) = k(i) * tau(i)
      tau_over_mu(i) = tau(i)/mu0(i)
    END DO 
    
    exp_minus_tau_over_mu0(:) = inv_expon(MIN(tau_over_mu, 500._wp), kproma)
    exp_ktau              (:) =     expon(MIN(k_tau,       500._WP), kproma)
        
!IBM* ASSERT(NODEPS)
    DO i = 1, kproma
      k_mu    (i) = k(i) * mu0(i)
      k_gamma3(i) = k(i) * gamma3(i)
      k_gamma4(i) = k(i) * gamma4(i)

      exp_2ktau(i) = exp_ktau(i) * exp_ktau(i)
      one_minus_kmu_sqr(i) = 1.0_wp - k_mu(i) * k_mu(i)

      ! Equation 14, multiplying top and bottom by exp(k*tau) 
      RT_denom(i) = &
        & one_minus_kmu_sqr(i) * (k(i) + gamma1(i)) * exp_2ktau(i) + &
        & one_minus_kmu_sqr(i) * (k(i) - gamma1(i)) 
      
      Rdir(i) = w0(i) * (                                                      &
        & ((1._wp - k_mu(i)) * (alpha2(i) + k_gamma3(i)) * exp_2ktau(i) -     &
        &  (1._wp + k_mu(i)) * (alpha2(i) - k_gamma3(i)) -                 &
        &  2.0_wp * (k_gamma3(i) - alpha2(i) * k_mu(i)) * exp_minus_tau_over_mu0(i) * exp_ktau(i)) &
        & ) /  RT_denom(i)

      ! Equation 15, multiplying top and bottom by exp(k*tau), 
      !   and precomputing some factors for numerical stability
      T_coeff_1(i) = (1._wp + k_mu(i)) * (alpha1(i) + k_gamma4(i))
      T_coeff_2(i) = (1._wp - k_mu(i)) * (alpha1(i) - k_gamma4(i))
      T_coeff_3(i) = 2.0_wp * (k_gamma4(i) + alpha1(i) * k_mu(i))

      Tdir(i) = exp_minus_tau_over_mu0(i) - (w0(i) *                   &
        & ((T_coeff_1(i) * exp_minus_tau_over_mu0(i) * exp_2ktau(i)  - &
        &   T_coeff_2(i) * exp_minus_tau_over_mu0(i)  - &
        &   T_coeff_3(i) * exp_ktau (i) ))) / RT_denom(i)

      ! Diffuse reflection and transmission, multiplying top and bottom by exp(2*k*tau)
      RTdif_denom(i) = (k(i) + gamma1(i)) * exp_2ktau(i) + k(i) - gamma1(i)
         
      ! Equation 25
      Rdif(i) = (gamma2(i) * (exp_2ktau(i) - 1._wp)) / RTdif_denom(i)
      ! Equation 26
      Tdif(i) = (2._wp * k(i) * exp_ktau(i)) / RTdif_denom(i)

    END DO 
  END SUBROUTINE two_stream_general

! --------------------------------------------------------------------
  SUBROUTINE two_stream_ec(kproma, kbdim, klev, mu0, tau, w0, g, Rdir, Rdif, Tdir, Tdif, update)
    INTEGER,  INTENT( IN) :: kproma, kbdim, klev            !< number of columns, max. number of col., number of levels
    REAL(WP), INTENT( IN) :: mu0(kbdim),               & !< Cosine of solar zenith angle
                             tau(kbdim,klev),             & !< Optical depth
                             w0(kbdim,klev),              & !< Single scattering albedo
                             g(kbdim,klev)                  !< Asymmetery parameter
    REAL(WP), INTENT(INOUT) :: Rdir(kbdim,klev), Tdir(kbdim,klev), & !< Reflectance and transmittance for direct beam
                             Rdif(kbdim,klev), Tdif(kbdim,klev)    !< Reflectance and transmittance for diffuse illumination
    LOGICAL,  INTENT( IN), OPTIONAL :: update(kbdim,klev)   !< Update this cell? 
    
    REAL(WP)             :: zrmu0(kbdim)
              
    zrmu0(1:kproma)=1._wp/mu0(1:kproma)
    IF(PRESENT(update)) THEN 
      CALL srtm_reftra_ec(kproma, kbdim, klev,  &
         &   g, mu0, zrmu0, tau, w0, &
         &   Rdir  , Rdif, Tdir , Tdif, update) 
    ELSE
      CALL srtm_reftra_ec(kproma, kbdim, klev,  &
         &   g, mu0, zrmu0, tau, w0, &
         &   Rdir  , Rdif, Tdir , Tdif) 
    END IF

  END SUBROUTINE two_stream_ec

   ! -----------------
  SUBROUTINE srtm_reftra_ec &
       & ( icount, kbdim, klev  ,            &
       &   pgg   , prmuz, prmuzi, ptau , pw, &
       &   pref  , prefd, ptra , ptrad , &
       &   ldrtchk &
       & )  

    !**** *SRTM_REFTRA* - REFLECTIVITY AND TRANSMISSIVITY

    !     PURPOSE.
    !     --------
    !           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLEAR OR 
    !     CLOUDY LAYER USING A CHOICE OF VARIOUS APPROXIMATIONS.

    !**   INTERFACE.
    !     ----------
    !          *SRTM_REFTRA* IS CALLED BY *SRTM_SPCVRT*

    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! INPUTS
    ! ------ 
    !      KMODTS  = 1 EDDINGTON (JOSEPH ET AL., 1976)
    !              = 2 PIFM (ZDUNKOWSKI ET AL., 1980)
    !              = 3 DISCRETE ORDINATES (LIOU, 1973)
    !      LDRTCHK = .T. IF CLOUDY
    !              = .F. IF CLEAR-SKY
    !      PGG     = ASSYMETRY FACTOR
    !      PRMUZ   = COSINE SOLAR ZENITH ANGLE
    !      PRMUZI  = INVERSE COSINE SOLAR ZENITH ANGLE
    !      PTAU    = OPTICAL THICKNESS
    !      PW      = SINGLE SCATTERING ALBEDO

    ! OUTPUTS
    ! -------
    !      PREF    : COLLIMATED BEAM REFLECTIVITY
    !      PREFD   : DIFFUSE BEAM REFLECTIVITY 
    !      PTRA    : COLLIMATED BEAM TRANSMISSIVITY
    !      PTRAD   : DIFFUSE BEAM TRANSMISSIVITY

    !     METHOD.
    !     -------
    !          STANDARD DELTA-EDDINGTON, P.I.F.M., OR D.O.M. LAYER CALCULATIONS.

    !     EXTERNALS.
    !     ----------
    !          NONE

    !     REFERENCE.
    !     ----------

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 03-02-27
    !        M.Hamrud   01-Oct-2003      CY28 Cleaning
    !        Mike Iacono, AER, Mar 2004: bug fix 
    !        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
    !        M. Puetz   20-Apr-2010 Gather/Scatter for better pipelining on scalar cpus
    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------


    INTEGER,INTENT(in)  :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev 
    REAL(wp)   ,INTENT(in)    :: pgg(kbdim,klev) 
    REAL(wp)   ,INTENT(in)    :: prmuz(kbdim) 
    REAL(wp)   ,INTENT(in)    :: prmuzi(kbdim) 
    REAL(wp)   ,INTENT(in)    :: ptau(kbdim,klev) 
    REAL(wp)   ,INTENT(in)    :: pw(kbdim,klev) 
    REAL(wp)   ,INTENT(inout)   :: pref(kbdim,klev) 
    REAL(wp)   ,INTENT(inout)   :: prefd(kbdim,klev) 
    REAL(wp)   ,INTENT(inout)   :: ptra(kbdim,klev) 
    REAL(wp)   ,INTENT(inout)   :: ptrad(kbdim,klev) 
    LOGICAL,INTENT(in),optional :: ldrtchk(kbdim,klev) 
    !     ------------------------------------------------------------------
    INTEGER :: idxt(kbdim), idxf(kbdim), idxc(kbdim), idxn(kbdim)
    INTEGER :: jk, ic, jc, ict, icf, icc, icn

    REAL(wp) :: zgamma1(kbdim),zgamma2(kbdim),zgamma3(kbdim),zcrit(kbdim)
    REAL(wp) :: zrk(kbdim), zem1(kbdim), zem2(kbdim), zep1(kbdim), zep2(kbdim)
    REAL(wp) :: za, za1, za2, zemm
    REAL(wp) :: zbeta, zdend, zdenr, zdent
    REAL(wp) :: zg, zg3, zgamma4, zgt
    REAL(wp) :: zr1, zr2, zr3, zr4, zr5, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
    REAL(wp) :: zsr3, zt1, zt2, zt3, zt4, zt5
    REAL(wp) :: zw, zwcrit
    REAL(wp) :: ztemp, zzz
    REAL(wp), PARAMETER :: replog=1e-12     ! epsilon for lapace transformation
    INTEGER,  PARAMETER :: kmodts=2

    !     ------------------------------------------------------------------

    zsr3=SQRT(3._wp)
    zwcrit=w0Min

    IF (.NOT.PRESENT(ldrtchk)) THEN
      ict = icount
      icf = 0
      DO ic=1,icount
        idxt(ic) = ic
      END DO
    END IF

    DO jk=1,klev
      IF (PRESENT(ldrtchk)) THEN
        ict = 0
        icf = 0
        DO ic=1,icount
          IF (ldrtchk(ic,jk)) THEN
            ict = ict + 1
            idxt(ict) = ic
          ELSE
!            icf = icf + 1
!            idxf(icf) = ic
          END IF
        END DO
      END IF
      !-- GENERAL TWO-STREAM EXPRESSIONS

      IF (kmodts == 1) THEN
!IBM* ASSERT(NODEPS)
        DO jc = 1,ict
          ic = idxt(jc)

          zw  =pw(ic,jk)
          zg  =pgg(ic,jk)  

          zg3= 3._wp * zg
          zgamma1(ic)= (7._wp - zw * (4._wp + zg3)) * 0.25_wp
          zgamma2(ic)=-(1._wp - zw * (4._wp - zg3)) * 0.25_wp
          zgamma3(ic)= (2._wp - zg3 * prmuz(ic) ) * 0.25_wp
          zzz=(1._wp - zg)**2
          zcrit(jc) = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ELSEIF (kmodts == 2) THEN  
!IBM* ASSERT(NODEPS)
        DO jc = 1,ict
          ic = idxt(jc)

          zw  =pw(ic,jk)
          zg  =pgg(ic,jk)  

          zg3= 3._wp * zg
          zgamma1(ic)= (8._wp - zw * (5._wp + zg3)) * 0.25_wp
          zgamma2(ic)=  3._wp *(zw * (1._wp - zg )) * 0.25_wp
          zgamma3(ic)= (2._wp - zg3 * prmuz(ic) ) * 0.25_wp
          zzz=(1._wp - zg)**2
          zcrit(jc) = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ELSEIF (kmodts == 3) THEN  
!IBM* ASSERT(NODEPS)
        DO jc = 1,ict
          ic = idxt(jc)

          zw  =pw(ic,jk)
          zg  =pgg(ic,jk)  

          zg3= 3._wp * zg
          zgamma1(ic)= zsr3 * (2._wp - zw * (1._wp + zg)) * 0.5_wp
          zgamma2(ic)= zsr3 * zw * (1._wp - zg ) * 0.5_wp
          zgamma3(ic)= (1._wp - zsr3 * zg * prmuz(ic) ) * 0.5_wp
          zzz=(1._wp - zg)**2
          zcrit(jc) = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ENDIF

      icc = 0
      icn = 0
!IBM* ASSERT(NODEPS)
      DO jc = 1,ict
        ic = idxt(jc)

        !-- RECOMPUTE ORIGINAL S.S.A. TO TEST FOR CONSERVATIVE SOLUTION
        !   ZTEMP=(1._wp - ZG)**2
        !   ZWO= ZW*ZTEMP/ (ZTEMP - (1._wp - ZW)*(ZG **2))

        !       ZWO= ZW / (1._wp - (1._wp - ZW) * (ZG / (1._wp - ZG))**2)
        !       IF (ZWO >= ZWCRIT) THEN

        IF (zcrit(jc) >= 0._wp) THEN
          icc = icc + 1
          idxc(icc) = ic
        ELSE
          icn = icn + 1
          idxn(icn) = ic
        END IF
      END DO

      !-- conservative scattering
      !
      !-- Homogeneous reflectance and transmittance

!IBM* ASSERT(NODEPS)
      DO jc = 1,icc
        ic = idxc(jc)
        zem2(jc) = -MIN(ptau(ic,jk) * prmuzi(ic),500._wp)
      END DO

      zem2(1:icc) = EXP(zem2(1:icc))


!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
      DO jc = 1,icc
        ic = idxc(jc)

        za  = zgamma1(ic) * prmuz(ic) 
        za1 = za - zgamma3(ic)
        zgt = zgamma1(ic) * ptau(ic,jk)

        ! collimated beam

        ztemp=1.0_wp/(1._wp + zgt)
        pref(ic,jk) = (zgt - za1 * (1._wp - zem2(jc))) *ztemp
        ptra(ic,jk) = 1._wp - pref(ic,jk)

        ! isotropic incidence

        prefd(ic,jk) = zgt *ztemp
        ptrad(ic,jk) = 1._wp - prefd(ic,jk)        
      END DO

      !-- non-conservative scattering

      !-- Homogeneous reflectance and transmittance

!IBM* ASSERT(NODEPS)
      DO jc = 1,icn
        ic = idxn(jc)
        zzz = zgamma1(ic)**2 - zgamma2(ic)**2
        zrk(jc) = SQRT(MAX(zzz,replog))

        zep1(jc) = MIN(zrk(jc) * ptau(ic,jk), 500._wp)
        zep2(jc) = MIN(ptau(ic,jk) * prmuzi(ic),500._wp)
      END DO

      zep1(1:icn) = EXP(zep1(1:icn))
      zep2(1:icn) = EXP(zep2(1:icn))
      zem1(1:icn) = 1.0_wp/zep1(1:icn)
      zem2(1:icn) = 1.0_wp/zep2(1:icn)

!IBM* ASSERT(NODEPS)
      DO jc = 1,icn
        ic = idxn(jc)
        zw  =pw(ic,jk)

        zgamma4 = 1._wp - zgamma3(ic)

        za1 = zgamma1(ic) * zgamma4 + zgamma2(ic) * zgamma3(ic)
        za2 = zgamma1(ic) * zgamma3(ic) + zgamma2(ic) * zgamma4

        zrp = zrk(jc) * prmuz(ic)               
        zrp1 = 1._wp + zrp
        zrm1 = 1._wp - zrp
        zrk2 = 2._wp * zrk(jc)
        zrpp = 1._wp - zrp*zrp
        zrkg = zrk(jc) + zgamma1(ic)
        zr1  = zrm1 * (za2 + zrk(jc) * zgamma3(ic))
        zr2  = zrp1 * (za2 - zrk(jc) * zgamma3(ic))
        zr3  = zrk2 * (zgamma3(ic) - za2 * prmuz(ic) )
        zr4  = zrpp * zrkg
        zr5  = zrpp * (zrk(jc) - zgamma1(ic))
        zt1  = zrp1 * (za1 + zrk(jc) * zgamma4)
        zt2  = zrm1 * (za1 - zrk(jc) * zgamma4)
        zt3  = zrk2 * (zgamma4 + za1 * prmuz(ic) )
        zt4  = zr4
        zt5  = zr5
        zbeta = - zr5 / zr4

        ! collimated beam

        zdenr = zr4*zep1(jc) + zr5*zem1(jc)
        pref(ic,jk) = zw  * (zr1*zep1(jc) - zr2*zem1(jc) - zr3*zem2(jc)) / zdenr

        zdent = zt4*zep1(jc) + zt5*zem1(jc)
        ptra(ic,jk) = zem2(jc) * (1._wp - zw  * (zt1*zep1(jc) - zt2*zem1(jc) - zt3*zep2(jc)) / zdent)

        ! diffuse beam

        zemm = zem1(jc)*zem1(jc)
        zdend = 1._wp / ( (1._wp - zbeta*zemm ) * zrkg)
        prefd(ic,jk) =  zgamma2(ic) * (1._wp - zemm) * zdend
        ptrad(ic,jk) =  zrk2*zem1(jc)*zdend

      END DO

!IBM* ASSERT(NODEPS)
      DO jc = 1,icf
        ic=idxf(jc)
        pref(ic,jk) =0.0_wp
        ptra(ic,jk) =1.0_wp
        prefd(ic,jk)=0.0_wp
        ptrad(ic,jk)=1.0_wp
      END DO
    ENDDO

    !     ------------------------------------------------------------------
  END SUBROUTINE srtm_reftra_ec
   ! -----------------
!>
!! @brief Delta-scale optical properties
!!
  SUBROUTINE delta_scale_1d(kproma, kbdim, tau, g, w0)
    INTEGER,  INTENT(IN   ) :: kproma                   ! number of columns
    INTEGER,  INTENT(IN   ) :: kbdim                    ! max. number of col.
    REAL(WP), INTENT(INOUT) :: tau(kbdim), g(kbdim), w0(kbdim)
    
    REAL(WP) :: wf(kbdim)
    INTEGER  :: i 
    
!IBM* ASSERT(NODEPS)
    DO i = 1, kproma
      wf(i) = w0(i) * g(i) * g(i) ! f = g**2
    
      tau(i) = (1.0_wp - wf(i)) * tau(i)
      w0 (i) = (w0(i) - wf(i)) / (1.0_wp - wf(i))
      ! Special case for f = g**2; generally g = (g - f) / (1 - f) 
      g  (i) = g(i)  / (1.0_wp + g(i))
    END DO 
  END SUBROUTINE delta_scale_1d
   ! -----------------
  SUBROUTINE delta_scale_2d(kproma, kbdim, klev, tau, g, w0, update)
    INTEGER,  INTENT(IN   ) :: kproma, klev          ! number of columns
    INTEGER,  INTENT(IN   ) :: kbdim                 ! max number of col.
    REAL(WP), INTENT(INOUT) :: tau(kbdim,klev), g(kbdim,klev), w0(kbdim,klev)
    LOGICAL,  INTENT(IN   ), OPTIONAL :: update(kbdim,klev)
    
    REAL(WP) :: wf(kbdim,klev)
    LOGICAL  :: do_this(kbdim, klev) 
        
    do_this(:,:) = .true. 
    IF(PRESENT(update)) do_this(1:kproma,1:klev) = update(1:kproma,1:klev)  
    
    WHERE(do_this(1:kproma,1:klev)) 
      wf(1:kproma,1:klev) = w0(1:kproma,1:klev) * g(1:kproma,1:klev) * g(1:kproma,1:klev) ! f = g**2

      tau(1:kproma,1:klev) = (1.0_wp - wf(1:kproma,1:klev)) * tau(1:kproma,1:klev)
      w0 (1:kproma,1:klev) = (w0(1:kproma,1:klev) - wf(1:kproma,1:klev)) / (1.0_wp - wf(1:kproma,1:klev))
      ! Special case for f = g**2; generally g = (g - f) / (1 - f) 
      g  (1:kproma,1:klev) = g(1:kproma,1:klev)  / (1.0_wp + g(1:kproma,1:klev))
    END WHERE
  END SUBROUTINE delta_scale_2d
   ! -----------------
END MODULE mo_srtm_solver
