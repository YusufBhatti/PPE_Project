!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_aero_settling
!! 
!! \brief
!!  This routine vertically redistributes condensed phase HNO3.
!!
!! \author Douglas E. Kinnison (NCAR) 
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de 
!!
!! \revision_history
!!  - DK: original version (1999-11-08)
!!  - DK: Modified for WACCM2, removed ICE (2004-09-03)
!!  - DK: Modified for WACCM3, added NAT 7um mode (2004-11-08)
!!
!! \limitations
!!  none
!!
!! \details
!!  For each aerosol type the terminal velocity is calculated. This
!!  quantity is dependent on: 1) the mass density of the aerosol;
!!  2) the radius; 3) the dynamic viscosity; 4) shape; and the Cunningham
!!  correction factor for spherical particles. See Fuchs, The Mechanics
!!  of Aerosols Oxford, Pergmann Press, pp 27-31, 1964 and Kasten,
!!  Falling Speed of Aerosol Particles, J. of Appl. Met., 7, 944-947,
!!  1968 for details. For aerosol with a radius of 3 microns (e.g., NAT)
!!  and 10 microns (e.g., ICE) the terminal velocity (cm sec-1) is  0.2
!!  and 1.7 respectively.
!!
!!  The flux of condensed phase HNO3 is then derived using the
!!  following equation:
!!         Flux (molec cm-2 s-1) = V * C * exp (8*ln^2 sigma).
!!
!!         where: V is terminal velocity (cm sec-1)
!!                C is condensed phase conc (molecules cm-3)
!!                sigma is the width of the log normal distribution.
!!
!!  The approach of settling the entire aerosol size distribution is
!!  based on work of Considine et al., JGR, 1999.
!!
!!  The routine is a straighforward approach, starting at the top zone
!!  (highest altitude) and progressing towards the surface. The gross
!!  settling of condensed phase HNO3 is simply the quanity:
!!  Flux * dt/dz.
!!
!! \bibliographic_references
!!  - D.E. Kinnison et al., JGR 112, D20302, doi:10.1029/2006JD007879, 2007
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ARGUMENTS
!
!      All of the following components are grid center quanitities.
!
!   INPUT:
!      ad            Air Density (molecules cm-3)
!      press         Pressure (Pascals)
!      timestep      Gross chemistry timestep (in seconds)
!      temp          Temperature (K)
!      hno3_cond     Condensed phase HNO3 (mole fraction)
!      radius_nat    Mean radius of NAT (cm)
!      zstar         log pressure altitude coordinate (km)
!
!   OUTPUT:
!      hno3_cond     vertically modified
!
!----------------------------------------------------------------------------------

      module mo_moz_aero_settling

      USE mo_kind,     ONLY : dp

      private
      public  :: strat_aer_settling

      contains

      subroutine strat_aer_settling( ad, press, timestep, zstar, temp, &
                                     hno3_cond, radius_nat, ncol, lchnk)
!                                    hno3_cond, radius_nat, ncol, lchnk, aero_ndx )

!++mgs: replaced mo_grid and chem_mods with mo_moz_mods
      USE mo_moz_mods,  ONLY : pcols => plonl, pver => plev
!!baustelle!! Use adv_mass for molecular masses ?
!!      use mo_moz_mods,  only : adv_mass
!--mgs
      USE mo_physical_constants, ONLY : grav 
!     use history,      only : outfld

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)     ::  ncol                                    ! columns in chunk
      integer, intent(in)     ::  lchnk                                   ! chunk id
!     integer, intent(in)     ::  aero_ndx                                ! aerosol index
      real(dp), intent(in)    ::  timestep                                ! model time step (s)
      real(dp), intent(in)    ::  ad(ncol,pver)                           ! Air density (molecules cm-3)
      real(dp), intent(in)    ::  radius_nat(ncol,pver)                   ! Mean radius of NAT (cm)
      real(dp), intent(in)    ::  zstar(ncol,pver)                        ! altitude (km)
      real(dp), intent(in)    ::  press(pcols,pver)                       ! Pressure (Pa)
      real(dp), intent(in)    ::  temp(pcols,pver)                        ! temperature (K)
      real(dp), intent(inout) ::  hno3_cond(ncol,pver)                    ! Condensed Phase HNO3 (VMR)

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      real(dp), parameter :: avo_num       = 6.022e23_dp, &    ! molecules/mole
                             MW_air        = 28.8_dp, &        ! grams/mole air
                             nat_dens      = 1.6_dp, &         ! g/cm^3
                             shape_fac_nat = 1.0_dp, &         ! TBD
                             sigma_nat     = 1.6_dp, &         ! Width of distribution
                             av_const      = 2.117265e4_dp, &  ! (8*8.31448*1000 / PI)
                             km2cm         = 1.e5_dp, &        ! km to cm
                             m2cm          = 1.e2_dp, &        ! m to cm
                             c1            = 2._dp/9._dp, &
                             gravity       = grav*m2cm       ! (cm/s^2)

      integer  :: i, k, kp1
      real(dp) :: Cc_nat                    ! Cunningham Correction Factor
      real(dp) :: dt_dz                     ! dt / dz, sec cm-1
      real(dp) :: flux_nat                  ! aerosol flux, molec cm-2 sec-1
      real(dp) :: mean_vel                  ! mean velocity, cm sec-1
      real(dp) :: mfp                       ! Mean Free Path
      real(dp) :: vel_nat                   ! Terminal velocity (cm sec-1)
      real(dp) :: rad_nat                   ! wrk radius nat (cm)
      real(dp) :: visc                      ! Dynamic viscosity of air
      real(dp) :: depos                     ! molecules/cm**3 deposited
      real(dp) :: atm_dens, atm_densa       ! total atm density and inverse
      real(dp) :: cond_hno3                 ! wrk variables
      real(dp) :: const_nat                 ! wrk variables
      real(dp) :: t                         ! working temperatue
      real(dp) :: velnat(ncol,pver)         ! holding variable for output
      logical  :: lon_mask(ncol)            ! longitude logic mask

      const_nat = exp( 8._dp*(log(sigma_nat))**2 )
      do k = 1,pver
	 velnat(:,k) = 0._dp
      end do


!----------------------------------------------------------------------
!     ... derive Aerosol Settling (explicit approach)
!----------------------------------------------------------------------
Level_loop : &
      do k = 1,pver-1
!----------------------------------------------------------------------
!     ... operate between 2.0hPa and 300hPa and only where nat exist
!----------------------------------------------------------------------
	 kp1 = k + 1
	 lon_mask(:) = press(:ncol,k) >= 2.e2_dp .and. press(:ncol,k) <= 300.e2_dp &
                       .and. (radius_nat(:,k) > 0._dp )
	 if( any( lon_mask(:) ) ) then
Column_loop : &
            do i = 1,ncol
	       if( lon_mask(i) ) then
                  t = temp(i,k)
!----------------------------------------------------------------------
!     ... General Setup for NAT
!         Calculate the settling of the NAT Aerosol
!         NOTE: Index "k" is the box that is being calculated
!               Index "k-1" is flux from above into the k box
!               A positive NET {flux*dt/dz} adds to the "k" box
!----------------------------------------------------------------------
!     ... mean Molecular Velocity, cm sec-1
!----------------------------------------------------------------------
                  mean_vel  = sqrt( av_const*t/MW_air )*100._dp
!----------------------------------------------------------------------
!     ... dynamic Viscosity, g cm-1 sec-1
!----------------------------------------------------------------------
	          visc = (t*1.458e-6_dp)**1.5_dp /(t + 110.4_dp)*10000._dp
!----------------------------------------------------------------------
!     ... mean Free Path, cm
!----------------------------------------------------------------------
		  atm_dens  = ad(i,k)
		  atm_densa = 1._dp/atm_dens
                  mfp = visc* avo_num / (.499_dp*atm_dens*MW_air*mean_vel)
!----------------------------------------------------------------------
!     ... dt / dz, sec cm-1; NOTE: zstar is in km
!----------------------------------------------------------------------
                  dt_dz = timestep / ((zstar(i,k-1) - zstar(i,k))*km2cm)
!----------------------------------------------------------------------
!     ... calculate NAT Aerosol Settling
!----------------------------------------------------------------------
		  rad_nat   = radius_nat(i,k)
                  cond_hno3 = hno3_cond(i,k)*atm_dens
!----------------------------------------------------------------------
!     ... Cunningham Correction Factor, Unitless
!----------------------------------------------------------------------
                  Cc_nat = 1._dp + (mfp/rad_nat)*(1.246_dp + .42_dp*exp( -.87_dp*rad_nat/mfp ))
!----------------------------------------------------------------------
!     ... terminal Velocity of Aerosol, cm sec-1
!----------------------------------------------------------------------
	          vel_nat = c1*rad_nat**2 * nat_dens*gravity*Cc_nat/visc*shape_fac_nat
		  velnat(i,k) = vel_nat
!----------------------------------------------------------------------
!     ... aerosol Flux, Cond-phase molecules cm-2 sec-1
!----------------------------------------------------------------------
         	  flux_nat = cond_hno3*vel_nat* const_nat
!----------------------------------------------------------------------
!     ... calculate NAT Aerosol Settling (i.e., HNO3 redistribution)
!----------------------------------------------------------------------
		  depos = min( cond_hno3,flux_nat*dt_dz )
!----------------------------------------------------------------------
!     ... modify the HNO3_cond in level "k"
!----------------------------------------------------------------------
	          hno3_cond(i,k) = (cond_hno3 - depos)*atm_densa
!----------------------------------------------------------------------
!     ... modify the HNO3_cond in level "k+1"
!----------------------------------------------------------------------
		  hno3_cond(i,kp1) = hno3_cond(i,kp1) + depos/ad(i,kp1)
               end if
            end do Column_loop
         end if
      end do Level_loop

!----------------------------------------------------------------------
!	... output nat velocity
!----------------------------------------------------------------------
!     if( aero_ndx == 1 ) then
!        call outfld( 'VEL_NAT1', velnat, ncol, lchnk )
!     else if( aero_ndx == 2 ) then
!        call outfld( 'VEL_NAT2', velnat, ncol, lchnk )
!     end if

      end subroutine strat_aer_settling

      end module mo_moz_aero_settling
