!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_strato_sad
!!
!! \brief
!!  This module handles the stratospheric aerosol processes in MOZART
!!
!! \author Douglas E. Kinnison (NCAR)
!! \author Stacy Walters (NCAR)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - DK: original version (1999-10-14)
!!
!! \limitations
!!  none
!!
!! \details
!!  This routine has the logic to derive the surface area density for
!!  three types of aerosols: Sulfate; Nitric Acid Trihydrate (NAT);
!!  and ICE. The surface area density is stored in sad_strat(3). The
!!  first, second, and third dimensions are SULFATE, NAT, and ICE SAD
!!  respectively. The effective radius of each particle is also stored
!!  in radius_strat(3).
!!
!!  NOTE1: For Sulfate and H2O ICE
!!  The Surface Area and Volume Densities are calculated from the
!!  second and third moment of the LOG Normal Distribution. For an example
!!  see Binkowski et al., JGR, 100, 26191-26209, 1995. The Volume Density
!!  is substituted into the SAD equation so that the SAD is dependent on
!!  the # of particles cm-3, the width of the distribution (sigma), and
!!  the volume density of the aerosol. For the ternary solution calculation
!!  the total sulfate mass is derived from the SAGEII SAD data. This approach
!!  has been previously used in Considine et al., JGR, 1999. The thermodynamic
!!  models used in this routine are from A. Tabazedeh.
!!
!!  NOTE2: For NAT, there are two modes centered at: 1) 0.5um radius; 2) 6.5um
!!  radius. This is based on ER2 measurements during NASA SOLVE, see Fahey et al., 
!!  Science, 291, 1026-1031, 2001. . The number density of the large mode is set to 2.3e-4
!!  Particle cm-3 and the condensed phase HNO3 mass assigned to this mode to 
!!  produce particles with a mode radius of 6.5 microns. Any additional
!!  HNO3 is then assigned to the smaller mode. Denitrification occurs on the
!!  larger mode. The HNO3 in the smaller mode is assigned to STS aerosols.
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


module mo_moz_strato_sad

      USE mo_math_constants,      ONLY : pi
      use mo_moz_mods,            ONLY : pcols => plonl, pver => plev
      use mo_moz_strato_sad_data, ONLY : a, b
      USE mo_exception,           ONLY: finish, message, message_text, em_info, em_debug, &
                                        em_error
      USE mo_kind,                ONLY: dp

      implicit none

      private
      public  :: sad_inti
      public  :: sad_strat_calc
      public  :: sad_top

      save

      real(dp), parameter :: one_thrd = 1._dp/3._dp
      real(dp), parameter :: two_thrd = 2._dp/3._dp
      real(dp), parameter :: four_pi  = 4._dp*pi

      integer :: sad_top
      integer :: sad_topp

      contains

      subroutine sad_inti
!----------------------------------------------------------------------
!     ... initialize the sad module
!----------------------------------------------------------------------

      USE mo_control,   ONLY: nvclev, vct

      implicit none

!---------------------------------------------------------------------- 
!	... Local variables
!---------------------------------------------------------------------- 
      integer  ::  k
      real(dp) ::  pmid
      character(10) :: num_string

!---------------------------------------------------------------------- 
!	... find level where etamids are all > 1 hPa
!---------------------------------------------------------------------- 
      sad_top = 0
      do k = pver,1,-1
         pmid=0.5*(vct(k)+100000._dp*vct(nvclev+k)+ &
                      vct(k+1)+100000._dp*vct(nvclev+k+1))
         if( pmid < 100._dp ) then
             sad_top = k
             exit
         end if
      end do
      sad_topp = sad_top + 1

      write(num_string,'(i4)')sad_top
      CALL message('','sad_inti: sad capped at level'//num_string)
      write(num_string,'(es10.2)')pmid*.01
      CALL message('',' whose midpoint is '//num_string//' hPa')

      end subroutine sad_inti
!===============================================================================
! ROUTINE
!   sad_strat_calc
!
!   Date...
!     14 October 1999
!
!   Programmed by...
!     Douglas E. Kinnison
!   Modified by
!     Stacy Walters
!     1 November 1999
!   Modified by 
!     Doug Kinnison
!     1 September 2004; Condensed phase H2O passed in from CAM
!     2 November 2004; New treatment of denitrification (NAT)
!    14 November 2004; STS mode of operation.
!
! ARGUMENTS
!   INPUT:
!     hno3_gas		Nitric Acid gas phase abundance (mole fraction)
!     hno3_cond(3)	Nitric Acid condensed phase abundance (mole fraction)
!                       (1) in STS; (2) in 0.5um NAT; (3) 6.5um NAT
!     h2o_cond		Water condensed phase           (mole fraction)
!     h2o_gas           Water gas-phase abundance       (mole fraction)

!     sage_sad		SAGEII surface area density     (cm2 cm-3)
!     m                 Airdensity                      (molecules cm-3)
!     press             Pressure                        (hPa, i.e., mb)

!
!   OUTPUT:
!
!	hno3_gas     = Gas-phase HNO3... Used in chemical solver. 
!       hno3_cond(1) = Condensed HNO3 from STS... Not used in mo_moz_aero_settling
!       hno3_cond(2) = Condensed HNO3 from NAT, sm mode... Not used in mo_moz_aero_settling
!       hno3_cond(3) = Condensed HNO3 from NAT, lg mode... USED in mo_moz_aero_settling
!
!       SAD_strat(1) = Sulfate Aerosol... Used in heterogeneous rate module 
!       SAD_strat(2) = NAT small mode.... Used in heterogeneous rate module
!       SAD_strat(3) = NAT large mode.... Not used in heterogeneous rate module.
!       SAD_strat(4) = Water-Ice......... Used in heterogeneous rate module
!
!       RAD_strat(1) = Sulfate Aerosol... Used in heterogeneous rate module
!       RAD_strat(2) = NAT small mode.... Not used in mo_moz_aero_settling (0.5um)
!       RAD_strat(3) = NAT large mode.... Used in mo_moz_aero_settling (6.5um)
!       RAD_strat(4) = Water-Ice......... Not used in heterogeneous rate module
!
!   NOTE1: The sum of hno3_cond(1-3) will be added to hno3_gas for advection of HNO3 in
!         WACCM3.
!
!   NOTE2: This routine does not partition H2O.
!
!
! ROUTINES CallED (in and below this routine):
!
!	sad_strat_calc
!         nat_sat_temp ...... derives the NAT saturation temp
!
!         sulfate_sad_calc .. Calculates the sulfate SAD;  HNO3, H2O cond. phase
!         calc_radius_lbs ... Calculates the radius for a H2SO4 Binary soln. (T>200K)
!         sad_to_h2so4 ...... Derives H2SO4 abundance (micrograms m-3)
!                             from SAGEII SAD.
!         density............ A. Tabazedeh binary thermo model
!         equil.............. A. Tabazedeh ternary thermo. model
!
!         nat_sad_calc....... Calculates the NAT SAD; HNO3, H2O cond. phase
!         nat_cond........... Derives the NAT HNO3 gas/cond partitioning
!
!         ice_sad_calc....... derives the ICE SAD and H2O gas/cond partitioning
!===============================================================================

      subroutine sad_strat_calc( lchnk, m, press, temper, hno3_gas, &
                                 hno3_cond, h2o_gas, h2o_cond, sad_sage, radius_strat, &
                                 sad_strat, ncol)

      implicit none

!-------------------------------------------------------------------------------
!	... dummy arguments
!-------------------------------------------------------------------------------
      integer, intent(in)     ::  lchnk                      ! chnk id
      integer, intent(in)     ::  ncol                       ! columns in chunk
      real(dp), intent(in)    ::  m(ncol,pver)               ! Air density (molecules cm-3)
      real(dp), intent(in)    ::  sad_sage(ncol,pver)        ! SAGEII surface area density (cm2 aer. cm-3 air)
      real(dp), intent(in)    ::  press(ncol,pver)           ! Pressure, hPa
      real(dp), intent(in)    ::  temper(pcols,pver)         ! Temperature (K)
      real(dp), intent(in)    ::  h2o_gas(ncol,pver)         ! H2O gas-phase (mole fraction)
      real(dp), intent(in)    ::  h2o_cond(ncol,pver)        ! H2O condensed phase (mole fraction)
      real(dp), intent(inout) ::  hno3_gas(ncol,pver)        ! HNO3 condensed phase (mole fraction)
      real(dp), intent(inout) ::  hno3_cond(ncol,pver,3)     ! HNO3 condensed phase (mole fraction)
      real(dp), intent(out)   ::  radius_strat(ncol,pver,4)  ! Radius of Sulfate, NAT_sm, NAT_lg, and ICE (cm)
      real(dp), intent(out)   ::  sad_strat(ncol,pver,4)     ! Surface area density of Sulfate, NAT_sm, NAT_lg, ICE (cm2 cm-3)

!-------------------------------------------------------------------------------
!	... local variables
!-------------------------------------------------------------------------------
      real(dp), parameter :: temp_floor = 0._dp

      integer  ::  i, k, n
      integer  ::  dims(1)
      real(dp) ::  hno3_total(ncol,pver)         ! HNO3 total = gas-phase + condensed
      real(dp) ::  h2o_total(ncol,pver)          ! H2O total  = gas-phase + condensed
      real(dp) ::  radius_lbs(ncol,pver)         ! Radius of Liquid Binary Sulfate (cm)
      real(dp) ::  radius_sulfate(ncol,pver)     ! Radius of Sulfate aerosol (cm)
      real(dp) ::  radius_nat_sm(ncol,pver)      ! Radius of NAT aerosol at 0.5um  (cm)
      real(dp) ::  radius_nat_lg(ncol,pver)      ! Radius of NAT aerosol at 6.5um  (cm)
      real(dp) ::  radius_ice(ncol,pver)         ! Radius of ICE aerosol     (cm)
      real(dp) ::  sad_nat_sm(ncol,pver)         ! SAD of NAT aerosol at 0.5um (cm2 cm-3)
      real(dp) ::  sad_nat_lg(ncol,pver)         ! SAD of NAT aerosol at 6.5um (cm2 cm-3)
      real(dp) ::  sad_sulfate(ncol,pver)        ! SAD of Sulfate aerosol    (cm2 cm-3)
      real(dp) ::  sad_ice(ncol,pver)            ! SAD of ICE aerosol        (cm2 cm-3)
      real(dp) ::  tsat_nat(ncol,pver)           ! Temperature for NAT saturation
      real(dp) ::  h2o_avail(ncol,pver)          ! H2O temporary arrays
      real(dp) ::  hno3_avail(ncol,pver)         ! HNO3 temporary array
      real(dp) ::  hno3_gas_nat(ncol,pver)       ! HNO3 after call to NAT routines
      real(dp) ::  hno3_gas_sulf(ncol,pver)      ! HNO3 after call to STS routines
      real(dp) ::  hno3_cond_nat_sm(ncol,pver)   ! HNO3 condensed after call to NAT, small mode
      real(dp) ::  hno3_cond_nat_lg(ncol,pver)   ! HNO3 condensed after call to NAT, large mode
      real(dp) ::  hno3_cond_sulf(ncol,pver)     ! HNO3 condensed after call to STS routines
      real(dp) ::  temp(pcols,pver)              ! wrk temperature array

      logical  ::  z_val(ncol)
      logical  ::  mask(ncol,pver)
      logical  ::  mask1(ncol,pver)
      logical  ::  mask2(ncol,pver)

!----------------------------------------------------------------------
!     ... initialize to zero
!----------------------------------------------------------------------
      do n = 1,4
         do k = 1,pver
            radius_strat(:,k,n) = 0._dp
            sad_strat(:,k,n)    = 0._dp
         end do
      end do

!----------------------------------------------------------------------
!     ... limit temperature
!----------------------------------------------------------------------
      do k = 1,pver
         temp(:ncol,k) = max( temp_floor,temper(:ncol,k) )
      end do

!----------------------------------------------------------------------
!     ... total HNO3 and H2O gas and condensed phases
!----------------------------------------------------------------------
      do k = sad_topp,pver
         hno3_total(:,k) = hno3_gas(:,k) + hno3_cond(:,k,1) + hno3_cond(:,k,2) + hno3_cond(:,k,3)
         h2o_total(:,k)  = h2o_gas(:,k)  + h2o_cond(:,k)
      end do

!----------------------------------------------------------------------
!     ... Derive NAT Saturation Temperature
!     ... using total H2O - assuming NAT forms before ICE
!----------------------------------------------------------------------
      call nat_sat_temp( ncol, hno3_total, h2o_total, press, tsat_nat )

!----------------------------------------------------------------------
!     ... Set SAD to SAGEII if Temperature is GT 200K and Return
!     ... mask2 = true  .... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!     ... mask1 = false .... T <= TSAT_NAT
!     ... mask  = false .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
         mask2(:,k) = temp(:ncol,k) > 200._dp .or. sad_sage(:,k) <= 1.e-15_dp &
                      .or. press(:ncol,k) < 2._dp .or. press(:ncol,k) > 300._dp
      end do

sage_sad : &
      if( any( mask2(:,sad_topp:pver) ) ) then
         do k = sad_topp,pver
            where( mask2(:,k) )
               sad_strat(:,k,1)    = sad_sage(:,k)
               sad_strat(:,k,2)    = 0._dp
               sad_strat(:,k,3)    = 0._dp
               sad_strat(:,k,4)    = 0._dp
            endwhere
         end do
!----------------------------------------------------------------------
!     ... Calculate Liquid Binary Sulfate Radius
!----------------------------------------------------------------------
         call calc_radius_lbs( ncol, mask2, sad_sage, radius_lbs )
         do k = sad_topp,pver
            where( mask2(:,k) )
               radius_strat(:,k,1)     = radius_lbs(:,k)
               radius_strat(:,k,2)     = 0._dp
               radius_strat(:,k,3)     = 0._dp
               radius_strat(:,k,4)     = 0._dp
               hno3_gas(:,k)           = hno3_total(:,k)
               hno3_cond(:,k,1)        = 0._dp
               hno3_cond(:,k,2)        = 0._dp
               hno3_cond(:,k,3)        = 0._dp
            endwhere
         end do
         if( all( mask2(:,sad_topp:pver) ) ) then
            return
         end if
      end if sage_sad
!----------------------------------------------------------------------
!     ... Logic for deriving sad for sulfate, nat, and ice aerosol
!         Ice formation occurs here if condensed phase H2O exists (from CAM).
!         mask2 = false.... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!         mask1 = true .... T <= TSAT_NAT
!         mask  = true .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
         do i = 1,ncol
            if( .not. mask2(i,k) ) then
              mask(i,k) = (h2o_cond(i,k) > 0._dp)
            else
               mask(i,k) = .false.
            end if
         end do
      end do

all_sad : &
      if( any( mask(:,sad_topp:pver) ) ) then
         do k = sad_topp,pver
            where( mask(:,k) )
               h2o_avail(:,k) = h2o_cond(:,k)
            endwhere
         end do

         call ice_sad_calc( ncol, press, temp, m, h2o_avail, &
                            sad_ice, radius_ice, mask )

         do k = sad_topp,pver
            where( mask(:,k) )
               h2o_avail(:,k)  = h2o_gas(:,k)
               hno3_avail(:,k) = hno3_total(:,k)
            endwhere
            if( any(mask(:,k)) ) then
               where( mask(:,k) )
                  z_val(:) = hno3_avail(:,k) == 0._dp
               elsewhere
                  z_val(:) = .false.
               endwhere
               if( any( z_val(:) ) ) then
                  WRITE (message_text,'(a,i0,1x,i0)') 'Before NAT_SAD_CALC_1 has zero hno3_avail at lchnk,k = ',lchnk,k
                  CALL message('sad_strat_calc', message_text, level=em_debug)
               end if
            end if
         end do

         call nat_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, m, &
                            hno3_gas_nat, hno3_cond_nat_sm, hno3_cond_nat_lg, &
                            sad_nat_sm, sad_nat_lg, &
                            radius_nat_sm, radius_nat_lg, mask )

!----------------------------------------------------------------------
!  	... H2O Gas available for forming STS
!           The HNO3 avail for the STS aerosol is the amount in the gas
!           phase, plus any additional HNO3 not contained in the large
!           NAT mode at 6.5um.
!----------------------------------------------------------------------
         do k = sad_topp,pver
            where( mask(:,k) )
               h2o_avail(:,k)   = h2o_gas(:,k)
               hno3_avail(:,k)  = hno3_gas_nat(:,k)+hno3_cond_nat_sm(:,k)
            endwhere
            if( any(mask(:,k)) ) then
               where( mask(:,k) )
                  z_val(:) = hno3_avail(:,k) == 0._dp
               elsewhere
                  z_val(:) = .false.
               endwhere
               if( any( z_val(:) ) ) then
                  WRITE (message_text,'(a,i0,1x,i0)') 'After NAT_SAD_CALC_1 has zero hno3_avail at lchnk,k = ',lchnk,k
                  CALL message('sad_strat_calc', message_text, level=em_debug) 
              end if
            end if
         end do

         call sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, sad_sage, m, &
                                hno3_gas_sulf, hno3_cond_sulf, &
                                sad_sulfate, radius_sulfate, mask, lchnk, 1 )
         do k = sad_topp,pver
            where( mask(:,k) )
               sad_strat(:,k,1)    = sad_sulfate(:,k)
               sad_strat(:,k,2)    = 0.0_dp
               sad_strat(:,k,3)    = sad_nat_lg(:,k)
               sad_strat(:,k,4)    = sad_ice(:,k)
               radius_strat(:,k,1) = radius_sulfate(:,k)
               radius_strat(:,k,2) = 0.0_dp
               radius_strat(:,k,3) = radius_nat_lg(:,k)
               radius_strat(:,k,4) = radius_ice(:,k)
               hno3_gas(:,k)       = hno3_gas_sulf(:,k)
               hno3_cond(:,k,1)    = hno3_cond_sulf(:,k)  ! de-NOx the atmosphere
               hno3_cond(:,k,2)    = 0.0_dp
               hno3_cond(:,k,3)    = hno3_cond_nat_lg(:,k) ! de-NOY the atmosphere
            endwhere
         end do
      end if all_sad
!----------------------------------------------------------------------
!  	... Saturation of NAT occurs here it temp LE NAT saturation temp.
!           There is no condensed phase H2O (mask is false)
!           mask2 = false .... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!           mask1 = true  .... T <= TSAT_NAT
!           mask  = false .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
         do i = 1,ncol
            if( .not. mask2(i,k) .and. .not. mask(i,k) ) then
               mask1(i,k)   = temp(i,k) <= tsat_nat(i,k)
            else
               mask1(i,k) = .false.
            end if
         end do
      end do

nat_sad : &
      if( any( mask1(:,sad_topp:pver) ) ) then

         do k = sad_topp,pver
            where( mask1(:,k) )
               h2o_avail(:,k)  = h2o_gas(:,k)
               hno3_avail(:,k) = hno3_total(:,k)
            endwhere
            if( any(mask1(:,k)) ) then
               where( mask1(:,k) )
                  z_val(:) = hno3_avail(:,k) == 0._dp
               elsewhere
                  z_val(:) = .false.
               endwhere
               if( any( z_val(:) ) ) then
                  WRITE (message_text,'(a,i0,1x,i0)') 'Before NAT_SAD_CALC_2 has zero hno3_avail at lchnk,k = ',lchnk,k
                  CALL message('sad_strat_calc', message_text, level=em_debug) 
               end if
            end if
         end do

         call nat_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, m, &
                            hno3_gas_nat, hno3_cond_nat_sm, hno3_cond_nat_lg, &
                            sad_nat_sm, sad_nat_lg, &
                            radius_nat_sm, radius_nat_lg, mask1 )
!----------------------------------------------------------------------
!  	... H2O Gas available for forming STS
!           The HNO3 avail for the STS aerosol is the amount in the gas
!           phase, plus any additional HNO3 not contained in the large
!           NAT mode at 6.5um.
!----------------------------------------------------------------------
         do k = sad_topp,pver
            where( mask1(:,k) )
               h2o_avail(:,k)    = h2o_gas(:,k)
               hno3_avail(:,k)   = hno3_gas_nat(:,k)+hno3_cond_nat_sm(:,k)
            endwhere
            if( any(mask1(:,k)) ) then
               where( mask1(:,k) )
                  z_val(:) = hno3_avail(:,k) == 0._dp
               elsewhere
                  z_val(:) = .false.
               endwhere
               if( any( z_val(:) ) ) then
                  WRITE (message_text,'(a,i0,1x,i0)') 'After NAT_SAD_CALC_2 has zero hno3_avail at lchnk,k = ',lchnk,k
                  CALL message('sad_strat_calc', message_text, level=em_debug) 
               end if
            end if
         end do
         call sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, sad_sage, m, &
                                hno3_gas_sulf, hno3_cond_sulf, &
                                sad_sulfate, radius_sulfate, mask1, lchnk, 2 )
         do k = sad_topp,pver
            where( mask1(:,k) )
               sad_strat(:,k,1)    = sad_sulfate(:,k)
               sad_strat(:,k,2)    = 0._dp
               sad_strat(:,k,3)    = sad_nat_lg(:,k)
               sad_strat(:,k,4)    = 0._dp
               radius_strat(:,k,1) = radius_sulfate(:,k)
               radius_strat(:,k,2) = 0._dp
               radius_strat(:,k,3) = radius_nat_lg(:,k)
               radius_strat(:,k,4) = 0._dp
               hno3_gas(:,k)       = hno3_gas_sulf(:,k)
               hno3_cond(:,k,1)    = hno3_cond_sulf(:,k)  !de-NOx atmosphere
               hno3_cond(:,k,2)    = 0.0_dp
               hno3_cond(:,k,3)    = hno3_cond_nat_lg(:,k) !de-NOy atmosphere
            endwhere
         end do
      end if nat_sad
!----------------------------------------------------------------------
!          ... Sulfate Ternary (only)
!           No condensed phase H2O (mask = false)
!           No condensed phase NAT (mask1 = false)
!           mask2 = false .... T > 200K or SAD_SULF < 1e-15 or P >2hPa or P <300hPa
!           mask1 = false .... T <= TSAT_NAT
!           mask  = false .... H2O_COND > 0.0
!----------------------------------------------------------------------
      do k = sad_topp,pver
         mask(:,k) = .not. mask(:,k) .and. .not. mask1(:,k) .and. .not. mask2(:,k)
      end do

sulfate_sad : &
      if( any( mask(:,sad_topp:pver) ) ) then
         do k = sad_topp,pver
            where( mask(:,k) )
               h2o_avail(:,k)  = h2o_gas(:,k)
               hno3_avail(:,k) = hno3_total(:,k)
            endwhere
            if( any(mask(:,k)) ) then
               where( mask(:,k) )
                  z_val(:) = hno3_avail(:,k) == 0._dp
               elsewhere
                  z_val(:) = .false.
               endwhere
               if( any( z_val(:) ) ) then
                  WRITE (message_text,'(a,i0,1x,i0)') 'Before NAT_SAD_CALC_3 has zero hno3_avail at lchnk,k = ',lchnk,k
                  CALL message('sad_strat_calc', message_text, level=em_debug) 
               end if
            end if
         end do

         call sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, sad_sage, m, & 
                                hno3_gas_sulf, hno3_cond_sulf, &
                                sad_sulfate, radius_sulfate, mask, lchnk, 3 )
         do k = sad_topp,pver
            where( mask(:,k) )
               sad_strat(:,k,1)    = sad_sulfate(:,k)
               sad_strat(:,k,2)    = 0._dp
               sad_strat(:,k,3)    = 0._dp
               sad_strat(:,k,4)    = 0._dp
               radius_strat(:,k,1) = radius_sulfate(:,k)
               radius_strat(:,k,2) = 0._dp
               radius_strat(:,k,3) = 0._dp
               radius_strat(:,k,4) = 0._dp
               hno3_gas(:,k)       = hno3_gas_sulf(:,k)
               hno3_cond(:,k,1)    = hno3_cond_sulf(:,k)
               hno3_cond(:,k,2)    = 0._dp
               hno3_cond(:,k,3)    = 0._dp
            endwhere
         end do
      end if sulfate_sad

      end subroutine sad_strat_calc

      subroutine nat_sat_temp( ncol, hno3_total, h2o_avail, press, tsat_nat )

      implicit none

!----------------------------------------------------------------------
!        ... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(dp), intent(in)  :: press(ncol,pver)
      real(dp), intent(in)  :: h2o_avail(ncol,pver)
      real(dp), intent(in)  :: hno3_total(ncol,pver)
      real(dp), intent(out) :: tsat_nat(ncol,pver)

!----------------------------------------------------------------------
!        ... local variables
!----------------------------------------------------------------------
      real(dp), parameter :: ssrNAT = 10.0_dp
      real(dp), parameter :: ssrNATi = .1_dp
      real(dp), parameter :: aa = -2.7836_dp, &
                             bb = -0.00088_dp, &
                             cc = 38.9855_dp, &
                             dd = -11397.0_dp, &
                             ee = 0.009179_dp
      integer  :: k
      real(dp) :: bbb(ncol)                     ! temporary variable
      real(dp) :: ccc(ncol)                     ! temporary variable
      real(dp) :: wrk(ncol)                     ! temporary variable
      real(dp) :: tmp(ncol)                     ! temporary variable
      real(dp) :: phno3(ncol)                   ! hno3 partial pressure
      real(dp) :: ph2o(ncol)                    ! h2o  partial pressure

!----------------------------------------------------------------------
!             ... Derive HNO3 and H2O partial pressure (torr)
!          where: 0.7501 = 760/1013.
!----------------------------------------------------------------------
      do k = sad_topp,pver
         bbb(:)   = press(:,k) * .7501_dp
         phno3(:) = hno3_total(:,k) * bbb(:)
         ph2o(:)  = h2o_avail(:,k)  * bbb(:)
         where( (phno3(:) > 0._dp) .AND. (ph2o(:) > 0._dp) )
!----------------------------------------------------------------------
!             ... Calculate NAT Saturation Threshold Temperature
!           Hanson and Mauersberger: GRL, Vol.15, 8, p855-858, 1988.
!           Substitute m(T) and b(T) into Equation (1). Rearrange and
!           solve quadratic eqation. 
!----------------------------------------------------------------------
            tmp(:) = log10( ph2o(:) )
            wrk(:) = 1._dp / (ee + bb*tmp(:))
            bbb(:) = (aa*tmp(:) - log10( phno3(:)*ssrNATi ) + cc) * wrk(:)
            ccc(:) = dd *wrk(:)
            tsat_nat(:,k) = .5_dp * (-bbb(:) + sqrt( bbb(:)*bbb(:) - 4._dp*ccc(:) ))
         elsewhere
            tsat_nat(:,k) = 0._dp
         endwhere
      end do

      end subroutine nat_sat_temp

      subroutine calc_radius_lbs( ncol, mask, sad_sage, radius_lbs )

      implicit none

!----------------------------------------------------------------------
!        ... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(dp), intent(in)  :: sad_sage(ncol,pver)
      real(dp), intent(out) :: radius_lbs(ncol,pver)
      logical, intent(in)   :: mask(ncol,pver)

!----------------------------------------------------------------------
!        ... local variables
!----------------------------------------------------------------------
      integer  :: k
      real(dp) :: lbs_vol_dens(ncol)       ! Vol Density (cm3 aer / cm3 air)

!----------------------------------------------------------------------
!             ... parameters
!----------------------------------------------------------------------
      real(dp), parameter :: lbs_part_dens = 10._dp, sigma_lbs = 1.6_dp

!----------------------------------------------------------------------
!             ... calculate the volume density (cm3 aerosol / cm3 air)
!                 calculate the mean radius for binary soln
!----------------------------------------------------------------------
      do k = sad_topp,pver
         where( mask(:,k) )
            lbs_vol_dens(:) = ((sad_sage(:,k)**1.5_dp)/3._dp)/sqrt( four_pi*lbs_part_dens ) &
                              *exp( 1.5_dp*(log( sigma_lbs ))**2 )
            radius_lbs(:,k) = (3._dp*lbs_vol_dens(:)/(four_pi*lbs_part_dens))**one_thrd &
                              *exp( -1.5_dp*(log( sigma_lbs ))**2 )
         endwhere
      end do

      end subroutine calc_radius_lbs

      subroutine ice_sad_calc( ncol, press, temp, m, h2o_avail, &
                               sad_ice, radius_ice, mask )

      implicit none

!----------------------------------------------------------------------
!        ... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(dp), intent(in)  :: press(ncol,pver)
      real(dp), intent(in)  :: temp(pcols,pver)
      real(dp), intent(in)  :: m(ncol,pver)
      real(dp), intent(in)  :: h2o_avail(ncol,pver)
      real(dp), intent(out) :: sad_ice(ncol,pver)
      real(dp), intent(out) :: radius_ice(ncol,pver)
      logical, intent(in)   :: mask(ncol,pver)

!----------------------------------------------------------------------
!        ... local variables
!----------------------------------------------------------------------
      real(dp), parameter :: &
                 avo_num = 6.02214e23_dp, &
                 aconst  = -2663.5_dp, &
                 bconst  = 12.537_dp, &
                 ice_mass_dens = 1._dp, &
                 ice_part_dens = 1.e-3_dp, &
                 mwh2o         = 18._dp, &
                 sigma_ice     = 1.6_dp, &
                 ice_dens_aer  = ice_mass_dens / (mwh2o/avo_num), &
                 ice_dens_aeri = 1._dp/ice_dens_aer

      integer  :: k
      real(dp) :: h2o_cond_ice(ncol)              ! Condensed phase H2O (from CAM)
      real(dp) :: voldens_ice(ncol)               ! Volume Density, um3 cm-3

      do k = sad_topp,pver
         where( mask(:,k) )
!----------------------------------------------------------------------
!     .... Convert condensed phase to molecules cm-3 units
!----------------------------------------------------------------------
           h2o_cond_ice(:) = h2o_avail(:,k) * m(:,k)
!----------------------------------------------------------------------
!     .... ICE volume density .....
!----------------------------------------------------------------------
           voldens_ice(:) = h2o_cond_ice(:)*ice_dens_aeri
!----------------------------------------------------------------------
!     .... Calculate the SAD from log normal distribution .....
!----------------------------------------------------------------------
           sad_ice(:,k) = (four_pi*ice_part_dens)**one_thrd &
                         *(3._dp*voldens_ice(:))**two_thrd &
                         *exp( -(log( sigma_ice ))**2 )
!----------------------------------------------------------------------
!    .... Calculate the radius from log normal distribution .....
!----------------------------------------------------------------------
           radius_ice(:,k) = (3._dp*h2o_cond_ice(:) &
                              /(ice_dens_aer*four_pi*ice_part_dens))**one_thrd &
                             *exp( -1.5_dp*(log( sigma_ice ))**2 )
         endwhere
      end do

      end subroutine ice_sad_calc

      subroutine sulfate_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, sad_sage, m, &
                                   hno3_gas, hno3_cond, sad_sulfate, radius_sulfate, mask, &
                                   lchnk, flag )

      implicit none

!----------------------------------------------------------------------
!        ... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      integer, intent(in)   :: lchnk, flag
      real(dp), intent(in)  :: temp(pcols,pver)
      real(dp), intent(in)  :: press(ncol,pver)
      real(dp), intent(in)  :: m(ncol,pver)
      real(dp), intent(in)  :: h2o_avail(ncol,pver)
      real(dp), intent(in)  :: hno3_avail(ncol,pver)
      real(dp), intent(in)  :: sad_sage(ncol,pver)
      real(dp), intent(out) :: hno3_gas(ncol,pver)             ! Gas-phase HNO3, mole fraction
      real(dp), intent(out) :: hno3_cond(ncol,pver)            ! Condensed phase HNO3, mole fraction
      real(dp), intent(out) :: sad_sulfate(ncol,pver)   
      real(dp), intent(out) :: radius_sulfate(ncol,pver)
      logical, intent(in)   :: mask(ncol,pver)

!----------------------------------------------------------------------
!        ... local variables
!----------------------------------------------------------------------
      real(dp), parameter :: avo_num = 6.02214e23_dp, &
                             mwh2so4 = 98.076_dp, &
                             sigma_sulfate     = 1.6_dp, &
                             sulfate_part_dens = 10._dp

      integer  :: i, k
      real(dp) :: h2so4m(ncol,pver)               ! mass per volume, micro grams m-3
      real(dp) :: h2so4_aer_dens(ncol,pver)       ! grams cm-3 solution
      real(dp) :: h2so4_cond(ncol,pver)           ! Condensed H2SO4 (moles cm-3 air)
      real(dp) :: sulfate_vol_dens(ncol,pver)     ! Volume Density, um3 cm-3
      real(dp) :: wtf(ncol,pver)                  ! weight fraction of H2SO4 in ternary soln
      real(dp) :: wts(ncol,pver)                  ! weight percent of ternary solution

!----------------------------------------------------------------------
!          ... derive H2SO4 (micro grams / m3) from SAGEII SAD
!----------------------------------------------------------------------
      call sad2h2so4( h2o_avail, press, sad_sage, temp, sulfate_vol_dens, & 
                      h2so4_aer_dens, h2so4m, mask, ncol )

!----------------------------------------------------------------------
!          ... limit h2so4m
!----------------------------------------------------------------------
      do k = sad_topp,pver
         do i = 1,ncol
            if( mask(i,k) ) then
!              if( h2so4m(i,k) <= 0._dp ) then   ! changed by H. Schmidt, MPI, Nov. 2005
               if( h2so4m(i,k) <= 1.e-4_dp ) then
                  h2so4m(i,k) = 1.e-4_dp
               end if
            end if
         end do
      end do

!----------------------------------------------------------------------
!     .... Calculate the ternary soln volume density
!----------------------------------------------------------------------

      call equil( temp, h2so4m, hno3_avail, h2o_avail, press, & 
                  hno3_cond, h2so4_cond, wts, mask, &
                  ncol, lchnk, flag )

      do k = sad_topp,pver
         where( mask(:,k) )
!----------------------------------------------------------------------
!     .... convert h2o, hno3 from moles cm-3 air to molecules cm-3 air
!----------------------------------------------------------------------
            hno3_cond(:,k) = min( hno3_cond(:,k),hno3_avail(:,k) )
            hno3_gas(:,k)  = hno3_avail(:,k) - hno3_cond(:,k)
!----------------------------------------------------------------------
!     .... Derive ternary volume density (cm3 aer / cm3 air)
!----------------------------------------------------------------------
            wtf(:,k) = .01_dp* wts(:,k)
            sulfate_vol_dens(:,k) = h2so4_cond(:,k)*mwh2so4/(wtf(:,k)*h2so4_aer_dens(:,k))
!----------------------------------------------------------------------
!     .... Calculate the SAD (assuming ternary solution)
!----------------------------------------------------------------------
            sad_sulfate(:,k) = (four_pi*sulfate_part_dens)**one_thrd &
                               *(3._dp*sulfate_vol_dens(:,k))**two_thrd &
                               *exp( -(log( sigma_sulfate ))**2 )
!----------------------------------------------------------------------
!     .... Calculate the radius (assuming ternary solution)
!----------------------------------------------------------------------
            radius_sulfate(:,k) = (3._dp*sulfate_vol_dens(:,k) &
                                   /(four_pi*sulfate_part_dens))**one_thrd &
                                  *exp( -1.5_dp*(log( sigma_sulfate ))**2 )
         endwhere
      end do

      end subroutine sulfate_sad_calc


      subroutine nat_sad_calc( ncol, press, temp, h2o_avail, hno3_avail, m, &
                               hno3_gas, hno3_cond_sm, hno3_cond_lg, sad_nat_sm, &
                               sad_nat_lg, radius_nat_sm, radius_nat_lg, mask )

      implicit none

!----------------------------------------------------------------------
!        ... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(dp), intent(in)  :: press(ncol,pver)
      real(dp), intent(in)  :: m(ncol,pver)
      real(dp), intent(in)  :: temp(pcols,pver)
      real(dp), intent(in)  :: h2o_avail(ncol,pver)
      real(dp), intent(in)  :: hno3_avail(ncol,pver)
      real(dp), intent(out) :: hno3_cond_sm(ncol,pver)             ! HNO3 in condensed phase (mole fraction)
      real(dp), intent(out) :: hno3_cond_lg(ncol,pver)             ! HNO3 in condensed phase (mole fraction)
      real(dp), intent(out) :: hno3_gas(ncol,pver)                 ! HNO3 in gas-phase (mole fraction)
      real(dp), intent(out) :: sad_nat_sm(ncol,pver)   
      real(dp), intent(out) :: sad_nat_lg(ncol,pver)   
      real(dp), intent(out) :: radius_nat_sm(ncol,pver)   
      real(dp), intent(out) :: radius_nat_lg(ncol,pver)
      logical, intent(in)   :: mask(ncol,pver)                     ! grid mask

!----------------------------------------------------------------------
!        ... local variables
!----------------------------------------------------------------------
      integer  :: k, i
      real(dp) :: mi
      real(dp) :: wrk
      real(dp) :: nat_dens_condphase                  ! Condensed phase NAT, molec cm-3
      real(dp) :: nat_dens_cp_sm                      ! Condensed phase NAT, sm mode, molec cm-3
      real(dp) :: nat_dens_cp_lg                      ! Condensed phase NAT, lg mode, molec cm-3
      real(dp) :: voldens_lg                          ! Volume Density, sm mode, um3 cm-3
      real(dp) :: voldens_nat_lg                      ! Volume Density, sm mode, um3 cm-3
      real(dp) :: voldens_nat_sm                      ! Volume Density, sm mode, um3 cm-3
      real(dp) :: hno3_cond_total(ncol,pver)          ! Total Condensed phase HNO3 (both modes)

!----------------------------------------------------------------------
!     ... parameters
!----------------------------------------------------------------------
      real(dp), parameter :: avo_num          = 6.02214e23_dp, &
                             nat_mass_dens    = 1.6_dp, &
                             nat_part_dens_lg = 2.3e-4_dp, &
                             mwnat            = 117._dp, &
                             sigma_nat        = 1.3_dp, &
                             radius_sm        = 5e-5_dp, &
                             radius_lg        = 6.5e-4_dp, &
                             nat_dens_aer     = nat_mass_dens / (mwnat/avo_num), &
                             nat_dens_aeri    = 1._dp/nat_dens_aer

!----------------------------------------------------------------------
!     ... Derive HNO3 paritioning (call A. Tabazedeh routine for NAT)
!----------------------------------------------------------------------
      call nat_cond( ncol, press, temp, h2o_avail, hno3_avail, &
                     hno3_gas, hno3_cond_total, mask )

      nat_dens_cp_lg = (four_pi*nat_dens_aer*nat_part_dens_lg)*(radius_lg**3) &
                       /(3._dp*(exp( -1.5_dp*(log(sigma_nat))**2._dp ))**3 )

      voldens_nat_lg = nat_dens_cp_lg * nat_dens_aeri

      do k = sad_topp,pver
         do i = 1,ncol
masked :   if( mask(i,k) ) then
!----------------------------------------------------------------------
!     ... Calculated Condensed Phase NAT (i.e. HNO3) in molecules cm-3 of air units.
!         This is total (sm and lg modes)
!----------------------------------------------------------------------
              nat_dens_condphase = hno3_cond_total(i,k) * m(i,k)

              if( nat_dens_condphase > nat_dens_cp_lg ) then
                 nat_dens_cp_sm      = nat_dens_condphase - nat_dens_cp_lg
                 mi                  = 1._dp/m(i,k)
                 hno3_cond_sm(i,k)   = nat_dens_cp_sm * mi
                 hno3_cond_lg(i,k)   = nat_dens_cp_lg * mi
              else
                 nat_dens_cp_sm      = 0._dp
                 hno3_cond_lg(i,k)   = hno3_cond_total(i,k)
                 hno3_cond_sm(i,k)   = 0._dp
              end if
!----------------------------------------------------------------------
!     ... Calculate the Volume Density (for SAD, reff)
!----------------------------------------------------------------------
              if( hno3_cond_lg(i,k) > 0._dp ) then
                 voldens_lg = voldens_nat_lg
              else
                 voldens_lg = 0._dp
              end if
              voldens_nat_sm = nat_dens_cp_sm * nat_dens_aeri
!----------------------------------------------------------------------
!     ... Calculate the SAD from log normal distribution
!         Assuming sigma, nat_part_dens (# particles per cm3 of air),
!         and the mean radius. 
!----------------------------------------------------------------------
              wrk             = exp( -2.5_dp*(log( sigma_nat )**2 ))
              sad_nat_lg(i,k) = 3._dp/radius_lg * voldens_lg * wrk
              sad_nat_sm(i,k) = 0.0_dp
!----------------------------------------------------------------------
!     ... Calculate the radius of NAT from log normal distribution
!         Assuming sigma and nat_part_dens (# particles per cm3 of air)
!----------------------------------------------------------------------
              radius_nat_lg(i,k) = radius_lg
              radius_nat_sm(i,k) = 0.0_dp

           end if masked
         end do
      end do

      end subroutine nat_sad_calc

      subroutine nat_cond( ncol, press, temp, h2o_avail, hno3_avail, &
                           hno3_gas, hno3_cond, mask )

      implicit none

!----------------------------------------------------------------------
!         ... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol
      real(dp), intent(in)  :: press(ncol,pver)
      real(dp), intent(in)  :: temp(pcols,pver)
      real(dp), intent(in)  :: h2o_avail(ncol,pver)
      real(dp), intent(in)  :: hno3_avail(ncol,pver)
      real(dp), intent(out) :: hno3_gas(ncol,pver)
      real(dp), intent(out) :: hno3_cond(ncol,pver)
      logical, intent(in)   :: mask(ncol,pver)

!----------------------------------------------------------------------
!         ... local variables
!----------------------------------------------------------------------
      real(dp), parameter ::  aa = -2.7836_dp,  &
                              bb = -0.00088_dp, &
                              cc = 38.9855_dp,  &
                              dd = -11397.0_dp, &
                              ee = 0.009179_dp

      integer  :: i, k
      real(dp) :: bt                                 ! temporary variable
      real(dp) :: mt                                 ! temporary variable
      real(dp) :: t                                  ! temporary variable
      real(dp) :: logPhno3                           ! temporary variable
      real(dp) :: phno3                              ! hno3 partial pressure
      real(dp) :: ph2o                               ! h2o  partial pressure
      real(dp) :: phno3_eq                           ! partial pressure above NAT
      real(dp) :: wrk      
      
      do k = sad_topp,pver
         do i = 1,ncol
!----------------------------------------------------------------------
!     .... Derive HNO3 and H2O partial pressure (torr)
!          where: 0.7501 = 760/1013.
!----------------------------------------------------------------------        
            if( mask(i,k) ) then
               wrk   = press(i,k) * .7501_dp
               phno3 = hno3_avail(i,k) * wrk
               ph2o  = h2o_avail(i,k)  * wrk
!----------------------------------------------------------------------
!     Calculating the temperature coefficients for the variation of HNO3
!     and H2O vapor pressure (torr) over a trihydrate solution of HNO3/H2O
!     The coefficients are taken from Hanson and Mauersberger: 
!     GRL, Vol.15, 8, p855-858, 1988.
!----------------------------------------------------------------------
               t   = temp(i,k)
               bt  = cc + dd/t + ee*t
               mt  = aa + bb*t
  
               logphno3 = mt*log10( ph2o ) + bt
               phno3_eq = 10._dp**logphno3

               if( phno3 > phno3_eq ) then
                  wrk            = 1._dp / wrk
                  hno3_cond(i,k) = (phno3 - phno3_eq) * wrk
                  hno3_gas(i,k)  = phno3_eq * wrk
               else
                  hno3_cond(i,k) = 0._dp
                  hno3_gas(i,k)  = hno3_avail(i,k)
               end if
            end if
         end do
      end do

      end subroutine nat_cond

      subroutine sad2h2so4( h2o, press, sad_sage, temp, lbs_vol_dens, &
                            h2so4_aer_dens, h2so4m, mask, ncol )

      implicit none

!----------------------------------------------------------------------
!        ... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: ncol                                    ! columns in chunk
      real(dp), intent(in)  :: h2o(ncol,pver)                          ! h2o mole fraction
      real(dp), intent(in)  :: sad_sage(ncol,pver)                     ! sad from SAGEII cm2 aer, cm-3 air
      real(dp), intent(in)  :: press(ncol,pver)                        ! pressure (hPa)
      real(dp), intent(in)  :: temp(pcols,pver)                        ! temperature (K)
      real(dp), intent(out) :: h2so4m(ncol,pver)                       ! microgram/m**3 of air,
      real(dp), intent(out) :: h2so4_aer_dens(ncol,pver)               ! units: grams / cm3-aerosol
      real(dp), intent(out) :: lbs_vol_dens(ncol,pver)                 ! cm3 aer / cm3 air
      logical, intent(in)   :: mask(ncol,pver)                         ! activation mask

!----------------------------------------------------------------------
!        ... local variables
!----------------------------------------------------------------------
      real(dp), parameter :: lbs_part_dens = 10._dp
      real(dp), parameter :: sigma_lbs     = 1.6_dp
!     real(dp), parameter :: t_floor       = 180._dp
      real(dp), parameter :: t_floor       = 183._dp  ! changed by H. Schmidt, MPI, Nov. 2005.: temperatures below
                                                      ! 182K have led to molh2so4=0 and molhno3=0 and a subsequent
                                                      ! division by zero error

      integer  :: i, k, l
      real(dp) :: wts0   
      real(dp) :: p                         ! pressure, torr
      real(dp) :: tr                        ! inverse temperature
      real(dp) :: c(6)

!----------------------------------------------------------------------
!             ... Calculate the volume density (cm3 aerosol / cm3 air)
!----------------------------------------------------------------------
      do k = sad_topp,pver
         do i = 1,ncol
            if( mask(i,k) ) then
               lbs_vol_dens(i,k) = ((sad_sage(i,k)**1.5_dp)/3._dp)/sqrt( four_pi*lbs_part_dens ) &
                                   *exp( 1.5_dp*(log( sigma_lbs ))**2 )
!----------------------------------------------------------------------
!             ... calculate Molality from Tabazadeh EQUIL routine (binary soln)
!----------------------------------------------------------------------
!             ... DEK, added a minimum to temperature
!----------------------------------------------------------------------
               p   = h2o(i,k) * press(i,k) * .7501_dp
               tr  = 1._dp / max( t_floor,temp(i,k) )
             
               do l = 1,6
                  c(l) = exp( a(1,l) + tr*(a(2,l) + tr*(a(3,l) + tr*(a(4,l) + tr*a(5,l)))) )
               end do
!----------------------------------------------------------------------
!             ... H2SO4/H2O pure weight percent and molality (moles gram-1)
!----------------------------------------------------------------------
               wts0 = max( 0._dp,c(1) + p*(-c(2) + p*(c(3) + p*(-c(4) + p*(c(5) - p*c(6))))) )
!----------------------------------------------------------------------
!             ... derive weight fraction for density routine
!----------------------------------------------------------------------
               wts0 = .01_dp *wts0
!----------------------------------------------------------------------
!             ... calculate the binary solution density, grams / cm3-aerosol
!----------------------------------------------------------------------
               h2so4_aer_dens(i,k) = max( 0._dp,density( temp(i,k), wts0 ) )
!----------------------------------------------------------------------
!             ... calculate the H2SO4 micrograms m-3 abundance for binary soln
!----------------------------------------------------------------------
               h2so4m(i,k) = lbs_vol_dens(i,k)*h2so4_aer_dens(i,k)*wts0*1.e12_dp
            end if
         end do
      end do
   
      end subroutine sad2h2so4

!======================================================================
!
! ROUTINE
!   EQUIL
!
!   Date...
!     7 October 1999
!
!   Programmed by...
!     A. Tabazadeh
!
! DESCRIPTION
!        Ternary solution routine
!
! ARGUMENTS
!
!....  INPUT:
!
!        H2SO4m    = microgram/m**3 of air
!        HNO3r     = mole fraction
!        H2Or      = mole fraction
!        PTOTAL    = hPa
!
!....  Output
!
!        Cwater    = Total moles of liguid water / cm**3 of air
!        hno3_cond = HNO3 Condensed phase (mole fraction)
!        CH2SO4    = Total H2SO4 moles / cm**3 of air
!        WTS       = Weight percent of H2SO4 in the ternary aerosol
!
!======================================================================

      subroutine equil( temper, h2so4m, hno3_avail, h2o_avail, press, &
                        hno3_cond, ch2so4, wts, mask, ncol, lchnk, flag )
!----------------------------------------------------------------------
!     *******************************************************************
!     ******************Written by Azadeh Tabazadeh (1993)***************
!     ***** (modified from EQUISOLV -- M. Z. Jacobson -- see below) *****
!     **** NASA Ames Research Center , Tel. (415) 604 - 1096 ******* ****
!     *******************************************************************
!
!     *******************************************************************
!     * This program solves the equilibrium composition for the ternary *
!     * system of H2SO4/HNO3/H2O under typical stratospheric conditions.*
!     * The formulation of this work is described by Tabazadeh, A.,     *
!     * Turco, R. P., and Jacobson, M. Z. (1994), "A model for studying *
!     * the composition and chemical effects of stratospheric aerosols,"*
!     * J. Geophys. Res., 99, 12,897, 1994.        *
!     *******************************************************************
!
!     *******************************************************************
!     * The solution mechanism for the equilibrium equations is des-    *
!     * cribed by Jacobson, M. Z., Turco, R. P., and Tabazadeh, A.      *
!     * (1994), "Simulating Equilibrium within aerosols and non-equil-  *
!     * ibrium between gases and aerosols," J. Geophys. Res., in review.*
!     * The mechanism is also codified in the fortran program, EQUISOLV,*
!     * by M.Z. Jacobson (1991-3). EQUISOLV solves any number of        *
!     * gas / liquid / solid / ionic equilibrium equations simultan-    *
!     * eously and includes treatment of the water equations and act-   *
!     * ivity coefficients. The activity coeffients currently in        *
!     * EQUISOLV are valid for tropospheric temperatures. The acitiv-   *
!     * ities listed here are valid for stratospheric temperatures only.*
!     *******************************************************************
!
!        *******************************************************************
!
!        DEFINING PARAMETERS
!
!       *NOTE*          Solver parameters including, F, Z, QN, QD, and deltaX
!                 are described in Jacobson et al.
!
!        PTOTAL        = Total atmospheric pressure in mb
!        H2SO4m        = Total mass of H2SO4 (microgram/m**3 of air)
!        HNO3r        = HNO3 mixing ratio
!        H2Or        = H2O mixing ratio
!        P            = Partial pressure of water in units of torr
!        pures        = molality for a pure H2SO4/H2O system
!        puren        = molality for a pure HNO3/H2O sytem
!        WTS0        = weight percent of H2SO4 in a pure H2SO4/H2O system
!        WTN0        = weight percent of HNO3 in a pure HNO3/H2O system
!        WTS            = weight percent of H2SO4 in the ternary aerosol
!        WTN         = weight percent of HNO3 in the ternary aerosol
!        PHNO3        = HNO3 vapor pressure over the ternary system in atm
!        HNO3        = HNO3 vapor concentration over the ternary system (#/cm3)
!        CH2SO4        = Total H2SO4 moles / cm**3 of air
!        CHNO3        = Total HNO3 moles / cm**3 of air
!        CHplus        = Total H+ moles / cm**3 0f air
!        CPHNO3        = Total moles of HNO3 gas / cm**3 of air
!        CNO3        = Total moles of NO3- / cm**3 0f air
!        Cwater        = Total moles of liguid water / cm**3 of air
!        KS          = Solubility constant for dissolution of HNO3 in
!                      water ( HNO3(gas) === H+(aq) + NO3- (aq) )
!        nm          = HNO3 molality at the STREN of the ternary solution
!        sm          = H2SO4 molality at the STREN of the ternary solution
!        molHNO3        = Equilibrium molality of HNO3 in the ternary solution
!        molH2SO4= Equilibrium molality of H2SO4 in the ternary solution
!     STREN   = ionic strenght for the ternary solutin, which in
!                      this case is = 3 * molH2SO4 + molHNO3
!        acts        = Pure mean binary activity coefficient for the H2SO4/
!                      H2O system evaluated at the STREN of the ternary system
!        actn        = Pure mean binary activity coefficient for the HNO3/
!                      H2O system evaluated at the STREN of the ternary system
!     ymix    = Mixed binary activity coefficient for the HNO3/H2O in
!                      the ternary solution
!----------------------------------------------------------------------

      implicit none

!----------------------------------------------------------------------
!        ... dummy arguments
!----------------------------------------------------------------------
      integer, intent(in)   :: lchnk
      integer, intent(in)   :: flag
      integer, intent(in)   :: ncol                         ! columns in chunk
      real(dp), intent(in)  :: h2so4m(ncol,pver)    
      real(dp), intent(in)  :: hno3_avail(ncol,pver)    
      real(dp), intent(in)  :: h2o_avail(ncol,pver)    
      real(dp), intent(in)  :: press(ncol,pver)
      real(dp), intent(in)  :: temper(pcols,pver)
      real(dp), intent(out) :: hno3_cond(ncol,pver)    
      real(dp), intent(out) :: ch2so4(ncol,pver)    
      real(dp), intent(out) :: wts(ncol,pver)
      logical, intent(in)   :: mask(ncol,pver)              ! activation mask

!----------------------------------------------------------------------
!        ... local variables
!----------------------------------------------------------------------
      integer, parameter  :: itermax = 50
      real(dp), parameter :: t0  = 298.15_dp
      real(dp), parameter :: ks0 = 2.45e6_dp

      integer  :: i, iter, k, l, nstep
      real(dp) :: h2o_cond(ncol,pver)
      real(dp) :: p
      real(dp) :: tr
      real(dp) :: wts0
      real(dp) :: wtn0
      real(dp) :: pures
      real(dp) :: puren
      real(dp) :: chno3
      real(dp) :: chplus
      real(dp) :: cphno3
      real(dp) :: cno3
      real(dp) :: qd
      real(dp) :: qn
      real(dp) :: z, num, den
      real(dp) :: deltax
      real(dp) :: chplusnew
      real(dp) :: cno3new
      real(dp) :: stren
      real(dp) :: sm
      real(dp) :: actn
      real(dp) :: acts
      real(dp) :: nm
      real(dp) :: ks
      real(dp) :: lnks
      real(dp) :: lnks0
      real(dp) :: mixyln
      real(dp) :: molhno3
      real(dp) :: molh2so4
      real(dp) :: wrk_h2so4
      real(dp) :: cphno3new
      real(dp) :: phno3
      real(dp) :: t, t1, t2, f, f1, f2, ymix, hplus, wtotal, wtn, ratio 
      real(dp) :: c(12)
      real(dp) :: d(13:22)
      logical  :: converged

      lnks0 = log( ks0 )
Level_loop : &
      do k = sad_topp,pver
Column_loop : &
         do i = 1,ncol
            if( mask(i,k) ) then
               p = h2o_avail(i,k) * press(i,k) * .7501_dp
!----------------------------------------------------------------------
!        Calculating the molality for pure binary systems of H2SO4/H2O
!        and HNO3/H2O at a given temperature and water vapor pressure
!        profile (relative humiditiy). Water activities were used to
!        calculate the molalities as described in Tabazadeh et al. (1994).
!----------------------------------------------------------------------
               t  = max( 180._dp,temper(i,k) )
               tr = 1._dp/t
               do l = 1,12
                  c(l) = exp( a(1,l) + tr*(a(2,l) + tr*(a(3,l) + tr*(a(4,l) + tr*a(5,l)))) )
               end do
!----------------------------------------------------------------------
!        ... H2SO4/H2O pure weight percent and molality
!----------------------------------------------------------------------
               wts0  = max( 0._dp,c(1) + p*(-c(2) + p*(c(3) + p*(-c(4) + p*(c(5) - p*c(6))))) )
               pures = (wts0 * 1000._dp)/(100._dp - wts0)
               pures = pures / 98._dp
!----------------------------------------------------------------------
!        ... HNO3/H2O pure weight percent and molality
!----------------------------------------------------------------------
               puren = c(7) + p*(-c(8) + p*(c(9) + p*(-c(10) + p*(c(11) - p*c(12)))))
!              wtn0 = (puren * 6300._dp) /(puren * 63._dp + 1000._dp)
!----------------------------------------------------------------------
!        The solving scheme is described both in Jacobson et al. and Tabazadeh
!        et al.. Assumptions:
!        (1) H2SO4 is present only in the aqueous-phase
!        (2) H2SO4 and HNO3 in solution are fully dissocated into H+
!            SO42- and NO3-
!        (3) PHNO3 + NO3- = constant
!----------------------------------------------------------------------
               ch2so4(i,k) = (h2so4m(i,k)*1.e-12_dp) / 98._dp
               if( pures > 0._dp ) then
                  wrk_h2so4 = (1000._dp*ch2so4(i,k))/(pures*18._dp)
               else
                  wrk_h2so4 = 0._dp
               end if
               chno3 = 1.2029e-5_dp * press(i,k) * tr * hno3_avail(i,k)
!----------------------------------------------------------------------
!        Setting up initial guesses for the equations above.  Note that
!        for the initial choices the mass and the charge must be conserved.
!----------------------------------------------------------------------
!              chplus        = 2._dp * ch2so4(i,k) + .5_dp * chno3
               cphno3        = .5_dp * chno3
               cno3        = .5_dp * chno3
               qd        = cphno3
               qn        = cno3
               z        = .5_dp * (qd + qn)
               deltax        = qd - z
!              chplusnew = chplus + deltax
               cno3new   = cno3 + deltax
               cphno3new = cphno3 - deltax
               do l = 13,22
                  d(l) = b(1,l) + t*(b(2,l) + t*(b(3,l) + t*(b(4,l) + t*b(5,l))))
               end do
!----------------------------------------------------------------------
!        Note that KS depends only on the temperature
!----------------------------------------------------------------------
               t1        = (t - t0)/(t*t0)
               t2        = t0/t - 1._dp - log( t0/t )
               lnks     = lnks0 - 8792.3984_dp * t1  - 16.8439_dp * t2
               ks        = exp( lnks )

               converged = .false.
Iter_loop : &
               do iter = 1,itermax
!----------------------------------------------------------------------
!        Cwater is the water equation as described in Tabazadeh et
!        al. and Jacobson et al.
!----------------------------------------------------------------------
                  if( puren > 0._dp ) then
                     t1 = (1000._dp*cno3new)/(puren*18._dp)
                  else
                     t1 = 0._dp
                  end if
                  h2o_cond(i,k) = t1 + wrk_h2so4
                  if( h2o_cond(i,k) > 0._dp ) then
                     molhno3  = (1000._dp * cno3new) / (18._dp * h2o_cond(i,k))
                     molh2so4 = (1000._dp * ch2so4(i,k)) / (18._dp * h2o_cond(i,k))
                  else
!                    molhno3  = 0._dp
!                    molh2so4 = 0._dp
!hotfix: prevent division by zero! #569
!problem: mask derived from h2o_cond outside of this routine, new calculation in
!         equil can lead to h2o_cond = 0
                     molhno3  = 1.e-25_dp
                     molh2so4 = 1.e-25_dp
!!++mgs20140331: DEBUG ####
!                    call message("EQUIL", "verbotener IF Zweig!", level=em_warn, all_print=.true.)
!                    write(message_text, '(a,5i5)') 'equil>> i,k,lchnk,iter = ',i,k,lchnk,iter
!                    CALL message('equil', message_text, level=em_warn, all_print=.true.)
!                    write(message_text,'(a,5e15.5)') ' h2o_cond, t1, wrk_h2so4,puren,cno3new = ', &
!                          h2o_cond(i,k),t1,wrk_h2so4,puren,cno3new
!                    CALL message('equil', message_text, level=em_warn, all_print=.true.)
!                    write(message_text,'(a,5e15.5)') ' pures,ch2so4,h2so4m,temper,h2o_avail = ', &
!                          pures,ch2so4(i,k),h2so4m(i,k),temper(i,k),h2o_avail(i,k)
!                    CALL message('equil', message_text, level=em_warn, all_print=.true.)
!                    write(message_text,'(a,4e15.5)') ' h2o_avail,hno3_cond,press,wts = ', &
!                           h2o_avail(i,k),hno3_cond(i,k),press(i,k),wts(i,k)
!                    CALL message('equil', message_text, level=em_warn, all_print=.true.)
                  end if
                  stren        = molhno3 + 3._dp * molh2so4
                  stren=max(stren,1.e-10_dp)
                  phno3        = 1000._dp * cphno3new * .0820578_dp * t
!----------------------------------------------------------------------
!        (1) Calculating the activity of H2SO4 at a given STREN
!----------------------------------------------------------------------
                  sm        = stren/3._dp
                  acts = d(13) + sm*(d(14) + sm*(d(15) + sm*d(16)))
!----------------------------------------------------------------------
!        (2) Calculating the activity for HNO3 at a given STREN
!----------------------------------------------------------------------
                  nm        = stren
                  actn         = d(17) + nm*(d(18) + nm*(d(19) + nm*(d(20) + nm*(d(21) + nm*d(22)))))
!----------------------------------------------------------------------
!        (3) Calculating the mixed activity coefficient for HNO3 at STREN
!            as described by Tabazadeh et al.
!----------------------------------------------------------------------
                  f1        = 2._dp * (molh2so4 + molhno3) * actn
                  f2        = 2.25_dp * molh2so4 * acts
                  mixyln = (f1 + f2) / (2._dp * stren)
                  ymix         = exp( mixyln )
                  hplus        = 2._dp * molh2so4 + molhno3
!----------------------------------------------------------------------
!       Calculating the ratio F and resetting the deltaX (see Jacobson et al.)
!----------------------------------------------------------------------
                  num = ymix**2 * hplus * molhno3
                  den = phno3 * ks
                  if( den == 0._dp ) then
!                     write(*,*) 'equil: i,k,iter = ',i,k,iter
!                     write(*,*) 'equil: num,den,ymix,hplus,molhno3,phno3,ks,temp'
!                     write(*,*) num,den,ymix,hplus,molhno3,phno3,ks,temper(i,k)
                     if( num == 0._dp ) then
                        converged = .true.
                        exit
                     else
                        CALL finish('equil', 'Not converged (num/=0).')
                     end if
                  end if
                  f = num / den
                  z = .5_dp * z
                  if( abs( f - 1._dp ) <= .00005_dp ) then
                     converged = .true.
                     exit
                  else
                     if( f > 1._dp ) then
                        deltax = -z
                     else
                        if( f < 1.e-2_dp .and. iter > 9 ) then
                           cno3new   = chno3
                           cphno3new = 0._dp
                           if( puren > 0._dp ) then
                              t1 = (1000._dp*cno3new)/(puren*18._dp)
                           else
                              t1 = 0._dp
                           end if
                           h2o_cond(i,k) = t1 + wrk_h2so4
                           if( h2o_cond(i,k) > 0._dp ) then
                              molhno3  = (1000._dp * cno3new) / (18._dp * h2o_cond(i,k))
                              molh2so4 = (1000._dp * ch2so4(i,k)) / (18._dp * h2o_cond(i,k))
                           else
!                             molhno3  = 0._dp
!                             molh2so4 = 0._dp
!hotfix: prevent division by zero! #569
!problem: mask derived from h2o_cond outside of this routine, new calculation in
!         equil can lead to h2o_cond = 0
                              molhno3  = 1.e-25_dp
                              molh2so4 = 1.e-25_dp
                           end if
                           converged = .true.
!                           write(*,*) 'Equil : f = ',f
                           exit
                        end if
                        deltax = z
                     end if
                     cno3new   = cno3new   + deltax
                     cphno3new = cphno3new - deltax
                     cycle
                  end if
               end do Iter_loop
               wtotal   = molhno3 * 63._dp + molh2so4 * 98._dp + 1000._dp
               wts(i,k) = (molh2so4 * 9800._dp) / wtotal
               if( cno3new /= 0._dp .or. cphno3new /= 0._dp ) then
                  ratio        = max( 0._dp,min( 1._dp,cno3new/(cphno3new + cno3new) ) )
                  hno3_cond(i,k) = ratio*hno3_avail(i,k)
               else
                  hno3_cond(i,k) = 0._dp
               end if
               if( .not. converged ) then
                  WRITE(message_text,'(a,4(i0,1x),e25.15)')                                      &
                                     'Failed to converge @ flag,lchnk,i,k,f = ',flag,lchnk,i,k,f
                  CALL message('sad:equil', message_text, level=em_error, all_print=.true.)
                  WRITE(message_text,'(a,4(e25.15))')                                      &
                                     'h2o_avail,hno3_avail,p,t = ',h2o_avail(i,k),hno3_avail(i,k),press(i,k),temper(i,k)
                  CALL message('sad:equil', message_text, level=em_error, all_print=.true.)
                  WRITE(message_text,'(a,4(e25.15))')                                      &
                                     'molhno3,molh2so4,h2o_cond,hno3_cond = ',molhno3,molh2so4,h2o_cond(i,k),hno3_cond(i,k)
                  CALL message('equil', message_text, level=em_error, all_print=.true.)
                  CALL finish('sad:equil', 'Failed to converge!')
               end if
            end if
         end do Column_loop
      end do Level_loop

      end subroutine equil

!======================================================================
!
!
! ROUTINE
!   DENSITY
!
!   Date...
!     7 October 1999
!
!   Programmed by...
!     A. Tabazadeh
!
! DESCRIPTION
!     Calculates the density (g cm-3) of a binary sulfate solution.
!
! ARGUMENTS
!   INPUT
!      T           Temperature
!      w           Weight fraction
!
!   OUTPUT
!        den       Density of the Binary Solution (g cm-3)
!
!======================================================================
       
      function density( temp, w )

      implicit none

!----------------------------------------------------------------------
!        ... Dummy arguments
!----------------------------------------------------------------------
      real(dp), intent(in) :: temp, w

!----------------------------------------------------------------------
!        ... Function declarations
!----------------------------------------------------------------------
      real(dp) :: density

!----------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------
      real(dp), parameter :: a9 = -268.2616e4_dp, a10 = 576.4288e3_dp

      real(dp) :: a0, a1, a2, a3, a4, a5, a6, a7 ,a8
      real(dp) :: c1, c2, c3, c4

!----------------------------------------------------------------------
!        ... Temperature variables
!----------------------------------------------------------------------
      c1 = temp - 273.15_dp
      c2 = c1**2
      c3 = c1*c2
      c4 = c1*c3
!----------------------------------------------------------------------
!        Polynomial Coefficients
!----------------------------------------------------------------------
      a0 = 999.8426_dp + 334.5402e-4_dp*c1 - 569.1304e-5_dp*c2
      a1 = 547.2659_dp - 530.0445e-2_dp*c1 + 118.7671e-4_dp*c2 + 599.0008e-6_dp*c3
      a2 = 526.295e1_dp + 372.0445e-1_dp*c1 + 120.1909e-3_dp*c2 - 414.8594e-5_dp*c3 + 119.7973e-7_dp*c4
      a3 = -621.3958e2_dp - 287.7670_dp*c1 - 406.4638e-3_dp*c2 + 111.9488e-4_dp*c3 + 360.7768e-7_dp*c4
      a4 = 409.0293e3_dp + 127.0854e1_dp*c1 + 326.9710e-3_dp*c2 - 137.7435e-4_dp*c3 - 263.3585e-7_dp*c4
      a5 = -159.6989e4_dp - 306.2836e1_dp*c1 + 136.6499e-3_dp*c2 + 637.3031e-5_dp*c3
      a6 = 385.7411e4_dp + 408.3717e1_dp*c1 - 192.7785e-3_dp*c2
      a7 = -580.8064e4_dp - 284.4401e1_dp*c1
      a8 = 530.1976e4_dp + 809.1053_dp*c1
!----------------------------------------------------------------------
!        ... Summation
!----------------------------------------------------------------------
      density = .001_dp*(a0 + w*(a1 + w*(a2 + w*(a3 + w*(a4 + w*(a5 + w*(a6 + w*(a7 + w*(a8 + w*(a9 + w*a10))))))))))

      end function density

end module mo_moz_strato_sad
