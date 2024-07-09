!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_eslookup
!!
!! \brief
!!  This module collects various utilities to calculate saturation vapor pressure
!!  of water
!!
!! \author Stacy Walters (NCAR)
!! \author Douglas E. Kinnison (NCAR)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - SW: original version (pre 2004)
!!
!! \limitations
!!  none
!!
!! \details
!!  - ESINTI sets up a lookup table for saturation vapor pressure
!!  - AQSAT  utility procedure to look up and return saturation vapor
!!           pressure from precomputed table, calculate and return
!!           saturation specific humidity (g/g),for input arrays of
!!           temperature and pressure (dimensioned ii,kk)
!!           This routine is useful for evaluating only a selected
!!           region in the vertical.
!!  - AQSATD as above, but also computes gamma (l/cp)*(d(qsat)/dT).
!!  - GESTBL builds saturation pressure table
!!  - VQSATD the same function as qsatd, but operates on vectors of 
!!           temperature and pressure
!!  - GFFGCH computes saturation vapor pressure over water and/or over 
!!           ice using Goff & Gratch (1946) relationships
!!
!! \bibliographic_references
!!  - Goff, J. A., and S. Gratch (1946) Low-pressure properties of water from .160 to 212 °F, 
!!  in Transactions of the American Society of Heating and Ventilating Engineers, pp 95.122, 
!!  presented at the 52nd annual meeting of the American Society of Heating and Ventilating 
!!  Engineers, New York, 1946.
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

module MO_MOZ_ESLOOKUP

      use mo_kind,      only: dp
      use mo_exception, only: finish
      
      implicit none

      save

      integer, parameter :: plenest = 250  ! length of saturation vapor pressure table
!-----------------------------------------------------------------------
!         ... Table of saturation vapor pressure values es from tmin degrees
!           to tmax+1 degrees k in one degree increments.  ttrice defines the
!           transition region where es is a combination of ice & water values
!-----------------------------------------------------------------------

!++vector
      real(dp) :: &
           ttrice, &             ! transition range from es over water to es over ice
           epsqs, &              ! Ratio of h2o to dry air molecular weights 
           rgasv, &              ! Gas constant for water vapor
           hlatf, &              ! Latent heat of vaporization
           hlatv, &              ! Latent heat of fusion
           cp                    ! specific heat of dry air
      real(dp) :: &
              pcf(6)             ! polynomial coeffs -> es transition water to ice
      logical :: &
           icephs                ! false => saturation vapor pressure over water only

!     estblf1 is external now, i.e. not part of module mo_eslookup (for inlining) 
      interface
         function ESTBLF( td )
            use mo_kind,    only: dp
            real(dp), intent(in) :: td(:)
            real(dp), dimension( SIZE(td) ) :: ESTBLF
         end function ESTBLF
         real function ESTBLF1( td )
            use mo_kind,    only: dp
            real(dp), intent(in) :: td
         end function ESTBLF1
      end interface
!--vector

      private :: GESTBL, GFFGCH

      CONTAINS

      subroutine ESINTI( epslon, latvap, latice, rh2o, cpair, xip )
!-----------------------------------------------------------------------
!         ... Initialize es lookup tables 
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!         ... Dummy arguments
!-----------------------------------------------------------------------
      real(dp), intent(in) :: &
           epslon, &       ! Ratio of h2o to dry air molecular weights 
           latvap, &       ! Latent heat of vaporization
           latice, &       ! Latent heat of fusion
           rh2o, &         ! Gas constant for water vapor
           cpair           ! Specific heat of dry air
      logical, intent(in) :: &
           xip             ! Ice phase (true or false)

!-----------------------------------------------------------------------
!         ... Local variables
!-----------------------------------------------------------------------
      real(dp) :: &
           tmn, &          ! Minimum temperature entry in table
           tmx, &          ! Maximum temperature entry in table
           trice           ! Trans range from es over h2o to es over ice
      logical :: ip        ! Ice phase (true or false)

!-----------------------------------------------------------------------
!         ... Specify control parameters first
!-----------------------------------------------------------------------
      tmn   = 173.16_dp
      tmx   = 375.16_dp
      trice =  20.00_dp
      ip    = xip

!-----------------------------------------------------------------------
!         ... Call gestbl to build saturation vapor pressure table.
!-----------------------------------------------------------------------
      call GESTBL( tmn, tmx, trice, ip, epslon, latvap, latice, rh2o, cpair )

      end subroutine ESINTI

      subroutine AQSAT( t, p, es, qs, ii, ilen, kk, kstart, kend )
!-----------------------------------------------------------------------
!         ... Utility procedure to look up and return saturation vapor
!           pressure from precomputed table, calculate and return
!           saturation specific humidity (g/g),for input arrays of
!           temperature and pressure (dimensioned ii,kk)
!           This routine is useful for evaluating only a selected
!           region in the vertical.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!         ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
              ii, &          ! I dimension of arrays t, p, es, qs
              kk, &          ! K dimension of arrays t, p, es, qs 
              ilen, &        ! length of vectors in I direction which are assumed to start at 1
              kstart, &      ! starting location in K direction
              kend           ! ending location in K direction
      real(dp), dimension(ii,kk), intent(in) :: &
           t, &              ! Temperature
           p                 ! Pressure
      real(dp), dimension(ii,kk), intent(out) :: &
           es, &             ! Saturation vapor pressure
           qs                ! Saturation specific humidity

!-----------------------------------------------------------------------
!         ... Local variables
!-----------------------------------------------------------------------
      real(dp)    :: omeps          ! 1 - 0.622
      integer :: i, &           ! i index
                 k              ! k index

      omeps = 1._dp - epsqs
      do k = kstart,kend
!++vector
         do i = 1,ilen
            es(i,k) = ESTBLF1( t(i,k) )
!--vector
!-----------------------------------------------------------------------
!         ... Saturation specific humidity
!-----------------------------------------------------------------------
            qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
!-----------------------------------------------------------------------
!         ... The following check is to avoid the generation of negative values
!           that can occur in the upper stratosphere and mesosphere
!-----------------------------------------------------------------------
            qs(i,k) = MIN( 1._dp,qs(i,k) )
!++jsr changed qs(:,:) < 0._dp to qs(:,:) <= 0, otherwise 0 can occur. Later,
!      a quantity is divided by qs.
            if( qs(i,k) <= 0._dp ) then
               qs(i,k) = 1._dp
               es(i,k) = p(i,k)
            end if
!--jsr
         end do
      end do

      end subroutine AQSAT

      subroutine AQSATD( t, p, es, qs, gam, ii, ilen, kk, kstart, kend )
!-----------------------------------------------------------------------
!         ... Utility procedure to look up and return saturation vapor pressure from 
!           precomputed table, calculate and return saturation specific humidity 
!           (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT) for input 
!           arrays of temperature and pressure (dimensioned ii,kk).
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!         ... Input arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
              ii, &         ! I dimension of arrays t, p, es, qs
              kk, &         ! K dimension of arrays t, p, es, qs 
              ilen, &       ! vector length in I direction
              kstart, &     ! starting location in K direction
              kend          ! ending location in K direction
      real(dp), dimension(ii,kk), intent(in) :: &
           t, &             ! Temperature
           p                ! Pressure

      real(dp), dimension(ii,kk), intent(out) :: &
           es, &            ! Saturation vapor pressure
           qs, &            ! Saturation specific humidity
           gam              ! (l/cp)*(d(qs)/dt)

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      logical :: lflg          ! true if in temperature transition region
      integer :: i, &          ! i index for vector calculations
                 k             ! k index 
      real(dp) :: &
           omeps, &         ! 1. - 0.622
           trinv, &         ! reciprocal of ttrice (transition range)
           tc, &            ! temperature (in degrees C)
           weight, &        ! weight for es transition from water to ice
           hltalt, &        ! appropriately modified hlat for T derivatives
           hlatsb, &        ! hlat weighted in transition region
           hlatvp, &        ! hlat modified for t changes above 273.16
           tterm, &         ! account for d(es)/dT in transition region
           desdt            ! d(es)/dT

      omeps = 1._dp - epsqs
      do k = kstart,kend
!++vector
         do i = 1,ilen
            es(i,k) = ESTBLF1( t(i,k) )
!--vector
!-----------------------------------------------------------------------
!         ... Saturation specific humidity
!-----------------------------------------------------------------------
            qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
!-----------------------------------------------------------------------
!         ... The following check is to avoid the generation of negative
!           values that can occur in the upper stratosphere and mesosphere
!-----------------------------------------------------------------------
            qs(i,k) = MIN( 1._dp,qs(i,k) )
!++jsr changed qs < 0 in qs <=0 for having non-zero qs-values guaranteed
            if( qs(i,k) < 0._dp ) then
               qs(i,k) = 1._dp
               es(i,k) = p(i,k)
            end if
!--jsr
         end do
      end do
!-----------------------------------------------------------------------
!         ... "generalized" analytic expression for t derivative of es
!            accurate to within 1 percent for 173.16 < t < 373.16
!-----------------------------------------------------------------------
      trinv = 0._dp
      if( icephs .and. ttrice /= 0._dp ) then
         trinv = 1._dp/ttrice
         do k = kstart,kend
            do i = 1,ilen
               if( qs(i,k) == 1._dp ) then
                  gam(i,k) = 0._dp
               else
!-----------------------------------------------------------------------
!         ... Weighting of hlat accounts for transition from water to ice
!           polynomial expression approximates difference between es over
!           water and es over ice from 0 to -ttrice (C) (min of ttrice is
!           -40): required for accurate estimate of es derivative in transition 
!           range from ice to water also accounting for change of hlatv with t 
!           above 273.16 where constant slope is given by -2369 j/(kg c) =cpv - cw
!-----------------------------------------------------------------------
                  tc     = t(i,k) - 273.16_dp
                  lflg   = tc >= ttrice .and. tc < 0._dp
                  weight = MIN( -tc*trinv,1._dp )
                  hlatsb = hlatv + weight*hlatf
                  hlatvp = hlatv - 2369._dp*tc
                  if( t(i,k) < 273.16_dp ) then
                     hltalt = hlatsb
                  else
                     hltalt = hlatvp
                  end if
                  if( lflg ) then
                     tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
                  else
                     tterm = 0._dp
                  end if
                  desdt  = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k)) + tterm*trinv
                  gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt &
                             /(cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
               end if
            end do
         end do
      else
!-----------------------------------------------------------------------
!         ... No icephs or water to ice transition
!-----------------------------------------------------------------------
         do k = kstart,kend
            do i = 1,ilen
               if( qs(i,k) == 1._dp ) then
                  gam(i,k) = 0._dp
               else
!-----------------------------------------------------------------------
!         ... Account for change of hlatv with t above 273.16 where
!           constant slope is given by -2369 j/(kg c) = cpv - cw
!-----------------------------------------------------------------------
                  hlatvp = hlatv - 2369._dp*(t(i,k) - 273.16_dp)
                  if( icephs ) then
                     hlatsb = hlatv + hlatf
                  else
                     hlatsb = hlatv
                  end if
                  if( t(i,k) < 273.16_dp ) then
                     hltalt = hlatsb
                  else
                     hltalt = hlatvp
                  end if
                  desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k))
                  gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt &
                             /(cp*es(i,k)*(p(i,k) - omeps*es(i,k)))
               end if
            end do
         end do
      end if

      end subroutine AQSATD

      subroutine GESTBL( tmn, tmx, trice, ip, epsil, latvap, latice, rh2o, cpair )
!-----------------------------------------------------------------------
!         ... Builds saturation vapor pressure table for later lookup procedure.
!           Uses Goff & Gratch (1946) relationships to generate the table
!           according to a set of free parameters defined below.  Auxiliary
!           routines are also included for making rapid estimates (well with 1%)
!           of both es and d(es)/dt for the particular table configuration.
!-----------------------------------------------------------------------

      USE mo_mpi,        ONLY: p_parallel_io
      USE mo_exception,  ONLY: finish, message, message_text, em_debug

      implicit none

!-----------------------------------------------------------------------
!         ... Dummy arguments
!-----------------------------------------------------------------------
      real(dp), intent(in) :: &
           tmn, &        ! Minimum temperature entry in es lookup table
           tmx, &        ! Maximum temperature entry in es lookup table
           epsil, &      ! Ratio of h2o to dry air molecular weights
           trice, &      ! Transition range from es over range to es over ice
           latvap, &     ! Latent heat of vaporization
           latice, &     ! Latent heat of fusion
           rh2o, &       ! Gas constant for water vapor
           cpair         ! Specific heat of dry air
      logical, intent(in) :: &
           ip       ! Ice phase logical flag

!++vector
! --- define common
      real(dp) :: &
           tmin, &               ! min temperature (K) for table
           tmax                  ! max temperature (K) for table
      real(dp) :: estbl(plenest)     ! table values of saturation vapor pressure
      common /kk_com/ tmin,tmax,estbl
!--vector

!-----------------------------------------------------------------------
!         ... Local variables
!-----------------------------------------------------------------------
      real(dp) :: t             ! Temperature
      integer :: &
              n, &       ! Increment counter
              lentbl, &  ! Calculated length of lookup table
              itype      ! Ice phase: 0 -> no ice phase
                         !            1 -> ice phase, no transition
                         !           -x -> ice phase, x degree transition

!-----------------------------------------------------------------------
!         ... Set es table parameters
!-----------------------------------------------------------------------
      tmin   = tmn       ! Minimum temperature entry in table
      tmax   = tmx       ! Maximum temperature entry in table
      lentbl = INT( tmax - tmin + 2.000001_dp )
      if( lentbl > plenest ) then
         write(message_text,9000) tmax, tmin, plenest
         call finish('gestbl (MOZ)', message_text)
      end if
      ttrice = trice     ! Trans. range from es over h2o to es over ice
      icephs = ip        ! Ice phase (true or false)

!-----------------------------------------------------------------------
!         ... Set physical constants required for es calculation
!-----------------------------------------------------------------------
      epsqs  = epsil
      hlatv  = latvap
      hlatf  = latice
      rgasv  = rh2o
      cp     = cpair

!-----------------------------------------------------------------------
!         ... Begin building es table.
!           Check whether ice phase requested.
!           If so, set appropriate transition range for temperature
!-----------------------------------------------------------------------
      if( icephs ) then
         if( ttrice /= 0._dp ) then
            itype = INT( -ttrice )
         else
            itype = 1
         end if
      else
         itype = 0
      end if

      t = tmin - 1._dp
      do n = 1,lentbl
         t = t + 1._dp
         estbl(n) = GFFGCH( t, itype )
      end do

      do n = lentbl+1,plenest
         estbl(n) = -99999.0_dp
      end do

!-----------------------------------------------------------------------
!         ... Table complete -- Set coefficients for polynomial approximation of
!         ... difference between saturation vapor press over water and saturation
!         ... pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
!         ... is valid in the range -40 < t < 0 (degrees C).
!-----------------------------------------------------------------------
!                  --- Fifth order approximation ---
!-----------------------------------------------------------------------
      pcf(1) =  5.04469588506e-01_dp
      pcf(2) = -5.47288442819e+00_dp
      pcf(3) = -3.67471858735e-01_dp
      pcf(4) = -8.95963532403e-03_dp
      pcf(5) = -7.78053686625e-05_dp

      WRITE (message_text, '(a,a,i0,a,i0)') 'Saturation vapor pressure table completed.',  &
                              ' lentbl = ',lentbl, ', n = ', n
      CALL message('gestbl (MOZ)', message_text, level=em_debug)

 9000 format('GESTBL: FATAL ERROR *********************************',/, &
           ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH', &
           ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/, &
           ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)

      end subroutine GESTBL

      subroutine VQSATD( t, p, es, qs, gam, len )
!-----------------------------------------------------------------------
!         ... Utility procedure to look up and return saturation vapor pressure from 
!           precomputed table, calculate and return saturation specific humidity 
!           (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
!           function as qsatd, but operates on vectors of temperature and pressure
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!         ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
           len       ! vector length
      real(dp), dimension(len), intent(in) :: &
           t, &      ! temperature
           p         ! pressure

      real(dp), dimension(len), intent(out) :: &
           es, &     ! saturation vapor pressure
           qs, &     ! saturation specific humidity
           gam       ! (l/cp)*(d(qs)/dt)

!-----------------------------------------------------------------------
!         ... Local variables
!-----------------------------------------------------------------------
      logical :: &
           lflg      ! true if in temperature transition region
      integer :: &
           i         ! index for vector calculations
      real(dp) :: &
           omeps, &  ! 1. - 0.622
           trinv, &  ! reciprocal of ttrice (transition range)
           tc, &     ! temperature (in degrees C)
           weight, & ! weight for es transition from water to ice
           hltalt, & ! appropriately modified hlat for T derivatives  
           hlatsb, & ! hlat weighted in transition region
           hlatvp, & ! hlat modified for t changes above 273.16
           tterm, &  ! account for d(es)/dT in transition region
           desdt     ! d(es)/dT

      omeps     = 1._dp - epsqs
!++vector
      do i = 1,len
         es(i) = ESTBLF1( t(i) )
!--vector
!-----------------------------------------------------------------------
!         ... Saturation specific humidity
!-----------------------------------------------------------------------
         qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
!-----------------------------------------------------------------------
!         ... The following check is to avoid the generation of negative
!           values that can occur in the upper stratosphere and mesosphere
!-----------------------------------------------------------------------
         qs(i) = MIN( 1._dp,qs(i) )
         if( qs(i) < 0._dp ) then
            qs(i) = 1._dp
            es(i) = p(i)
         end if
      end do

!-----------------------------------------------------------------------
!         ... "generalized" analytic expression for t derivative of es
!            accurate to within 1 percent for 173.16 < t < 373.16
!-----------------------------------------------------------------------
      trinv = 0._dp
      if( icephs .and. ttrice /= 0._dp ) then
         trinv = 1._dp/ttrice
         do i = 1,len
            if( qs(i) == 1._dp ) then
               gam(i) = 0._dp
            else
!-----------------------------------------------------------------------
!         ... Weighting of hlat accounts for transition from water to ice
!           polynomial expression approximates difference between es over
!           water and es over ice from 0 to -ttrice (C) (min of ttrice is
!           -40): required for accurate estimate of es derivative in transition 
!           range from ice to water also accounting for change of hlatv with t 
!           above 273.16 where const slope is given by -2369 j/(kg c) = cpv - cw
!-----------------------------------------------------------------------
               tc     = t(i) - 273.16_dp
               lflg   = tc >= ttrice .and. tc < 0._dp
               weight = MIN( -tc*trinv,1._dp )
               hlatsb = hlatv + weight*hlatf
               hlatvp = hlatv - 2369._dp*tc
               if( t(i) < 273.16_dp ) then
                  hltalt = hlatsb
               else
                  hltalt = hlatvp
               end if
               if( lflg ) then
                  tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
               else
                  tterm = 0._dp
               end if
               desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
               gam(i) = hltalt*qs(i)*p(i)*desdt &
                        /(cp*es(i)*(p(i) - omeps*es(i)))
            end if
         end do
      else
!-----------------------------------------------------------------------
!         ... No icephs or water to ice transition
!-----------------------------------------------------------------------
         do i = 1,len
            if( qs(i) == 1._dp ) then
               gam(i) = 0._dp
            else
!-----------------------------------------------------------------------
!         ... Account for change of hlatv with t above 273.16 where
!           constant slope is given by -2369 j/(kg c) = cpv - cw
!-----------------------------------------------------------------------
               hlatvp = hlatv - 2369.0_dp*(t(i)-273.16_dp)
               if( icephs ) then
                  hlatsb = hlatv + hlatf
               else
                  hlatsb = hlatv
               end if
               if( t(i) < 273.16_dp ) then
                  hltalt = hlatsb
               else
                  hltalt = hlatvp
               end if
               desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
               gam(i) = hltalt*qs(i)*p(i)*desdt &
                        /(cp*es(i)*(p(i) - omeps*es(i)))
            end if
         end do
      end if

      end subroutine VQSATD

      real(dp) function GFFGCH( t, itype )
!-----------------------------------------------------------------------
!         ... Computes saturation vapor pressure over water and/or over ice using
!           Goff & Gratch (1946) relationships.  T (temperature), and itype are
!           input parameters, while es (saturation vapor pressure) is an output
!           parameter.  The input parameter itype serves two purposes: a value of
!           zero indicates that saturation vapor pressures over water are to be
!           returned (regardless of temperature), while a value of one indicates
!           that saturation vapor pressures over ice should be returned when t is
!           less than 273.16 degrees k.  If itype is negative, its absolute value
!           is interpreted to define a temperature transition region below 273.16
!           degrees k in which the returned saturation vapor pressure is a
!           weighted average of the respective ice and water value.  That is, in
!           the temperature range 0 => -itype degrees c, the saturation vapor
!           pressures are assumed to be a weighted average of the vapor pressure
!           over supercooled water and ice (all water at 0 c; all ice at -itype
!           c).  Maximum transition range => 40 c
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!         ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: &
                          itype      ! Flag for ice phase and associated transition
      real(dp), intent(in) :: t          ! Temperature

!-----------------------------------------------------------------------
!         ... Local variables
!-----------------------------------------------------------------------
      real(dp), parameter :: t0 = 273.16_dp
      real(dp), parameter :: ps = 1013.246_dp
      real(dp), parameter :: ts = 373.16_dp
      integer :: ltype
      real(dp) :: &
           e1, &      ! Intermediate scratch variable for es over water
           e2, &      ! Intermediate scratch variable for es over water
           eswtr, &   ! Saturation vapor pressure over water
           es, &      ! Saturation vapor pressure
           f, &       ! Intermediate scratch variable for es over water
           f1, &      ! Intermediate scratch variable for es over water
           f2, &      ! Intermediate scratch variable for es over water
           f3, &      ! Intermediate scratch variable for es over water
           f4, &      ! Intermediate scratch variable for es over water
           f5, &      ! Intermediate scratch variable for es over water
           term1, &   ! Intermediate scratch variable for es over ice
           term2, &   ! Intermediate scratch variable for es over ice
           term3, &   ! Intermediate scratch variable for es over ice
           tr, &      ! Transition range for es over water to es over ice
           tcrit, &
           temp, &
           weight     ! Intermediate scratch variable for es transition
      integer :: &
           itypo      ! Intermediate scratch variable for holding itype

!-----------------------------------------------------------------------
!         ... Check on whether there is to be a transition region for es
!-----------------------------------------------------------------------
      ltype = itype
      if( ltype < 0 ) then
         tr    = ABS(REAL(ltype,dp))
         ltype = 1
      else
         tr    = 0._dp
      end if
      if( tr > 40._dp ) then
         write(*,900) tr
        CALL finish('GFFGCH', 'tr > 40')
      end if

      tcrit = t0 - tr
      if( ltype == 1 .and. t < tcrit ) then
!-----------------------------------------------------------------------
!         ... Ice
!-----------------------------------------------------------------------
         temp  = t0/t
         term1 = 2.01889049_dp/temp
         term2 = 3.56654_dp*LOG( temp )
         term3 = 20.947031_dp*temp
         es    = 575.185606e10_dp*EXP(-(term1 + term2 + term3))
      else
!-----------------------------------------------------------------------
!         ... Water
!-----------------------------------------------------------------------
         temp = ts/t
         e1 = 11.344_dp*(1._dp - t/ts)
         e2 = -3.49149_dp*(temp - 1._dp)
         f1 = -7.90298_dp*(temp - 1._dp)
         f2 = 5.02808_dp*LOG10( temp )
         f3 = -1.3816_dp*(10._dp**e1 - 1._dp)/10000000._dp
         f4 = 8.1328_dp*(10._dp**e2 - 1._dp)/1000._dp
         f5 = LOG10( ps )
         f  = f1 + f2 + f3 + f4 + f5
         es = (10._dp**f)*100._dp
         eswtr = es
         if( ltype /= 0 .and. t < t0 ) then
            temp = t0/t
            term1 = 2.01889049_dp/temp
            term2 = 3.56654_dp*LOG( temp )
            term3 = 20.947031_dp*temp
            es    = 575.185606e10_dp*EXP( -(term1 + term2 + term3) )
            if( t >= tcrit ) then
!-----------------------------------------------------------------------
!         ... Weighted transition between water and ice
!-----------------------------------------------------------------------
               weight = MIN( (t0 - t)/tr,1._dp )
               es = weight*es + (1._dp - weight)*eswtr
            end if
         end if
      end if

      GFFGCH = es

  900 format('GFFGCH: FATAL ERROR ******************************',/, &
             'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', &
             ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', &
             ' 40.0 DEGREES C',/, ' TR = ',f7.2)

      end function GFFGCH

end module MO_MOZ_ESLOOKUP

      function ESTBLF( td )
!-----------------------------------------------------------------------
!         ... Form saturation specific humidity
!-----------------------------------------------------------------------
      use mo_kind,    only: dp

      implicit none

!-----------------------------------------------------------------------
!         ... Dummy arguments
!-----------------------------------------------------------------------
      real(dp), intent(in) :: td(:)

!++vector
     integer, parameter :: plenest = 250  ! length of saturation vapor pressure table
      real(dp) :: &
           tmin, &               ! min temperature (K) for table
           tmax                  ! max temperature (K) for table
      real(dp) :: estbl(plenest)     ! table values of saturation vapor pressure
      common /kk_com/ tmin,tmax,estbl
!--vector

!-----------------------------------------------------------------------
!         ... Local variables
!-----------------------------------------------------------------------
      integer :: i, k
      real(dp)    :: t1, t2, temp

      real(dp), dimension( SIZE(td) ) :: ESTBLF

      do k = 1,SIZE( td )
         t1        = MAX( MIN( td(k),tmax ),tmin )
         t2        = AINT( t1 - tmin )
         i         = INT( t1 - tmin )
         temp      = tmin + t2 - t1
         ESTBLF(k) = (temp + 1._dp)*estbl(i+1) - temp*estbl(i+2)
      end do

      end function ESTBLF

!++vector
      real function ESTBLF1( td )
!-----------------------------------------------------------------------
!       ... Form saturation specific humidity
!-----------------------------------------------------------------------
      use mo_kind,    only: dp

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      real(dp), intent(in) :: td

      integer, parameter :: plenest = 250  ! length of saturation vapor pressure table
      real(dp) :: &
           tmin, &               ! min temperature (K) for table
           tmax                  ! max temperature (K) for table
      real(dp) :: estbl(plenest)     ! table values of saturation vapor pressure
      common /kk_com/ tmin,tmax,estbl

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer :: i, k
      real(dp)    :: t1, t2, temp

      t1        = MAX( MIN( td,tmax ),tmin )
      t2        = AINT( t1 - tmin )
      i         = INT( t1 - tmin )
      temp      = tmin + t2 - t1
      ESTBLF1 = (temp + 1._dp)*estbl(i+1) - temp*estbl(i+2)

      end function ESTBLF1
!--vector
