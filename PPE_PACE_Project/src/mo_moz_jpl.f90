!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_jpl
!!
!! \brief
!!  This routine calculates third order reaction rate constants according to the Troe 
!!  formula
!!
!! \author Stacy Walters (NCAR)
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
!!  Even though this module contains only one small function it should be preserved 
!!  because jpl() is called from preprocessor-generated code.
!!
!! \bibliographic_references
!!  http://jpldataeval.jpl.nasa.gov/
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


      module mo_moz_jpl

      use mo_kind,  only : dp

      contains

      subroutine jpl( rate, m, factor, ko, kinf, plnplv )
!-----------------------------------------------------------------
!        ... Calculate JPL troe rate
!-----------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------
      integer , intent(in)  ::   plnplv   ! number of lons * number of levs
      real(dp), intent(in)  ::   factor
      real(dp), intent(in)  ::   ko(plnplv)
      real(dp), intent(in)  ::   kinf(plnplv)
      real(dp), intent(in)  ::   m(plnplv)
      real(dp), intent(out) ::   rate(plnplv)

!-----------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------
      real(dp)  ::     xpo(plnplv)
      

      xpo(:)  = ko(:) * m(:) / kinf(:)
      rate(:) = ko(:) / (1._dp + xpo(:))
      xpo(:)  = log10( xpo(:) )
      xpo(:)  = 1._dp / (1._dp + xpo(:)*xpo(:))
      rate(:) = rate(:) * factor**xpo(:)

      end subroutine jpl

      end module mo_moz_jpl
