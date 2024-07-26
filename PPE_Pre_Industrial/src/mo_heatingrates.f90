!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_heatingrates
!-------------------------------------------------------------------------
!
!    Sebastian Rast, MPI Met, Hamburg, January 2010
!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: mo_heatingrates
!
! !DESCRIPTION:
! Diagnostic output of heating rates for solar and longwave radiation
!
! !REVISION HISTORY:
! original source by J.S.Rast, (2010-01-15)
!
! !USES:

  USE mo_kind,              ONLY: wp

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC              :: set_heatingrates

! !LOCAL VARIABLES

CONTAINS
!EOP
!-------------------------------------------------------------------------
!BOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_heatingrates
!
! !SUBROUTINE INTERFACE:
SUBROUTINE set_heatingrates( &
                          kproma       ,kbdim           ,klev             ,&
                          klevp1                        ,pconvfact        ,&
                          pflxs        ,pflxt           ,pdtdts           ,&
                          pdtdtt                                           )
!
! !DESCRIPTION:
! sets heating rates for output
!
! !REVISION HISTORY:
! original source by J.S. Rast (2010-01-15)
!
! !USES:
! !PARAMETER
  INTEGER, INTENT(in)        :: kproma, kbdim, klev, klevp1
  REAL(wp), INTENT(in)       :: pconvfact(kbdim,klev)
  REAL(wp), INTENT(in)       :: pflxs(kbdim,klevp1), pflxt(kbdim,klevp1)
!< pdtdtt, pdtdts: heating rates thermal and solar spectrum in K/day 
  REAL(wp), INTENT(out)      :: pdtdtt(kbdim,klev),pdtdts(kbdim,klev) 
! !LOCAL VARIABLES
  INTEGER                    :: jk

  DO jk = 1, klev
     pdtdts(1:kproma,jk)= pconvfact(1:kproma,jk) * &
          (pflxs(1:kproma,jk+1)-pflxs(1:kproma,jk))
     pdtdtt(1:kproma,jk)= pconvfact(1:kproma,jk) * &
          (pflxt(1:kproma,jk+1)-pflxt(1:kproma,jk))
  END DO


END SUBROUTINE set_heatingrates
!EOP
!-------------------------------------------------------------------------
!EOC
!-------------------------------------------------------------------------
END MODULE mo_heatingrates
