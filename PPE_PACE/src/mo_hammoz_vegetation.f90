!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_vegetation.f90
!!
!! \brief
!! This is a module to temporarily replace the former mo_vegetation.f90 from
!! echam6.1, which has disapeared in echam6.3, because the variables
!! it used to define are now available as dynamically computed quantities
!! in jsbach. In waiting for implementing the usage of dynamic jsbach
!! evapotranspiration variables (see also #77), it is necessary to
!! re-introduce these as constants.
!!
!! \author S. Ferrachat (ETH Zurich)
!!
!! \responsible_coder
!! S. Ferrachat (sylvaine.ferrachat@env.ethz.ch)
!!
!! \revision_history
!!   -# S. Ferrachat (ETH Zurich) - original code (2014-11-19)
!!
!! \limitations
!! None
!!
!! \details
!! None
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
MODULE mo_hammoz_vegetation
  
  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! _______________________________________________________________________
  !
  ! module *mo_hammoz_vegetation* variables for computation of evapotranspiration
  !
  !
  !      Version 1  C.A Blondin  2/12/86  ECMWF
  !
  ! _______________________________________________________________________

  REAL(dp), PARAMETER :: cva = 5000._dp  !< constant to define the stomatal resistance
  REAL(dp), PARAMETER :: cvb = 10._dp    !< constant to define the stomatal resistance
  REAL(dp), PARAMETER :: cvc = 100._dp   !< minimum stomatal resistance
  REAL(dp), PARAMETER :: cvk = .9_dp     !< 
  REAL(dp), PARAMETER :: cvbc = 1000._dp !< cvb*cvc
  REAL(dp), PARAMETER :: cvkc = 90._dp   !< cvk*cvc
  REAL(dp), PARAMETER :: cvabc = 60._dp  !< (cva+cvbc)/cvc
  REAL(dp), PARAMETER :: cvrad = 0.55_dp !< fraction of the net s.w radiation contributing to p.a.r

END MODULE mo_hammoz_vegetation
