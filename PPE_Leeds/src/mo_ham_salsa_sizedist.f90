!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_salsa_sizedist.f90
!!
!! \brief
!! This module defines the sizedistribution used in SALSA.
!!
!! \author Harri Kokkola, harri.kokkola@fmi.fi
!!
!! \responsible_coder
!! Harri Kokkola, harri.kokkola@fmi.fi
!!
!! \revision_history
!!   -# Harri Kokkola (FMI) - Implementation of SALSA aerosol microphysics model (2014)
!!   -# Tommi Bergman (FMI) - Comment cleanup (November 2014)
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
MODULE mo_ham_salsa_sizedist

CONTAINS 

  SUBROUTINE size_distribution(kproma, kbdim, klev, &
       n, dpg, sigmag, naero)

    USE mo_ham_salsactl, ONLY :     &
         nreg,                      &
         vhilim,                    &
         vlolim,                    &
         dpmid,                     &
         in1a,                      &
         fn2b

    USE mo_kind, ONLY : dp

    USE mo_math_constants, ONLY : pi, &
         pi_6

    IMPLICIT NONE

    INTEGER, PARAMETER :: nmod = 7

    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid points 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical levels 

    REAL(dp), INTENT(IN) ::         &
         n(nmod)                  , & ! total concentration of a mode
         dpg(nmod)                , & ! geometric-mean diameter of a mode
         sigmag(nmod)                 ! standard deviation of a mode

    REAL(dp), INTENT(OUT) ::        &
         naero(kbdim,klev,fn2b)      ! number concentration  [#/m3]

    !-- local variables
    REAL(dp) ::                     &
         deltadp                      ! bin width [m]

    INTEGER :: ii, jj, kk

    naero = 0.

    DO jj = 1,klev    ! vertical grid
       DO ii = 1,kproma ! horizontal grid

          DO kk = in1a, fn2b

             deltadp = (vhilim(kk)**(1._dp/3._dp)-vlolim(kk)**(1._dp/3._dp))/   &
                  pi_6**(1._dp/3._dp)

             !-- size distribution
             !   ntot = total number, total area, or total volume concentration
             !   dpg = geometric-mean number, area, or volume diameter
             !   n(kk) = number, area, or volume concentration in a bin
             naero(ii,jj,kk) = sum(n*deltadp/                        &
                  (dpmid(kk)*sqrt(2._dp*pi)*log(sigmag))*                 &
                  exp(-log(dpmid(kk)/dpg)**2/(2._dp*log(sigmag)**2)))

          END DO

       END DO

    END DO

  END SUBROUTINE size_distribution

END MODULE mo_ham_salsa_sizedist
