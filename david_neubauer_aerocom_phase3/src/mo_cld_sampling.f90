!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form. 
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module to provide Monte Carlo samples based on cloud fraction
!!
!! @remarks
!!   This module draws Monte Carlo samples from profiles of the atmospheric
!!   state and on overlap assumption ( 1=max-ran, 2=maximum, 3=random). 
!!   Users provide profiles of cloud fraction and in-cloud optical thickness. 
!!   A single sample is drawn; cloud fraction is either 1 or 0 and optical 
!!   thickness is set to 0 where cloud fraction is 0. 
!!
!! @author Robert Pincus, U. Colorado, while visiting MPI-M, Hamburg (2010-08)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Written by Robert Pincus; simplifed from code originally written for the 
!!   GFDL atmospheric model AM2. 
!!
!! @par Copyright
!!   The AER copyright
!!   on the original code is as follows: Copyright 2002-2009, Atmospheric and
!!   Environmental Research, Inc. (AER). This software may be used, copied, or
!!   redistributed as long as it is not sold and this copyright notice is
!!   reproduced on each copy made.  This model is provided as is without any
!!   express or implied warranties. (http://www.rtweb.aer.com/)               
!! 

MODULE mo_cld_sampling 

  USE mo_kind, ONLY          : wp
  USE mo_exception, ONLY     : finish
  USE mo_random_numbers, ONLY: get_random

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sample_cld_state

CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! @brief Returns a sample of the cloud state
  !!
  !! @remarks  
  !
  SUBROUTINE sample_cld_state(kproma, kbdim, klev, ksamps, rnseeds, i_overlap, cld_frac, cldy)

    INTEGER,  INTENT(IN)    :: kproma, kbdim, klev, ksamps !< numbers of columns, levels, samples 
    INTEGER,  INTENT(INOUT) :: rnseeds(:, :)           !< Seeds for random number generator (kbdim, :) 
    INTEGER,  INTENT(IN)    :: i_overlap               !< 1=max-ran, 2=maximum, 3=random 
    REAL(wp), INTENT(IN)    :: cld_frac(kbdim,klev)    !< cloud fraction
    LOGICAL, INTENT(OUT)    :: cldy(kbdim,klev,ksamps) !< Logical: cloud present? 

    REAL(wp) :: rank(kbdim,klev,ksamps)
    INTEGER  :: jk, js 

    ! Here cldy(:,:,1) indicates whether any cloud is present 
    !
    cldy(1:kproma,1:klev,1) = cld_frac(1:kproma,1:klev) > 0._wp
    SELECT CASE(i_overlap) 
    CASE(1) 
      ! Maximum-random overlap
      DO js = 1, ksamps
        DO jk = 1, klev
          ! mask means we compute random numbers only when cloud is present 
          CALL get_random(kproma, kbdim, rnseeds, cldy(:,jk,1), rank(:,jk,js))
        END DO 
      END DO 
      ! There may be a better way to structure this calculation...
      DO jk = klev-1, 1, -1
        DO js = 1, ksamps
          rank(1:kproma,jk,js) = MERGE(rank(1:kproma,jk+1,js),                                 & 
                                     ! Max overlap... 
                                     rank(1:kproma,jk,js) * (1._wp - cld_frac(1:kproma,jk+1)), & 
                                     ! ... or random overlap in the clear sky portion,  
                                     ! depending on whether or not you have cloud in the layer above 
                                     rank(1:kproma,jk+1,js) > 1._wp - cld_frac(1:kproma,jk+1) ) 
        END DO
      END DO  
    CASE(2) 
      !
      !  Max overlap means every cell in a column is identical 
      ! 
      DO js = 1, ksamps
        CALL get_random(kproma, kbdim, rnseeds, rank(:, 1, js))
        rank(1:kproma,2:klev,js) = SPREAD(rank(1:kproma,1,js), DIM=2, NCOPIES=(klev-1))
      END DO 
    CASE(3) 
      !
      !  Random overlap means every cell is independent
      ! 
      DO js = 1, ksamps
        DO jk = 1, klev
          ! mask means we compute random numbers only when cloud is present 
          CALL get_random(kproma, kbdim, rnseeds, cldy(:,jk,1), rank(:,jk,js))
        END DO 
      END DO 
    CASE DEFAULT
      CALL finish('In sample_cld_state: unknown overlap assumption') 
    END SELECT
    

    ! Now cldy indicates whether the sample (ks) is cloudy or not.    
    DO js = 1, ksamps
      cldy(1:kproma,1:klev,js) = rank(1:kproma,1:klev,js) > (1. - cld_frac(1:kproma,1:klev))
    END DO 
  
  END SUBROUTINE sample_cld_state

END MODULE mo_cld_sampling
