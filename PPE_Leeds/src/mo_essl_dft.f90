!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
#ifdef ESSL_DFT
MODULE mo_essl_dft

  USE mo_kind,      ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_essl_dft
  PUBLIC :: cleanup_essl_dft
  PUBLIC :: essl_dft

  INTEGER :: naux1, naux2
  INTEGER :: naux3 = 1
  
  REAL(dp), ALLOCATABLE :: aux1_r2c(:), aux1_c2r(:), aux2(:)
  REAL(dp) :: aux3(1)

CONTAINS

  SUBROUTINE init_essl_dft(nlon, nofs, nseq)
    
    INTEGER, INTENT(in) :: nlon, nofs, nseq
    REAL(dp) :: x(1),y(1) ! dummy arrays for initialization of ESSL
    
    ! use processor independent formulas for allocating ESSL work arrays
    
    IF (nlon <= 4096) THEN
      naux1 = 32000
      naux2 = 22000
    ELSE
      naux1 = 30000 + 1.64_dp * nlon
      naux2 = 20000 + 1.14_dp * nlon
    END IF
    
    IF (ALLOCATED(aux1_r2c)) DEALLOCATE(aux1_r2c)
    IF (ALLOCATED(aux1_c2r)) DEALLOCATE(aux1_c2r)
    IF (ALLOCATED(aux2))     DEALLOCATE(aux2)
    
    ALLOCATE(aux1_r2c(naux1))
    ALLOCATE(aux1_c2r(naux1))
    ALLOCATE(aux2(naux2))
    
    ! initialize forward transform
    
    CALL DCRFT(1, x, nofs/2, y, nofs, nlon, nseq, -1, 1.0_dp, &
         aux1_c2r, naux1, aux2, naux2, aux3, naux3)
    
    ! initialize backward transform
    
    CALL DRCFT(1, x, nofs, y, nofs/2, nlon, nseq, 1, 1.0_dp/nlon, &
         aux1_r2c, naux1, aux2, naux2, aux3, naux3)
    
  END SUBROUTINE init_essl_dft

  SUBROUTINE cleanup_essl_dft
    
    IF (ALLOCATED(aux1_r2c)) DEALLOCATE(aux1_r2c)
    IF (ALLOCATED(aux1_c2r)) DEALLOCATE(aux1_c2r)
    IF (ALLOCATED(aux2))     DEALLOCATE(aux2)
    
  END SUBROUTINE cleanup_essl_dft
  
  SUBROUTINE essl_dft(a, nlon, nofs, nseq, isign)

    ! .. Scalar Arguments ..
    INTEGER, INTENT(in) :: isign, nlon, nofs, nseq
    ! ..
    ! .. Array Arguments ..
    REAL(dp), INTENT(inout) :: a(*)

    INTEGER :: i

    IF (isign > 0) THEN

      ! complex to real transform

      CALL DCRFT(0, a, nofs/2, a, nofs, nlon, nseq, -1, 1.0_dp, &
           aux1_c2r, naux1, aux2, naux2, aux3, naux3)
              
    ELSE

      ! complex to real transform

      CALL DRCFT(0, a, nofs, a, nofs/2, nlon, nseq, 1, 1.0_dp/nlon, &
           aux1_r2c, naux1, aux2, naux2, aux3, naux3)
       
    END IF

  END SUBROUTINE essl_dft

END MODULE mo_essl_dft

#else

MODULE mo_essl_dft
END MODULE mo_essl_dft

#endif
