!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_truncation

  USE mo_kind,      ONLY: dp
  USE mo_control,   ONLY: nlev, nhgl
  USE mo_exception, ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nmp, nnp, mcrit, am, ntrn, ntrm, ntrk
  PUBLIC :: init_truncation
  PUBLIC :: cleanup_truncation

  ! ---------------------------------------------------------------
  !
  ! module *mo_truncation* - quantities related to the spectral truncation.
  !
  ! ---------------------------------------------------------------

  INTEGER, ALLOCATABLE :: ntrm(:)       ! max zonal wave number.
  INTEGER, ALLOCATABLE :: ntrk(:)       ! max meridional wave number.
  INTEGER, ALLOCATABLE :: ntrn(:)       ! max meridional wave number for m=0.
  INTEGER, ALLOCATABLE :: nmp(:)        ! displacement of the first point of 
                                        ! columns for computations in 
                                        ! spectral space.
  INTEGER, ALLOCATABLE :: nnp(:)        ! number of points on each column.
  INTEGER, ALLOCATABLE :: mcrit(:)      ! critical zonal wave number depending
                                        ! on the latitude line,beyond which
                                        ! Fourier components are ignored.
  REAL(dp), ALLOCATABLE :: am(:)        ! REAL(m).i.e jm-1 in the jm loops.

CONTAINS

  SUBROUTINE init_truncation (nm, nn, nk)

    INTEGER ,INTENT(in) :: nm ! max zonal wave number
    INTEGER ,INTENT(in) :: nn ! max meridional wave number for m=0
    INTEGER ,INTENT(in) :: nk ! max meridional wave number

    ! Description:
    !
    ! Computes parameters used for computations in spectral space.
    !
    ! Method:
    !
    ! This subroutine computes some parameters related to the
    ! truncation and used for the computations in spectral space and
    ! for the *Legendre transforms.
    !
    ! *scpar* is called from *initialise*
    !
    ! Results:
    ! The results are stored in arrays in module *mo_truncation*
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, March 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! 
    ! for more details see file AUTHORS
    !

!   USE mo_control,    only: nkp1, nmp1, nn
!   USE mo_truncation, only: am, nmp, nnp

    !  Local scalars: 
    INTEGER :: jm, jmax
    CHARACTER(len=6) :: yfmt

    !  Intrinsic functions 
    INTRINSIC MIN

    !  Executable statements 

    IF (.NOT. ALLOCATED (ntrm)) THEN
      ALLOCATE (ntrm(nlev))
      ALLOCATE (ntrk(nlev))
      ALLOCATE (ntrn(nlev))
      ALLOCATE (mcrit(nhgl))
      ALLOCATE (nmp(nm+1))
      ALLOCATE (nnp(nm+1))
      ALLOCATE (am(nm+1))
    END IF

!-- 0. These parameters may be overwritten later in 'setdyn'

    ntrm (:) = nm
    ntrn (:) = nn
    ntrk (:) = nk
    mcrit(:) = nm+1

!-- 1. Preliminary computations

!-- 2. Compute parameters

    DO jm = 1, nm+1
      nnp(jm) = MIN(nk+1-jm,nn) + 1
    END DO

    nmp(1) = 0
    DO jm = 2, nm+1
      nmp(jm) = nmp(jm-1) + nnp(jm-1)
    END DO

!-- 3. Fill arrays *am* and *annp1*

    DO jm = 1, nm+1
      am(jm) = jm - 1.0_dp
    END DO

    WRITE (message_text, '(a)') ' Number of points on each column (NNP): '  
    CALL message('', message_text)
    DO jm = 1, nm+1, 11
      jmax = MIN(jm+10, nm+1)
      WRITE (yfmt,'(a1,i2,a3)') '(', jmax-jm+1, 'i7)' 
      WRITE (message_text, yfmt) nnp(jm:jmax)
      CALL message('', message_text)
    ENDDO
    WRITE (message_text, '(a)') ' Displacement of the first point of columns (NMP): '  
    CALL message('', message_text)
    DO jm = 1, nm+1, 11
      jmax = MIN(jm+10, nm+1)
      WRITE (yfmt,'(a1,i2,a3)') '(', jmax-jm+1, 'i7)' 
      WRITE (message_text, yfmt) nmp(jm:jmax)                       
      CALL message('', message_text)
    ENDDO

  END SUBROUTINE init_truncation

  SUBROUTINE cleanup_truncation

    IF (ALLOCATED (ntrm)) THEN
      DEALLOCATE (ntrm)
      DEALLOCATE (ntrk)
      DEALLOCATE (ntrn)
      DEALLOCATE (mcrit)
      DEALLOCATE (nmp)
      DEALLOCATE (nnp)
      DEALLOCATE (am)
    END IF

  END SUBROUTINE cleanup_truncation

END MODULE mo_truncation
