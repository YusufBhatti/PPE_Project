!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
      MODULE mo_string_utls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! utilities for easier string handling
!
! Authors:
!
! J.S. Rast: original source
! M. Schultz, FZ Juelich  (2009-10-02)  - changed rteurn value of st1_in_st2_proof to boolean
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PRIVATE
      PUBLIC        :: st1_in_st2_proof, st1_in_st2_idx

INTERFACE st1_in_st2_proof
   MODULE PROCEDURE st1_in_st2_proof_0d
   MODULE PROCEDURE st1_in_st2_proof_1d
END INTERFACE st1_in_st2_proof


      CONTAINS

      FUNCTION st1_in_st2_proof_0d(st1, st2, ierr)   RESULT (lexist)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! verifies that string st1 is contained in the vector st2.
! st1_in_st2_proof_0d is 
!   true:  if string st1 exists in st2 
!          or st1 only contains the zero string
!   false: if string st1 does not exist in st2
! The optional ierr parameter will return the following values:
!       0: if string st1 exists in st2
!      -1: if string st1 is empty   (function result is still true!)
!       1: if string st1 does not exist in st2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

      CHARACTER (len=*), INTENT(in)    :: st1,st2(:)
      INTEGER, OPTIONAL, INTENT(out)   :: ierr
      INTEGER                          :: jj,n2
      LOGICAL                          :: lexist

      lexist = .false.
      IF (PRESENT(ierr)) ierr = -999   ! ### should not occur!
      n2=SIZE(st2)
      IF (TRIM(st1) /= '') THEN
         DO jj=1,n2
            IF (TRIM(st1) == TRIM(st2(jj))) THEN
               lexist=.true.
               IF (PRESENT(ierr)) ierr = 0
               EXIT
            END IF
         END DO
         IF (.not. lexist) THEN
            IF (PRESENT(ierr)) ierr = 1
         END IF
      ELSE
         IF (PRESENT(ierr)) ierr = -1
      END IF
      END FUNCTION st1_in_st2_proof_0d

      FUNCTION st1_in_st2_proof_1d(st1, st2, ierr)   RESULT(lexist)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! verifies that all string elements of st1 are also contained in st2.
! st1_in_st2_proof_1d is 
!   true:  if all string elements in st1 also lexist in st2 
!          or st1 only contains zero strings
!   false: if any string element in st1 is not contained in the vector st2.
! The optional ierr parameter will return the following values:
!       0: if all strings in st1 exist in st2
!      -1: if any string in st1 is empty   (function result may still be true!)
!       n: if the function result is false with n the index of first string 
!          element in st1 which does not match any element of st2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

! local variables

      CHARACTER (len=*), INTENT(in)    :: st1(:),st2(:)
      INTEGER, OPTIONAL, INTENT(out)   :: ierr
      INTEGER                          :: ii,jj,n1,n2
      LOGICAL                          :: lexist

      lexist = .false.
      IF (PRESENT(ierr)) ierr = 0
      n1=SIZE(st1)
      n2=SIZE(st2)
      DO ii=1,n1
         IF (TRIM(st1(ii)) /= '') THEN
            lexist = .false.
            DO jj=1,n2
               IF (TRIM(st1(ii)) == TRIM(st2(jj))) THEN
                  lexist=.true.
                  EXIT
               END IF
            END DO
            IF (.not. lexist) THEN
               IF (PRESENT(ierr)) ierr = ii
               EXIT
            END IF
         ELSE
            IF (PRESENT(ierr)) THEN
               IF (ierr == 0) ierr = -1
            END IF
         END IF
      END DO
      END FUNCTION st1_in_st2_proof_1d

      INTEGER FUNCTION st1_in_st2_idx(st1,st2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gives index of st1 in vector of strings st2, 0 if st1 is not in st2
! st1_in_st2_idx is 
!   i: if string st1 lexists as i-th element in st2
!   0: if string st1 does not lexists in st2 or st1 is the zero string
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

! local variables

      CHARACTER (len=*), INTENT(in)    :: st1, st2(:)
      INTEGER                          :: jj, n2
!      LOGICAL                          :: lexist

      st1_in_st2_idx=0
      n2=SIZE(st2)
      IF ( TRIM(st1) == '' ) THEN
         st1_in_st2_idx = 0 
      ELSE
         DO jj=1,n2
            IF (TRIM(st1) == TRIM(st2(jj))) THEN
               st1_in_st2_idx = jj 
               EXIT
            END IF
         END DO
      END IF
      END FUNCTION st1_in_st2_idx

      END MODULE mo_string_utls
