!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_spectral
  ! I. Kirchner, MPI, May 2002

  ! simple spectral transform

  USE mo_kind,       ONLY: dp
  USE mo_exception,  ONLY: finish, message

  IMPLICIT NONE

  PRIVATE

  PUBLIC gp2sp
  PUBLIC sp2gp
  PUBLIC sp2sp

  PUBLIC corrsp
  INTERFACE corrsp
     MODULE PROCEDURE corrsp3d
     MODULE PROCEDURE corrsp2d
  END INTERFACE

  PUBLIC expand
  INTERFACE expand
     MODULE PROCEDURE expand_2d
     MODULE PROCEDURE expand_3d
     MODULE PROCEDURE expand_4d
  END INTERFACE

  PUBLIC compress
  INTERFACE compress
     MODULE PROCEDURE compress_2d
     MODULE PROCEDURE compress_3d
     MODULE PROCEDURE compress_4d
  END INTERFACE

  CHARACTER(len=256) :: my_mess

CONTAINS

  !--------------------------------------------------------------------
  ! correlation between two 3-dim spectral fields
  FUNCTION corrsp3d(sp1,sp2) RESULT(cc)

    REAL(kind=dp), INTENT(in)  :: sp1(:,:,:), sp2(:,:,:)  ! (lev,2,nsp)
    REAL(kind=dp)              :: cc           ! spatial correlation coefficient

    REAL(kind=dp) :: sum1, sum2, sum12, val1, val2, ww
    INTEGER       :: jl, jk, js, nsp, nlev, nnp1
    nlev = SIZE(sp1,1)
    IF (SIZE(sp2,1) /= nlev) CALL finish('mo_spectral:corrsp3d',&
         'No. of levels mismatch')
    nsp  = SIZE(sp1,3)
    IF (SIZE(sp2,3) /= nsp)  CALL finish('mo_spectral:corrsp3d',&
         'No. of spectral components mismatch')
    nnp1  = get_nnp1(nsp)
    cc    = 0.0_dp
    sum12 = 0.0_dp
    sum1  = 0.0_dp
    sum2  = 0.0_dp

    DO js=1,nsp
       IF (js > nnp1) THEN
          ww = 2.0_dp
       ELSE IF (js > 1) THEN
          ww = 1.0_dp
       ELSE
          ww = 0.0_dp
       END IF
       DO jk=1,2
          DO jl=1,nlev
             val1  = sp1(jl,jk,js)
             val2  = sp2(jl,jk,js)
             sum12 = sum12  + ww*val1*val2
             sum1  = sum1   + ww*val1*val1
             sum2  = sum2   + ww*val2*val2
          END DO
       END DO
    END DO
    IF (sum1>0.0_dp .AND. sum2>0.0_dp) then
       cc = sum12/(SQRT(sum1)*SQRT(sum2))
    ELSE
       WRITE(my_mess,*) 'SUM1= ',sum1,' SUM2= ',sum2,' EPSILON= ',EPSILON(1._dp)
       CALL message('mo_spectral:corrsp3d',my_mess)
    END IF

  END FUNCTION corrsp3d

  ! correlation between 2-dim spectral fields
  FUNCTION corrsp2d(sp1,sp2) RESULT(cc)
    REAL(kind=dp), INTENT(in)  :: sp1(:,:), sp2(:,:)  ! 2,nsp
    REAL(kind=dp) :: cc

    REAL(kind=dp) :: sum1, sum2, sum12, val1, val2, ww
    INTEGER       :: jk, js
    INTEGER       :: nsp, nnp1

    nsp  = SIZE(sp1,2)
    IF (SIZE(sp2,2) /= nsp)  CALL finish('mo_spectral:corrsp3d',&
         'No. of spectral components mismatch')
    nnp1  = get_nnp1(nsp)
    cc    = 0.0_dp
    sum12 = 0.0_dp
    sum1  = 0.0_dp
    sum2  = 0.0_dp

    DO js=1,nsp
       IF (js > nnp1) THEN
          ww = 2.0_dp
       ELSE IF (js > 1) THEN
          ww = 1.0_dp
       ELSE
          ww = 0.0_dp
       END IF
       DO jk=1,2
          val1  = sp1(jk,js)
          val2  = sp2(jk,js)
          sum12 = sum12  + ww*val1*val2
          sum1  = sum1   + ww*val1*val1
          sum2  = sum2   + ww*val2*val2
       END DO
    END DO
    IF (sum1>0.0_dp .AND. sum2>0.0_dp) then
       cc = sum12/(SQRT(sum1)*SQRT(sum2))
    ELSE
       WRITE(my_mess,*) 'SUM1= ',sum1,' SUM2= ',sum2,' EPSILON= ',EPSILON(1.0_dp)
       CALL message('mo_spectral:corrsp2d',my_mess)
    END IF

  END FUNCTION corrsp2d

  !--------------------------------------------------------------------
  ! routines only for scalar fields

!  USE mo_truncation, ONLY: nmp, nnp
!  USE mo_control,    ONLY: nsp, ngl, nmp1, nhgl, nlon 

  SUBROUTINE gp2sp(grid,spec)
    REAL(kind=dp), INTENT(in)  :: grid(:,:,:)  ! lon, no, lat
    REAL(kind=dp), INTENT(out) :: spec(:,:,:)  ! no, 2, nsp

    INTEGER :: no

    no  = SIZE(grid,2)
    IF(SIZE(spec,1)/=no) CALL finish('mo_spectral:gp2sp',&
         'No. of levels mismatch')
    IF(SIZE(spec,2)/=2) CALL finish('mo_spectral:gp2sp',&
         'second dimension of spectral array not 2')

    spec(:,:,:) = 0.0_dp

    ! fftd
    ! sym1
    ! ldt

  END SUBROUTINE gp2sp

  SUBROUTINE sp2gp(spec,grid)
    REAL(kind=dp), INTENT(in)  :: spec(:,:,:)  ! no, 2, nsp
    REAL(kind=dp), INTENT(out) :: grid(:,:,:)  ! lon, no, lat

    INTEGER :: no

    no  = SIZE(grid,2)
    IF(SIZE(spec,1)/=no) CALL finish('mo_spectral:sp2gp',&
         'No. of levels mismatch')
    IF(SIZE(spec,2)/=2) CALL finish('mo_spectral:sp2gp',&
         'second dimension of spectral array not 2')

    grid(:,:,:) = 0.0_dp

    ! lti
    ! sym2
    ! fft

  END SUBROUTINE sp2gp

  ! change truncation
  SUBROUTINE sp2sp(sp1,sp2)
    REAL(kind=dp), INTENT(in)  :: sp1(:,:,:)  ! no, 2, nsp
    REAL(kind=dp), INTENT(out) :: sp2(:,:,:)  ! no, 2, nsp
    INTEGER :: nsp1, nsp2, no, tr1, tr2, k
    INTEGER, ALLOCATABLE :: nnp1(:), nmp1(:), nnp2(:), nmp2(:)

    no = SIZE(sp1,1)
    IF (SIZE(sp2,1)/=no) CALL finish('mo_spectral:sp2sp',&
         'No. of fields not equal')
    IF(SIZE(sp1,2)/=2) CALL finish('mo_spectral:sp2sp',&
         'second dimension of first spectral array not 2')
    IF(SIZE(sp2,2)/=2) CALL finish('mo_spectral:sp2sp',&
         'second dimension of second spectral array not 2')
    nsp1 = SIZE(sp1,3)
    nsp2 = SIZE(sp2,3)

    IF (nsp1 == nsp2) THEN  ! nothing to do
       sp2 = sp1
    ELSE
       ! prepare index fields
       tr1 = find_trunc(nsp1)
       ALLOCATE(nnp1(tr1+1),nmp1(tr1+1))
       CALL get_mp(tr1,nmp1,nnp1)
       tr2 = find_trunc(nsp2)
       ALLOCATE(nnp2(tr2+1),nmp2(tr2+1))
       CALL get_mp(tr2,nmp2,nnp2)
       sp2(:,:,:) = 0.0_dp
       IF (nsp1 < nsp2) THEN ! expand field
          DO k = 1, tr1+1
             sp2(:,:,nmp2(k)+1:nmp2(k)+nnp1(k)) = sp1(:,:,nmp1(k)+1:nmp1(k)+nnp1(k))
          END DO
       ELSE ! truncate field
          DO k = 1, tr2+1
             sp2(:,:,nmp2(k)+1:nmp2(k)+nnp2(k)) = sp1(:,:,nmp1(k)+1:nmp1(k)+nnp2(k))
          END DO
       END IF
       DEALLOCATE(nnp1,nmp1,nnp2,nmp2)
    END IF

  END SUBROUTINE sp2sp

  !--------------------------------------------------------------------
  ! simple functions only for triangular truncation
  FUNCTION get_nnp1(nsp) RESULT(nnp1)
    INTEGER, INTENT(in) :: nsp  ! number of spectral coefficients
    INTEGER             :: nnp1 ! wave 0 zonal index
    ! works only for triangular truncation
    INTEGER :: idx, count
    idx   = 1
    count = nsp
    DO
       nnp1  = count
       count = count - idx
       IF (count < 1) EXIT 
       idx   = idx + 1
    END DO
  END FUNCTION get_nnp1

  FUNCTION find_trunc(nsp) RESULT(tr)
    ! calculates the truncation from the number of
    ! spectral coefficients
    INTEGER, INTENT(in) :: nsp
    INTEGER :: tr, count
    tr    = 1
    count = nsp
    DO
       count = count - tr
       IF (count < 1) EXIT
       tr = tr + 1
    END DO
    tr = tr - 1
  END FUNCTION find_trunc

  SUBROUTINE get_mp(tr,nmp,nnp)
    INTEGER, INTENT(in) :: tr  ! truncation
    INTEGER, INTENT(out) :: nmp(:), nnp(:)
    INTEGER :: i, idx, nmp1

    nmp1 = tr+1
    IF (SIZE(nmp)/=nmp1) CALL finish('mo_spectral:get_mp',&
         'pointer array NMP too small')
    IF (SIZE(nnp)/=nmp1) CALL finish('mo_spectral:get_mp',&
         'pointer array NNP too small')
    idx = 1
    DO i=nmp1,1,-1
       nnp(i) = idx
       idx = idx + 1
    END DO
    idx = 0
    DO i=1,nmp1
       nmp(i) = idx
       idx = idx + nnp(i)
    END DO
  END SUBROUTINE get_mp

  !--------------------------------------------------------------------
  ! compress/expand fields for FFT
  SUBROUTINE expand_2d(fout,f2d,offset)
    REAL(dp), INTENT(out) :: fout(:)
    REAL(dp), INTENT(in)  :: f2d(:,:)
    INTEGER,  INTENT(in)  :: offset
    INTEGER :: i, j1, d1, nlon
    i = 0
    nlon = SIZE(f2d,1)
    d1   = SIZE(f2d,2)
    DO j1=1,d1
       fout(i+1:i+nlon) = f2d(:,j1)
       i = i + nlon + offset
    END DO
  END SUBROUTINE expand_2d

  SUBROUTINE expand_3d(fout,f3d,offset)
    REAL(dp), INTENT(out) :: fout(:)
    REAL(dp), INTENT(in)  :: f3d(:,:,:)
    INTEGER,  INTENT(in)  :: offset
    INTEGER :: i, j1, j2, d1, d2, nlon
    nlon = SIZE(f3d,1)
    d1   = SIZE(f3d,2)
    d2   = SIZE(f3d,3)
    i = 0
    DO j2=1,d2
       DO j1=1,d1
          fout(i+1:i+nlon) = f3d(:,j1,j2)
          i = i + nlon + offset
       END DO
    END DO
  END SUBROUTINE expand_3d

  SUBROUTINE expand_4d(fout,f4d,offset)
      REAL(dp), INTENT(out) :: fout(:)
      REAL(dp), INTENT(in)  :: f4d(:,:,:,:)
      INTEGER,  INTENT(in)  :: offset
      INTEGER :: i, j1, j2, j3, d1, d2, d3, nlon

      nlon = SIZE(f4d,1)
      d1   = SIZE(f4d,2)
      d2   = SIZE(f4d,3)
      d3   = SIZE(f4d,4)
      i = 0
      DO j3=1,d3
         DO j2=1,d2
            DO j1=1,d1
               fout(i+1:i+nlon) = f4d(:,j1,j2,j3)
               i = i + nlon + offset
            END DO
         END DO
      END DO
    END SUBROUTINE expand_4d

    SUBROUTINE compress_2d(finp,f2d,offset)
      REAL(dp), INTENT(in)  :: finp(:)
      REAL(dp), INTENT(out) :: f2d(:,:)
      INTEGER,  INTENT(in)  :: offset
      INTEGER :: i, j1, d1, nlon

      i = 0
      nlon = SIZE(f2d,1)
      d1   = SIZE(f2d,2)
      DO j1=1,d1
         f2d(:,j1) = finp(i+1:i+nlon)
         i = i + nlon + offset
      END DO
    END SUBROUTINE compress_2d

    SUBROUTINE compress_3d(finp,f3d,offset)
      REAL(dp), INTENT(in)  :: finp(:)
      REAL(dp), INTENT(out) :: f3d(:,:,:)
      INTEGER,  INTENT(in)  :: offset
      INTEGER :: i, j1, j2, d1, d2, nlon

      i = 0
      nlon = SIZE(f3d,1)
      d1   = SIZE(f3d,2)
      d2   = SIZE(f3d,3)
      DO j2=1,d2
         DO j1=1,d1
            f3d(:,j1,j2) = finp(i+1:i+nlon)
            i = i + nlon + offset
         END DO
      END DO
    END SUBROUTINE compress_3d

    SUBROUTINE compress_4d(finp,f4d,offset)
      REAL(dp), INTENT(in)  :: finp(:)
      REAL(dp), INTENT(out) :: f4d(:,:,:,:)
      INTEGER,  INTENT(in)  :: offset
      INTEGER :: i, j1, j2, j3, d1, d2, d3, nlon

      i = 0
      nlon = SIZE(f4d,1)
      d1   = SIZE(f4d,2)
      d2   = SIZE(f4d,3)
      d3   = SIZE(f4d,4)
      DO j3=1,d3
         DO j2=1,d2
            DO j1=1,d1
               f4d(:,j1,j2,j3) = finp(i+1:i+nlon)
               i = i + nlon + offset
            END DO
         END DO
      END DO
    END SUBROUTINE compress_4d

END MODULE mo_spectral
