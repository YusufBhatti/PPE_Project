!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_fft992

  USE mo_kind,      ONLY: dp
  USE mo_exception, ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_fft992
  PUBLIC :: cleanup_fft992
  PUBLIC :: fft992

  REAL(dp), PARAMETER :: pi    = 3.14159265358979323846_dp
  REAL(dp), PARAMETER :: sin36 = 0.58778525229247312917_dp
  REAL(dp), PARAMETER :: sin45 = 0.70710678118654752440_dp
  REAL(dp), PARAMETER :: sin60 = 0.86602540378443864676_dp
  REAL(dp), PARAMETER :: sin72 = 0.95105651629515572116_dp
  REAL(dp), PARAMETER :: qrt5  = 0.55901699437494742410_dp

  REAL(dp), PARAMETER :: qqrt5  = 2.0_dp*qrt5
  REAL(dp), PARAMETER :: ssin36 = 2.0_dp*sin36
  REAL(dp), PARAMETER :: ssin45 = 2.0_dp*sin45
  REAL(dp), PARAMETER :: ssin60 = 2.0_dp*sin60
  REAL(dp), PARAMETER :: ssin72 = 2.0_dp*sin72

  INTEGER :: ifax(10)
  REAL(dp), ALLOCATABLE :: trigs(:)

CONTAINS

  SUBROUTINE init_fft992(nlon)

    INTEGER, INTENT(in) :: nlon

    !  Executable statements

    IF (.NOT. ALLOCATED(trigs)) ALLOCATE(trigs(nlon))

    ! Set up trigonometric tables

    CALL set992(nlon)

  END SUBROUTINE init_fft992

  SUBROUTINE cleanup_fft992

    IF (ALLOCATED(trigs)) DEALLOCATE(trigs)

  ENDSUBROUTINE cleanup_fft992

  !===========================================================================
  !
  !     'fft992' - multiple fast REAL periodic transform
  !
  !     author: clive temperton, january 1998
  !
  !     this routine is a modernized and enhanced version of fft991
  !         - cray directives and ancient fortran constructs removed
  !         - "vector chopping" removed
  !         - work array is now dynamically allocated
  !         - stride in work array is now always 1
  !
  !     REAL transform of length n performed by removing redundant
  !     operations from complex transform of length n
  !
  !     a is the array containing input & output data
  !     trigs is a previously prepared list of trig function values
  !     ifax is a previously prepared list of factors of n
  !     inc is the increment within each data 'vector'
  !         (e.g. inc=1 for consecutively stored data)
  !     jump is the increment between the start of each data vector
  !     n is the length of the data vectors
  !     lot is the number of data vectors
  !     isign = +1 for transform from spectral to gridpoint
  !           = -1 for transform from gridpoint to spectral
  !
  !     ordering of coefficients:
  !         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
  !         where b(0)=b(n/2)=0; (n+2) locations required
  !
  !     ordering of data:
  !         x(0),x(1),x(2),...,x(n-1), 0 , 0 ; (n+2) locations required
  !
  !     vectorization is achieved by doing the transforms in parallel
  !
  !     n must be composed of factors 2,3 & 5 but does not have to be even
  !
  !     definition of transforms:
  !     -------------------------
  !
  !     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
  !         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
  !
  !     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
  !               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
  !
  !===========================================================================

  SUBROUTINE fft992(a,inc,jump,n,lot,isign)

    ! .. Scalar Arguments ..
    INTEGER, INTENT(in) :: inc, isign, jump, lot, n
    ! ..
    ! .. Array Arguments ..
    REAL(dp), INTENT(inout) :: a(*)
    ! ..
    ! .. Local Scalars ..
    INTEGER :: i, ia, ibase, ierr, ifac, igo, ii, inca, ix, j, jbase, jj, &
         jumpa, k, la, nfax
    LOGICAL :: lipl
    ! ..
    ! .. Local Arrays ..
    REAL(dp) :: work(n*lot+1)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC mod
    ! ..

    nfax = ifax(1)
    IF (isign==+1) THEN

      !     isign=+1, spectral to gridpoint transform
      !     -----------------------------------------

      i = 1
!DIR$ IVDEP
!CDIR NODEP
      DO j = 1, lot
        a(i+inc) = 0.5_dp*a(i)
        i = i + jump
      END DO
      IF (mod(n,2)==0) THEN
        i = n*inc + 1
!DIR$ IVDEP
!CDIR NODEP
        DO j = 1, lot
          a(i) = 0.5_dp*a(i)
          i = i + jump
        END DO
      END IF

      ia = inc + 1
      la = 1
      igo = + 1

      DO k = 1, nfax
        ifac = ifax(k+1)
        ierr = -1
        IF (k==nfax .AND. nfax>2 .AND. igo==+1) THEN
          lipl = .TRUE.
          ELSE
            lipl = .FALSE.
          END IF
          IF (inc==1 .AND. jump<(2*n) .AND. k>1 .AND. k<(nfax-mod(nfax,2))) THEN
            inca = lot
            jumpa = 1
          ELSE
            inca = inc
            jumpa = jump
          END IF
          IF (igo==+1) THEN
            CALL rpassf(a(ia),a(ia+la*inca),work(1),work(ifac*la*lot+1), &
              inca,lot,jumpa,1,lot,n,ifac,la,ierr,lipl)
          ELSE
            CALL rpassf(work(1),work(la*lot+1),a(ia),a(ia+ifac*la*inca), &
              lot,inca,1,jumpa,lot,n,ifac,la,ierr,lipl)
          END IF
          IF (ierr/=0) THEN
            IF (ierr==2) THEN
              WRITE (message_text,'(a,i3,a)') &
                   'factor = ', ifac, ' not catered for'
              CALL message('fft992',TRIM(message_text))
            END IF
            IF (ierr==3) THEN
              WRITE (message_text,'(a,i3,a)') &
                   'factor = ', ifac, ' only catered for if la*ifac=n'
              CALL message('fft992',TRIM(message_text))
            END IF
            RETURN
          END IF
          la = ifac*la
          igo = -igo
          ia = 1
        END DO

!     if necessary, copy results back to a
!     ------------------------------------
        IF (nfax==1) THEN
          ibase = 1
          jbase = 1
          DO jj = 1, n
            i = ibase
            j = jbase
            DO ii = 1, lot
              a(j) = work(i)
              i = i + 1
              j = j + jump
            END DO
            ibase = ibase + lot
            jbase = jbase + inc
          END DO
        END IF

!     fill in zeros at end
!     --------------------
        ix = n*inc + 1
!DIR$ IVDEP
!CDIR NODEP
        DO j = 1, lot
          a(ix) = 0.0_dp
          a(ix+inc) = 0.0_dp
          ix = ix + jump
        END DO

      ELSE IF (isign==-1) THEN

!     isign=-1, gridpoint to spectral transform
!     -----------------------------------------
        ia = 1
        la = n
        igo = + 1

        DO k = 1, nfax
          ifac = ifax(nfax+2-k)
          la = la/ifac
          ierr = -1
          IF (k==1 .AND. nfax>2 .AND. mod(nfax,2)==1) THEN
            lipl = .TRUE.
          ELSE
            lipl = .FALSE.
          END IF
          IF (inc==1 .AND. jump<(2*n) .AND. k>(1+mod(nfax, &
              2)) .AND. k<nfax) THEN
            inca = lot
            jumpa = 1
          ELSE
            inca = inc
            jumpa = jump
          END IF
          IF (igo==+1) THEN
            CALL qpassf(a(ia),a(ia+ifac*la*inca),work(1),work(la*lot+1), &
              inca,lot,jumpa,1,lot,n,ifac,la,ierr,lipl)
          ELSE
            CALL qpassf(work(1),work(ifac*la*lot+1),a(ia),a(ia+la*inca), &
              lot,inca,1,jumpa,lot,n,ifac,la,ierr,lipl)
          END IF
          IF (ierr/=0) THEN
            IF (ierr==2) THEN
              WRITE (message_text,'(a,i3,a)') &
                   'factor = ', ifac, ' not catered for'
              CALL message('fft992',TRIM(message_text))
            END IF
            IF (ierr==3) THEN
              WRITE (message_text,'(a,i3,a)') &
                   'factor = ', ifac, ' only catered for if la*ifac=n'
              CALL message('fft992',TRIM(message_text))
            END IF
            RETURN
          END IF
          IF (lipl) THEN
            ia = 1
          ELSE
            igo = -igo
            ia = inc + 1
          END IF
        END DO

!     if necessary, copy results back to a
!     ------------------------------------
        IF (nfax==1) THEN
          ibase = 1
          jbase = inc + 1
          DO jj = 1, n
            i = ibase
            j = jbase
            DO ii = 1, lot
              a(j) = work(i)
              i = i + 1
              j = j + jump
            END DO
            ibase = ibase + lot
            jbase = jbase + inc
          END DO
        END IF

!     shift a(0) & fill in zero imag parts
!     ------------------------------------
        ix = 1
!DIR$ IVDEP
!CDIR NODEP
!DIR$ VECTOR ALWAYS
        DO j = 1, lot
          a(ix) = a(ix+inc)
          a(ix+inc) = 0.0_dp
          ix = ix + jump
        END DO
        IF (mod(n,2)==0) THEN
          ix = (n+1)*inc + 1
          DO j = 1, lot
            a(ix) = 0.0_dp
            ix = ix + jump
          END DO
        END IF

      END IF

      RETURN
    END SUBROUTINE fft992

!=============================================================================
!
!     subroutine 'set99' - computes factors of n & trigonometric
!     functions required by fft99 & fft991
!
!=============================================================================
    SUBROUTINE set992(n)

! .. Scalar Arguments ..
      INTEGER, INTENT(in) :: n
! ..
! .. Local Scalars ..
      REAL(dp) :: angle, del
      INTEGER :: i, ifac, il, k, nfax, nhl, nil, nu
! ..
! .. Local Arrays ..
      INTEGER :: jfax(10)
      INTEGER, PARAMETER :: nlfax(7) =(/ 6, 8, 5, 4, 3, 2, 1 /)
! ..
!      del = 4.0E0_dp*asin(1.0E0_dp)/REAL(n,dp)
      del = 2.0_dp*pi/REAL(n,dp)
      nil = 0
      nhl = (n/2) - 1
      DO k = nil, nhl
        angle = REAL(k,dp)*del
        trigs(2*k+1) = cos(angle)
        trigs(2*k+2) = sin(angle)
      ENDDO

!     find factors of n (8,6,5,4,3,2; only one 8 allowed)
!     look for sixes first, store factors in descending order
      nu = n
      ifac = 6
      k = 0
      il = 1
20    CONTINUE
      IF (mod(nu,ifac)/=0) GO TO 30
      k = k + 1
      jfax(k) = ifac
      IF (ifac/=8) GO TO 25
      IF (k==1) GO TO 25
      jfax(1) = 8
      jfax(k) = 6
25    CONTINUE
      nu = nu/ifac
      IF (nu==1) GO TO 50
      IF (ifac/=8) GO TO 20
30    CONTINUE
      il = il + 1
      ifac = nlfax(il)
      IF (ifac>1) GO TO 20

      WRITE (message_text,'(a,i4,a)') &
           'n = ', n, ' contains illegal factors'
      CALL message('set992',TRIM(message_text))
      RETURN

!     now reverse order of factors
50    CONTINUE
      nfax = k
      ifax(1) = nfax
      DO i = 1, nfax
        ifax(nfax+2-i) = jfax(i)
      ENDDO
      ifax(10) = n
      RETURN
    END SUBROUTINE set992

!=============================================================================
!
!     subroutine 'qpassf' - performs one pass through data as part
!     of multiple REAL fft (fourier analysis) routine
!
!     a is first REAL input vector
!         equivalence b(1) with a(ifac*la*inc1+1)
!     c is first REAL output vector
!         equivalence d(1) with c(la*inc2+1)
!     trigs is a precalculated list of sines & cosines
!     inc1 is the addressing increment for a
!     inc2 is the addressing increment for c
!     inc3 is the increment between input vectors a
!     inc4 is the increment between output vectors c
!     lot is the number of vectors
!     n is the length of the vectors
!     ifac is the current factor of n
!     la = n/(product of factors used so far)
!     ierr is an error indicator:
!              0 - pass completed without error
!              1 - lot greater than 64
!              2 - ifac not catered for
!              3 - ifac only catered for if la=n/ifac
!     lipl=.t. => results are returned to input array
!              (only valid if la=n/ifac, i.e. on first pass)
!
!=============================================================================

    SUBROUTINE qpassf(a,b,c,d,inc1,inc2,inc3,inc4,lot,n,ifac,la,ierr, &
        lipl)


! .. Scalar Arguments ..
      INTEGER :: ierr, ifac, inc1, inc2, inc3, inc4, la, lot, n
      LOGICAL :: lipl
! ..
! .. Array Arguments ..
      REAL(dp) :: a(*), b(*), c(*), d(*)
! ..
! .. Local Scalars ..
      REAL(dp) :: a0, a1, a10, a11, a2, a20, a21, a3, a4, a5, a6, b0, b1, &
        b10, b11, b2, b20, b21, b3, b4, b5, b6, c1, c2, c3, c4, c5, s1, &
        s2, s3, s4, s5, t1, t2, t3, t4, t5, t6, &
        t7, z, zqrt5, zsin36, zsin45, zsin60, zsin72
      INTEGER :: i, ia, ib, ibad, ibase, ic, id, ie, if, ig, ih, iink, ijk, &
        ijump, ila, ilot, inc11, j, ja, jb, jbase, jc, jd, je, jf, jink, k, &
        kb, kc, kd, ke, kf, kstop, l, m

! ..
      m = n/ifac
      iink = la*inc1
      jink = la*inc2
      ijump = (ifac-1)*iink
      kstop = (n-ifac)/(2*ifac)

      ibase = 0
      jbase = 0
      ibad = 0

!     increase the vector length by fusing the loops if the
!     data layout is appropriate:
      IF (inc1==lot .AND. inc2==lot .AND. inc3==1 .AND. inc4==1) THEN
        ila = 1
        ilot = la*lot
        inc11 = la*lot
      ELSE
        ila = la
        ilot = lot
        inc11 = inc1
      END IF


      IF (ifac==2) THEN

!     coding for factor 2
!     -------------------
! 200     CONTINUE
        ia = 1
        ib = ia + iink
        ja = 1
        jb = ja + (2*m-la)*inc2

        IF (la/=m) THEN

          DO l = 1, ila
            i = ibase
            j = jbase
!DIR$ IVDEP
!CDIR NODEP
!DIR$ VECTOR ALWAYS
            DO ijk = 1, ilot
              c(ja+j) = a(ia+i) + a(ib+i)
              c(jb+j) = a(ia+i) - a(ib+i)
              i = i + inc3
              j = j + inc4
            END DO
            ibase = ibase + inc11
            jbase = jbase + inc2
          END DO
          ja = ja + jink
          jink = 2*jink
          jb = jb - jink
          ibase = ibase + ijump
          ijump = 2*ijump + iink

          IF (ja<jb) THEN
            DO k = la, kstop, la
              kb = k + k
              c1 = trigs(kb+1)
              s1 = trigs(kb+2)
              jbase = 0
              DO l = 1, ila
                i = ibase
                j = jbase
!DIR$ IVDEP
!CDIR NODEP
                DO ijk = 1, ilot
                  c(ja+j) = a(ia+i) + (c1*a(ib+i)+s1*b(ib+i))
                  c(jb+j) = a(ia+i) - (c1*a(ib+i)+s1*b(ib+i))
                  d(ja+j) = (c1*b(ib+i)-s1*a(ib+i)) + b(ia+i)
                  d(jb+j) = (c1*b(ib+i)-s1*a(ib+i)) - b(ia+i)
                  i = i + inc3
                  j = j + inc4
                END DO
                ibase = ibase + inc11
                jbase = jbase + inc2
              END DO
              ibase = ibase + ijump
              ja = ja + jink
              jb = jb - jink
            END DO
          END IF

          IF (ja==jb) THEN
            jbase = 0
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
!DIR$ VECTOR ALWAYS
              DO ijk = 1, ilot
                c(ja+j) = a(ia+i)
                d(ja+j) = -a(ib+i)
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF

!!! case la=m
        ELSE
          z = 1.0_dp/REAL(n,dp)
          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                t1 = z*(a(ia+i)-a(ib+i))
                a(ia+i) = z*(a(ia+i)+a(ib+i))
                a(ib+i) = t1
                i = i + inc3
              END DO
              ibase = ibase + inc11
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = z*(a(ia+i)+a(ib+i))
                c(jb+j) = z*(a(ia+i)-a(ib+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF
        END IF

      ELSE IF (ifac==3) THEN

!     coding for factor 3
!     -------------------
! 300     CONTINUE
        ia = 1
        ib = ia + iink
        ic = ib + iink
        ja = 1
        jb = ja + (2*m-la)*inc2
        jc = jb

        IF (la/=m) THEN

          DO l = 1, ila
            i = ibase
            j = jbase
!DIR$ IVDEP
!CDIR NODEP
            DO ijk = 1, ilot
              c(ja+j) = a(ia+i) + (a(ib+i)+a(ic+i))
              c(jb+j) = a(ia+i) - 0.5_dp*(a(ib+i)+a(ic+i))
              d(jb+j) = sin60*(a(ic+i)-a(ib+i))
              i = i + inc3
              j = j + inc4
            END DO
            ibase = ibase + inc11
            jbase = jbase + inc2
          END DO
          ja = ja + jink
          jink = 2*jink
          jb = jb + jink
          jc = jc - jink
          ibase = ibase + ijump
          ijump = 2*ijump + iink

          IF (ja<jc) THEN
            DO k = la, kstop, la
              kb = k + k
              kc = kb + kb
              c1 = trigs(kb+1)
              s1 = trigs(kb+2)
              c2 = trigs(kc+1)
              s2 = trigs(kc+2)
              jbase = 0
              DO l = 1, ila
                i = ibase
                j = jbase
!DIR$ IVDEP
!CDIR NODEP
                DO ijk = 1, ilot
                  a1 = (c1*a(ib+i)+s1*b(ib+i)) + (c2*a(ic+i)+s2*b(ic+i))
                  b1 = (c1*b(ib+i)-s1*a(ib+i)) + (c2*b(ic+i)-s2*a(ic+i))
                  a2 = a(ia+i) - 0.5_dp*a1
                  b2 = b(ia+i) - 0.5_dp*b1
                  a3 = sin60*((c1*a(ib+i)+s1*b(ib+i))-(c2*a(ic+i)+s2*b(ic+i)))
                  b3 = sin60*((c1*b(ib+i)-s1*a(ib+i))-(c2*b(ic+i)-s2*a(ic+i)))
                  c(ja+j) = a(ia+i) + a1
                  d(ja+j) = b(ia+i) + b1
                  c(jb+j) = a2 + b3
                  d(jb+j) = b2 - a3
                  c(jc+j) = a2 - b3
                  d(jc+j) = -(b2+a3)
                  i = i + inc3
                  j = j + inc4
                END DO
                ibase = ibase + inc11
                jbase = jbase + inc2
              END DO
              ibase = ibase + ijump
              ja = ja + jink
              jb = jb + jink
              jc = jc - jink
            END DO
          END IF

          IF (ja==jc) THEN
            jbase = 0
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = a(ia+i) + 0.5_dp*(a(ib+i)-a(ic+i))
                d(ja+j) = -sin60*(a(ib+i)+a(ic+i))
                c(jb+j) = a(ia+i) - (a(ib+i)-a(ic+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF

!!! case la=m
        ELSE
          z = 1.0_dp/REAL(n,dp)
          zsin60 = z*sin60
          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                t1 = z*(a(ia+i)-0.5_dp*(a(ib+i)+a(ic+i)))
                t2 = zsin60*(a(ic+i)-a(ib+i))
                a(ia+i) = z*(a(ia+i)+(a(ib+i)+a(ic+i)))
                a(ib+i) = t1
                a(ic+i) = t2
                i = i + inc3
              END DO
              ibase = ibase + inc11
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = z*(a(ia+i)+(a(ib+i)+a(ic+i)))
                c(jb+j) = z*(a(ia+i)-0.5_dp*(a(ib+i)+a(ic+i)))
                d(jb+j) = zsin60*(a(ic+i)-a(ib+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF
        END IF

      ELSE IF (ifac==4) THEN

!     coding for factor 4
!     -------------------
! 400     CONTINUE
        ia = 1
        ib = ia + iink
        ic = ib + iink
        id = ic + iink
        ja = 1
        jb = ja + (2*m-la)*inc2
        jc = jb + 2*m*inc2
        jd = jb

        IF (la/=m) THEN

          DO l = 1, ila
            i = ibase
            j = jbase
!DIR$ IVDEP
!CDIR NODEP
!DIR$ VECTOR ALWAYS
            DO ijk = 1, ilot
              c(ja+j) = (a(ia+i)+a(ic+i)) + (a(ib+i)+a(id+i))
              c(jc+j) = (a(ia+i)+a(ic+i)) - (a(ib+i)+a(id+i))
              c(jb+j) = a(ia+i) - a(ic+i)
              d(jb+j) = a(id+i) - a(ib+i)
              i = i + inc3
              j = j + inc4
            END DO
            ibase = ibase + inc11
            jbase = jbase + inc2
          END DO
          ja = ja + jink
          jink = 2*jink
          jb = jb + jink
          jc = jc - jink
          jd = jd - jink
          ibase = ibase + ijump
          ijump = 2*ijump + iink

          IF (jb<jc) THEN
            DO k = la, kstop, la
              kb = k + k
              kc = kb + kb
              kd = kc + kb
              c1 = trigs(kb+1)
              s1 = trigs(kb+2)
              c2 = trigs(kc+1)
              s2 = trigs(kc+2)
              c3 = trigs(kd+1)
              s3 = trigs(kd+2)
              jbase = 0
              DO l = 1, ila
                i = ibase
                j = jbase
!DIR$ IVDEP
!CDIR NODEP
                DO ijk = 1, ilot
                  a0 = a(ia+i) + (c2*a(ic+i)+s2*b(ic+i))
                  a2 = a(ia+i) - (c2*a(ic+i)+s2*b(ic+i))
                  a1 = (c1*a(ib+i)+s1*b(ib+i)) + (c3*a(id+i)+s3*b(id+i))
                  a3 = (c1*a(ib+i)+s1*b(ib+i)) - (c3*a(id+i)+s3*b(id+i))
                  b0 = b(ia+i) + (c2*b(ic+i)-s2*a(ic+i))
                  b2 = b(ia+i) - (c2*b(ic+i)-s2*a(ic+i))
                  b1 = (c1*b(ib+i)-s1*a(ib+i)) + (c3*b(id+i)-s3*a(id+i))
                  b3 = (c1*b(ib+i)-s1*a(ib+i)) - (c3*b(id+i)-s3*a(id+i))
                  c(ja+j) = a0 + a1
                  c(jc+j) = a0 - a1
                  d(ja+j) = b0 + b1
                  d(jc+j) = b1 - b0
                  c(jb+j) = a2 + b3
                  c(jd+j) = a2 - b3
                  d(jb+j) = b2 - a3
                  d(jd+j) = -(b2+a3)
                  i = i + inc3
                  j = j + inc4
                END DO
                ibase = ibase + inc11
                jbase = jbase + inc2
              END DO
              ibase = ibase + ijump
              ja = ja + jink
              jb = jb + jink
              jc = jc - jink
              jd = jd - jink
            END DO
          END IF

          IF (jb==jc) THEN
            jbase = 0
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = a(ia+i) + sin45*(a(ib+i)-a(id+i))
                c(jb+j) = a(ia+i) - sin45*(a(ib+i)-a(id+i))
                d(ja+j) = -a(ic+i) - sin45*(a(ib+i)+a(id+i))
                d(jb+j) = a(ic+i) - sin45*(a(ib+i)+a(id+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF

!!! case la=m
        ELSE
          z = 1.0_dp/REAL(n,dp)
          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                t1 = z*(a(ia+i)-a(ic+i))
                t3 = z*(a(id+i)-a(ib+i))
                t2 = z*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
                a(ia+i) = z*((a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i)))
                a(ib+i) = t1
                a(ic+i) = t2
                a(id+i) = t3
                i = i + inc3
              END DO
              ibase = ibase + inc11
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = z*((a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i)))
                c(jc+j) = z*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
                c(jb+j) = z*(a(ia+i)-a(ic+i))
                d(jb+j) = z*(a(id+i)-a(ib+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF
        END IF

      ELSE IF (ifac==5) THEN

!     coding for factor 5
!     -------------------
! 500     CONTINUE
        ia = 1
        ib = ia + iink
        ic = ib + iink
        id = ic + iink
        ie = id + iink
        ja = 1
        jb = ja + (2*m-la)*inc2
        jc = jb + 2*m*inc2
        jd = jc
        je = jb

        IF (la/=m) THEN

          DO l = 1, ila
            i = ibase
            j = jbase
!DIR$ IVDEP
!CDIR NODEP
!DIR$ VECTOR ALWAYS
            DO ijk = 1, ilot
              a1 = a(ib+i) + a(ie+i)
              a3 = a(ib+i) - a(ie+i)
              a2 = a(ic+i) + a(id+i)
              a4 = a(ic+i) - a(id+i)
              a5 = a(ia+i) - 0.25_dp*(a1+a2)
              a6 = qrt5*(a1-a2)
              c(ja+j) = a(ia+i) + (a1+a2)
              c(jb+j) = a5 + a6
              c(jc+j) = a5 - a6
              d(jb+j) = -sin72*a3 - sin36*a4
              d(jc+j) = -sin36*a3 + sin72*a4
              i = i + inc3
              j = j + inc4
            END DO
            ibase = ibase + inc11
            jbase = jbase + inc2
          END DO
          ja = ja + jink
          jink = 2*jink
          jb = jb + jink
          jc = jc + jink
          jd = jd - jink
          je = je - jink
          ibase = ibase + ijump
          ijump = 2*ijump + iink

          IF (jb<jd) THEN
            DO k = la, kstop, la
              kb = k + k
              kc = kb + kb
              kd = kc + kb
              ke = kd + kb
              c1 = trigs(kb+1)
              s1 = trigs(kb+2)
              c2 = trigs(kc+1)
              s2 = trigs(kc+2)
              c3 = trigs(kd+1)
              s3 = trigs(kd+2)
              c4 = trigs(ke+1)
              s4 = trigs(ke+2)
              jbase = 0
              DO l = 1, ila
                i = ibase
                j = jbase
!DIR$ IVDEP
!CDIR NODEP
                DO ijk = 1, ilot
                  a1 = (c1*a(ib+i)+s1*b(ib+i)) + (c4*a(ie+i)+s4*b(ie+i))
                  a3 = (c1*a(ib+i)+s1*b(ib+i)) - (c4*a(ie+i)+s4*b(ie+i))
                  a2 = (c2*a(ic+i)+s2*b(ic+i)) + (c3*a(id+i)+s3*b(id+i))
                  a4 = (c2*a(ic+i)+s2*b(ic+i)) - (c3*a(id+i)+s3*b(id+i))
                  b1 = (c1*b(ib+i)-s1*a(ib+i)) + (c4*b(ie+i)-s4*a(ie+i))
                  b3 = (c1*b(ib+i)-s1*a(ib+i)) - (c4*b(ie+i)-s4*a(ie+i))
                  b2 = (c2*b(ic+i)-s2*a(ic+i)) + (c3*b(id+i)-s3*a(id+i))
                  b4 = (c2*b(ic+i)-s2*a(ic+i)) - (c3*b(id+i)-s3*a(id+i))
                  a5 = a(ia+i) - 0.25_dp*(a1+a2)
                  a6 = qrt5*(a1-a2)
                  b5 = b(ia+i) - 0.25_dp*(b1+b2)
                  b6 = qrt5*(b1-b2)
                  a10 = a5 + a6
                  a20 = a5 - a6
                  b10 = b5 + b6
                  b20 = b5 - b6
                  a11 = sin72*b3 + sin36*b4
                  a21 = sin36*b3 - sin72*b4
                  b11 = sin72*a3 + sin36*a4
                  b21 = sin36*a3 - sin72*a4
                  c(ja+j) = a(ia+i) + (a1+a2)
                  c(jb+j) = a10 + a11
                  c(je+j) = a10 - a11
                  c(jc+j) = a20 + a21
                  c(jd+j) = a20 - a21
                  d(ja+j) = b(ia+i) + (b1+b2)
                  d(jb+j) = b10 - b11
                  d(je+j) = -(b10+b11)
                  d(jc+j) = b20 - b21
                  d(jd+j) = -(b20+b21)
                  i = i + inc3
                  j = j + inc4
                END DO
                ibase = ibase + inc11
                jbase = jbase + inc2
              END DO
              ibase = ibase + ijump
              ja = ja + jink
              jb = jb + jink
              jc = jc + jink
              jd = jd - jink
              je = je - jink
            END DO
          END IF

          IF (jb==jd) THEN
            jbase = 0
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                a1 = a(ib+i) + a(ie+i)
                a3 = a(ib+i) - a(ie+i)
                a2 = a(ic+i) + a(id+i)
                a4 = a(ic+i) - a(id+i)
                a5 = a(ia+i) + 0.25_dp*(a3-a4)
                a6 = qrt5*(a3+a4)
                c(ja+j) = a5 + a6
                c(jb+j) = a5 - a6
                c(jc+j) = a(ia+i) - (a3-a4)
                d(ja+j) = -sin36*a1 - sin72*a2
                d(jb+j) = -sin72*a1 + sin36*a2
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF

!!! case la=m
        ELSE
          z = 1.0_dp/REAL(n,dp)
          zqrt5 = z*qrt5
          zsin36 = z*sin36
          zsin72 = z*sin72
          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                a1 = a(ib+i) + a(ie+i)
                a3 = a(ib+i) - a(ie+i)
                a2 = a(ic+i) + a(id+i)
                a4 = a(ic+i) - a(id+i)
                a5 = z*(a(ia+i)-0.25_dp*(a1+a2))
                a6 = zqrt5*(a1-a2)
                a(ia+i) = z*(a(ia+i)+(a1+a2))
                a(ib+i) = a5 + a6
                a(id+i) = a5 - a6
                a(ic+i) = -zsin72*a3 - zsin36*a4
                a(ie+i) = -zsin36*a3 + zsin72*a4
                i = i + inc3
              END DO
              ibase = ibase + inc11
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                a1 = a(ib+i) + a(ie+i)
                a3 = a(ib+i) - a(ie+i)
                a2 = a(ic+i) + a(id+i)
                a4 = a(ic+i) - a(id+i)
                a5 = z*(a(ia+i)-0.25_dp*(a1+a2))
                a6 = zqrt5*(a1-a2)
                c(ja+j) = z*(a(ia+i)+(a1+a2))
                c(jb+j) = a5 + a6
                c(jc+j) = a5 - a6
                d(jb+j) = -zsin72*a3 - zsin36*a4
                d(jc+j) = -zsin36*a3 + zsin72*a4
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF
        END IF

      ELSE IF (ifac==6) THEN

!     coding for factor 6
!     -------------------
! 600     CONTINUE
        ia = 1
        ib = ia + iink
        ic = ib + iink
        id = ic + iink
        ie = id + iink
        if = ie + iink
        ja = 1
        jb = ja + (2*m-la)*inc2
        jc = jb + 2*m*inc2
        jd = jc + 2*m*inc2
        je = jc
        jf = jb

        IF (la/=m) THEN

          DO l = 1, ila
            i = ibase
            j = jbase
!DIR$ IVDEP
!CDIR NODEP
            DO ijk = 1, ilot
              a11 = (a(ic+i)+a(if+i)) + (a(ib+i)+a(ie+i))
              c(ja+j) = (a(ia+i)+a(id+i)) + a11
              c(jc+j) = (a(ia+i)+a(id+i)-0.5_dp*a11)
              d(jc+j) = sin60*((a(ic+i)+a(if+i))-(a(ib+i)+a(ie+i)))
              a11 = (a(ic+i)-a(if+i)) + (a(ie+i)-a(ib+i))
              c(jb+j) = (a(ia+i)-a(id+i)) - 0.5_dp*a11
              d(jb+j) = sin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(if+i)))
              c(jd+j) = (a(ia+i)-a(id+i)) + a11
              i = i + inc3
              j = j + inc4
            END DO
            ibase = ibase + inc11
            jbase = jbase + inc2
          END DO
          ja = ja + jink
          jink = 2*jink
          jb = jb + jink
          jc = jc + jink
          jd = jd - jink
          je = je - jink
          jf = jf - jink
          ibase = ibase + ijump
          ijump = 2*ijump + iink

          IF (jc<jd) THEN
            DO k = la, kstop, la
              kb = k + k
              kc = kb + kb
              kd = kc + kb
              ke = kd + kb
              kf = ke + kb
              c1 = trigs(kb+1)
              s1 = trigs(kb+2)
              c2 = trigs(kc+1)
              s2 = trigs(kc+2)
              c3 = trigs(kd+1)
              s3 = trigs(kd+2)
              c4 = trigs(ke+1)
              s4 = trigs(ke+2)
              c5 = trigs(kf+1)
              s5 = trigs(kf+2)
              jbase = 0
              DO l = 1, ila
                i = ibase
                j = jbase
!DIR$ IVDEP
!CDIR NODEP
                DO ijk = 1, ilot
                  a1 = c1*a(ib+i) + s1*b(ib+i)
                  b1 = c1*b(ib+i) - s1*a(ib+i)
                  a2 = c2*a(ic+i) + s2*b(ic+i)
                  b2 = c2*b(ic+i) - s2*a(ic+i)
                  a3 = c3*a(id+i) + s3*b(id+i)
                  b3 = c3*b(id+i) - s3*a(id+i)
                  a4 = c4*a(ie+i) + s4*b(ie+i)
                  b4 = c4*b(ie+i) - s4*a(ie+i)
                  a5 = c5*a(if+i) + s5*b(if+i)
                  b5 = c5*b(if+i) - s5*a(if+i)
                  a11 = (a2+a5) + (a1+a4)
                  a20 = (a(ia+i)+a3) - 0.5_dp*a11
                  a21 = sin60*((a2+a5)-(a1+a4))
                  b11 = (b2+b5) + (b1+b4)
                  b20 = (b(ia+i)+b3) - 0.5_dp*b11
                  b21 = sin60*((b2+b5)-(b1+b4))
                  c(ja+j) = (a(ia+i)+a3) + a11
                  d(ja+j) = (b(ia+i)+b3) + b11
                  c(jc+j) = a20 - b21
                  d(jc+j) = a21 + b20
                  c(je+j) = a20 + b21
                  d(je+j) = a21 - b20
                  a11 = (a2-a5) + (a4-a1)
                  a20 = (a(ia+i)-a3) - 0.5_dp*a11
                  a21 = sin60*((a4-a1)-(a2-a5))
                  b11 = (b5-b2) - (b4-b1)
                  b20 = (b3-b(ia+i)) - 0.5_dp*b11
                  b21 = sin60*((b5-b2)+(b4-b1))
                  c(jb+j) = a20 - b21
                  d(jb+j) = a21 - b20
                  c(jd+j) = a11 + (a(ia+i)-a3)
                  d(jd+j) = b11 + (b3-b(ia+i))
                  c(jf+j) = a20 + b21
                  d(jf+j) = a21 + b20
                  i = i + inc3
                  j = j + inc4
                END DO
                ibase = ibase + inc11
                jbase = jbase + inc2
              END DO
              ibase = ibase + ijump
              ja = ja + jink
              jb = jb + jink
              jc = jc + jink
              jd = jd - jink
              je = je - jink
              jf = jf - jink
            END DO
          END IF

          IF (jc==jd) THEN
            jbase = 0
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = (a(ia+i)+0.5_dp*(a(ic+i)-a(ie+i))) + &
                  sin60*(a(ib+i)-a(if+i))
                d(ja+j) = -(a(id+i)+0.5_dp*(a(ib+i)+a(if+i))) - &
                  sin60*(a(ic+i)+a(ie+i))
                c(jb+j) = a(ia+i) - (a(ic+i)-a(ie+i))
                d(jb+j) = a(id+i) - (a(ib+i)+a(if+i))
                c(jc+j) = (a(ia+i)+0.5_dp*(a(ic+i)-a(ie+i))) - &
                  sin60*(a(ib+i)-a(if+i))
                d(jc+j) = -(a(id+i)+0.5_dp*(a(ib+i)+a(if+i))) + &
                  sin60*(a(ic+i)+a(ie+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF

!!! case la=m
        ELSE
          z = 1.0_dp/REAL(n,dp)
          zsin60 = z*sin60
          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                a11 = (a(ic+i)-a(if+i)) + (a(ie+i)-a(ib+i))
                t1 = z*((a(ia+i)-a(id+i))-0.5_dp*a11)
                t5 = z*((a(ia+i)-a(id+i))+a11)
                t2 = zsin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(if+i)))
                t4 = zsin60*((a(ic+i)+a(if+i))-(a(ib+i)+a(ie+i)))
                a11 = (a(ic+i)+a(if+i)) + (a(ib+i)+a(ie+i))
                t3 = z*((a(ia+i)+a(id+i))-0.5_dp*a11)
                a(ia+i) = z*((a(ia+i)+a(id+i))+a11)
                a(ib+i) = t1
                a(ic+i) = t2
                a(id+i) = t3
                a(ie+i) = t4
                a(if+i) = t5
                i = i + inc3
              END DO
              ibase = ibase + inc11
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                a11 = (a(ic+i)+a(if+i)) + (a(ib+i)+a(ie+i))
                c(ja+j) = z*((a(ia+i)+a(id+i))+a11)
                c(jc+j) = z*((a(ia+i)+a(id+i))-0.5_dp*a11)
                d(jc+j) = zsin60*((a(ic+i)+a(if+i))-(a(ib+i)+a(ie+i)))
                a11 = (a(ic+i)-a(if+i)) + (a(ie+i)-a(ib+i))
                c(jb+j) = z*((a(ia+i)-a(id+i))-0.5_dp*a11)
                d(jb+j) = zsin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(if+i)))
                c(jd+j) = z*((a(ia+i)-a(id+i))+a11)
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF
        END IF

      ELSE IF (ifac==8) THEN

!     coding for factor 8
!     -------------------
! 800     CONTINUE
        IF (la/=m) THEN
          ibad = 3
        ELSE
          ia = 1
          ib = ia + iink
          ic = ib + iink
          id = ic + iink
          ie = id + iink
          if = ie + iink
          ig = if + iink
          ih = ig + iink
          ja = 1
          jb = ja + la*inc2
          jc = jb + 2*m*inc2
          jd = jc + 2*m*inc2
          je = jd + 2*m*inc2
          z = 1.0_dp/REAL(n,dp)
          zsin45 = z*sin45

          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                t3 = z*((a(ia+i)+a(ie+i))-(a(ic+i)+a(ig+i)))
                t4 = z*((a(id+i)+a(ih+i))-(a(ib+i)+a(if+i)))
                t1 = z*(a(ia+i)-a(ie+i)) + zsin45*((a(ih+i)-a(id+i))-(a(if+i)-a(ib+i)))
                t5 = z*(a(ia+i)-a(ie+i)) - zsin45*((a(ih+i)-a(id+i))-(a(if+i)-a(ib+i)))
                t2 = zsin45*((a(ih+i)-a(id+i))+(a(if+i)-a(ib+i))) + z*(a(ig+i)-a(ic+i))
                t6 = zsin45*((a(ih+i)-a(id+i))+(a(if+i)-a(ib+i))) - z*(a(ig+i)-a(ic+i))
                t7 = z*(((a(ia+i)+a(ie+i))+(a(ic+i)+a(ig+i)))-((a(id+i)+ a(ih+i))+(a(ib+i)+a(if+i))))
                a(ia+i) = z*(((a(ia+i)+a(ie+i))+(a(ic+i)+a(ig+i)))+((a(id+i)+ a(ih+i))+(a(ib+i)+a(if+i))))
                a(ib+i) = t1
                a(ic+i) = t2
                a(id+i) = t3
                a(ie+i) = t4
                a(if+i) = t5
                a(ig+i) = t6
                a(ih+i) = t7
                i = i + inc3
              END DO
              ibase = ibase + inc11
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = z*(((a(ia+i)+a(ie+i))+(a(ic+i)+a(ig+i)))+((a(id+i)+ &
                  a(ih+i))+(a(ib+i)+a(if+i))))
                c(je+j) = z*(((a(ia+i)+a(ie+i))+(a(ic+i)+a(ig+i)))-((a(id+i)+ &
                  a(ih+i))+(a(ib+i)+a(if+i))))
                c(jc+j) = z*((a(ia+i)+a(ie+i))-(a(ic+i)+a(ig+i)))
                d(jc+j) = z*((a(id+i)+a(ih+i))-(a(ib+i)+a(if+i)))
                c(jb+j) = z*(a(ia+i)-a(ie+i)) + zsin45*((a(ih+i)-a(id+i))-(a(if+i)-a(ib+i)))
                c(jd+j) = z*(a(ia+i)-a(ie+i)) - zsin45*((a(ih+i)-a(id+i))-(a(if+i)-a(ib+i)))
                d(jb+j) = zsin45*((a(ih+i)-a(id+i))+(a(if+i)-a(ib+i))) + z*(a(ig+i)-a(ic+i))
                d(jd+j) = zsin45*((a(ih+i)-a(id+i))+(a(if+i)-a(ib+i))) - z*(a(ig+i)-a(ic+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc11
              jbase = jbase + inc2
            END DO
          END IF

        END IF

      ELSE


        ibad = 2
!!! illegal factor
      END IF

!     return
!     ------
! 900   CONTINUE
      ierr = ibad
      RETURN
    END SUBROUTINE qpassf
!     subroutine 'rpassf' - performs one pass through data as part
!     of multiple REAL fft (fourier synthesis) routine

!     a is first REAL input vector
!         equivalence b(1) with a (la*inc1+1)
!     c is first REAL output vector
!         equivalence d(1) with c(ifac*la*inc2+1)
!     trigs is a precalculated list of sines & cosines
!     inc1 is the addressing increment for a
!     inc2 is the addressing increment for c
!     inc3 is the increment between input vectors a
!     inc4 is the increment between output vectors c
!     lot is the number of vectors
!     n is the length of the vectors
!     ifac is the current factor of n
!     la is the product of previous factors
!     ierr is an error indicator:
!              0 - pass completed without error
!              1 - lot greater than 64
!              2 - ifac not catered for
!              3 - ifac only catered for if la=n/ifac
!     lipl=.t. => results are returned to input array
!              (only valid if la=n/ifac, i.e. on last pass)

!-----------------------------------------------------------------------

    SUBROUTINE rpassf(a,b,c,d,inc1,inc2,inc3,inc4,lot,n,ifac,la,ierr, &
        lipl)

! .. Scalar Arguments ..
      INTEGER :: ierr, ifac, inc1, inc2, inc3, inc4, la, lot, n
      LOGICAL :: lipl
! ..
! .. Array Arguments ..
      REAL(dp) :: a(*), b(*), c(*), d(*)
! ..
! .. Local Scalars ..
      REAL(dp) :: a10, a11, a20, a21, b10, b11, b20, b21, c1, c2, c3, c4, c5, &
        s1, s2, s3, s4, s5, t1, t2, t3, t4, t5, t6, t7
      INTEGER :: i, ia, ib, ibad, ibase, ic, id, ie, if, iink, ijk, ila, ilot, &
        inc21, j, ja, jb, jbase, jc, jd, je, jf, jg, jh, jink, jump, k, kb, &
        kc, kd, ke, kf, kstop, l, m

! ..
      m = n/ifac
      iink = la*inc1
      jink = la*inc2
      jump = (ifac-1)*jink
      kstop = (n-ifac)/(2*ifac)

      ibase = 0
      jbase = 0
      ibad = 0

!     increase the vector length by fusing the loops if the
!     data layout is appropriate:
      IF (inc1==lot .AND. inc2==lot .AND. inc3==1 .AND. inc4==1) THEN
        ila = 1
        ilot = la*lot
        inc21 = la*lot
      ELSE
        ila = la
        ilot = lot
        inc21 = inc2
      END IF

      IF (ifac==2) THEN

!     coding for factor 2
!     -------------------

        ia = 1
        ib = ia + (2*m-la)*inc1
        ja = 1
        jb = ja + jink

        IF (la/=m) THEN

          DO l = 1, ila
            i = ibase
            j = jbase
!DIR$ IVDEP
!CDIR NODEP
!DIR$ VECTOR ALWAYS
            DO ijk = 1, ilot
              c(ja+j) = a(ia+i) + a(ib+i)
              c(jb+j) = a(ia+i) - a(ib+i)
              i = i + inc3
              j = j + inc4
            END DO
            ibase = ibase + inc1
            jbase = jbase + inc21
          END DO
          ia = ia + iink
          iink = 2*iink
          ib = ib - iink
          ibase = 0
          jbase = jbase + jump
          jump = 2*jump + jink

          IF (ia<ib) THEN
            DO k = la, kstop, la
              kb = k + k
              c1 = trigs(kb+1)
              s1 = trigs(kb+2)
              ibase = 0
              DO l = 1, ila
                i = ibase
                j = jbase
!DIR$ IVDEP
!CDIR NODEP
                DO ijk = 1, ilot
                  c(ja+j) = a(ia+i) + a(ib+i)
                  d(ja+j) = b(ia+i) - b(ib+i)
                  c(jb+j) = c1*(a(ia+i)-a(ib+i)) - s1*(b(ia+i)+b(ib+i))
                  d(jb+j) = s1*(a(ia+i)-a(ib+i)) + c1*(b(ia+i)+b(ib+i))
                  i = i + inc3
                  j = j + inc4
                END DO
                ibase = ibase + inc1
                jbase = jbase + inc21
              END DO
              ia = ia + iink
              ib = ib - iink
              jbase = jbase + jump
            END DO
          END IF

          IF (ia==ib) THEN
            ibase = 0
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
!DIR$ VECTOR ALWAYS
              DO ijk = 1, ilot
                c(ja+j) = a(ia+i)
                c(jb+j) = -b(ia+i)
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF

!!! case la=m
        ELSE
          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                t1 = 2.0_dp*(a(ia+i)-a(ib+i))
                a(ia+i) = 2.0_dp*(a(ia+i)+a(ib+i))
                a(ib+i) = t1
                i = i + inc3
              END DO
              ibase = ibase + inc1
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = 2.0_dp*(a(ia+i)+a(ib+i))
                c(jb+j) = 2.0_dp*(a(ia+i)-a(ib+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF
        END IF

      ELSE IF (ifac==3) THEN

!     coding for factor 3
!     -------------------

        ia = 1
        ib = ia + (2*m-la)*inc1
        ic = ib
        ja = 1
        jb = ja + jink
        jc = jb + jink

        IF (la/=m) THEN

          DO l = 1, ila
            i = ibase
            j = jbase
!DIR$ IVDEP
!CDIR NODEP
            DO ijk = 1, ilot
              c(ja+j) = a(ia+i) + a(ib+i)
              c(jb+j) = (a(ia+i)-0.5_dp*a(ib+i)) - (sin60*(b(ib+i)))
              c(jc+j) = (a(ia+i)-0.5_dp*a(ib+i)) + (sin60*(b(ib+i)))
              i = i + inc3
              j = j + inc4
            END DO
            ibase = ibase + inc1
            jbase = jbase + inc21
          END DO
          ia = ia + iink
          iink = 2*iink
          ib = ib + iink
          ic = ic - iink
          jbase = jbase + jump
          jump = 2*jump + jink

          IF (ia<ic) THEN
            DO k = la, kstop, la
              kb = k + k
              kc = kb + kb
              c1 = trigs(kb+1)
              s1 = trigs(kb+2)
              c2 = trigs(kc+1)
              s2 = trigs(kc+2)
              ibase = 0
              DO l = 1, ila
                i = ibase
                j = jbase
!DIR$ IVDEP
!CDIR NODEP
                DO ijk = 1, ilot
                  c(ja+j) = a(ia+i) + (a(ib+i)+a(ic+i))
                  d(ja+j) = b(ia+i) + (b(ib+i)-b(ic+i))
                  c(jb+j) = c1*((a(ia+i)-0.5_dp*(a(ib+i)+ &
                    a(ic+i)))-(sin60*(b(ib+i)+b(ic+i)))) - s1*((b(ia+ &
                    i)-0.5_dp*(b(ib+i)-b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
                  d(jb+j) = s1*((a(ia+i)-0.5_dp*(a(ib+i)+ &
                    a(ic+i)))-(sin60*(b(ib+i)+b(ic+i)))) + c1*((b(ia+ &
                    i)-0.5_dp*(b(ib+i)-b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
                  c(jc+j) = c2*((a(ia+i)-0.5_dp*(a(ib+i)+ &
                    a(ic+i)))+(sin60*(b(ib+i)+b(ic+i)))) - s2*((b(ia+ &
                    i)-0.5_dp*(b(ib+i)-b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
                  d(jc+j) = s2*((a(ia+i)-0.5_dp*(a(ib+i)+ &
                    a(ic+i)))+(sin60*(b(ib+i)+b(ic+i)))) + c2*((b(ia+ &
                    i)-0.5_dp*(b(ib+i)-b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
                  i = i + inc3
                  j = j + inc4
                END DO
                ibase = ibase + inc1
                jbase = jbase + inc21
              END DO
              ia = ia + iink
              ib = ib + iink
              ic = ic - iink
              jbase = jbase + jump
            END DO
          END IF

          IF (ia==ic) THEN
            ibase = 0
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = a(ia+i) + a(ib+i)
                c(jb+j) = (0.5_dp*a(ia+i)-a(ib+i)) - (sin60*b(ia+i))
                c(jc+j) = -(0.5_dp*a(ia+i)-a(ib+i)) - (sin60*b(ia+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF

!!! case la=m
        ELSE
          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                t1 = (2.0_dp*a(ia+i)-a(ib+i)) - (ssin60*b(ib+i))
                t2 = (2.0_dp*a(ia+i)-a(ib+i)) + (ssin60*b(ib+i))
                a(ia+i) = 2.0_dp*(a(ia+i)+a(ib+i))
                a(ib+i) = t1
                b(ib+i) = t2
                i = i + inc3
              END DO
              ibase = ibase + inc1
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = 2.0_dp*(a(ia+i)+a(ib+i))
                c(jb+j) = (2.0_dp*a(ia+i)-a(ib+i)) - (ssin60*b(ib+i))
                c(jc+j) = (2.0_dp*a(ia+i)-a(ib+i)) + (ssin60*b(ib+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF
        END IF

      ELSE IF (ifac==4) THEN

!     coding for factor 4
!     -------------------

        ia = 1
        ib = ia + (2*m-la)*inc1
        ic = ib + 2*m*inc1
        id = ib
        ja = 1
        jb = ja + jink
        jc = jb + jink
        jd = jc + jink

        IF (la/=m) THEN

          DO l = 1, ila
            i = ibase
            j = jbase
!DIR$ IVDEP
!CDIR NODEP
!DIR$ VECTOR ALWAYS
            DO ijk = 1, ilot
              c(ja+j) = (a(ia+i)+a(ic+i)) + a(ib+i)
              c(jb+j) = (a(ia+i)-a(ic+i)) - b(ib+i)
              c(jc+j) = (a(ia+i)+a(ic+i)) - a(ib+i)
              c(jd+j) = (a(ia+i)-a(ic+i)) + b(ib+i)
              i = i + inc3
              j = j + inc4
            END DO
            ibase = ibase + inc1
            jbase = jbase + inc21
          END DO
          ia = ia + iink
          iink = 2*iink
          ib = ib + iink
          ic = ic - iink
          id = id - iink
          jbase = jbase + jump
          jump = 2*jump + jink

          IF (ib<ic) THEN
            DO k = la, kstop, la
              kb = k + k
              kc = kb + kb
              kd = kc + kb
              c1 = trigs(kb+1)
              s1 = trigs(kb+2)
              c2 = trigs(kc+1)
              s2 = trigs(kc+2)
              c3 = trigs(kd+1)
              s3 = trigs(kd+2)
              ibase = 0
              DO l = 1, ila
                i = ibase
                j = jbase
!DIR$ IVDEP
!CDIR NODEP
                DO ijk = 1, ilot
                  c(ja+j) = (a(ia+i)+a(ic+i)) + (a(ib+i)+a(id+i))
                  d(ja+j) = (b(ia+i)-b(ic+i)) + (b(ib+i)-b(id+i))
                  c(jc+j) = c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) - &
                    s2*((b(ia+i)-b(ic+i))-(b(ib+i)-b(id+i)))
                  d(jc+j) = s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) + &
                    c2*((b(ia+i)-b(ic+i))-(b(ib+i)-b(id+i)))
                  c(jb+j) = c1*((a(ia+i)-a(ic+i))-(b(ib+i)+b(id+i))) - &
                    s1*((b(ia+i)+b(ic+i))+(a(ib+i)-a(id+i)))
                  d(jb+j) = s1*((a(ia+i)-a(ic+i))-(b(ib+i)+b(id+i))) + &
                    c1*((b(ia+i)+b(ic+i))+(a(ib+i)-a(id+i)))
                  c(jd+j) = c3*((a(ia+i)-a(ic+i))+(b(ib+i)+b(id+i))) - &
                    s3*((b(ia+i)+b(ic+i))-(a(ib+i)-a(id+i)))
                  d(jd+j) = s3*((a(ia+i)-a(ic+i))+(b(ib+i)+b(id+i))) + &
                    c3*((b(ia+i)+b(ic+i))-(a(ib+i)-a(id+i)))
                  i = i + inc3
                  j = j + inc4
                END DO
                ibase = ibase + inc1
                jbase = jbase + inc21
              END DO
              ia = ia + iink
              ib = ib + iink
              ic = ic - iink
              id = id - iink
              jbase = jbase + jump
            END DO
          END IF

          IF (ib==ic) THEN
            ibase = 0
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = a(ia+i) + a(ib+i)
                c(jb+j) = sin45*((a(ia+i)-a(ib+i))-(b(ia+i)+b(ib+i)))
                c(jc+j) = b(ib+i) - b(ia+i)
                c(jd+j) = -sin45*((a(ia+i)-a(ib+i))+(b(ia+i)+b(ib+i)))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF

!!! case la=m
        ELSE
          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                t1 = 2.0_dp*((a(ia+i)-a(ic+i))-b(ib+i))
                t2 = 2.0_dp*((a(ia+i)+a(ic+i))-a(ib+i))
                t3 = 2.0_dp*((a(ia+i)-a(ic+i))+b(ib+i))
                a(ia+i) = 2.0_dp*((a(ia+i)+a(ic+i))+a(ib+i))
                a(ib+i) = t1
                b(ib+i) = t2
                a(ic+i) = t3
                i = i + inc3
              END DO
              ibase = ibase + inc1
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = 2.0_dp*((a(ia+i)+a(ic+i))+a(ib+i))
                c(jb+j) = 2.0_dp*((a(ia+i)-a(ic+i))-b(ib+i))
                c(jc+j) = 2.0_dp*((a(ia+i)+a(ic+i))-a(ib+i))
                c(jd+j) = 2.0_dp*((a(ia+i)-a(ic+i))+b(ib+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF
        END IF

      ELSE IF (ifac==5) THEN

!     coding for factor 5
!     -------------------

        ia = 1
        ib = ia + (2*m-la)*inc1
        ic = ib + 2*m*inc1
        id = ic
        ie = ib
        ja = 1
        jb = ja + jink
        jc = jb + jink
        jd = jc + jink
        je = jd + jink

        IF (la/=m) THEN

          DO l = 1, ila
            i = ibase
            j = jbase
!DIR$ IVDEP
!CDIR NODEP
            DO ijk = 1, ilot
              c(ja+j) = a(ia+i) + (a(ib+i)+a(ic+i))
              c(jb+j) = ((a(ia+i)-0.25_dp*(a(ib+i)+a(ic+i)))+qrt5*(a(ib+ &
                i)-a(ic+i))) - (sin72*b(ib+i)+sin36*b(ic+i))
              c(jc+j) = ((a(ia+i)-0.25_dp*(a(ib+i)+a(ic+i)))-qrt5*(a(ib+ &
                i)-a(ic+i))) - (sin36*b(ib+i)-sin72*b(ic+i))
              c(jd+j) = ((a(ia+i)-0.25_dp*(a(ib+i)+a(ic+i)))-qrt5*(a(ib+ &
                i)-a(ic+i))) + (sin36*b(ib+i)-sin72*b(ic+i))
              c(je+j) = ((a(ia+i)-0.25_dp*(a(ib+i)+a(ic+i)))+qrt5*(a(ib+ &
                i)-a(ic+i))) + (sin72*b(ib+i)+sin36*b(ic+i))
              i = i + inc3
              j = j + inc4
            END DO
            ibase = ibase + inc1
            jbase = jbase + inc21
          END DO
          ia = ia + iink
          iink = 2*iink
          ib = ib + iink
          ic = ic + iink
          id = id - iink
          ie = ie - iink
          jbase = jbase + jump
          jump = 2*jump + jink

          IF (ib<id) THEN
            DO k = la, kstop, la
              kb = k + k
              kc = kb + kb
              kd = kc + kb
              ke = kd + kb
              c1 = trigs(kb+1)
              s1 = trigs(kb+2)
              c2 = trigs(kc+1)
              s2 = trigs(kc+2)
              c3 = trigs(kd+1)
              s3 = trigs(kd+2)
              c4 = trigs(ke+1)
              s4 = trigs(ke+2)
              ibase = 0
              DO l = 1, ila
                i = ibase
                j = jbase
!DIR$ IVDEP
!CDIR NODEP
                DO ijk = 1, ilot

                  a10 = (a(ia+i)-0.25_dp*((a(ib+i)+a(ie+i))+(a(ic+i)+ &
                    a(id+i)))) + qrt5*((a(ib+i)+a(ie+i))-(a(ic+i)+a(id+i)))
                  a20 = (a(ia+i)-0.25_dp*((a(ib+i)+a(ie+i))+(a(ic+i)+ &
                    a(id+i)))) - qrt5*((a(ib+i)+a(ie+i))-(a(ic+i)+a(id+i)))
                  b10 = (b(ia+i)-0.25_dp*((b(ib+i)-b(ie+i))+(b(ic+i)- &
                    b(id+i)))) + qrt5*((b(ib+i)-b(ie+i))-(b(ic+i)-b(id+i)))
                  b20 = (b(ia+i)-0.25_dp*((b(ib+i)-b(ie+i))+(b(ic+i)- &
                    b(id+i)))) - qrt5*((b(ib+i)-b(ie+i))-(b(ic+i)-b(id+i)))
                  a11 = sin72*(b(ib+i)+b(ie+i)) + sin36*(b(ic+i)+b(id+i))
                  a21 = sin36*(b(ib+i)+b(ie+i)) - sin72*(b(ic+i)+b(id+i))
                  b11 = sin72*(a(ib+i)-a(ie+i)) + sin36*(a(ic+i)-a(id+i))
                  b21 = sin36*(a(ib+i)-a(ie+i)) - sin72*(a(ic+i)-a(id+i))

                  c(ja+j) = a(ia+i) + ((a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i)))
                  d(ja+j) = b(ia+i) + ((b(ib+i)-b(ie+i))+(b(ic+i)-b(id+i)))
                  c(jb+j) = c1*(a10-a11) - s1*(b10+b11)
                  d(jb+j) = s1*(a10-a11) + c1*(b10+b11)
                  c(je+j) = c4*(a10+a11) - s4*(b10-b11)
                  d(je+j) = s4*(a10+a11) + c4*(b10-b11)
                  c(jc+j) = c2*(a20-a21) - s2*(b20+b21)
                  d(jc+j) = s2*(a20-a21) + c2*(b20+b21)
                  c(jd+j) = c3*(a20+a21) - s3*(b20-b21)
                  d(jd+j) = s3*(a20+a21) + c3*(b20-b21)

                  i = i + inc3
                  j = j + inc4
                END DO
                ibase = ibase + inc1
                jbase = jbase + inc21
              END DO
              ia = ia + iink
              ib = ib + iink
              ic = ic + iink
              id = id - iink
              ie = ie - iink
              jbase = jbase + jump
            END DO
          END IF

          IF (ib==id) THEN
            ibase = 0
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = (a(ia+i)+a(ib+i)) + a(ic+i)
                c(jb+j) = (qrt5*(a(ia+i)-a(ib+i))+(0.25_dp*(a(ia+i)+ &
                  a(ib+i))-a(ic+i))) - (sin36*b(ia+i)+sin72*b(ib+i))
                c(je+j) = -(qrt5*(a(ia+i)-a(ib+i))+(0.25_dp*(a(ia+i)+ &
                  a(ib+i))-a(ic+i))) - (sin36*b(ia+i)+sin72*b(ib+i))
                c(jc+j) = (qrt5*(a(ia+i)-a(ib+i))-(0.25_dp*(a(ia+i)+ &
                  a(ib+i))-a(ic+i))) - (sin72*b(ia+i)-sin36*b(ib+i))
                c(jd+j) = -(qrt5*(a(ia+i)-a(ib+i))-(0.25_dp*(a(ia+i)+ &
                  a(ib+i))-a(ic+i))) - (sin72*b(ia+i)-sin36*b(ib+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF

!!! case la=m
        ELSE
          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                t1 = (2.0_dp*(a(ia+i)-0.25_dp*(a(ib+i)+a(ic+i)))+qqrt5*(a(ib+ &
                  i)-a(ic+i))) - (ssin72*b(ib+i)+ssin36*b(ic+i))
                t2 = (2.0_dp*(a(ia+i)-0.25_dp*(a(ib+i)+a(ic+i)))-qqrt5*(a(ib+ &
                  i)-a(ic+i))) - (ssin36*b(ib+i)-ssin72*b(ic+i))
                t3 = (2.0_dp*(a(ia+i)-0.25_dp*(a(ib+i)+a(ic+i)))-qqrt5*(a(ib+ &
                  i)-a(ic+i))) + (ssin36*b(ib+i)-ssin72*b(ic+i))
                t4 = (2.0_dp*(a(ia+i)-0.25_dp*(a(ib+i)+a(ic+i)))+qqrt5*(a(ib+ &
                  i)-a(ic+i))) + (ssin72*b(ib+i)+ssin36*b(ic+i))
                a(ia+i) = 2.0_dp*(a(ia+i)+(a(ib+i)+a(ic+i)))
                a(ib+i) = t1
                b(ib+i) = t2
                a(ic+i) = t3
                b(ic+i) = t4
                i = i + inc3
              END DO
              ibase = ibase + inc1
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = 2.0_dp*(a(ia+i)+(a(ib+i)+a(ic+i)))
                c(jb+j) = (2.0_dp*(a(ia+i)-0.25_dp*(a(ib+i)+ &
                  a(ic+i)))+qqrt5*(a(ib+i)-a(ic+i))) - &
                  (ssin72*b(ib+i)+ssin36*b(ic+i))
                c(jc+j) = (2.0_dp*(a(ia+i)-0.25_dp*(a(ib+i)+ &
                  a(ic+i)))-qqrt5*(a(ib+i)-a(ic+i))) - &
                  (ssin36*b(ib+i)-ssin72*b(ic+i))
                c(jd+j) = (2.0_dp*(a(ia+i)-0.25_dp*(a(ib+i)+ &
                  a(ic+i)))-qqrt5*(a(ib+i)-a(ic+i))) + &
                  (ssin36*b(ib+i)-ssin72*b(ic+i))
                c(je+j) = (2.0_dp*(a(ia+i)-0.25_dp*(a(ib+i)+ &
                  a(ic+i)))+qqrt5*(a(ib+i)-a(ic+i))) + &
                  (ssin72*b(ib+i)+ssin36*b(ic+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF
        END IF

      ELSE IF (ifac==6) THEN

!     coding for factor 6
!     -------------------

        ia = 1
        ib = ia + (2*m-la)*inc1
        ic = ib + 2*m*inc1
        id = ic + 2*m*inc1
        ie = ic
        if = ib
        ja = 1
        jb = ja + jink
        jc = jb + jink
        jd = jc + jink
        je = jd + jink
        jf = je + jink

        IF (la/=m) THEN

          DO l = 1, ila
            i = ibase
            j = jbase
!DIR$ IVDEP
!CDIR NODEP
            DO ijk = 1, ilot
              c(ja+j) = (a(ia+i)+a(id+i)) + (a(ib+i)+a(ic+i))
              c(jd+j) = (a(ia+i)-a(id+i)) - (a(ib+i)-a(ic+i))
              c(jb+j) = ((a(ia+i)-a(id+i))+0.5_dp*(a(ib+i)-a(ic+i))) - &
                (sin60*(b(ib+i)+b(ic+i)))
              c(jf+j) = ((a(ia+i)-a(id+i))+0.5_dp*(a(ib+i)-a(ic+i))) + &
                (sin60*(b(ib+i)+b(ic+i)))
              c(jc+j) = ((a(ia+i)+a(id+i))-0.5_dp*(a(ib+i)+a(ic+i))) - &
                (sin60*(b(ib+i)-b(ic+i)))
              c(je+j) = ((a(ia+i)+a(id+i))-0.5_dp*(a(ib+i)+a(ic+i))) + &
                (sin60*(b(ib+i)-b(ic+i)))
              i = i + inc3
              j = j + inc4
            END DO
            ibase = ibase + inc1
            jbase = jbase + inc21
          END DO
          ia = ia + iink
          iink = 2*iink
          ib = ib + iink
          ic = ic + iink
          id = id - iink
          ie = ie - iink
          if = if - iink
          jbase = jbase + jump
          jump = 2*jump + jink

          IF (ic<id) THEN
            DO k = la, kstop, la
              kb = k + k
              kc = kb + kb
              kd = kc + kb
              ke = kd + kb
              kf = ke + kb
              c1 = trigs(kb+1)
              s1 = trigs(kb+2)
              c2 = trigs(kc+1)
              s2 = trigs(kc+2)
              c3 = trigs(kd+1)
              s3 = trigs(kd+2)
              c4 = trigs(ke+1)
              s4 = trigs(ke+2)
              c5 = trigs(kf+1)
              s5 = trigs(kf+2)
              ibase = 0
              DO l = 1, ila
                i = ibase
                j = jbase
!DIR$ IVDEP
!CDIR NODEP
                DO ijk = 1, ilot

                  a11 = (a(ie+i)+a(ib+i)) + (a(ic+i)+a(if+i))
                  a20 = (a(ia+i)+a(id+i)) - 0.5_dp*a11
                  a21 = sin60*((a(ie+i)+a(ib+i))-(a(ic+i)+a(if+i)))
                  b11 = (b(ib+i)-b(ie+i)) + (b(ic+i)-b(if+i))
                  b20 = (b(ia+i)-b(id+i)) - 0.5_dp*b11
                  b21 = sin60*((b(ib+i)-b(ie+i))-(b(ic+i)-b(if+i)))

                  c(ja+j) = (a(ia+i)+a(id+i)) + a11
                  d(ja+j) = (b(ia+i)-b(id+i)) + b11
                  c(jc+j) = c2*(a20-b21) - s2*(b20+a21)
                  d(jc+j) = s2*(a20-b21) + c2*(b20+a21)
                  c(je+j) = c4*(a20+b21) - s4*(b20-a21)
                  d(je+j) = s4*(a20+b21) + c4*(b20-a21)

                  a11 = (a(ie+i)-a(ib+i)) + (a(ic+i)-a(if+i))
                  b11 = (b(ie+i)+b(ib+i)) - (b(ic+i)+b(if+i))
                  a20 = (a(ia+i)-a(id+i)) - 0.5_dp*a11
                  a21 = sin60*((a(ie+i)-a(ib+i))-(a(ic+i)-a(if+i)))
                  b20 = (b(ia+i)+b(id+i)) + 0.5_dp*b11
                  b21 = sin60*((b(ie+i)+b(ib+i))+(b(ic+i)+b(if+i)))

                  c(jd+j) = c3*((a(ia+i)-a(id+i))+a11) - s3*((b(ia+i)+b(id+ &
                    i))-b11)
                  d(jd+j) = s3*((a(ia+i)-a(id+i))+a11) + c3*((b(ia+i)+b(id+ &
                    i))-b11)
                  c(jb+j) = c1*(a20-b21) - s1*(b20-a21)
                  d(jb+j) = s1*(a20-b21) + c1*(b20-a21)
                  c(jf+j) = c5*(a20+b21) - s5*(b20+a21)
                  d(jf+j) = s5*(a20+b21) + c5*(b20+a21)

                  i = i + inc3
                  j = j + inc4
                END DO
                ibase = ibase + inc1
                jbase = jbase + inc21
              END DO
              ia = ia + iink
              ib = ib + iink
              ic = ic + iink
              id = id - iink
              ie = ie - iink
              if = if - iink
              jbase = jbase + jump
            END DO
          END IF

          IF (ic==id) THEN
            ibase = 0
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = a(ib+i) + (a(ia+i)+a(ic+i))
                c(jd+j) = b(ib+i) - (b(ia+i)+b(ic+i))
                c(jb+j) = (sin60*(a(ia+i)-a(ic+i))) - (0.5_dp*(b(ia+i)+b(ic+ &
                  i))+b(ib+i))
                c(jf+j) = -(sin60*(a(ia+i)-a(ic+i))) - (0.5_dp*(b(ia+i)+b(ic+ &
                  i))+b(ib+i))
                c(jc+j) = sin60*(b(ic+i)-b(ia+i)) + (0.5_dp*(a(ia+i)+a(ic+ &
                  i))-a(ib+i))
                c(je+j) = sin60*(b(ic+i)-b(ia+i)) - (0.5_dp*(a(ia+i)+a(ic+ &
                  i))-a(ib+i))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF

!!! case la=m
        ELSE
          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                t1 = (2.0_dp*(a(ia+i)-a(id+i))+(a(ib+i)-a(ic+i))) - &
                  (ssin60*(b(ib+i)+b(ic+i)))
                t5 = (2.0_dp*(a(ia+i)-a(id+i))+(a(ib+i)-a(ic+i))) + &
                  (ssin60*(b(ib+i)+b(ic+i)))
                t2 = (2.0_dp*(a(ia+i)+a(id+i))-(a(ib+i)+a(ic+i))) - &
                  (ssin60*(b(ib+i)-b(ic+i)))
                t4 = (2.0_dp*(a(ia+i)+a(id+i))-(a(ib+i)+a(ic+i))) + &
                  (ssin60*(b(ib+i)-b(ic+i)))
                t3 = (2.0_dp*(a(ia+i)-a(id+i))) - (2.0_dp*(a(ib+i)-a(ic+i)))
                a(ia+i) = (2.0_dp*(a(ia+i)+a(id+i))) + (2.0_dp*(a(ib+i)+a(ic+ &
                  i)))
                a(ib+i) = t1
                b(ib+i) = t2
                a(ic+i) = t3
                b(ic+i) = t4
                a(id+i) = t5
                i = i + inc3
              END DO
              ibase = ibase + inc1
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = (2.0_dp*(a(ia+i)+a(id+i))) + (2.0_dp*(a(ib+i)+a(ic+ &
                  i)))
                c(jd+j) = (2.0_dp*(a(ia+i)-a(id+i))) - (2.0_dp*(a(ib+i)-a(ic+ &
                  i)))
                c(jb+j) = (2.0_dp*(a(ia+i)-a(id+i))+(a(ib+i)-a(ic+i))) - &
                  (ssin60*(b(ib+i)+b(ic+i)))
                c(jf+j) = (2.0_dp*(a(ia+i)-a(id+i))+(a(ib+i)-a(ic+i))) + &
                  (ssin60*(b(ib+i)+b(ic+i)))
                c(jc+j) = (2.0_dp*(a(ia+i)+a(id+i))-(a(ib+i)+a(ic+i))) - &
                  (ssin60*(b(ib+i)-b(ic+i)))
                c(je+j) = (2.0_dp*(a(ia+i)+a(id+i))-(a(ib+i)+a(ic+i))) + &
                  (ssin60*(b(ib+i)-b(ic+i)))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF
        END IF

      ELSE IF (ifac==8) THEN

!     coding for factor 8
!     -------------------

        IF (la/=m) THEN
          ibad = 3
        ELSE
          ia = 1
          ib = ia + la*inc1
          ic = ib + 2*la*inc1
          id = ic + 2*la*inc1
          ie = id + 2*la*inc1
          ja = 1
          jb = ja + jink
          jc = jb + jink
          jd = jc + jink
          je = jd + jink
          jf = je + jink
          jg = jf + jink
          jh = jg + jink

          IF (lipl) THEN
            DO l = 1, ila
              i = ibase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                t2 = 2.0_dp*(((a(ia+i)+a(ie+i))-a(ic+i))-(b(ib+i)-b(id+i)))
                t6 = 2.0_dp*(((a(ia+i)+a(ie+i))-a(ic+i))+(b(ib+i)-b(id+i)))
                t1 = 2.0_dp*((a(ia+i)-a(ie+i))-b(ic+i)) + ssin45*((a(ib+ &
                  i)-a(id+i))-(b(ib+i)+b(id+i)))
                t5 = 2.0_dp*((a(ia+i)-a(ie+i))-b(ic+i)) - ssin45*((a(ib+ &
                  i)-a(id+i))-(b(ib+i)+b(id+i)))
                t3 = 2.0_dp*((a(ia+i)-a(ie+i))+b(ic+i)) - ssin45*((a(ib+ &
                  i)-a(id+i))+(b(ib+i)+b(id+i)))
                t7 = 2.0_dp*((a(ia+i)-a(ie+i))+b(ic+i)) + ssin45*((a(ib+ &
                  i)-a(id+i))+(b(ib+i)+b(id+i)))
                t4 = 2.0_dp*(((a(ia+i)+a(ie+i))+a(ic+i))-(a(ib+i)+a(id+i)))
                a(ia+i) = 2.0_dp*(((a(ia+i)+a(ie+i))+a(ic+i))+(a(ib+i)+a(id+ &
                  i)))
                a(ib+i) = t1
                b(ib+i) = t2
                a(ic+i) = t3
                b(ic+i) = t4
                a(id+i) = t5
                b(id+i) = t6
                a(ie+i) = t7
                i = i + inc3
              END DO
              ibase = ibase + inc1
            END DO
          ELSE
            DO l = 1, ila
              i = ibase
              j = jbase
!DIR$ IVDEP
!CDIR NODEP
              DO ijk = 1, ilot
                c(ja+j) = 2.0_dp*(((a(ia+i)+a(ie+i))+a(ic+i))+(a(ib+i)+a(id+ &
                  i)))
                c(je+j) = 2.0_dp*(((a(ia+i)+a(ie+i))+a(ic+i))-(a(ib+i)+a(id+ &
                  i)))
                c(jc+j) = 2.0_dp*(((a(ia+i)+a(ie+i))-a(ic+i))-(b(ib+i)-b(id+ &
                  i)))
                c(jg+j) = 2.0_dp*(((a(ia+i)+a(ie+i))-a(ic+i))+(b(ib+i)-b(id+ &
                  i)))
                c(jb+j) = 2.0_dp*((a(ia+i)-a(ie+i))-b(ic+i)) + ssin45*((a(ib+ &
                  i)-a(id+i))-(b(ib+i)+b(id+i)))
                c(jf+j) = 2.0_dp*((a(ia+i)-a(ie+i))-b(ic+i)) - ssin45*((a(ib+ &
                  i)-a(id+i))-(b(ib+i)+b(id+i)))
                c(jd+j) = 2.0_dp*((a(ia+i)-a(ie+i))+b(ic+i)) - ssin45*((a(ib+ &
                  i)-a(id+i))+(b(ib+i)+b(id+i)))
                c(jh+j) = 2.0_dp*((a(ia+i)-a(ie+i))+b(ic+i)) + ssin45*((a(ib+ &
                  i)-a(id+i))+(b(ib+i)+b(id+i)))
                i = i + inc3
                j = j + inc4
              END DO
              ibase = ibase + inc1
              jbase = jbase + inc21
            END DO
          END IF

        END IF

      ELSE


        ibad = 2
!!! illegal factor
      END IF

!     return
!     ------

      ierr = ibad
      RETURN
    END SUBROUTINE rpassf

END MODULE mo_fft992
