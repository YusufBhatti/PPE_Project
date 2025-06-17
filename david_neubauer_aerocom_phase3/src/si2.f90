#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE si2

  ! Description:
  !
  ! 2nd part of the semi-implicit scheme (done in fourier space).
  !
  ! Method:
  !
  ! This subroutine computes the contribution in
  ! *Fourier space to the semi-implicit scheme (mainly for
  ! the vorticity and humidity equations).
  !
  ! *si2* is called from *fcc1*.
  !
  ! Reference:
  ! 1-appendix *b1:organisation of the spectral model.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! I. Kirchner, MPI, December 1998, vorticity semi-implicit diagnosticsterm modified
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! I. Kirchner, MPI, May 2002, bugfix semi-implicit part vorticity
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp 
  USE mo_control,       ONLY: ltimer
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_si2
  USE mo_gaussgrid,     ONLY: gl_racst, gl_twomu
  USE mo_truncation,    ONLY: am
  USE mo_semi_impl,     ONLY: betazq
  USE mo_physical_constants,     ONLY: earth_radius
  USE mo_time_control,  ONLY: time_step_len
  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE mo_buffer_fft,    ONLY: fdm1, fvol, fvom, fu0, fdu0, fvom1

  IMPLICIT NONE

  !  Local scalars: 
  REAL(dp):: z1, z2, z3, zdt, zbmi, zbmr, zcmi, zcmr, zdl, zdtma2,     &
             zdtmda2, zdu0, zrcst, zu0, zvoli, zvolr, zvomi, zvomr
  INTEGER :: jlev, jm, jrow, jglat
  INTEGER :: nflev, nmp1, nflat

  !  Executable statements 

  ! loop bounds on this PE

  nflev = dc% nflev  ! local (Fourier space) : number of levels 
  nflat = dc% nflat  ! local (Fourier space) : number of latitudes  
  nmp1  = dc% nm+1   ! global                : number of coefficients m

!-- 1. Locate and allocate storage

!-- 1.1 Fourier components

  ! latitude loop indices

!$OMP PARALLEL
  IF (ltimer) CALL timer_start(timer_si2)
!$OMP DO PRIVATE(jrow,jglat,zdt,z1,z2,z3,zrcst,jlev,zu0,zdu0,jm,zdtma2,zdtmda2,   &
!$OMP            zdl,zbmr,zbmi,zcmr,zcmi,zvolr,zvoli,zvomr,zvomi)
  lat_loop: DO jrow = 1, nflat

    jglat = dc% glat(jrow)                 ! global continuous north -> south

!-- 2. Skip over *si2* during initialisation iterations
!      or compute temporary quantities.

!-- 2.1 Compute temporary quantities

    zdt   = 0.5_dp*time_step_len
    z1    = betazq*zdt*gl_racst(jglat)
    z2    = z1*earth_radius
    z3    = gl_twomu(jglat)*gl_racst(jglat)
    zrcst = gl_racst(jglat)*earth_radius

!-- 3. Semi implicit modifications

!-- 3.1 Initiate scans

    level_loop: DO jlev = 1, nflev

      zu0  = z1*fu0(jlev,jrow)
      zdu0 = z2*(fdu0(jlev,jrow)+z3*fu0(jlev,jrow))

!DIR$ IVDEP
!OCL NOVREC
      spec_loop: DO jm = 1, nmp1

        zdtma2 = am(jm)*zu0
        zdtmda2 = am(jm)*zdu0
        
        zdl = am(jm)*zrcst

        zbmr = 1._dp/(1._dp+zdtma2**2)
        zbmi = -zdtma2*zbmr

        zcmr = -2._dp*zdtmda2*    zdtma2    *zbmr**2
        zcmi =    -zdtmda2*(1._dp-zdtma2**2)*zbmr**2

!-- 3.2 Divergence equation

        fdm1(2*jm-1,jlev,jrow) =fdm1(2*jm-1,jlev,jrow)+zdl*fvom(2*jm,  jlev,jrow)
        fdm1(2*jm,  jlev,jrow) =fdm1(2*jm,  jlev,jrow)-zdl*fvom(2*jm-1,jlev,jrow)

!-- 3.3 Vorticity equation
        zvolr = zbmr*(fvom1(2*jm-1,jlev,jrow) - am(jm)*fvol(2*jm,  jlev,jrow))  &
              - zbmi*(fvom1(2*jm,  jlev,jrow) + am(jm)*fvol(2*jm-1,jlev,jrow))  &
              - zcmr* fvom (2*jm-1,jlev,jrow) &
              + zcmi* fvom (2*jm,  jlev,jrow)

        zvoli = zbmi* (fvom1(2*jm-1,jlev,jrow) - am(jm)*fvol(2*jm,  jlev,jrow)) &
              + zbmr* (fvom1(2*jm,  jlev,jrow) + am(jm)*fvol(2*jm-1,jlev,jrow)) &
              - zcmi* fvom (2*jm-1,jlev,jrow) &
              - zcmr* fvom (2*jm,  jlev,jrow)

        zvomr = zbmr*fvom(2*jm-1,jlev,jrow) - zbmi*fvom(2*jm,jlev,jrow)
        zvomi = zbmi*fvom(2*jm-1,jlev,jrow) + zbmr*fvom(2*jm,jlev,jrow)

!!$        IF (ltdiag) THEN
!!$          ! store terms without semi-implicit part
!!$
!!$          ! L-Term before adjustment
!!$          zdigslr = fvom1(2*jm-1,jlev,jrow) - am(jm)*fvol(2*jm  ,jlev,jrow)
!!$          zdigsli = fvom1(2*jm  ,jlev,jrow) + am(jm)*fvol(2*jm-1,jlev,jrow)
!!$
!!$          ! M-Term before adjustment
!!$          zdigsmr = fvom(2*jm-1,jlev,jrow)
!!$          zdigsmi = fvom(2*jm  ,jlev,jrow)
!!$
!!$        ENDIF

        fvol(2*jm-1,jlev,jrow) = zvolr
        fvol(2*jm,  jlev,jrow) = zvoli
        fvom(2*jm-1,jlev,jrow) = zvomr
        fvom(2*jm,  jlev,jrow) = zvomi

!!$        IF (ltdiag) THEN
!!$          ! d/dl of explicit part (accounted in SI1)
!!$          zdigsr = - am(jm)*pdigs(2*jm  ,jlev,1,jrow)
!!$          zdigsi =   am(jm)*pdigs(2*jm-1,jlev,1,jrow)
!!$
!!$          ! adjustment of semi-implicit part of vorticity, L-Term
!!$          pdigs(2*jm-1,jlev,1,jrow) = zdigsr + (zvolr - zdigslr)
!!$          pdigs(2*jm,  jlev,1,jrow) = zdigsi + (zvoli - zdigsli)
!!$
!!$          ! adjustment of semi-implicit part of vorticity, M-Term
!!$          pdigs(2*jm-1,jlev,4,jrow) = zvomr - zdigsmr
!!$          pdigs(2*jm,  jlev,4,jrow) = zvomi - zdigsmi
!!$
!!$        ENDIF

      END DO spec_loop

!!$      IF (ltdiag) THEN
!!$        IF(l_putdata(dio_index)) THEN
!!$
!!$          DO jm = 1,nmp1
!!$            ! correct zonal derivatives with factor m/(1-mu**2)
!!$            zdl = am(jm)*zrcst
!!$            pdigb(2*jm-1,jlev,:,jrow) = - zdl * pdiga(2*jm,  jlev,1:10,jrow)
!!$            pdigb(2*jm,  jlev,:,jrow) =   zdl * pdiga(2*jm-1,jlev,1:10,jrow)
!!$          ENDDO
!!$
!!$        END IF
!!$      END IF

    END DO level_loop

  END DO lat_loop
!$OMP END DO
  IF (ltimer) CALL timer_stop(timer_si2)
!$OMP END PARALLEL

END SUBROUTINE si2
