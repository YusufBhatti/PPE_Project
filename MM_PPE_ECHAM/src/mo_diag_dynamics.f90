!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_diag_dynamics

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_diag_dynamics
  PUBLIC :: diag_dynamics
  PUBLIC :: print_diag_dynamics

  ! global mean statistics

  REAL(dp), ALLOCATABLE :: voh(:)  ! horizontal RMS of vorticity (nlev)
  REAL(dp), ALLOCATABLE :: dh(:)   ! horizontal RMS of divergence (nlev)

  REAL(dp), ALLOCATABLE :: th(:)   ! horizontal average of temperature (nlev)
  REAL(dp), ALLOCATABLE :: qh(:)   ! horizontal average of humidity (nlev)
  REAL(dp), ALLOCATABLE :: qlh(:)  ! horizontal average of liquid water (nlev)
  REAL(dp), ALLOCATABLE :: qih(:)  ! horizontal average of ice water (nlev)

  REAL(dp) :: gvo         ! global RMS of vorticity.
  REAL(dp) :: gd          ! global RMS of divergence.

  REAL(dp) :: gt          ! global average of temperature.
  REAL(dp) :: gq          ! global average of humidity.
  REAL(dp) :: gql         ! global average of liquid water
  REAL(dp) :: gqi         ! global average of ice water

  REAL(dp) :: gps         ! mean surface pressure.

  REAL(dp) :: gke         ! global kinetic energy.
  REAL(dp) :: gpe         ! global potential energy.
  REAL(dp) :: gte         ! global total energy (dry).
  REAL(dp) :: glq         ! global latent heat energy.
  REAL(dp) :: gtpe        ! global total energy (wet, gte+glq).

!  REAL(dp) :: gqm         ! equivalent water content.
!
!  REAL(dp) :: gts         ! land top layer energy.
!  REAL(dp) :: gtd         ! land deep layer energy.
!  REAL(dp) :: gws         ! land top layer water content.
!  REAL(dp) :: gwd         ! land deep layer water content.
!  REAL(dp) :: gsn         ! land snow equivalent depth.

  ! zonal mean statistics

  REAL(dp), ALLOCATABLE :: delph(:)   ! nlev

  REAL(dp), ALLOCATABLE :: voz(:)     ! zonal RMS of vorticity (nlat)
  REAL(dp), ALLOCATABLE :: dz(:)      ! zonal RMS of divergence (nlat)
  REAL(dp), ALLOCATABLE :: qz(:)      ! zonal average of humidity (nlat)
  REAL(dp), ALLOCATABLE :: qlz(:)     ! zonal average of liquid water (nlat)
  REAL(dp), ALLOCATABLE :: qiz(:)     ! zonal average of ice water (nlat)
  REAL(dp), ALLOCATABLE :: tz(:)      ! zonal average of temperature (nlat)
  REAL(dp), ALLOCATABLE :: psz(:)     ! zonal mean surface pressure (nlat)

  REAL(dp), ALLOCATABLE :: gkez(:)    ! zonal sum of kinetic energy (nlat)
  REAL(dp), ALLOCATABLE :: gpez(:)    ! zonal sum of potential energy (nlat)
  REAL(dp), ALLOCATABLE :: glqz(:)    ! zonal sum of latent heat energy (nlat)
  REAL(dp), ALLOCATABLE :: gpsz(:)    ! zonal sum of surface pressure (nlat)

  REAL(dp), ALLOCATABLE :: gtsz(:)    ! zonal sum of land top layer energy (nlat)
  REAL(dp), ALLOCATABLE :: gtdz(:)    ! zonal sum of land deep layer energy (nlat)
  REAL(dp), ALLOCATABLE :: gwsz(:)    ! zonal sum of land top layer water content (nlat)
  REAL(dp), ALLOCATABLE :: gwdz(:)    ! zonal sum of land deep layer water content (nlat)
  REAL(dp), ALLOCATABLE :: gsnz(:)    ! zonal sum of land snow equivalent depth (nlat

CONTAINS

  SUBROUTINE diag_dynamics

    ! Description:
    !
    ! Computes some statistics for the prognostic dynamical variables.
    !
    ! Method:
    !
    ! diag_dynamics is called from scan1
    !
    ! Results:
    !
    ! diag_dynamics computes and print with a chosen frequency the RMS of
    ! vorticity, divergence and means of temperature, surface pressure,
    ! humidity and energy(kinetic,potential and total).
    ! All these quantities are accumulated for each level
    ! in arrays.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, September 1982, original source
    ! U. Schlese, MPI,  1990,  zonal diagnostics added
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! U. Schlese, DKRZ, November 1999, cloud ice added
    ! L. Kornblueh, MPI, October 2001, parallelized and packed in a module
    !
    ! for more details see file AUTHORS
    !

    USE mo_control,       ONLY: nlev, nlevp1, nlon, ngl
    USE mo_gaussgrid,     ONLY: gl_budw, gl_rcsth
    USE mo_physical_constants, ONLY: grav, als, alv, tmelt, cpd, vtmpc2
    USE mo_scan_buffer,   ONLY: d, t, u, v, vo
    USE mo_memory_gl,     ONLY: q, xl, xi
    USE mo_memory_g3a,    ONLY: geospm
    USE mo_memory_g3b,    ONLY: aps
    USE mo_mpi,           ONLY: p_pe, p_io
    USE mo_decomposition, ONLY: gl_dc => global_decomposition
    USE mo_transpose,     ONLY: gather_gp

    !  Local scalars:
    REAL(dp) :: zgpe, zgps, zke, zlq, zpe, zqnlon, zqpsz, zw, zwke, zwpe, &
                zwsg, zzdelp
    INTEGER :: jlev, jlon, jglat

    !  Local arrays:
    REAL(dp) :: zdelp(nlev), zdh(nlev), zth(nlev), zvoh(nlev), &
                zqh(nlev), zqlh(nlev), zqih(nlev)
    REAL(dp) :: zph(nlon,nlevp1)

    REAL(dp), POINTER :: zaps(:,:), zgeo(:,:)
    REAL(dp), POINTER :: zvo(:,:,:), zd(:,:,:), zq(:,:,:),  &
                         zxl(:,:,:), zxi(:,:,:), zt(:,:,:), &
                         zu(:,:,:), zv(:,:,:)

    !  Intrinsic functions
    INTRINSIC DOT_PRODUCT, SQRT, SUM

    !  External functions
    EXTERNAL pres

    !  Executable statements

    ! gather fields
    
    IF (p_pe == p_io) THEN
      ALLOCATE (zaps(nlon,       ngl))
      ALLOCATE (zgeo(nlon,       ngl))
      ALLOCATE (zvo (nlon, nlev, ngl))
      ALLOCATE (zd  (nlon, nlev, ngl))
      ALLOCATE (zq  (nlon, nlev, ngl))
      ALLOCATE (zxl (nlon, nlev, ngl))
      ALLOCATE (zxi (nlon, nlev, ngl))
      ALLOCATE (zt  (nlon, nlev, ngl))
      ALLOCATE (zu  (nlon, nlev, ngl))
      ALLOCATE (zv  (nlon, nlev, ngl))
    END IF
    
    CALL gather_gp (zaps, aps,    gl_dc)
    CALL gather_gp (zgeo, geospm, gl_dc)
    CALL gather_gp (zvo,  vo,     gl_dc)
    CALL gather_gp (zd,   d,      gl_dc) 
    CALL gather_gp (zq,   q,      gl_dc)     
    CALL gather_gp (zxl,  xl,     gl_dc)    
    CALL gather_gp (zxi,  xi,     gl_dc)    
    CALL gather_gp (zt,   t,      gl_dc)    
    CALL gather_gp (zu,   u,      gl_dc)    
    CALL gather_gp (zv,   v,      gl_dc)

    IF (p_pe == p_io) THEN

      ! Compute statistics

      DO jglat = 1, ngl

        ! Set up weighting constants
        
        zw = gl_budw(jglat)
        zwsg = zw/grav
        zwke = zwsg*gl_rcsth(jglat)
        zwpe = zwsg*cpd
        zqnlon = 1.0_dp/REAL(nlon,dp)
        
        ! Accumulate horizontal statistics
        
        voz(jglat) = 0.0_dp
        dz(jglat)  = 0.0_dp
        qz(jglat)  = 0.0_dp
        qlz(jglat) = 0.0_dp
        qiz(jglat) = 0.0_dp
        tz(jglat)  = 0.0_dp
        
        zph(:,nlevp1) = zaps(1:nlon,jglat)
        call pres (zph,nlon,zph(:,nlevp1),nlon)
      
        psz(jglat) = SUM(zph(:,nlevp1))*zqnlon
      
        DO jlev = 1, nlev
          zdelp(jlev) = 0.0_dp
          zvoh(jlev)  = 0.0_dp
          zdh(jlev)   = 0.0_dp
          zqh(jlev)   = 0.0_dp
          zqlh(jlev)  = 0.0_dp
          zqih(jlev)  = 0.0_dp
          zth(jlev)   = 0.0_dp
          
          DO jlon = 1, nlon
            zzdelp = zph(jlon,jlev+1) - zph(jlon,jlev)
            zdelp(jlev) = zdelp(jlev) + zzdelp
            zvoh (jlev) = zvoh (jlev) + zzdelp* zvo(jlon,jlev,jglat)**2
            zdh  (jlev) = zdh  (jlev) + zzdelp* zd (jlon,jlev,jglat)**2
            zqh  (jlev) = zqh  (jlev) + zzdelp* zq (jlon,jlev,jglat)
            zqlh (jlev) = zqlh (jlev) + zzdelp* zxl(jlon,jlev,jglat)
            zqih (jlev) = zqih (jlev) + zzdelp* zxi(jlon,jlev,jglat)
            zth  (jlev) = zth  (jlev) + zzdelp* zt (jlon,jlev,jglat)
          END DO
          
          voz(jglat) = voz(jglat) + zvoh(jlev)*zqnlon
          dz(jglat)  = dz(jglat)  + zdh(jlev)*zqnlon
          qz(jglat)  = qz(jglat)  + zqh(jlev)*zqnlon
          qlz(jglat) = qlz(jglat) + zqlh(jlev)*zqnlon
          qiz(jglat) = qiz(jglat) + zqih(jlev)*zqnlon
          tz(jglat)  = tz(jglat)  + zth(jlev)*zqnlon
        END DO
    
      delph(:) = delph(:) + zw*zdelp(:)
      voh(:)   = voh(:)   + zw*zvoh(:)
      dh(:)    = dh(:)    + zw*zdh(:)
      qh(:)    = qh(:)    + zw*zqh(:)
      qlh(:)   = qlh(:)   + zw*zqlh(:)
      qih(:)   = qih(:)   + zw*zqih(:)
      th(:)    = th(:)    + zw*zth(:)
      
      zqpsz = 1.0_dp/psz(jglat)
      voz(jglat) = SQRT(voz(jglat)*zqpsz)
      dz(jglat)  = SQRT(dz(jglat)*zqpsz)
      qz(jglat)  = qz(jglat)*zqpsz
      qlz(jglat) = qlz(jglat)*zqpsz
      qiz(jglat) = qiz(jglat)*zqpsz
      tz(jglat)  = tz(jglat)*zqpsz
      
      zgps = zw*psz(jglat)*nlon
      zgpe = zwsg*DOT_PRODUCT(zph(1:nlon,nlevp1),zgeo(1:nlon,jglat))
      
      ! Accumulate global statistics.
      
      zke = 0.0_dp
      zpe = 0.0_dp
      zlq = 0.0_dp
      
      DO jlev = 1, nlev
        DO jlon = 1, nlon
          zzdelp = zph(jlon,jlev+1) - zph(jlon,jlev)
          zke = zke+zzdelp*(zu(jlon,jlev,jglat)**2+zv(jlon,jlev,jglat)**2)
          zpe = zpe+zzdelp*(1._dp+vtmpc2*zq(jlon,jlev,jglat))*zt(jlon,jlev,jglat)
          IF (zt(jlon,jlev,jglat) > tmelt) THEN
            zlq = zlq+zzdelp*zq(jlon,jlev,jglat)*alv
          ELSE
            zlq = zlq+zzdelp*zq(jlon,jlev,jglat)*als
          END IF
        END DO
      END DO

      gkez(jglat) = zwke*zke
      gpez(jglat) = zwpe*zpe + zgpe
      glqz(jglat) = zwsg*zlq
      gpsz(jglat) = zgps
  
    END DO

    DEALLOCATE (zaps)
    DEALLOCATE (zgeo)
    DEALLOCATE (zvo)
    DEALLOCATE (zd)
    DEALLOCATE (zq)
    DEALLOCATE (zxl)
    DEALLOCATE (zxi)
    DEALLOCATE (zt)
    DEALLOCATE (zu)
    DEALLOCATE (zv)

  END IF
    
  END SUBROUTINE diag_dynamics

  SUBROUTINE init_diag_dynamics

    ! Description:
    !
    ! Blank statistics at beginning of time step
    !
    ! Method:
    !
    ! init_diag_dynamics is called from scan1
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, June 1983, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! U. Schlese, DKRZ, March 2000, Physics diagnostics removed
    ! I. Kirchner, MPI, December 2000, time control
    ! L. Kornblueh, MPI, October 2001, parallelization and module packing
    !
    ! for more details see file AUTHORS
    !

    USE mo_control, ONLY: nlev, ngl

    LOGICAL, SAVE :: not_used = .TRUE.

    IF (not_used) THEN
      ALLOCATE (voh(nlev))
      ALLOCATE (dh (nlev))
      ALLOCATE (th (nlev))
      ALLOCATE (qh (nlev))
      ALLOCATE (qlh(nlev))
      ALLOCATE (qih(nlev))

      ALLOCATE (delph(nlev))

      ALLOCATE (voz(ngl))
      ALLOCATE (dz(ngl))
      ALLOCATE (qz(ngl))
      ALLOCATE (qlz(ngl))
      ALLOCATE (qiz(ngl))
      ALLOCATE (tz(ngl))
      ALLOCATE (psz(ngl))

      ALLOCATE (gkez(ngl))
      ALLOCATE (gpez(ngl))
      ALLOCATE (glqz(ngl))
      ALLOCATE (gpsz(ngl))

      ALLOCATE (gtsz(ngl))
      ALLOCATE (gtdz(ngl))
      ALLOCATE (gwsz(ngl))
      ALLOCATE (gwdz(ngl))
      ALLOCATE (gsnz(ngl))

      not_used = .FALSE.
    END IF


    !  Executable statements

    ! Blank statistics at beginning of time step

    voh(:)   = 0.0_dp
    dh(:)    = 0.0_dp
    qh(:)    = 0.0_dp
    qlh(:)   = 0.0_dp
    qih(:)   = 0.0_dp
    th(:)    = 0.0_dp
    delph(:) = 0.0_dp
    gvo      = 0.0_dp
    gd       = 0.0_dp
    gq       = 0.0_dp
    gql      = 0.0_dp
    gqi      = 0.0_dp
    gt       = 0.0_dp
    gps      = 0.0_dp
    gke      = 0.0_dp
    gpe      = 0.0_dp
    glq      = 0.0_dp
    gte      = 0.0_dp
    gtpe     = 0.0_dp

  END SUBROUTINE init_diag_dynamics

  SUBROUTINE print_diag_dynamics

    ! Description:
    !
    ! Complete and print statistics for dynamics
    !
    ! Method:
    !
    ! print_diag_dynamics is called from scan1
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, June 1983, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! I. Kirchner, MPI, November 2000, date/time control
    ! L. Kornblueh, MPI, October 2001, parallelized and packed in a module
    !
    ! for more details see file AUTHORS
    !

    USE mo_io_units,            ONLY: nout
    USE mo_physical_constants,  ONLY: grav
    USE mo_control,           ONLY: nlev, ngl
    USE mo_gaussgrid,         ONLY: philat
    USE mo_semi_impl,         ONLY: ulm, uvmax, nulev, ulat
    USE mo_forecast_switches, ONLY: lumax, lzondia
    USE mo_time_control,      ONLY: get_time_step, l_diagvert
    USE mo_mpi,               ONLY: p_pe, p_io

    INTEGER :: istep

    !  Local scalars:
    REAL(dp) :: gqp, gqlp, gqip, qhd, qlhd, qihd, zlat, zqdelp
    INTEGER  :: jglat, jlev

    !  Intrinsic functions
    INTRINSIC SQRT, SUM


    !  Executable statements

    istep = get_time_step()

    ! Complete statistics

    IF (p_pe == p_io) THEN

      ! Complete global sums
      
      gke = SUM(gkez)
      gpe = SUM(gpez)
      glq = SUM(glqz)
      gps = SUM(gpsz)
      
      gvo = SUM(voh)
      gd  = SUM(dh)
      gq  = SUM(qh)
      gql = SUM(qlh)
      gqi = SUM(qih)
      gt  = SUM(th)
      
      IF (l_diagvert) THEN
        DO jlev = 1, nlev
          zqdelp    = 1.0_dp/delph(jlev)
          voh(jlev) = SQRT(voh(jlev)*zqdelp)
          dh(jlev)  = SQRT(dh(jlev)*zqdelp)
          qh(jlev)  = qh(jlev)*zqdelp
          qlh(jlev) = qlh(jlev)*zqdelp
          qih(jlev) = qih(jlev)*zqdelp
          th(jlev)  = th(jlev)*zqdelp
        END DO
      END IF
      
      gvo  = SQRT(gvo/gps)
      gd   = SQRT(gd/gps)
      gq   = gq/gps
      gql  = gql/gps
      gqi  = gqi/gps
      gt   = gt/gps
      gte  = gke+gpe
      gtpe = gte+glq
      
      ! Print statistics
      
      WRITE (nout, &
           '(/,a,i6,/,2(a,e11.5,/),a,f10.2,/,3(a,e11.5,/),a,f10.2,/,5(a,e11.5,/))') &
           ' Step: ', istep,                                &
           '  std VO [1/s]                : ', gvo,         &
           '  std D  [1/s]                : ', gd,          &
           '  T [K]                       : ', gt,          &
           '  q                           : ', gq,          &
           '  ql                          : ', gql,         &
           '  qi                          : ', gqi,         &
           '  ps [hPa]                    : ', gps*0.01_dp, &
           '  kinetic energy [J/m**2]     : ', gke,         &
           '  potential energy [J/m**2]   : ', gpe,         &
           '  total energy (dry) [J/m**2] : ', gte,         &
           '  latent energy [J/m**2]      : ', glq,         &
           '  total energy (wet) [J/m**2] : ', gtpe
      
      IF (l_diagvert) THEN
        gqp  = gq*gps/grav
        gqlp = gql*gps/grav
        gqip = gqi*gps/grav
        WRITE (nout,'(/,a,3(e11.4),/)') &
             'q*ps/g, ql*ps/g, qi*ps/g =', gqp, gqlp, gqip
        
        DO jlev = 1, nlev
          qhd  = qh(jlev)*delph(jlev)/grav
          qlhd = qlh(jlev)*delph(jlev)/grav
          qihd = qih(jlev)*delph(jlev)/grav
          WRITE (nout,'(4x,a,i2,a,2(1p,e11.4),0p,f10.2,1x,6(1p,e11.4))') &
               '(',jlev,')', voh(jlev), dh(jlev), th(jlev), &
               qh(jlev), qlh(jlev), qih(jlev), qhd, qlhd, qihd
        END DO
        
        WRITE (nout,*)
        
      END IF
      
      IF (lumax) THEN
        WRITE (nout,'(a,f7.2,a,i3,a,f7.1,a,f7.1)') &
             '=> lat = ',ulat, 'deg  level = ', nulev, &
             '  max sqrt(u**2+v**2) = ',uvmax, '  ul=', ulm
      END IF
      
      IF (lzondia) THEN
        WRITE (nout,'(a,i7)') 'Zonal statistics, nstep= ', istep
        DO jglat = 1, ngl
          zlat = philat(jglat)
          WRITE (nout,'(2x,f5.0,2(3(1p,e11.4),0p,f10.2,1x))') &
               zlat, voz(jglat), dz(jglat), tz(jglat), qz(jglat), &
               qlz(jglat), qiz(jglat), psz(jglat)*0.01_dp
        END DO
        
        WRITE (nout,*)
        
      END IF

    END IF

  END SUBROUTINE print_diag_dynamics

END MODULE mo_diag_dynamics
