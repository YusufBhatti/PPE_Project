!OPTION! -O nomovediv

#ifdef __xlC__
@PROCESS NOSTRICT
@PROCESS HOT
@PROCESS XLF90(nosignedzero)
#else
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#endif

#ifndef _STRIP_
#define _STRIP_ 256
#endif

#if (defined _UNICOSMP) || (defined __SX__) || (defined ES)
#define VECTOR 1
#else
#define __tuned__
#endif

#define _NOMODULO

#define _PACKED_TRANSPOSITION
#undef CHECK_TRANSPOSITION
#undef STANDARD_FORTRAN

#ifndef ECHAM_GPC_SCHEDULE
#define ECHAM_GPC_SCHEDULE SCHEDULE(STATIC)
#endif

!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

MODULE mo_tpcore

  ! S.-J. Lin, GSFC, unknown, Original source
  ! S.-J. Lin, GSFC, July 2001, Last modified
  ! L. Kornblueh, MPI, August 2001, modified for Echam 4/5
  ! L. Kornblueh, MPI, October 2001, removed bug for yms
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! A. Rhodin DWD/MPI, June    2002, transport only selected tracers
  ! L. Kornblueh, MPI, November 2002, optimizations for vector architecture
  ! L. Kornblueh, MPI, February 2004, more optimizations for vector 
  !                                   architecture and parallelization
  ! S. Shingu, NEC, March 2004,  more optimizations for vector 
  !                              architecture and parallelization
  !                              proposed more or less the same things
  !                              as L. Kornblueh above plus some
  !                              additional communication improvements 
  ! M. Rosenhauer, DKRZ, August 2004, bugfix in vertical transport
  ! L. Kornblueh, MPI, January 2008, added optimizations from 
  !                                  Stefan Borowski, NEC
  !===========================================================================
  !
  ! TransPort module for NASA Goddard Chemistry Transport Model
  !
  ! Purpose: perform the transport of  3-D mixing ratio fields using
  !          externally specified winds and surface pressure on the
  !          hybrid Eta-coordinate.
  !          One call to tpcore updates the 3-D mixing ratio
  !          fields for one time step (DT).
  !
  ! Schemes: Multi-dimensional Flux Form Semi-Lagrangian (FFSL) schemes
  !          (Lin and Rood 1996, MWR) with many unpublished modifications
  !
  ! Programer: S.-J. Lin
  ! Messaging passing library based on "Pilgrim" developed by W. Sawyer
  ! Last modified: July, 2001
  !
  ! Send comments/suggestions to
  !
  !                 S.-J. Lin
  !                 Code 910.3, NASA/GSFC, Greenbelt, MD 20771
  !                 E-mail: slin@dao.gsfc.nasa.gov
  !
  ! The algorithm is primarily based on the following papers:
  !
  ! 1. Lin, S.-J., and R. B. Rood, 1996: Multidimensional flux form semi-
  !    Lagrangian transport schemes. Mon. Wea. Rev., 124, 2046-2070.
  !
  ! 2. Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994: A class of
  !    the van Leer-type transport schemes and its applications to the moist-
  !    ure transport in a General Circulation Model. Mon. Wea. Rev., 122,
  !    1575-1593.

  USE mo_kind,          ONLY: wp
#if defined(DEBUG) || defined(CHECK_TRANSPOSITION)
  USE mo_exception,     ONLY: message, message_text, finish
#else
  USE mo_exception,     ONLY: finish
#endif
#ifndef NOMPI
  USE mo_decomposition, ONLY: gdc => global_decomposition, &
                              ldc => local_decomposition
#else
  USE mo_decomposition, ONLY: ldc => local_decomposition
#endif
  USE mo_tracdef,       ONLY: trlist
  USE mo_tr_gp_ffsl,    ONLY: pack_gp, unpack_ffsl, gp2ffsl, &
                              pack_ffsl, unpack_gp, ffsl2gp, &
                              init_tr_gp_ffsl
#ifdef _OPENMP
  USE omp_lib,          ONLY: omp_get_num_threads, omp_get_thread_num
#endif
#ifdef _PROFILE
  USE mo_profile,       ONLY: trace_start, trace_stop
#endif
  USE mo_control,       ONLY: lyaxt_transposition
  USE mo_echam_yaxt,    ONLY: yaxt_tr_gp_to_ffsl, yaxt_tr_ffsl_to_gp, &
       &                      yaxt_tr_tracer_gp_to_ffsl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_tpcore
  PUBLIC :: setup_tpcore
  PUBLIC :: cleanup_tpcore
  PUBLIC :: tpcore_transport
  PUBLIC :: tpcore_tendencies

  REAL(wp) ,ALLOCATABLE :: dtdx5(:), dtdx5b(:)
  REAL(wp) ,ALLOCATABLE :: dtdy5(:), dtdy5b(:)
  REAL(wp) ,ALLOCATABLE :: cosp(:)
  REAL(wp) ,ALLOCATABLE :: cose(:)
  REAL(wp) ,ALLOCATABLE ::  gw(:)
  REAL(wp) ,ALLOCATABLE :: rgw(:)

  REAL(wp) ,ALLOCATABLE :: clat(:)
  REAL(wp) ,ALLOCATABLE, TARGET :: ub(:,:,:), vb(:,:,:), fb(:,:,:,:), delp(:,:,:)
  REAL(wp) ,ALLOCATABLE, TARGET :: psb(:,:), ps1b(:,:), ps2b(:,:)
  REAL(wp) ,ALLOCATABLE :: as(:), bs(:)

  !-- fb in gridpoint space

  REAL(wp) ,ALLOCATABLE,TARGET :: ps_gp(:,:), delp_gp(:,:,:)
  REAL(wp) ,ALLOCATABLE,TARGET :: fb_gp(:,:,:,:)

  INTEGER :: plon, plat, plev, pcnst

  INTEGER :: kstart, kstep, klev

  INTEGER, PARAMETER :: kord = 7  ! Vertical mapping scheme
  INTEGER, PARAMETER :: iv = 1    ! Monotonicity constraints for top and bottom
                         ! Enforce strong constraint at top & bottom
                         ! May want to change to iv=0 if diffusion is a problem

                         ! monotonicity at top and bottom
                         ! iv=0 : weak constraint
                         ! iv=1 : strong constraint
                         ! iv=-1: for vector

                         ! iv =-1: winds
                         ! iv = 0: positive definite scalars
                         ! iv = 1: others

  ! possible gp<->ffsl transposition methods:
  INTEGER, PARAMETER :: original_transposition_method = 1
  INTEGER, PARAMETER :: packed_transposition_method   = 2
  INTEGER, PARAMETER :: yaxt_transposition_method     = 3

  ! set default transposition:
#ifdef _PACKED_TRANSPOSITION
  INTEGER, SAVE :: transposition_choice = packed_transposition_method
#else
  INTEGER, SAVE :: transposition_choice = original_transposition_method
#endif

  ! link xtm1(:,:,it,:) to fb(:,:,:,ib) using it=map_fb2ti(ib):
  INTEGER, ALLOCATABLE, SAVE :: map_fb2ti(:)

CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE init_tpcore
    !-------------------------------------------------------------------------

    USE mo_math_constants,     ONLY: pi
    USE mo_physical_constants, ONLY: earth_radius
    USE mo_gaussgrid,     ONLY: gauaw
    USE mo_control,       ONLY: vct, nvclev, ngl
    USE mo_tracdef,       ONLY: trlist
    USE mo_advection,     ONLY: jps, nalatd, nalond, nalev, nacnst, nalat, nalon

    !-----
    ! Local
    !-----

    REAL(wp) :: ae           ! Earth's radius (m)
    REAL(wp) :: zgw(ngl)     ! dummy for Gaussian weights
    REAL(wp) :: zgmu(ngl)    ! dummy for sin of Gaussian latitudes (N->S)

    REAL(wp) :: elat(ngl+1)  ! cell edge latitude in radian
    REAL(wp) :: sine(ngl+1)
    REAL(wp) :: dlat(ngl)    ! delta-latitude in radian
    REAL(wp) :: dlon

    INTEGER :: jfirst, jlast
    INTEGER :: j, ib, it

    jfirst = ldc%ffsl%lats
    jlast  = ldc%ffsl%latn

    plon  = ldc%nlon
    plat  = ldc%nglat
    plev  = ldc%nlev
    pcnst = jps + COUNT(trlist% ti(1:trlist% ntrac)% ntran /= 0)

    kstart = ldc%set_b
    kstep  = ldc%nprocb
    klev   = (plev-kstart) / kstep + 1

    nalev  = plev
    nacnst = pcnst
    nalat  = plat
    nalon  = plon
    nalatd = plat/2
    nalond = plon+1

    IF ( jlast - jfirst < 2) THEN
       CALL finish ('init_tpcore', 'Minimum size of subdomain is 3')
    ENDIF

    !----------------
    ! Allocate arrays
    !----------------

    ALLOCATE ( cosp(ldc%nlat) )
    ALLOCATE ( cose(ldc%nlat) )
    ALLOCATE (   gw(ldc%nlat) )
    ALLOCATE (  rgw(ldc%nlat) )

    ALLOCATE ( clat(ldc%nlat) )

    ALLOCATE (ub (plon, jfirst:jlast, klev))
    ALLOCATE (vb (plon, jfirst:jlast, klev))
    ALLOCATE (fb (plon, jfirst:jlast, klev, pcnst))
    ALLOCATE (delp (plon, jfirst:jlast, klev))

    ALLOCATE (ps_gp(ldc%nproma, ldc%ngpblks))
    ALLOCATE (fb_gp(ldc%nproma, plev, pcnst, ldc%ngpblks))
    ALLOCATE (delp_gp(ldc%nproma, plev, ldc%ngpblks))

    ALLOCATE (psb(plon,jfirst:jlast))
    ALLOCATE (ps1b(plon,jfirst:jlast))
    ALLOCATE (ps2b(plon,jfirst:jlast))

    ALLOCATE (as(nvclev) )
    ALLOCATE (bs(nvclev) )

    ALLOCATE ( dtdx5(ldc%nlat), dtdx5b(ldc%nlat) )
    ALLOCATE ( dtdy5(ldc%nlat), dtdy5b(ldc%nlat) )

    ALLOCATE(map_fb2ti(pcnst))
    IF (jps /= 3) CALL finish('init_tpcore', 'False assumption: jps==3')
    ib = jps
    map_fb2ti(1:ib) = 0 ! ib=1..jps does not map into trlist%ti
    DO it = 1, trlist%ntrac
      IF(trlist%ti(it)%ntran /= 0) THEN
        ib = ib + 1
        map_fb2ti(ib) = it
      ENDIF
    ENDDO
    IF (ib /= pcnst) CALL finish('init_tpcore', 'Internal error: (ib /= pcnst).')

    ae = earth_radius

    as(:) = vct(1:nvclev)
    bs(:) = vct(nvclev+1:2*nvclev)

    CALL gauaw (zgmu, zgw, ldc%nlat)

    DO j = 1, ldc%nlat
       clat(ldc%nlat+1-j) = ASIN(zgmu(j))
    ENDDO

    dlon = 2*pi / REAL(plon,wp)

    elat(1) = -0.5_wp*pi         ! S. Pole
    sine(1) = -1.0_wp
    cose(1) =  0.0_wp

    DO j = 2, ldc%nlat
       elat(j) = 0.5_wp*(clat(j-1) + clat(j))
       sine(j) = SIN(elat(j))
       cose(j) = COS(elat(j))
    ENDDO

    elat(ldc%nlat+1) = 0.5_wp*pi       ! N. Pole
    sine(ldc%nlat+1) = 1.0_wp

    dlat(1) = 2*(elat(2) - elat(1))  ! Polar cap
    DO j = 2, ldc%nlat-1
       dlat(j) = elat(j+1) - elat(j)
    ENDDO
    dlat(ldc%nlat) = 2*(elat(ldc%nlat+1) - elat(ldc%nlat))    ! Polar cap

    DO j = 1, ldc%nlat
       gw(j) = sine(j+1) - sine(j)  ! Gaussian weights                | gw
       cosp(j) = gw(j) / dlat(j)    ! SQRT(1 - sin(Gaussian lat)**2)  | sqcst
       rgw(j) =  1.0_wp / gw(j)     ! 1 / Gaussian weight             | 1/gw
       ! *dt put into setup_tpcore
       !    dtdx5(j) = 0.5*dt / (dlon*ae*cosp(j))
       !    dtdy5(j) = 0.5*dt / (ae*dlat(j))
       dtdx5b(j) = 0.5_wp / (dlon*ae*cosp(j))
       dtdy5b(j) = 0.5_wp / (ae*dlat(j))
    ENDDO

    IF (lyaxt_transposition) transposition_choice = yaxt_transposition_method

    SELECT CASE (transposition_choice)
    CASE (packed_transposition_method)
      ! ub, vb, fb(pcnst) -> 2+pcnst
      CALL init_tr_gp_ffsl(default_buffer_size=2+pcnst)
    CASE (yaxt_transposition_method, original_transposition_method)
      ! no initialization
    CASE DEFAULT
      CALL finish ('tpcore','unexpected transposition_choice')
    END SELECT

  END SUBROUTINE init_tpcore

  !-------------------------------------------------------------------------
  SUBROUTINE setup_tpcore
    !-------------------------------------------------------------------------
    USE mo_gaussgrid,     ONLY: gl_rsqcst
    USE mo_scan_buffer,   ONLY: u, v, alps
    USE mo_memory_g1a,    ONLY: qm1, xlm1, xim1, xtm1, alpsm1
    USE mo_time_control,  ONLY: time_step_len
    USE mo_transpose,     ONLY: tr_gp_ffsl

    INTEGER :: i, it, ib, j, jrow, k
    REAL(wp) :: ztmst

    INTEGER :: jfirst, jlast
    jfirst = ldc%ffsl%lats
    jlast  = ldc%ffsl%latn

!$OMP PARALLEL PRIVATE(ztmst,j,jrow,it,ib)

    ztmst = time_step_len

!$omp do
    DO j = 1, ldc%nlat
    dtdx5(j) = dtdx5b(j)*ztmst
    dtdy5(j) = dtdy5b(j)*ztmst
    ENDDO
!$omp end do nowait

!$OMP END PARALLEL

     SELECT CASE (transposition_choice)
     CASE (packed_transposition_method)
     CALL pack_gp(u)
     CALL pack_gp(v)
     CALL pack_gp(qm1)
     CALL pack_gp(xlm1)
     CALL pack_gp(xim1)

     CALL gp2ffsl
       
     CALL unpack_ffsl(ub)
     CALL unpack_ffsl(vb)
     CALL unpack_ffsl(fb(:,:,:,1))
     CALL unpack_ffsl(fb(:,:,:,2))
     CALL unpack_ffsl(fb(:,:,:,3))
     CASE (yaxt_transposition_method)
       CALL yaxt_tr_gp_to_ffsl(u,v,qm1,xlm1,xim1, &
            &                  alpsm1,alps, &
            &                  ub,vb,fb(:,:,:,1),fb(:,:,:,2),fb(:,:,:,3),  &
            &                  ps1b, ps2b)
     CASE (original_transposition_method)
      CALL tr_gp_ffsl (ldc ,1 ,u           ,ub          )
      CALL tr_gp_ffsl (ldc ,1 ,v           ,vb          )
      CALL tr_gp_ffsl (ldc ,1 ,qm1         ,fb(:,:,:,1 ))
      CALL tr_gp_ffsl (ldc ,1 ,xlm1        ,fb(:,:,:,2 ))
      CALL tr_gp_ffsl (ldc ,1 ,xim1        ,fb(:,:,:,3 ))
     END SELECT

#ifdef CHECK_TRANSPOSITION
     CALL check3d_gp2ffsl('u->ub', u, ub)
     CALL check3d_gp2ffsl('v->vb', v, vb)
     CALL check3d_gp2ffsl('qm1->fb1', qm1, fb(:,:,:,1))
     CALL check3d_gp2ffsl('xlm1->fb2', xlm1, fb(:,:,:,2))
     CALL check3d_gp2ffsl('xim1->fb3', xim1, fb(:,:,:,3))
#endif

     SELECT CASE (transposition_choice)
     CASE (yaxt_transposition_method)
       IF (trlist% ntrac > 0) THEN
         CALL yaxt_tr_tracer_gp_to_ffsl(xtm1, fb)
#ifdef CHECK_TRANSPOSITION
         ib = 3
         DO it = 1, trlist% ntrac
           IF(trlist% ti(it)% ntran /= 0) THEN
             ib = ib + 1
             CALL check3d_gp2ffsl('xtm1->fb', xtm1(:,:,it,:), fb(:,:,:,ib))
           ENDIF
         END DO
#endif
       ENDIF
     CASE DEFAULT
       DO ib = 4, pcnst
         it = map_fb2ti(ib)
         CALL tr_gp_ffsl (ldc ,1 ,xtm1(:,:,it,:), fb(:,:,:,ib))
       ENDDO
     END SELECT

     SELECT CASE (transposition_choice)
     CASE (packed_transposition_method)
      CALL pack_gp(alpsm1, with_exp=.FALSE.)
      CALL pack_gp(alps, with_exp=.FALSE.)
      CALL gp2ffsl
      CALL unpack_ffsl(ps1b)
      CALL unpack_ffsl(ps2b)
     CASE (yaxt_transposition_method)
       ! communication of alpsm1 and alps already done
     CASE (original_transposition_method)
      CALL tr_gp_ffsl (ldc ,1 ,alpsm1  ,ps1b)
      CALL tr_gp_ffsl (ldc ,1 ,alps    ,ps2b)
     END SELECT

#ifdef CHECK_TRANSPOSITION
     CALL check2d_gp2ffsl('alpsm1->ps1b', alpsm1, ps1b)
     CALL check2d_gp2ffsl('alps->ps2b', alps, ps2b)
#endif

!$OMP PARALLEL private(jrow,i,j,k)
!$OMP DO
       DO jrow = jfirst, jlast
         ps1b(:,jrow) = EXP(ps1b(:,jrow))
         ps2b(:,jrow) = EXP(ps2b(:,jrow))
       ENDDO
!$OMP END DO nowait

!$OMP DO
   DO k = 1, klev
    DO i = ldc%ffsl%lats, ldc%ffsl%latn
      ub(:,i,k) = ub(:,i,k) * gl_rsqcst(i)
      vb(:,i,k) = vb(:,i,k) * gl_rsqcst(i)
    END DO
  ENDDO
!$OMP END DO nowait

!$OMP DO
  DO j = jfirst, jlast
    psb (:,j) = 0.0_wp
  ENDDO
!$OMP END DO

!$OMP END PARALLEL

!=============================================================================
! Original copy changed into transpose above:
!
!    DO i = 1, plat
!      l = plat-i+1  ! S->N latitude index
!      jglat  = ldc%glat(l)   ! global index (continuous from north to south)
!      igprow = MIN(2*jglat-1,2*(plat+1-jglat))   ! global ping pong index
!      zrcst  = 1.0_wp/sqcst(igprow)                 ! and 1./cos(latitude)
!
!       ! arrays must have the format south -> north
!
!       DO jk = 1, plev
!          DO jl = 1, plon
!             ub(jl,i,jk)   = u(jl,jk,l)*zrcst
!             vb(jl,i,jk)   = v(jl,jk,l)*zrcst
!             fb(jl,i,jk,1) = qm1(jl,jk,l)
!             fb(jl,i,jk,2) = xlm1(jl,jk,l)
!             fb(jl,i,jk,3) = xim1(jl,jk,l)
!          END DO
!       END DO
!
!       DO jt = 1, ntrac
!          DO jk = 1, plev
!             DO jl = 1, plon
!                fb(jl,i,jk,jps+jt) = xtm1(jl,jk,jt,l)
!             END DO
!          END DO
!       END DO
!
!       DO jl = 1, plon
!          ps1b(jl,i) = EXP(alpsm1(jl,l))
!          ps2b(jl,i) = EXP(alpsm1(jl,l))+EXP(alpste(jl,l))*ztmst
!       END DO
!
!    END DO
!=============================================================================

  END SUBROUTINE setup_tpcore

  !-------------------------------------------------------------------------
  SUBROUTINE cleanup_tpcore
    !-------------------------------------------------------------------------

    DEALLOCATE ( dtdy5b, dtdy5 )
    DEALLOCATE ( dtdx5b, dtdx5 )

    DEALLOCATE (   bs )
    DEALLOCATE (   as )

    DEALLOCATE ( ps2b )
    DEALLOCATE ( ps1b )
    DEALLOCATE (  psb )

    DEALLOCATE (   fb )
    DEALLOCATE (   vb )
    DEALLOCATE (   ub )
    DEALLOCATE (   delp )

    DEALLOCATE (ps_gp )
    DEALLOCATE (fb_gp )
    DEALLOCATE (delp_gp )

    DEALLOCATE ( clat )

    DEALLOCATE (  rgw )
    DEALLOCATE (   gw )
    DEALLOCATE ( cose )
    DEALLOCATE ( cosp )

  END SUBROUTINE cleanup_tpcore

  !------------------------------------------------------------
  SUBROUTINE tpcore_transport
    !------------------------------------------------------------

    USE mo_physical_constants,     ONLY: earth_radius
    USE mo_time_control,  ONLY: time_step_len

    CALL tpcore (time_step_len, earth_radius, ldc%nlon, ldc%nlat, ldc%nlev, &
         ldc%ffsl%lats, ldc%ffsl%latn, &
         pcnst, as, bs, ub, vb, ps1b, ps2b, psb, fb, 4, 4, 0)

  END SUBROUTINE tpcore_transport

  !------------------------------------------------------------
  SUBROUTINE tpcore_tendencies(psm1cor, pscor)
    !------------------------------------------------------------

#ifdef DEBUG
    USE mo_gaussgrid,     ONLY: gl_budw
#endif
    USE mo_time_control,  ONLY: time_step_len
    USE mo_scan_buffer,   ONLY: qte, xlte, xite, xtte, alps
    USE mo_memory_g1a,    ONLY: qm1, xlm1, xim1, xtm1, alpsm1
#ifdef DEBUG
    USE mo_transpose,     ONLY: tr_gp_ffsl
#endif

    REAL(wp), INTENT(in) :: pscor
    REAL(wp), INTENT(in) :: psm1cor

#ifdef DEBUG
    REAL(wp) :: zpsz(plat), zps1z(plat), zps2z(plat)
    REAL(wp) :: hwps, hwps1, hwps2

    INTEGER  :: jlat
#endif

    INTEGER  :: it, ib
    REAL(wp) :: zrtmst, zlpsm1c, zlpsc
    INTEGER  :: jrow, jl, nproma

    !  Executable statements

!$OMP PARALLEL PRIVATE(zrtmst,zlpsm1c,zlpsc,it,ib,jrow,jl,nproma)

    zrtmst = 1.0_wp/time_step_len

    zlpsm1c = LOG(psm1cor)
    zlpsc = LOG(pscor)

!$OMP DO ECHAM_GPC_SCHEDULE
    DO jrow = 1, ldc%ngpblks
    qte(:,:,jrow)  = qte(:,:,jrow)  + (fb_gp(:,:,1,jrow) -  qm1(:,:,jrow)) * zrtmst
    xlte(:,:,jrow) = xlte(:,:,jrow) + (fb_gp(:,:,2,jrow) - xlm1(:,:,jrow)) * zrtmst
    xite(:,:,jrow) = xite(:,:,jrow) + (fb_gp(:,:,3,jrow) - xim1(:,:,jrow)) * zrtmst
    DO ib = 4, pcnst
      it = map_fb2ti(ib)
      xtte(:,:,it,jrow) = xtte(:,:,it,jrow) &
                        + (fb_gp(:,:,ib,jrow)-xtm1(:,:,it,jrow)) * zrtmst
    ENDDO

      IF ( jrow == ldc% ngpblks ) THEN
        nproma = ldc% npromz
      ELSE
        nproma = ldc% nproma
      END IF
      
      ! fix mass of air of spectral dynamics

      DO jl = 1, nproma
        alpsm1(jl,jrow) = alpsm1(jl,jrow)+zlpsm1c
        alps(jl,jrow) = alps(jl,jrow)+zlpsc
      END DO
      
    END DO
!$OMP END DO

!$OMP END PARALLEL

#ifdef DEBUG
    ! do checking of global mass

    CALL tr_gp_ffsl(ldc , 1, ps_gp, psb)

    DO jlat = 1, plat
       zps1z(jlat) = gl_budw(jlat)*SUM(ps1b(:,jlat))
       zps2z(jlat) = gl_budw(jlat)*SUM(ps2b(:,jlat))
       zpsz(jlat)  = gl_budw(jlat)*SUM(psb(:,jlat))
    END DO

    hwps1 = SUM(zps1z(:))
    hwps2 = SUM(zps1z(:))
    hwps  = SUM(zpsz(:))


    message_text = ''
    WRITE (message_text,'(3(a,f10.3))') &
         ' PS(spectral): ',     hwps2, &
         ' PS(spectral,old): ', hwps1, &
         ' PS(tpcore): ',       hwps
    CALL message('tpcore',TRIM(message_text))
#endif

  END SUBROUTINE tpcore_tendencies

  !------------------------------------------------------------
  SUBROUTINE tpcore(dt, ae, im, jm, km, jfirst, jlast, nq,   &
                    ak, bk, u, v, ps1, ps2, ps,  q,          &
                    iord, jord, n_adj)
    !------------------------------------------------------------

    USE mo_transpose,     ONLY: tr_gp_ffsl
!    USE mo_transpose,     ONLY: reorder
    USE mo_decomposition, ONLY: ldc => local_decomposition

    INTEGER, INTENT(in):: im         ! Global E-W dimension
    INTEGER, INTENT(in):: jm         ! Global N-S dimension
    INTEGER, INTENT(in):: km         ! Vertical dimension
    INTEGER, INTENT(in):: jfirst     ! Local first index for N-S
    INTEGER, INTENT(in):: jlast      ! Local last index for N-S
    INTEGER, INTENT(in):: nq         ! Ghosted latitudes (3 required by PPM)
    INTEGER, INTENT(in):: iord       ! E-W transport scheme
    INTEGER, INTENT(in):: jord       ! N-S transport scheme
    INTEGER, INTENT(in):: n_adj      ! Number of adjustemnt to air_mass_flux
                                     ! 0 --> no adjustment

    ! Recommended values : iord=jord=4, kord=7
    !  _ord:
    !---------------------------------------------------------------------------
    !        1: 1st order upstream scheme
    !        2: 2nd order van Leer (full monotonicity constraint;
    !           see Lin et al 1994, MWR)
    !        3: Standard monotonic PPM* (Collela & Woodward 1984)
    !        4: New & Improved monotonic PPM
    !        5: positive-definite PPM (constraint on the subgrid distribution is
    !           only strong enough to prevent generation of negative values;
    !           both overshoots & undershootes are possible).
    !        6: un-constrained PPM (nearly diffusion free; faster but
    !           positivity of the subgrid distribution is not quaranteed.
    !        7: Huynh/Van Leer/Lin full monotonicity constraint
    !---------------------------------------------------------------------------
    ! Only kord can be set to 7 to enable the use of Huynh's 2nd monotonicity
    ! constraint for piece-wise parabolic distribution.
    ! *PPM: Piece-wise Parabolic Method

    REAL(wp), INTENT(in):: ak(km+1)                ! See below
    REAL(wp), INTENT(in):: bk(km+1)                ! See below
    REAL(wp), INTENT(in):: u(im,jfirst:jlast,klev) ! u-wind (m/s) at mid-time-level (t=t+dt/2)
    REAL(wp), INTENT(in):: v(im,jfirst:jlast,klev) ! v-wind (m/s) at mid-time-level (t=t+dt/2)

    !------------------------------------------------------
    ! The hybrid ETA-coordinate:
    ! pressure at layer edges are defined as follows:
    !
    !        p(i,j,k) = ak(k) + bk(k)*ps(i,j)
    !------------------------------------------------------
    ! ak and bk are defined at layer edges.
    !
    !                  /////////////////////////////////
    !              / \ ------ Model top P=ak(1) --------- ak(1), bk(1)
    !               |
    !    delp(1)    |  ........... q(i,j,1) ............
    !               |
    !              \ / ---------------------------------  ak(2), bk(2)
    !
    !
    !
    !              / \ ---------------------------------  ak(k), bk(k)
    !               |
    !    delp(k)    |  ........... q(i,j,k) ............
    !               |
    !              \ / ---------------------------------  ak(k+1), bk(k+1)
    !
    !
    !
    !              / \ ---------------------------------  ak(km), bk(km)
    !               |
    !    delp(km)   |  ........... q(i,j,km) .........
    !               |
    !              \ / -----Earth's surface P=Psfc ------ ak(km+1), bk(km+1)
    !                 //////////////////////////////////
    !
    ! Note: surface pressure can be of any unit (e.g., pascal or mb) as long as it is
    ! consistent with the definition of (ak, bk) defined above
    ! Winds (u,v), ps, and q are assumed to be defined at the same points.
    ! The latitudes are given by clat, input to the initialization routine: init_tpcore.

    REAL(wp), INTENT(in):: ps1(im,jfirst:jlast)  ! surface pressure at current time
    REAL(wp), INTENT(in):: ps2(im,jfirst:jlast)  ! surface pressure at future time=t+dt
    REAL(wp), INTENT(in):: dt                    ! Transport time step in seconds
    REAL(wp), INTENT(in):: ae                    ! Earth's radius (m)

    REAL(wp), INTENT(inout):: q(im,jfirst:jlast,klev,nq)  ! Tracer "mixing ratios"
                                                          ! q could easily be re-dimensioned
    REAL(wp), INTENT(inout):: ps(im,jfirst:jlast)         ! "predicted" surface pressure

    !------
    ! Local
    !------

#ifndef NOMPI
    INTEGER, PARAMETER :: ng = 3      ! Primary ghost zones in N-S direction
                                      ! -- tracres and u-wind
    INTEGER, PARAMETER :: mg = 1      ! Secondary ghost zones in N-S direction
                                      ! Meridional flux and surface pressure
#else
    INTEGER, PARAMETER :: ng = 0      ! No ghosting needed!
    INTEGER, PARAMETER :: mg = 0      ! No ghosting needed!
#endif

    ! REAL(wp) :: delp(im,jfirst:jlast,klev)   ! Predicted thickness at future time (t=t+dt)

    REAL(wp) :: fx(im,jfirst:jlast,klev)     ! E-W air mass flux
    REAL(wp) :: va(im,jfirst:jlast,klev)     ! N-S CFL at cell center (scalar points)

    !-----------------------
    ! Ghosted local arrays:
    !-----------------------


    REAL(wp) :: delp1(im,jfirst-ng:jlast+ng,klev) ! Pressure thickness at current time (t)
    REAL(wp) :: fy(im,jfirst:jlast+mg,klev)       ! N-S air mass flux
    REAL(wp) :: cx(im,jfirst-ng:jlast+ng,klev)    ! E-W CFL number on C-grid
    REAL(wp) :: cy(im,jfirst:jlast+mg,klev)       ! N-S CFL number on C-grid
    REAL(wp) :: psm(im,jfirst-ng:jlast+ng)
    REAL(wp) :: psn(im,jfirst-mg:jlast+mg)
    REAL(wp) :: q2(im,jfirst-ng:jlast+ng,klev,nq) ! local 2D q array

    ! REAL(wp) :: delp_gp(ldc%nproma,km,ldc%ngpblks) ! delp in gridpoint space

    LOGICAL :: ffsl(jfirst-ng:jlast+ng,klev)      ! Flag to compute Integer fluxes

    REAL(wp) :: yms(im,jfirst:jlast+mg,klev)
    REAL(wp) :: v3(im,jfirst-mg:jlast+mg,klev)
    
    ! Local variables:

    INTEGER :: i,j,k,iq
    INTEGER :: iord_bg              ! E-W scheme for background mass flux
    INTEGER :: jord_bg              ! N-S scheme for background mass flux

#ifdef VECTOR
    REAL(wp) :: qtmp((jlast-jfirst+1)*klev*nq+1,im)

    INTEGER :: ipx
#endif

    iord_bg = 3
    jord_bg = 3

!$OMP PARALLEL PRIVATE(i,j,k,iq)

    ! Ensure inputs are single-valued at poles:

!$OMP DO
    DO j=jfirst,jlast
      DO i=1,im
        psm(i,j) = ps1(i,j)
        psn(i,j) = ps2(i,j)
      ENDDO
    ENDDO
!$OMP END DO

!$OMP MASTER
    IF ( jfirst == 1 ) THEN
      CALL xpavg(psm(1,1), im)
      CALL xpavg(psn(1,1), im)
    ENDIF

    IF ( jlast == jm ) THEN
      CALL xpavg(psm(1,jm), im)
      CALL xpavg(psn(1,jm), im)
    ENDIF

#ifndef NOMPI
    ! Ghost psm and psn north/south
    CALL ghost_update(psm, im, jm, 1, 1, jfirst, jlast, ng)
    CALL ghost_update(psn, im, jm, 1, 1, jfirst, jlast, mg)
#endif

!$OMP END MASTER
!$OMP BARRIER

    ! Average q at both poles
    DO iq=1,nq
!$OMP DO
      DO k=1,klev
        IF ( jfirst == 1 ) THEN
          CALL xpavg(q(1,1,k,iq), im)
        ENDIF
        IF ( jlast == jm ) THEN
          CALL xpavg(q(1,jm,k,iq), im)
        ENDIF
      ENDDO
!$OMP END DO NOWAIT
    ENDDO
!$OMP BARRIER

    !----------------------------------------------
    ! Compute background air mass fluxes
    !----------------------------------------------

    CALL air_mass_flux(im, jm, km, jfirst, jlast,      &
                       iord_bg, jord_bg,   ak, bk,     &
                       psm, psn, ps,  u, v,            &
                       cx, cy, va, fx, fy, ng, mg,     &
                       ffsl, delp1, delp, dt,     &
                       ae, n_adj, yms, v3)

#ifdef DEBUG
!$OMP MASTER
    IF (ANY(cy >= 1.0_wp)) THEN
      message_text = ''
      WRITE (message_text,*) &
           ' Courant number in North-South is exceeding 1 ... ', MAXVAL(cy)
      CALL message('tpcore',TRIM(message_text))
      CALL finish ('tpcore','No save transport possible')
    ELSE
      message_text = ''
      WRITE (message_text,'(a,f6.2,a,f6.2)') &
           ' Max Courant number in East-West ', MAXVAL(cx), &
           ' and North-South ', MAXVAL(cy)
      CALL message('tpcore',TRIM(message_text))
    END IF
!$OMP END MASTER
!$OMP BARRIER
#endif

    !---------------------------------------------------
    ! Do tracer transport
    !---------------------------------------------------

    DO iq=1,nq
!$OMP DO
      DO k=1,klev
        DO j=jfirst,jlast
          DO i=1,im
            q2(i,j,k,iq) = q(i,j,k,iq)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
    ENDDO
!$OMP BARRIER

#ifndef NOMPI
!$OMP MASTER
    CALL ghost_update(q2, im, jm, klev, nq, jfirst, jlast, ng)
!$OMP END MASTER
!$OMP BARRIER
#endif

    Multi_Tracer: DO iq=1,nq
!$OMP DO
    Vertical:  DO k=1,klev

          CALL tp2g(q2(1,jfirst-ng,k,iq),    va(1,jfirst,k),     &
               cx(1,jfirst-ng,k),  cy(1,jfirst,k),          &
               im,  jm,  iord,     jord,                    &
               ng,  mg,  fx(1,jfirst,k), fy(1,jfirst,k),    &
               ffsl(jfirst-ng,k),    jfirst,   jlast,       &
               delp1(1,jfirst-mg,k),    delp(1,jfirst,k) )
          
          DO j=jfirst,jlast
             DO i=1,im
                q(i,j,k,iq) = q2(i,j,k,iq) ! q2 has ghost cells, q hasn't
             END DO
          END DO
          
       ENDDO Vertical
!$OMP END DO
    ENDDO Multi_Tracer
   
#ifdef VECTOR
    
    !-------------------------------------------------------------------
    ! Apply a simple nearest neighbor flux correction to reduce negatives
    ! fct_x calculation in non-vector case in tp2g (original version)
    !-------------------------------------------------------------------
    IF ( iv /= -1 ) THEN
      CALL fct_x(q, im, jm, jfirst, jlast, klev, nq, ipx, qtmp)
    ENDIF
    
!$OMP DO
    DO k=1,klev
       DO iq=1,nq
          DO j=jfirst,jlast
             DO i=1,im
                q(i,j,k,iq) = q(i,j,k,iq) / delp(i,j,k)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!$OMP END DO
    
#endif
!$OMP END PARALLEL

    !---------------------------------------------------------------
    ! Perform Remapping back to the hybrid sigma-pressure coordinate
    ! Mass will be conserved if predicted ps2 == psn (data/model)
    !---------------------------------------------------------------

    !   Transform back to gridpoint space

    SELECT CASE (transposition_choice)
    CASE (packed_transposition_method)
        CALL pack_ffsl(q)
        CALL pack_ffsl(delp)
!lk check for ps_gp ...
        CALL ffsl2gp
        CALL unpack_gp(fb_gp)
        CALL unpack_gp(delp_gp)
    CASE (yaxt_transposition_method)
      CALL yaxt_tr_ffsl_to_gp(fb_gp, ps_gp, delp_gp, &
           &                  q, ps, delp)
    CASE (original_transposition_method)
      CALL tr_gp_ffsl(ldc ,-1 ,fb_gp, q)
      CALL tr_gp_ffsl(ldc ,-1 ,ps_gp, ps)
      CALL tr_gp_ffsl(ldc ,-1 ,delp_gp, delp)
    END SELECT

#ifdef CHECK_TRANSPOSITION
    CALL check4d_ffsl2gp('q->fb_gp', q, fb_gp)
    IF (transposition_choice == yaxt_transposition_method) THEN
      CALL check2d_ffsl2gp('ps->ps_gp', ps, ps_gp)
    ENDIF
    CALL check3d_ffsl2gp('delp->delp_gp', delp, delp_gp)
#endif

!$OMP PARALLEL
   IF ( km >= 5 ) THEN
!lk either reorder or use the commented CALL tr_gp_ffsl 6 lines above ...
!lk      CALL reorder(ps_gp, psb)
      CALL qmap_gp(ps_gp, delp_gp, fb_gp, ldc%nproma, ldc%npromz, ldc%ngpblks, &
                   km, nq, ak, bk)
    ENDIF
!$OMP END PARALLEL

  END SUBROUTINE tpcore

#ifdef CHECK_TRANSPOSITION
  SUBROUTINE check2d_gp2ffsl(label, gp, ffsl, with_exp)
    USE mo_transpose,     ONLY: tr_gp_ffsl
    CHARACTER(len=*), INTENT(in) :: label
    REAL(wp), INTENT(in) :: gp(:,:)
    REAL(wp), INTENT(in) :: ffsl(:,:)
    LOGICAL, OPTIONAL, INTENT(in) :: with_exp


    REAL(wp) :: ref_gp(SIZE(gp,1), SIZE(gp,2))
    REAL(wp) :: ref_ffsl(SIZE(ffsl,1), SIZE(ffsl,2))
    LOGICAL :: lexp

    IF (PRESENT(with_exp)) THEN
      lexp = with_exp
    ELSE
      lexp = .FALSE.
    ENDIF

    ref_gp = gp
    CALL tr_gp_ffsl(ldc , 1 ,ref_gp, ref_ffsl)
    IF (lexp) ref_ffsl = EXP(ref_ffsl)

    IF (ANY(ffsl /= ref_ffsl)) CALL finish('mo_tpcore::check2d_gp2ffsl', label)
    IF (ldc%pe == 0) CALL message('mo_tpcore::check2d_gp2ffsl:',label//' passed')
  END SUBROUTINE check2d_gp2ffsl

  SUBROUTINE check3d_gp2ffsl(label, gp, ffsl)
    USE mo_transpose,     ONLY: tr_gp_ffsl
    CHARACTER(len=*), INTENT(in) :: label
    REAL(wp), INTENT(in) :: gp(:,:,:)
    REAL(wp), INTENT(in) :: ffsl(:,:,:)

    REAL(wp) :: ref_gp(SIZE(gp,1), SIZE(gp,2), SIZE(gp,3))
    REAL(wp) :: ref_ffsl(SIZE(ffsl,1), SIZE(ffsl,2), SIZE(ffsl,3))

    ref_gp = gp
    CALL tr_gp_ffsl(ldc , 1 ,ref_gp, ref_ffsl)

    IF (ANY(ffsl /= ref_ffsl)) CALL finish('mo_tpcore::check3d_gp2ffsl', label)
    IF (ldc%pe == 0) CALL message('mo_tpcore::check3d_gp2ffsl:',label//' passed')
  END SUBROUTINE check3d_gp2ffsl

  SUBROUTINE check2d_ffsl2gp(label, ffsl, gp)
    USE mo_transpose,     ONLY: tr_gp_ffsl
    CHARACTER(len=*), INTENT(in) :: label
    REAL(wp), INTENT(in) :: ffsl(:,:)
    REAL(wp), INTENT(in) :: gp(:,:)

    REAL(wp) :: ref_ffsl(SIZE(ffsl,1), SIZE(ffsl,2))
    REAL(wp) :: ref_gp(SIZE(gp,1), SIZE(gp,2))

    ref_ffsl = ffsl
    CALL tr_gp_ffsl(ldc ,-1 ,ref_gp, ref_ffsl)

    IF (ANY(gp /= ref_gp)) CALL finish('mo_tpcore::check2d_ffsl2gp', label)
    IF (ldc%pe == 0) CALL message('mo_tpcore::check2d_ffsl2gp:',label//' passed')
  END SUBROUTINE check2d_ffsl2gp

  SUBROUTINE check3d_ffsl2gp(label, ffsl, gp)
    USE mo_transpose,     ONLY: tr_gp_ffsl
    CHARACTER(len=*), INTENT(in) :: label
    REAL(wp), INTENT(in) :: ffsl(:,:,:)
    REAL(wp), INTENT(in) :: gp(:,:,:)

    REAL(wp) :: ref_ffsl(SIZE(ffsl,1), SIZE(ffsl,2), SIZE(ffsl,3))
    REAL(wp) :: ref_gp(SIZE(gp,1), SIZE(gp,2), SIZE(gp,3))

    ref_ffsl = ffsl
    CALL tr_gp_ffsl(ldc ,-1 ,ref_gp, ref_ffsl)

    IF (ANY(gp /= ref_gp)) CALL finish('mo_tpcore::check3d_ffsl2gp', label)
    IF (ldc%pe == 0) CALL message('mo_tpcore::check3d_ffsl2gp:',label//' passed')
  END SUBROUTINE check3d_ffsl2gp

  SUBROUTINE check4d_ffsl2gp(label, ffsl, gp)
    USE mo_transpose,     ONLY: tr_gp_ffsl
    CHARACTER(len=*), INTENT(in) :: label
    REAL(wp), INTENT(in) :: ffsl(:,:,:,:)
    REAL(wp), INTENT(in) :: gp(:,:,:,:)

    REAL(wp) :: ref_ffsl(SIZE(ffsl,1), SIZE(ffsl,2), SIZE(ffsl,3), SIZE(ffsl,4))
    REAL(wp) :: ref_gp(SIZE(gp,1), SIZE(gp,2), SIZE(gp,3), SIZE(gp,4))

    ref_ffsl = ffsl
    CALL tr_gp_ffsl(ldc ,-1 ,ref_gp, ref_ffsl)

    IF (ANY(gp /= ref_gp)) CALL finish('mo_tpcore::check4d_ffsl2gp', label)
    IF (ldc%pe == 0) CALL message('mo_tpcore::check4d_ffsl2gp:', label//' passed')
  END SUBROUTINE check4d_ffsl2gp
#endif

  SUBROUTINE air_mass_flux(im, jm, km, jfirst, jlast, iord, jord,    &
                           ak, bk, ps1, ps2, ps, u, v, cx, cy, va,   &
                           fx, fy, ng,  mg,  ffsl, delp1,  delp,     &
                           dt, ae,  n_adj, yms, v3)
    !------------------------------------------------------
    ! The hybrid coordinate systems pressure at layer edges
    ! is defined as follows:
    !
    !        p(i,j,k) = ak(k) + bk(k)*ps(i,j)          (1)
    !------------------------------------------------------
    !
    ! Input from Data/Model:
    ! (u,v) is the time mean wind at Time=t+dt/2
    ! delp1 is the layer thickness at Time=t
    !
    ! Output:
    ! delp is the predicted thickness at Time=t+dt
    ! (fx,fy): background air mass flxues
    ! (cx,cy): CFL number

    INTEGER, INTENT(in):: im
    INTEGER, INTENT(in):: jm
    INTEGER, INTENT(in):: km
    INTEGER, INTENT(in):: jfirst
    INTEGER, INTENT(in):: jlast
    INTEGER, INTENT(in):: iord
    INTEGER, INTENT(in):: jord
    INTEGER, INTENT(in):: ng
    INTEGER, INTENT(in):: mg
    INTEGER, INTENT(in):: n_adj

    REAL(wp), INTENT(in):: dt
    REAL(wp), INTENT(in):: ae
    REAL(wp), INTENT(in):: ak(km+1)
    REAL(wp), INTENT(in):: bk(km+1)
    REAL(wp), INTENT(in):: ps1(im,jfirst-ng:jlast+ng)
    REAL(wp), INTENT(in):: ps2(im,jfirst-mg:jlast+mg)
    REAL(wp), INTENT(in):: u(im,jfirst:jlast,klev)
    REAL(wp), INTENT(in):: v(im,jfirst:jlast,klev)

    LOGICAL,INTENT(out):: ffsl(jfirst-ng:jlast+ng,klev)

    REAL(wp), INTENT(out):: cx(im,jfirst-ng:jlast+ng,klev)
    REAL(wp), INTENT(out):: delp (im,jfirst:jlast,klev)
    REAL(wp), INTENT(out):: fx(im,jfirst:jlast,klev)
    REAL(wp), INTENT(out):: cy(im,jfirst:jlast+mg,klev)
    REAL(wp), INTENT(out):: fy(im,jfirst:jlast+mg,klev)
    REAL(wp), INTENT(out):: va(im,jfirst:jlast,klev)

    REAL(wp), INTENT(out):: delp1(im,jfirst-ng:jlast+ng,klev)

    REAL(wp), INTENT(inout):: ps(im,jfirst:jlast)

    REAL(wp),INTENT(inout) :: yms(im,jfirst:jlast+mg,klev)
    REAL(wp),INTENT(inout) :: v3(im,jfirst-mg:jlast+mg,klev)

    ! Local:
    REAL(wp), PARAMETER  :: tiny = 1.0e-10_wp

    REAL(wp) :: dak, dbk
    REAL(wp) :: dtoa, vt, upd

    INTEGER :: i,j,k,k2
    INTEGER :: js2g0
    INTEGER :: jn2g0
    INTEGER :: js2gd, jn2gd

    js2g0  = MAX(2,jfirst)        ! No ghosting
    jn2g0  = MIN(jm-1,jlast)      ! No ghosting
    js2gd = MAX(2,  jfirst-ng)    ! NG latitudes on S (starting at 2)
    jn2gd = MIN(jm-1,jlast+ng)    ! NG latitudes on S (ending at jm-1)

    dtoa = 0.5_wp*dt/ae

    ! set to zero to prevent diagnostic problems - parts of cx and cy 
    ! are not set

!$OMP DO
    DO k=1,klev
       cx(:,:,k) = 0.0_wp
       cy(:,:,k) = 0.0_wp
       v3(:,jfirst:jlast,k) = v(:,:,k)
    END DO
!$OMP END DO

#ifndef NOMPI
!$OMP MASTER
    CALL ghost_update(v3, im, jm, klev, 1, jfirst, jlast, mg)
!$OMP END MASTER
!$OMP BARRIER
#endif

!$OMP DO
    DO k=1,klev
      DO j=MAX(jfirst,2), MIN(jlast+1,jm)
        DO i=1,im
          vt = v3(i,j,k) + v3(i,j-1,k)
          upd = FSEL(vt,dtdy5(j-1),dtdy5(j))
          cy(i,j,k)  = upd*vt
          yms(i,j,k) = dtoa*vt*cose(j)
        ENDDO
      ENDDO

      DO j=js2g0,jn2g0
        DO i=1,im
           upd = FSEL(cy(i,j,k),cy(i,j,k),cy(i,j+1,k))
           va(i,j,k) = FSEL(cy(i,j,k)*cy(i,j+1,k), upd, 0.0_wp)
        ENDDO
        cx(1,j,k) = dtdx5(j)*(u(1,j,k)+u(im,j,k))
        DO i=2,im
          cx(i,j,k) = dtdx5(j)*(u(i,j,k)+u(i-1,j,k))
        ENDDO
      ENDDO
    ENDDO ! DO k=1,klev
!$OMP END DO

#ifndef NOMPI
!$OMP MASTER
    CALL ghost_update(cx(1,jfirst-ng,1), im, jm, klev, 1, jfirst, jlast, ng)
!$OMP END MASTER
!$OMP BARRIER
#endif

!$OMP DO
    DO k = 1, klev
      DO j=js2gd,jn2gd                ! ffsl needed on N*ng S*ng
        ffsl(j,k) = .FALSE.
        DO i=1,im
          IF( ABS(cx(i,j,k)) > 1.0_wp ) THEN
            ffsl(j,k) = .TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDDO

      k2 = (k-1)*kstep+kstart
      dak = ak(k2+1) - ak(k2)
      dbk = bk(k2+1) - bk(k2)

      DO j=MAX(1,jfirst-ng),MIN(jm,jlast+ng)
        DO i=1,im
          delp1(i,j,k) = dak + dbk*ps1(i,j)
        ENDDO
      ENDDO

      CALL tp2d(va(1,jfirst,k), delp1(1,jfirst-ng,k), cx(1,jfirst-ng,k), &
           cy(1,jfirst,k), im, jm, iord, jord, ng,  mg,             &
           fx(1,jfirst,k), fy(1,jfirst,k),  ffsl(jfirst-ng,k),      &
           cx(1,jfirst,k), yms(1,jfirst,k), 0, jfirst, jlast)

!DIR$ CONCURRENT
      DO j=js2g0,jn2g0
        DO i=1,im-1
          delp(i,j,k) = delp1(i,j,k) + fx(i,j,k) - fx(i+1,j,k) +          &
               (fy(i,j,k)-fy(i,j+1,k))*rgw(j)
        ENDDO
        delp(im,j,k) = delp1(im,j,k) + fx(im,j,k) - fx(1,j,k) +         &
             (fy(im,j,k)-fy(im,j+1,k))*rgw(j)
      ENDDO

      IF ( jfirst ==  1 ) THEN
        DO i=1,im
          delp(i,1,k) = delp1(i,1,k) - fy(i,2,k)*rgw(1)
        ENDDO
        CALL xpavg(delp(1,1,k), im)
      ENDIF

      IF ( jlast == jm ) THEN
        DO i=1,im
          delp(i,jm,k) = delp1(i,jm,k) + fy(i,jm,k)*rgw(jm)
        ENDDO
        CALL xpavg(delp(1,jm,k), im)
      ENDIF

      IF ( n_adj == 0 ) THEN
        DO j=js2g0,jn2g0
          IF( ffsl(j,k) ) THEN
            DO i=1,im
              fx(i,j,k) = fx(i,j,k)/SIGN(MAX(ABS(cx(i,j,k)),tiny),cx(i,j,k))
            ENDDO
          ENDIF
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO

    !-- RJ: The final computation of ps is done in qmap_gp

    !--------------------------------------------------------------
    ! Apply mass_flux adjuster to nudge predicted ps towards "data"
    !--------------------------------------------------------------

    IF ( n_adj > 0 ) THEN

    ! Compute ps for adj_fx

#ifndef NOMPI
!$OMP MASTER
      CALL bcast_delp(delp)
!$OMP END MASTER
!$OMP BARRIER
#endif

!$OMP DO
      DO j=jfirst,jlast
        DO i=1,im
          ps(i,j) = ak(1)
        ENDDO

        DO k=1,klev
          DO i=1,im
            ps(i,j) = ps(i,j) + delp(i,j,k)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

      CALL adj_fx(im, jm, km, jfirst, jlast, ak, bk, ffsl,  &
           ps, ps2, delp, fx, cx, fy, ng, mg,    &
           tiny, n_adj)
    ENDIF

  END SUBROUTINE air_mass_flux

  SUBROUTINE tp2g(h,  va, crx, cry, im, jm,              &
                  iord, jord, ng, mg, xfx, yfx, ffsl,    &
                  jfirst, jlast, dp, dpp)

    INTEGER, INTENT(in):: im, jm             ! Dimensions
    INTEGER, INTENT(in):: jfirst, jlast      ! Latitude strip
    INTEGER, INTENT(in):: iord, jord         ! Interpolation order in x,y
    INTEGER, INTENT(in):: ng                 ! Max. NS dependencies
    INTEGER, INTENT(in):: mg                 ! Secondary ghosting zones
    LOGICAL, INTENT(in):: ffsl(jfirst-ng:jlast+ng)  ! Use flux-form semi-Lagrangian trans.?
    REAL(wp), INTENT(in):: va(im,jfirst:jlast)   ! CFL in y at cell center
    REAL(wp), INTENT(in):: dp(im,jfirst-mg:jlast+mg)
    REAL(wp), INTENT(in):: dpp(im,jfirst:jlast)

    REAL(wp), INTENT(in):: crx(im,jfirst-ng:jlast+ng) ! ( N*NG S*NG )
    REAL(wp), INTENT(in):: cry(im,jfirst:jlast+mg)    ! ( N like FY )

    REAL(wp), INTENT(in):: xfx(im,jfirst:jlast)       ! x-mass flux
    REAL(wp), INTENT(in):: yfx(im,jfirst:jlast+mg)     ! y-mass flux

    REAL(wp), INTENT(inout):: h(im,jfirst-ng:jlast+ng)

    ! Local
    REAL(wp) :: fx(im,jfirst:jlast)        ! tracer flux in x ( unghosted )
    REAL(wp) :: fy(im,jfirst:jlast+mg)     ! tracer flux in y ( N, see tp2c )

    INTEGER :: i, j, js2g0, jn2g0
    INTEGER :: jlast_odd, jsize

    js2g0  = MAX(2,jfirst)          !  No ghosting
    jn2g0  = MIN(jm-1,jlast)        !  No ghosting

    CALL tp2d(va, h(1,jfirst-ng), crx(1,jfirst-ng), cry, im, jm,      &
              iord, jord, ng, mg, fx, fy, ffsl(jfirst-ng),          &
              xfx, yfx, 1, jfirst, jlast)

    DO j=js2g0,jn2g0
      DO i=1,im-1
        h(i,j) = h(i,j)*dp(i,j)+fx(i,j)-fx(i+1,j)+(fy(i,j)-fy(i,j+1))*rgw(j)
      ENDDO
    ENDDO

    DO j=js2g0,jn2g0
      h(im,j) = h(im,j)*dp(im,j)+fx(im,j)-fx(1,j)+(fy(im,j)-fy(im,j+1))*rgw(j)
    ENDDO

    ! Poles
    IF ( jfirst == 1 ) THEN
      DO i=1,im
        h(i,1) = h(i,1)*dp(i,1) - fy(i,2)*rgw(1)
      ENDDO
      CALL xpavg(h(1, 1), im)
    ENDIF

    IF ( jlast == jm ) THEN
      DO i=1,im
        h(i,jm) = h(i,jm)*dp(i,jm) + fy(i,jm)*rgw(jm)
      ENDDO
      CALL xpavg(h(1,jm), im)
    ENDIF

#ifndef VECTOR
    !-------------------------------------------------------------------
    ! Apply a simple nearest neighbor flux correction to reduce negatives
    !-------------------------------------------------------------------
    IF ( iv /= -1 ) THEN
      jsize = (jlast-jfirst+1)
      IF ( MOD(jsize,2)== 0) jsize=jsize+1
      jlast_odd = jfirst + jsize -1
      CALL fct_x(h, im, jm, jfirst, jlast, jlast_odd, ng, i)
    ENDIF

    DO j=jfirst,jlast
      DO i=1,im
        h(i,j) = h(i,j) / dpp(i,j)
      ENDDO
    ENDDO
#endif

  END SUBROUTINE tp2g

  SUBROUTINE tp2d(va, q, crx, cry, im, jm, iord, jord, ng, mg, fx, fy,      &
                  ffsl, xfx, yfx, id, jfirst, jlast)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: im, jm         ! Dimensions
    INTEGER, INTENT(in) :: jfirst, jlast  ! Latitude strip

    INTEGER, INTENT(in) :: iord, jord       ! Interpolation order in x,y
    INTEGER, INTENT(in) :: ng               ! Max. NS dependencies
    INTEGER, INTENT(in) :: mg               !
    INTEGER, INTENT(in) :: id               ! density (0)  (mfx = C)
                                           ! mixing ratio (1) (mfx = mass flux)
    LOGICAL, INTENT(in) :: ffsl(jfirst-ng:jlast+ng) ! Use flux-form semi-Lagrangian trans.?
                                                   ! ghosted N*ng S*ng
    REAL(wp), INTENT(in) :: va(im,jfirst:jlast)        ! Courant  (unghosted)
    REAL(wp), INTENT(in) :: q(im,jfirst-ng:jlast+ng)   ! transported scalar ( N*NG S*NG )
    REAL(wp), INTENT(in) :: crx(im,jfirst-ng:jlast+ng) ! Ask S.-J. ( N*NG S*NG )
    REAL(wp), INTENT(in) :: cry(im,jfirst:jlast+mg)    ! Ask S.-J. ( N like FY )
    REAL(wp), INTENT(in) :: xfx(im,jfirst:jlast)       ! Ask S.-J. ( unghosted like FX )
    REAL(wp), INTENT(in) :: yfx(im,jfirst:jlast+mg)    ! Ask S.-J. ( N like FY )

    REAL(wp), INTENT(out) :: fx(im,jfirst:jlast)    ! Flux in x ( unghosted )
    REAL(wp), INTENT(out) :: fy(im,jfirst:jlast+mg) ! Flux in y ( N, see tp2c )

    ! Local:
    INTEGER :: i, j, jp, jt, js2g0, js2gng, jn2g0, jn2gng, jn1g1
    REAL(wp) :: adx(im,jfirst-ng:jlast+ng)
    REAL(wp) :: wk2(-1:im+2,jfirst-ng:jlast+ng)
    REAL(wp) :: dm(0:im,jfirst-ng:jlast+ng)
    REAL(wp) :: al(0:im,jfirst-ng:jlast+ng)
    REAL(wp) :: ar(0:im,jfirst-ng:jlast+ng)
    REAL(wp) :: a6(0:im,jfirst-ng:jlast+ng)
    REAL(wp) :: dm1(im,jfirst-ng:jlast+ng)
    REAL(wp) :: a1l(im,jfirst-ng:jlast+ng)
    REAL(wp) :: a1r(im,jfirst-ng:jlast+ng)
    REAL(wp) :: a16(im,jfirst-ng:jlast+ng)
    INTEGER :: ilow(im,jfirst-ng:jlast+ng), ihigh(im,jfirst-ng:jlast+ng), ilow0, ihigh0
    INTEGER :: isave1(im,jfirst-ng:jlast+ng), isave2(im,jfirst-ng:jlast+ng)
    INTEGER :: icrit(im),icrx(im),jtidx(im)
    REAL(wp):: qsum(im)

    INTEGER :: idelta(jfirst-ng:jlast+ng), ic, itmp, lmt
    REAL(wp) :: atmp, tmp, qmax, qmin, qmax1, qmin1, qmax2, qmin2, ru
    INTEGER :: jm1, im2, js2gng1, jn2gng1
    INTEGER, PARAMETER :: iv0 = 0             !  Scalar (==0) Vector (==1)
    INTEGER :: js1g1, js2g1, jn1g2, jn2g1

    REAL(wp), PARAMETER :: cos_upw = 0.25_wp  ! critical cosine for upwind, ~ at 87 deg
    REAL(wp), PARAMETER :: cos_van = 0.25_wp  ! critical cosine for van Leer, ~ at 75 deg
    REAL(wp), PARAMETER :: cos_ppm = 0.25_wp  ! critical cosine for ppm
    REAL(wp), PARAMETER :: r24 = 1.0_wp/24.0_wp
    REAL(wp), PARAMETER :: r3  = 1.0_wp/3.0_wp
    REAL(wp), PARAMETER :: r23 = 2.0_wp/3.0_wp

#ifdef _PROFILE
    CALL trace_start ('tp2d_xtp', 7)
#endif    

    ! Number of ghost latitudes
    js2g0  = MAX(2,jfirst)          !  No ghosting
    jn2g0  = MIN(jm-1,jlast)        !  No ghosting
    js2gng = MAX(2,jfirst-ng)       !  Number needed on S
    jn2gng = MIN(jm-1,jlast+ng)     !  Number needed on N
    jn1g1  = MIN(jm,jlast+1)        ! Ghost N*1
    js2gng1 = MAX(2,   jfirst-ng+1)     !  Number needed on S
    jn2gng1 = MIN(jm-1,jlast+ng-1)      !  Number needed on N
    js1g1  = MAX(1,jfirst-1)         ! Ghost S*1
    js2g1  = MAX(2,jfirst-1)         ! Ghost S*1
    jn1g2  = MIN(jm,jlast+2)         ! Ghost N*2
    jn2g1  = MIN(jm-1,jlast+1)       ! Ghost N*1
    jm1 = jm - 1
    im2 = im / 2

!CDIR NOLSTVAL
    DO j=js2gng,jn2gng

       ! mpuetz: we convert all velocities to integer values up front
       ! so we avoid mixed integer-float computations in the following loop

       DO i=1,im
          icrx(i)  = INT(crx(i,j))
          icrit(i) = INT(FSEL(-crx(i,j),0.0_wp,1.0_wp))
       END DO

        IF( ffsl(j) ) THEN

           idelta(j) = -1
           DO i=1,im
              IF(icrit(i) > 0) THEN
                 itmp = i - icrx(i) - 1
                 ilow0  = itmp + 1
                 ihigh0 = i - 1
                 ! mpuetz: periodic boundaries in west direction
#ifdef _NOMODULO
                 IF (itmp <= 0) itmp=itmp+im
#else
                 itmp=MODULO(itmp-1,im)+1
#endif
                 isave1(i,j) = itmp
                 ilow(i,j)   = ilow0
                 ihigh(i,j)  = ihigh0
                 idelta(j)   = MAX(idelta(j),ihigh0-ilow0)
              ELSE
                 itmp   = i - icrx(i)
                 ilow0  = i
                 ihigh0 = itmp - 1
                 ! mpuetz: periodic boundaries in east direction
#ifdef _NOMODULO
                 IF (itmp > im) itmp=itmp-im
#else
                 itmp=MODULO(itmp-1,im)+1
#endif
                 isave1(i,j) = itmp
                 ilow(i,j)   = ilow0
                 ihigh(i,j)  = ihigh0
                 idelta(j)   = MAX(idelta(j),ihigh0-ilow0)
              ENDIF
           ENDDO

           ! mpuetz: sum up all layer contributions
           qsum(:)=0.0_wp
           DO ic=0,idelta(j)
              DO i=1,im
                 IF (ic <= ihigh(i,j)-ilow(i,j) ) THEN
                    itmp=ilow(i,j)+ic
#ifndef _NOMODULO
                    itmp=MODULO(itmp-1,im)+1
#else
                    ! mpuetz: ilow can be 0, so we must check both directions
                    IF (itmp > im) itmp=itmp-im
                    IF (itmp <= 0) itmp=itmp+im
#endif
                    qsum(i)=qsum(i)+q(itmp,j)
                 ENDIF
              ENDDO
           ENDDO
           
           DO i=1,im
              tmp = SIGN(1.0_wp,crx(i,j)-1.0_wp)
              itmp = isave1(i,j)
              wk2(i,j) = (crx(i,j)-AINT(crx(i,j),wp)) * q(itmp,j) + tmp * qsum(i)
           ENDDO

      ELSE
        ! Regular PPM (Eulerian without FFSL extension)

         ! mpuetz: itmp can only be i-1 or i, since ABS(crx(i,j)) < 1
         
         DO i=1,im

            itmp = i - icrit(i)
#ifndef _NOMODULO
            itmp=MODULO(itmp-1,im)+1
#else
            IF (itmp <= 0) itmp=itmp+im
#endif
            isave2(i,j) = itmp
            wk2(i,j) = crx(i,j)*q(itmp,j)
         ENDDO
      ENDIF

    ENDDO

    DO j=js2gng,jn2gng               !  adx needed on N*ng S*ng
      DO i=1,im-1
        adx(i,j) = q(i,j) + 0.5_wp *                       &
             (wk2(i,j)-wk2(i+1,j) + q(i,j)*(crx(i+1,j)-crx(i,j)))
      ENDDO
      adx(im,j) = q(im,j) + 0.5_wp *                     &
           (wk2(im,j)-wk2(1,j) + q(im,j)*(crx(1,j)-crx(im,j)))
    ENDDO

    IF ( jfirst == 1 ) THEN
      DO i=1,im
        adx(i, 1) = q(i,1)
      ENDDO
    ENDIF
    IF ( jlast == jm ) THEN
      DO i=1,im
        adx(i,jm) = q(i,jm)
      ENDDO
    ENDIF

    IF(jord == 1) THEN

      DO j=js2g0,jn1g1

         jtidx(1:im) =  INT(REAL(j,wp) - cry(1:im,j))

         DO i=1,im
            jt = jtidx(i)
            fy(i,j) = adx(i,jt) * yfx(i,j)
         ENDDO
      ENDDO

    ELSE

      ! YMIST requires q on NS;  Only call to YMIST here

      DO j=js2gng1,jn2gng1
        DO i=1,im
          dm1(i,j) = 0.25_wp*(adx(i,j+1) - adx(i,j-1))
        ENDDO
      ENDDO

      IF( iv0 == 0 ) THEN

        IF ( jfirst == 1 ) THEN
          ! S pole
          DO i=1,im2
            tmp = 0.25_wp*(adx(i,2)-adx(i+im2,2))
            qmax = MAX(adx(i,2),adx(i,1), adx(i+im2,2)) - adx(i,1)
            qmin = adx(i,1) - MIN(adx(i,2),adx(i,1), adx(i+im2,2))
            dm1(i,1) = SIGN(MIN(ABS(tmp),qmax,qmin),tmp)
            dm1(i+im2,1) = -dm1(i,1)
          ENDDO
        ENDIF

        IF ( jlast == jm ) THEN
          ! N pole
          DO i=1,im2
            tmp = 0.25_wp*(adx(i+im2,jm1)-adx(i,jm1))
            qmax = MAX(adx(i+im2,jm1),adx(i,jm), adx(i,jm1)) - adx(i,jm)
            qmin = adx(i,jm) - MIN(adx(i+im2,jm1),adx(i,jm), adx(i,jm1))
            dm1(i,jm) = SIGN(MIN(ABS(tmp),qmax,qmin),tmp)
            dm1(i+im2,jm) = -dm1(i,jm)
          ENDDO
        ENDIF

      ELSE

        IF ( jfirst == 1 ) THEN
          ! South
          DO i=1,im2
            tmp  = 0.25_wp*(adx(i,2)+adx(i+im2,2))
            qmax = MAX(adx(i,2),adx(i,1), -adx(i+im2,2)) - adx(i,1)
            qmin = adx(i,1) - MIN(adx(i,2),adx(i,1),-adx(i+im2,2))
            dm1(i,1) = SIGN(MIN(ABS(tmp),qmax,qmin),tmp)
            dm1(i+im2,1) = -dm1(i,1)
          ENDDO
        ENDIF

        IF ( jlast == jm ) THEN
          ! North
          DO i=1,im2
            tmp  = -0.25_wp*(adx(i+im2,jm1)+adx(i,jm1))
            qmax = MAX(-adx(i+im2,jm1),adx(i,jm), adx(i,jm1)) - adx(i,jm)
            qmin = adx(i,jm) - MIN(-adx(i+im2,jm1),adx(i,jm), adx(i,jm1))
            dm1(i,jm) = SIGN(MIN(ABS(tmp),qmax,qmin),tmp)
            dm1(i+im2,jm) = -dm1(i,jm)
          ENDDO
        ENDIF

      ENDIF

      ! Applies monotonic slope constraint

      DO j=js2gng1,jn2gng1-MOD(jn2gng1-js2gng1+1,2),2
        DO i=1,im
          qmax = MAX(adx(i,j), adx(i,j+1))
          qmin = MIN(adx(i,j), adx(i,j+1))
          qmax1 = MAX(adx(i,j-1),qmax) - adx(i,j)
          qmin1 = adx(i,j) - MIN(adx(i,j-1),qmin)
          dm1(i,j) = SIGN(MIN(ABS(dm1(i,j)),qmin1,qmax1),dm1(i,j))
          qmax2 = MAX(adx(i,j+2),qmax) - adx(i,j+1)
          qmin2 = adx(i,j+1) - MIN(adx(i,j+2),qmin)
          dm1(i,j+1) = SIGN(MIN(ABS(dm1(i,j+1)),qmin2,qmax2),dm1(i,j+1))
        ENDDO
      ENDDO
      IF (MOD(jn2gng1-js2gng1+1,2) == 1) THEN
        j=jn2gng1
        DO i=1,im
          qmax = MAX(adx(i,j-1),adx(i,j),adx(i,j+1)) - adx(i,j)
          qmin = adx(i,j) - MIN(adx(i,j-1),adx(i,j),adx(i,j+1))
          dm1(i,j) = SIGN(MIN(ABS(dm1(i,j)),qmin,qmax),dm1(i,j))
        ENDDO
      ENDIF

      IF( jord == 2 ) THEN

        DO j=js2g0,jn1g1
           jtidx(1:im) =  INT(REAL(j,wp) - cry(1:im,j))

          DO i=1,im
             jt = jtidx(i)
             fy(i,j) = (adx(i,jt) + (SIGN(1.0_wp,cry(i,j))-cry(i,j))*dm1(i,jt)) * yfx(i,j)
          ENDDO
        ENDDO

      ELSEIF( jord >= 3 ) THEN

        DO j=js2g1,jn1g2                 ! AL needed N2S
          DO i=1,im                      ! P, dm ghosted N2S2 (at least)
            a1l(i,j) = 0.5_wp*(adx(i,j-1)+adx(i,j)) + r3*(dm1(i,j-1) - dm1(i,j))
          ENDDO
        ENDDO

        DO j=js1g1,jn2g1                 ! AR needed NS
          DO i=1,im
            a1r(i,j) = a1l(i,j+1)          ! AL ghosted N2S
          ENDDO
        ENDDO

        ! Poles:

        IF( iv0 == 0 ) THEN

          IF ( jfirst == 1 ) THEN
!CDIR NODEP
            DO i=1,im2
              a1l(i,    1) = a1l(i+im2,2)
              a1l(i+im2,1) = a1l(i,    2)
            ENDDO
          ENDIF

          IF ( jlast == jm ) THEN
!CDIR NODEP
            DO i=1,im2
              a1r(i,    jm) = a1r(i+im2,jm1)
              a1r(i+im2,jm) = a1r(i,    jm1)
            ENDDO
          ENDIF

        ELSE

          IF ( jfirst == 1 ) THEN
!CDIR NODEP
            DO i=1,im2
              a1l(i,    1) = -a1l(i+im2,2)
              a1l(i+im2,1) = -a1l(i,    2)
            ENDDO
          ENDIF

          IF ( jlast == jm ) THEN
!CDIR NODEP
            DO i=1,im2
              a1r(i,    jm) = -a1r(i+im2,jm1)
              a1r(i+im2,jm) = -a1r(i,    jm1)
            ENDDO
          ENDIF

        ENDIF

        IF( jord == 3 .OR. jord == 5 ) THEN

          DO j=js1g1,jn1g1               ! A6 needed NS
            DO i=1,im
              a16(i,j) = 3.0_wp*(2.0_wp*adx(i,j) - (a1l(i,j)+a1r(i,j)))
            ENDDO
          ENDDO

        ENDIF

        lmt = jord - 3

        CALL lmppm(dm1(1,js1g1), a16(1,js1g1), a1r(1,js1g1),               &
             a1l(1,js1g1),  adx(1,js1g1), im*(jn1g1-js1g1+1), lmt)

        DO j=js2g0,jn1g1                 ! flux needed N
          DO i=1,im
            IF(cry(i,j) > 0.0_wp) THEN
              fy(i,j) = (a1r(i,j-1) + 0.5_wp*cry(i,j)*(a1l(i,j-1) - a1r(i,j-1) +  &
                   a16(i,j-1)*(1.0_wp-r23*cry(i,j)))) * yfx(i,j)
            ELSE
              fy(i,j) = (a1l(i,j) - 0.5_wp*cry(i,j)*(a1r(i,j) - a1l(i,j) +        &
                   a16(i,j)*(1.0_wp+r23*cry(i,j)))) * yfx(i,j)
            ENDIF
          ENDDO
        ENDDO

      ENDIF

    ENDIF

    DO j=js2g0,jn2g0

       jtidx(1:im) = INT(REAL(j,wp) - va(1:im,j))

      DO i=1,im
        jp = jtidx(i)
        wk2(i,j) = q(i,j) +0.5_wp*va(i,j)*(q(i,jp)-q(i,jp+1))
      ENDDO

      wk2(  -1,j) = wk2(im-1,j)
      wk2(   0,j) = wk2(  im,j)
      wk2(im+1,j) = wk2(   1,j)
      wk2(im+2,j) = wk2(   2,j)

      IF( ffsl(j) ) THEN

        IF( iord == 1 .OR. cosp(j) < cos_upw) THEN

          DO i=1,im
            itmp = isave1(i,j)
            fx(i,j) = (crx(i,j)-AINT(crx(i,j),wp)) * wk2(itmp,j)
          ENDDO

        ELSE

          DO i=1,im
            ! 2nd order slope
            tmp = 0.25_wp*(wk2(i+1,j) - wk2(i-1,j))
            qmax = MAX(wk2(i-1,j), wk2(i,j), wk2(i+1,j)) - wk2(i,j)
            qmin = wk2(i,j) - MIN(wk2(i-1,j), wk2(i,j), wk2(i+1,j))
            dm(i,j) = SIGN(MIN(ABS(tmp),qmax,qmin), tmp)
          ENDDO
          dm(0,j) = dm(im,j)

          IF(iord >= 3 .AND. cosp(j) > cos_ppm) THEN

            al(1,j) = 0.5_wp*(wk2(0,j)+wk2(1,j)) + (dm(0,j) - dm(1,j))*r3
            DO i=2,im
              al(i,j) = 0.5_wp*(wk2(i-1,j)+wk2(i,j)) + (dm(i-1,j) - dm(i,j))*r3
              ar(i-1,j) = al(i,j)
            ENDDO
            ar(im,j) = al(1,j)

            IF(iord == 7) THEN
              CALL huynh(im, ar(1,j), al(1,j), wk2(1,j), a6(1,j), dm(1,j))
            ELSE
              IF(iord == 3 .OR. iord == 5) THEN
                DO i=1,im
                  a6(i,j) = 3.0_wp*(wk2(i,j)+wk2(i,j)  - (al(i,j)+ar(i,j)))
                ENDDO
              ENDIF
              lmt = iord - 3
              CALL lmppm( dm(1,j), a6(1,j), ar(1,j), al(1,j), wk2(1,j), im, lmt )
            ENDIF

            al(0,j) = al(im,j)
            ar(0,j) = ar(im,j)
            a6(0,j) = a6(im,j)

            DO i=1,im
              ru   = crx(i,j) - AINT(crx(i,j))
              itmp = isave1(i,j)

              atmp = FSEL(crx(i,j),ar(itmp,j),al(itmp,j))
              tmp  = FSEL(crx(i,j),1.0_wp,-1.0_wp)

              fx(i,j) = ru * (atmp + 0.5_wp * ru &
                                     * (al(itmp,j) - ar(itmp,j) &
                                     + a6(itmp,j) * (tmp - r23 * ru)))
            ENDDO

          ELSE

            DO i=1,im
              ru   = crx(i,j) - AINT(crx(i,j))
              itmp = isave1(i,j)
              tmp  = FSEL(cry(i,j),1.0_wp,-1.0_wp)
              fx(i,j) = ru * (wk2(itmp,j) + dm(itmp,j) * (tmp - ru))
            ENDDO

          ENDIF

       ENDIF

        qsum(:)=0.0_wp
        DO ic=0,idelta(j)
          DO i=1,im
            IF (ic <= ihigh(i,j)-ilow(i,j) ) THEN
              itmp=ilow(i,j)+ic
#ifndef _NOMODULO
              itmp=MODULO(itmp-1,im)+1
#else
              ! mpuetz: ilow can be 0, so we must check both directions
              IF (itmp > im) itmp=itmp-im
              IF (itmp <= 0) itmp=itmp+im
#endif
              qsum(i)=qsum(i)+wk2(itmp,j)
            ENDIF
          ENDDO
        ENDDO

        IF(id /= 0) THEN
          DO i=1,im
            tmp = SIGN(1.0_wp,crx(i,j)-1.0_wp)
            fx(i,j) = fx(i,j) + tmp * qsum(i)
            fx(i,j) =  fx(i,j) * xfx(i,j)
          ENDDO
        ELSE
          DO i=1,im
            tmp = SIGN(1.0_wp,crx(i,j)-1.0_wp)
            fx(i,j) = fx(i,j) + tmp * qsum(i)
          ENDDO
        ENDIF

      ELSE
        ! Regular PPM (Eulerian without FFSL extension)

        IF(iord == 1 .OR. cosp(j) < cos_upw) THEN

          DO i=1,im
            itmp = isave2(i,j)
            fx(i,j) = xfx(i,j)*wk2(itmp,j)
          ENDDO

        ELSE

          DO i=1,im
            dm(i,j) = r24*(8.0_wp*(wk2(i+1,j) - wk2(i-1,j)) + wk2(i-2,j) - wk2(i+2,j))
            ! Apply monotonicity constraint (Lin et al. 1994, MWR)
            qmax = MAX( wk2(i-1,j), wk2(i,j), wk2(i+1,j) ) - wk2(i,j)
            qmin = wk2(i,j) - MIN( wk2(i-1,j), wk2(i,j), wk2(i+1,j) )
            dm(i,j) = SIGN( MIN(ABS(dm(i,j)), qmax, qmin), dm(i,j) )
          ENDDO
          dm(0,j) = dm(im,j)

          IF( iord==2 .OR. cosp(j) < cos_van ) THEN
            DO i=1,im
              itmp = isave2(i,j)
              fx(i,j) =  xfx(i,j)*(wk2(itmp,j)+dm(itmp,j)*(SIGN(1.0_wp,crx(i,j))-crx(i,j)))
            ENDDO
          ELSE

            al(1,j) = 0.5_wp*(wk2(0,j)+wk2(1,j)) + (dm(0,j) - dm(1,j))*r3
            DO i=2,im
              al(i,j) = 0.5_wp*(wk2(i-1,j)+wk2(i,j)) + (dm(i-1,j) - dm(i,j))*r3
              ar(i-1,j) = al(i,j)
            ENDDO
            ar(im,j) = al(1,j)

            IF(iord == 7) THEN
              CALL huynh(im, ar(1,j), al(1,j), wk2(1,j), a6(1,j), dm(1,j))
            ELSE
              IF(iord == 3 .OR. iord == 5) THEN
                DO i=1,im
                  a6(i,j) = 3.0_wp*(wk2(i,j)+wk2(i,j)  - (al(i,j)+ar(i,j)))
                ENDDO
              ENDIF
              lmt = iord - 3
              CALL lmppm( dm(1,j), a6(1,j), ar(1,j), al(1,j), wk2(1,j), im, lmt )
            ENDIF

            al(0,j) = al(im,j)
            ar(0,j) = ar(im,j)
            a6(0,j) = a6(im,j)

            DO i=1,im
              IF(crx(i,j) > 0.0_wp) THEN
                fx(i,j) = (ar(i-1,j) &
                           + 0.5_wp * crx(i,j) & 
                             * (al(i-1,j) - ar(i-1,j) &
                                + a6(i-1,j) * (1.0_wp - r23 * crx(i,j)))) &
                          * xfx(i,j)
              ELSE
                fx(i,j) = (al(i,j) &
                           - 0.5_wp * crx(i,j) &
                             * (ar(i,j) - al(i,j) &
                                + a6(i,j) * (1.0_wp + r23 * crx(i,j)))) &
                          * xfx(i,j)
              ENDIF
            ENDDO

          ENDIF

        ENDIF

      ENDIF

   END DO

#ifdef _PROFILE
    CALL trace_stop ('tp2d_xtp', 7)
#endif

  END SUBROUTINE tp2d

  SUBROUTINE lmppm(dm, a6, ar, al, p, im, lmt)

    INTEGER, INTENT(in):: im   ! Total longitudes
    INTEGER, INTENT(in):: lmt  ! LMT = 0: full monotonicity
                               ! LMT = 1: Improved and simplified full monotonic constraint
                               ! LMT = 2: positive-definite constraint
                               ! LMT = 3: Quasi-monotone constraint
    REAL(wp), INTENT(in):: p(im)
    REAL(wp), INTENT(in):: dm(im)

    REAL(wp), INTENT(inout):: a6(im)
    REAL(wp), INTENT(inout):: ar(im)
    REAL(wp), INTENT(inout):: al(im)

    !Local:
    REAL(wp), PARAMETER :: r12 = 1.0_wp/12.0_wp

    REAL(wp) :: da1, da2, fmin, a6da
    REAL(wp) :: dr, dl

    INTEGER :: i

    ! LMT = 0: full monotonicity
    ! LMT = 1: Improved and simplified full monotonic constraint
    ! LMT = 2: positive-definite constraint
    ! LMT = 3: Quasi-monotone constraint

#ifdef _PROFILE
    CALL trace_start ('lmppm', 8)
#endif

    IF( lmt == 0 ) THEN

      ! Full constraint
      DO i=1,im
        IF(dm(i) == 0.0_wp) THEN
          ar(i) = p(i)
          al(i) = p(i)
          a6(i) = 0.0_wp
        ELSE
          da1  = ar(i) - al(i)
          da2  = da1**2
          a6da = a6(i)*da1
          IF(a6da < -da2) THEN
            a6(i) = 3.0_wp*(al(i)-p(i))
            ar(i) = al(i) - a6(i)
          ELSEIF(a6da > da2) THEN
            a6(i) = 3.0_wp*(ar(i)-p(i))
            al(i) = ar(i) - a6(i)
          ENDIF
        ENDIF
      ENDDO

    ELSEIF( lmt == 1 ) THEN

      ! Improved (Lin 200?) full constraint
      DO i=1,im
        da1 = dm(i) + dm(i)
        dl = SIGN(MIN(ABS(da1),ABS(al(i)-p(i))), da1)
        dr = SIGN(MIN(ABS(da1),ABS(ar(i)-p(i))), da1)
        ar(i) = p(i) + dr
        al(i) = p(i) - dl
        a6(i) = 3.0_wp*(dl-dr)
      ENDDO

    ELSEIF( lmt == 2 ) THEN
      ! Positive definite only constraint
      DO i=1,im
        IF(ABS(ar(i)-al(i)) >= -a6(i)) CYCLE
        fmin = p(i) + 0.25_wp*(ar(i)-al(i))**2/a6(i) + a6(i)*r12
        IF(fmin >= 0.0_wp) CYCLE
        IF(p(i) < ar(i) .AND. p(i) < al(i)) THEN
          ar(i) = p(i)
          al(i) = p(i)
          a6(i) = 0.0_wp
        ELSEIF(ar(i) > al(i)) THEN
          a6(i) = 3.0_wp*(al(i)-p(i))
          ar(i) = al(i) - a6(i)
        ELSE
          a6(i) = 3.0_wp*(ar(i)-p(i))
          al(i) = ar(i) - a6(i)
        ENDIF
      ENDDO

    ELSEIF(lmt == 3) THEN
      ! Quasi-monotone constraint
      DO i=1,im
        da1 = 4.0_wp*dm(i)
        dl = SIGN(MIN(ABS(da1),ABS(al(i)-p(i))), da1)
        dr = SIGN(MIN(ABS(da1),ABS(ar(i)-p(i))), da1)
        ar(i) = p(i) + dr
        al(i) = p(i) - dl
        a6(i) = 3.0_wp*(dl-dr)
      ENDDO
    ENDIF

#ifdef _PROFILE
    CALL trace_stop ('lmppm', 8)
#endif

  END SUBROUTINE lmppm

  SUBROUTINE huynh(im, ar, al, p, d2, d1)

    INTEGER, INTENT(in) :: im
    REAL(wp), INTENT(in) :: p(im)

    REAL(wp), INTENT(inout) :: ar(im)
    REAL(wp), INTENT(inout) :: al(im)
    REAL(wp), INTENT(out) :: d2(im)
    REAL(wp), INTENT(out) :: d1(im)

    ! Local:
    INTEGER :: i
    REAL(wp) :: pmp
    REAL(wp) :: lac
    REAL(wp) :: pmin
    REAL(wp) :: pmax

    ! Compute d1 and d2
    d1(1) = p(1) - p(im)
    DO i=2,im
      d1(i) = p(i) - p(i-1)
    ENDDO

    DO i=1,im-1
      d2(i) = d1(i+1) - d1(i)
    ENDDO
    d2(im) = d1(1) - d1(im)

    ! Constraint for AR
    !            i = 1
    pmp   = p(1) + 2.0_wp * d1(1)
    lac   = p(1) + 0.5_wp * (d1(1)+d2(im)) + d2(im)
    pmin  = MIN(p(1), pmp, lac)
    pmax  = MAX(p(1), pmp, lac)
    ar(1) = MIN(pmax, MAX(ar(1), pmin))

    DO i=2, im
      pmp   = p(i) + 2.0_wp*d1(i)
      lac   = p(i) + 0.5_wp*(d1(i)+d2(i-1)) + d2(i-1)
      pmin  = MIN(p(i), pmp, lac)
      pmax  = MAX(p(i), pmp, lac)
      ar(i) = MIN(pmax, MAX(ar(i), pmin))
    ENDDO

    ! Constraint for AL
    DO i=1, im-1
      pmp   = p(i) - 2.0_wp*d1(i+1)
      lac   = p(i) + 0.5_wp*(d2(i+1)-d1(i+1)) + d2(i+1)
      pmin  = MIN(p(i), pmp, lac)
      pmax  = MAX(p(i), pmp, lac)
      al(i) = MIN(pmax, MAX(al(i), pmin))
    ENDDO

    i = im
    pmp    = p(im) - 2.0_wp*d1(1)
    lac    = p(im) + 0.5_wp*(d2(1)-d1(1)) + d2(1)
    pmin   = MIN(p(im), pmp, lac)
    pmax   = MAX(p(im), pmp, lac)
    al(im) = MIN(pmax, MAX(al(im), pmin))

    ! compute A6 (d2)
    DO i=1, im
      d2(i) = 3.0_wp*(2*p(i)  - (al(i)+ar(i)))
    ENDDO

  END SUBROUTINE huynh

  SUBROUTINE xpavg(p, im)

    INTEGER, INTENT(in) :: im

    REAL(wp), INTENT(inout) :: p(im)

    ! Local
    INTEGER :: i
    REAL(wp) :: sum1

    sum1 = 0.0_wp
    DO i=1,im
      sum1 = sum1 + p(i)
    ENDDO
    sum1 = sum1 / REAL(im, wp)

    DO i=1,im
      p(i) = sum1
    ENDDO

  END SUBROUTINE xpavg



  SUBROUTINE qmap_gp(ps, delp, q, kproma, kpromz, jm, km, nq, ak, bk)

    INTEGER, INTENT(in) :: kproma, kpromz
    INTEGER, INTENT(in) :: jm, km                ! y, z dimensions
    INTEGER, INTENT(in) :: nq                    ! number of tracers
    REAL(wp), INTENT(in) :: ak(km+1)
    REAL(wp), INTENT(in) :: bk(km+1)
    REAL(wp), INTENT(in) :: delp(kproma,km,jm)

    REAL(wp), INTENT(inout) :: ps(kproma,jm)      ! surface pressure
    REAL(wp), INTENT(inout) :: q(kproma,km,nq,jm) ! tracers including specific humidity

    ! Local arrays:
    REAL(wp) :: pe1(kproma,km+1)
    REAL(wp) :: pe2(kproma,km+1)
#ifdef STANDARD_FORTRAN
    REAL(wp) :: qtmp(kproma,km)
#endif
    INTEGER :: i, j, k, iq, im

!$OMP DO ECHAM_GPC_SCHEDULE
   DO j=1,jm

      IF ( j == jm ) THEN
        im = kpromz
      ELSE
        im = kproma
      END IF

      DO i=1,im
        pe1(i,1) = ak(1)
      ENDDO 

      DO k=1,km
        DO i=1,im
          pe1(i,k+1) = pe1(i,k) + delp(i,k,j)
        ENDDO
      ENDDO

      DO i=1,im
        ps(i,j) = pe1(i,km+1)
      ENDDO


      ! k=1
      DO i=1,im
        pe2(i,1) = ak(1)
      ENDDO

      DO k=2,km
        DO i=1,im
          pe2(i,k) = ak(k) + bk(k)*ps(i,j)
        ENDDO
      ENDDO

      ! k=km+1
      DO i=1,im
        pe2(i,km+1) = ps(i,j)
      ENDDO

      DO iq=1,nq
#ifdef STANDARD_FORTRAN
        ! q is modified in map1_ppm_gp, so we can use it only once in the argument list
        DO k=1, km
          DO i=1,kproma
            qtmp(i,k) = q(i,k,iq,j)
          ENDDO
        ENDDO
        CALL map1_ppm_gp ( km, pe1,   qtmp(1,1),   &
                           km, pe2,   q(1,1,iq,j),   &
                           kproma, im)
#else
        CALL map1_ppm_gp ( km, pe1,   q(1,1,iq,j),   &
                           km, pe2,   q(1,1,iq,j),   &
                           kproma, im)
#endif
      ENDDO

    ENDDO
!$OMP END DO

  END SUBROUTINE qmap_gp

  SUBROUTINE map1_ppm_gp( km,   pe1,   q1,                         &
                          kn,   pe2,   q2,                         &
                          kbdim, im)

    INTEGER, INTENT(in) :: kbdim
    INTEGER, INTENT(in) :: im     ! E-W dimension
    INTEGER, INTENT(in) :: km     ! Original vertical dimension
    INTEGER, INTENT(in) :: kn     ! Target vertical dimension

    REAL(wp), INTENT(in) :: pe1(kbdim,km+1)  ! pressure at layer edges
                                          ! (from model top to bottom surface)
                                          ! in the original vertical coordinate
    REAL(wp), INTENT(in) :: pe2(kbdim,kn+1)  ! pressure at layer edges
                                          ! (from model top to bottom surface)
                                          ! in the new vertical coordinate

    REAL(wp), INTENT(inout) :: q1(kbdim,km)     ! Field input
    REAL(wp), INTENT(inout) :: q2(kbdim,kn)     ! Field output

    ! Local variables:
    REAL(wp), PARAMETER :: r3  = 1.0_wp/3.0_wp
    REAL(wp), PARAMETER :: r23 = 2.0_wp/3.0_wp

    REAL(wp) :: dp1(kbdim,km),dp1i(im,km)
    REAL(wp) :: q4(kbdim,km,4)
    INTEGER :: i, k, l, k0
    REAL(wp) :: pl, pr, esl

#if (defined _UNICOSMP) || (defined __SX__) || (defined ES)
    INTEGER, PARAMETER :: strip = _STRIP_
    INTEGER :: j, i0, l0, l1, l2, lmax, idx(im), idx1(im), icnt, icnt1, icnt2, icnt3, inew, k0new
!CDIR VREG(lv,qsum0)
    INTEGER :: lv(strip), lx(kbdim,kn+1), count
    REAL(wp) :: qsum0(strip), delp1, delp2, qsum
#elif (defined __tuned__)
    INTEGER :: idx(im),  idx1(im), lx(kbdim,kn+1)
    INTEGER :: j, inew, k0new, icnt, icnt1, icnt2, icnt3, l0, l1, l2, lmax
    REAL(wp) :: delp1, delp2
#else
    INTEGER :: ll
    REAL(wp) :: delp
#endif

    DO k=1,km
      DO i=1,im
        dp1(i,k) = pe1(i,k+1) - pe1(i,k)
        dp1i(i,k) = 1._wp/dp1(i,k) 
        q4(i,k,1) = q1(i,k)
      ENDDO
    ENDDO

    ! Compute vertical subgrid distribution
    CALL ppm2m( q4, dp1, km, kbdim, im )

#ifdef _PROFILE
    CALL trace_start ('map1_ppm_gp', 6)
#endif

#if (defined __tuned__)

    k0 = 1
    DO k = 1,kn+1
       DO i = 1,im
          idx(i) = i
       END DO
       icnt = im
       k0new = km+1
       DO l=k0,km
          inew = 0
!IBM* ASSERT(NODEPS)
          DO j = 1,icnt
             i = idx(j)
             idx1(j) = INT(FSEL(pe1(i,l+1) - pe2(i,k),0.0_wp,1.0_wp))
          END DO

          IF (icnt > 0) THEN
             DO j = 1,icnt
                i = idx(j)
                IF (idx1(j) > 0) THEN
                   inew = inew + 1
                   idx(inew) = i
                ELSE
                   lx(i,k) = l
                END IF
             END DO
          END IF

!             IF ((pe1(i,l+1) >= pe2(i,k)).AND.(pe2(i,k) >= pe1(i,l))) THEN
!             IF (pe1(i,l+1) < pe2(i,k)) THEN
!                inew = inew + 1
!                idx(inew) = i
!             ELSE
!                lx(i,k) = l
!             END IF

          IF (inew < icnt) k0new = MIN(k0new,l) 
          icnt = inew
          IF (inew == 0) EXIT
       END DO
       k0 = k0new

       DO j = 1,icnt
          lx(idx(j),k) = km+1
       END DO

       IF (k > 1) THEN

          icnt1=1
          icnt2=1
          icnt3=im 
          lmax = 0
          DO i=1,im
             l1=lx(i,k-1)
             l2=lx(i,k)
             IF (l1 < l2) THEN
                lmax=MAX(lmax,l2-l1)
                IF (l2 /= km+1) THEN
                   idx(icnt2) = i
                   icnt2 = icnt2 + 1
                ELSE
                   idx(icnt3) = i
                   icnt3 = icnt3 - 1
                END IF
             ELSE ! l1 == l2  (because l1 <= l2 always by construction)
                idx1(icnt1) = i
                icnt1 = icnt1 + 1
             END IF
          END DO

          ! completely internal case (l1 == l2)

!IBM* ASSERT(NODEPS)
          DO j=1,icnt1-1
             i=idx1(j)
             l1=lx(i,k-1)
             l2=lx(i,k)
             pl=(pe2(i,k-1)-pe1(i,l1))*dp1i(i,l1)
             pr=(pe2(i,k)-pe1(i,l1))*dp1i(i,l1)
             delp1=pe2(i,k)-pe2(i,k-1)
             q2(i,k-1) = delp1*(q4(i,l1,2)                   &
                  +0.5_wp*(q4(i,l1,4)+q4(i,l1,3)-q4(i,l1,2)) &
                  *(pr+pl)                                   &
                  -q4(i,l1,4)*r3*(pr*(pr+pl)+pl**2))
          END DO

!IBM* ASSERT(NODEPS)
          DO j=1,icnt2-1   ! (l1 < l2).AND.(l2 /= km+1) (this is the most likely case)
             i=idx(j)
             l1=lx(i,k-1)
             l2=lx(i,k)
             pl=(pe2(i,k-1)-pe1(i,l1))*dp1i(i,l1)
             pr=1.0_wp
             delp1=pe1(i,l1+1)-pe2(i,k-1)
             delp2=pe2(i,k)-pe1(i,l2)
             esl=delp2*dp1i(i,l2)
             q2(i,k-1)=delp2*(q4(i,l2,2)                       &
                  +0.5_wp*esl*(q4(i,l2,3)-q4(i,l2,2)         &
                  +q4(i,l2,4)*(1.0_wp-r23*esl)))             &
                  +delp1*(q4(i,l1,2)                         &
                  +0.5_wp*(q4(i,l1,4)+q4(i,l1,3)-q4(i,l1,2)) &
                  *(pr+pl)                                   &
                  -q4(i,l1,4)*r3*(pr*(pr+pl)+pl**2))
          ENDDO

!IBM* ASSERT(NODEPS)
          DO j=icnt3+1,im   ! (l1 < l2).AND.(l2 == km+1) (border case)
             i=idx(j)
             l1=lx(i,k-1)
             l2=lx(i,k)
             pl=(pe2(i,k-1)-pe1(i,l1))*dp1i(i,l1)
             pr=1.0_wp
             delp1=pe1(i,l1+1)-pe2(i,k-1)
             q2(i,k-1) = dp1(i,km)*q4(i,km,1)                  &
                  +delp1*(q4(i,l1,2)                         &
                  +0.5_wp*(q4(i,l1,4)+q4(i,l1,3)-q4(i,l1,2)) &
                  *(pr+pl)                                   &
                  -q4(i,l1,4)*r3*(pr*(pr+pl)+pl**2))
          ENDDO

          IF (lmax > 2) THEN ! this case is quite unlikely

             DO l0=1,lmax-1
                DO i=1,im
                   l1=lx(i,k-1)
                   l2=lx(i,k)
                   l=l0+l1
                   IF (l<l2) THEN
                      q2(i,k-1) = q2(i,k-1) + dp1(i,l)*q4(i,l,1)
                   ENDIF
                ENDDO
             ENDDO
             
          ELSE

             DO i=1,im
                l1=lx(i,k-1)
                l2=lx(i,k)
                IF (l2-l1 == 2) THEN
                   q2(i,k-1) = q2(i,k-1) + dp1(i,l1+1)*q4(i,l1+1,1)
                end if
             ENDDO
          
          ENDIF
!IBM* NOVECTOR
          DO i=1,im
             q2(i,k-1)=q2(i,k-1)/(pe2(i,k)-pe2(i,k-1))
          END DO
          
       END IF
    ENDDO

#elif (defined _UNICOSMP) || (defined __SX__) || (defined ES)

    DO i0=1,im,strip
      k0=1
      DO k=1,kn+1
!CDIR SHORTLOOP
        DO i=i0,min(i0+strip-1,im)
          lv(i-i0+1)=km+1
        ENDDO
        count=min(strip,im-i0+1)
        DO l=k0,km
!CDIR SHORTLOOP
          DO i=i0,min(i0+strip-1,im)
            IF (pe2(i,k)>=pe1(i,l).AND.pe2(i,k)<=pe1(i,l+1)) THEN
              lv(i-i0+1)=l
              count=count-1
            ENDIF
          ENDDO
          IF (count==0) THEN
            k0=km+1
!CDIR SHORTLOOP
            DO i=i0,min(i0+strip-1,im)
              k0=MIN(k0,lv(i-i0+1))
            ENDDO
            exit
          ENDIF
        ENDDO
!CDIR SHORTLOOP
        DO i=i0,min(i0+strip-1,im)
          lx(i,k)=lv(i-i0+1)
        ENDDO
      ENDDO
    ENDDO

    lmax=0
!CDIR OUTERUNROLL=8
    DO k=1,kn
      DO i=1,im
        l1=lx(i,k)
        l2=lx(i,k+1)
        lmax=MAX(lmax,l2-l1)
        qsum=0.0_wp
        pl=(pe2(i,k)-pe1(i,l1))*dp1i(i,l1)
        IF (l2-l1==0) THEN
          pr=(pe2(i,k+1)-pe1(i,l1))*dp1i(i,l1)
          delp1=pe2(i,k+1)-pe2(i,k)
        ELSE
          pr=1.0_wp
          delp1=pe1(i,l1+1)-pe2(i,k)
          IF (l2/=km+1) THEN
            delp2=pe2(i,k+1)-pe1(i,l2)
            esl=delp2*dp1i(i,l2)
            qsum=qsum                                      &
                 +delp2*(q4(i,l2,2)                        &
                        +0.5_wp*esl*(q4(i,l2,3)-q4(i,l2,2) &
                                    +q4(i,l2,4)*(1.0_wp-r23*esl)))
          ELSE
           qsum=qsum+dp1(i,km)*q4(i,km,1)
          ENDIF
        ENDIF
        q2(i,k)=qsum                                             &
               +delp1*(q4(i,l1,2)                                &
                      +0.5_wp*(q4(i,l1,4)+q4(i,l1,3)-q4(i,l1,2)) &
                             *(pr+pl)                            &
                      -q4(i,l1,4)*r3*(pr*(pr+pl)+pl**2))
      ENDDO
    ENDDO

    IF (lmax > 2) THEN

!CDIR NODEP
    DO i0=1,im,strip
      DO k=1,kn
!CDIR SHORTLOOP
        DO i=i0,min(i0+strip-1,im)
          l1=lx(i,k)
          l2=lx(i,k+1)
          qsum0(i-i0+1)=0.0_wp
        ENDDO
        DO l0=1,lmax-1
!CDIR SHORTLOOP
          DO i=i0,min(i0+strip-1,im)
           l1=lx(i,k)
           l2=lx(i,k+1)
           l=l0+l1
           IF (l<l2) THEN
            qsum0(i-i0+1)=qsum0(i-i0+1)+dp1(i,l)*q4(i,l,1)
           ENDIF
          ENDDO
        ENDDO
!CDIR SHORTLOOP
        DO i=i0,min(i0+strip-1,im)
          q2(i,k)= (q2(i,k)+qsum0(i-i0+1))/(pe2(i,k+1)-pe2(i,k))
        ENDDO
      ENDDO
    ENDDO

    ELSE

!CDIR OUTERUNROLL=8
    DO k=1,kn
      DO i=1,im
        l1=lx(i,k)
        l2=lx(i,k+1)
        qsum=0.0_wp
        IF (l2-l1 == 2) THEN
          qsum=dp1(i,l1+1)*q4(i,l1+1,1)
        end if
        q2(i,k)=(q2(i,k)+qsum)/(pe2(i,k+1)-pe2(i,k))
      ENDDO
    ENDDO

    ENDIF

#else

    ! mpuetz: old mapping, really slow

    DO i=1,im
      k0 = 1
      DO 122 k=1,kn
        DO l=k0,km
          ! locate the top edge: pe2(i,k)
          IF(pe2(i,k) >= pe1(i,l) .AND. pe2(i,k) <= pe1(i,l+1)) THEN
            pl = (pe2(i,k)-pe1(i,l)) * dp1i(i,l)
            IF(pe2(i,k+1) <= pe1(i,l+1)) THEN
              ! entire new grid is within the original grid
              pr = (pe2(i,k+1)-pe1(i,l)) * dp1i(i,l)
              q2(i,k) = q4(i,l,2) + 0.5_wp*(q4(i,l,4)+q4(i,l,3)-q4(i,l,2)) &
                   *(pr+pl)-q4(i,l,4)*r3*(pr*(pr+pl)+pl**2)
              k0 = l
              GOTO 122
            ELSE
              ! Fractional area...
              qsum = (pe1(i,l+1)-pe2(i,k))*(q4(i,l,2)+0.5_wp*(q4(i,l,4)+     &
                   q4(i,l,3)-q4(i,l,2))*(1.0_wp+pl)-q4(i,l,4)*             &
                   (r3*(1.0_wp+pl*(1.0_wp+pl))))
              DO ll=l+1,km
                ! locate the bottom edge: pe2(i,k+1)
                IF(pe2(i,k+1) > pe1(i,ll+1) ) THEN
                  ! Whole layer..
                  qsum = qsum + dp1(i,ll)*q4(i,ll,1)
                ELSE
                  delp = pe2(i,k+1)-pe1(i,ll)
                  esl = delp * dp1i(i,ll)
                  qsum = qsum + delp*(q4(i,ll,2)+0.5_wp*esl*               &
                       (q4(i,ll,3)-q4(i,ll,2)+q4(i,ll,4)*(1.0_wp-r23*esl)))
                  k0 = ll
                  GOTO 123
                ENDIF
              ENDDO
              GOTO 123
            ENDIF
          ENDIF
        ENDDO
123     q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
122   ENDDO
    ENDDO

#endif

#ifdef _PROFILE
    CALL trace_stop ('map1_ppm_gp', 6)
#endif

  END SUBROUTINE map1_ppm_gp

  SUBROUTINE ppm2m(a4, delp, km, kbdim, im)

    INTEGER, INTENT(in):: kbdim   ! real leading dimension
    INTEGER, INTENT(in):: im      ! actual leading dimension
    INTEGER, INTENT(in):: km      ! vertical dimension
    REAL(wp), INTENT(in):: delp(kbdim,km)     ! layer pressure thickness
    REAL(wp), INTENT(inout):: a4(kbdim,km,4)  ! Interpolated values

    ! local arrays.
    REAL(wp) :: dc(kbdim,km)
!    REAL(wp) :: h2(kbdim,km)
    REAL(wp) :: delq(kbdim,km)
    REAL(wp) :: df2(kbdim,km)
    REAL(wp) :: d4(kbdim,km)
    REAL(wp) :: d4i(kbdim,km)
    REAL(wp) :: d3(kbdim,km)
    REAL(wp) :: d3i(kbdim,km)

    REAL(wp) :: fac
    REAL(wp) :: a1, a2, c1, c2, c3, d1, d2
    REAL(wp) :: qmax, qmin, cmax, cmin
    REAL(wp) :: qm, dq, tmp
    REAL(wp) :: qmp, pmp
    REAL(wp) :: lac
    INTEGER :: lmt
    INTEGER :: i, k

    REAL(wp), PARAMETER :: r12 = 1.0_wp/12.0_wp
    REAL(wp) :: da1, da2, a6da
    REAL(wp) :: fmin

    ! AA <-- a4(1,i)
    ! AL <-- a4(2,i)
    ! AR <-- a4(3,i)
    ! A6 <-- a4(4,i)

#ifdef _PROFILE
    CALL trace_start ('ppm2m', 5)
#endif    

    DO i=1,im
       delq(i,1) =   a4(i,2,1) - a4(i,1,1)
       delq(i,2) =   a4(i,3,1) - a4(i,2,1)
       d4(i,2  ) = delp(i,1) + delp(i,2)
       d4(i,3  ) = delp(i,2) + delp(i,3)
       d3(i,2  ) = d4(i,2)+delp(i,3)

       d4i(i,2 ) = 1.0_wp/d4(i,2)
       d4i(i,3 ) = 1.0_wp/d4(i,3)
       d3i(i,2 ) = 1.0_wp/d3(i,2)
    END DO

    k=2
    DO i=1,im
       c1  = (delp(i,k-1)+0.5_wp*delp(i,k))*d4i(i,k+1)
       c2  = (delp(i,k+1)+0.5_wp*delp(i,k))*d4i(i,k)
       tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) * d3i(i,k)
       qmax = MAX(a4(i,k-1,1),a4(i,k,1),a4(i,k+1,1)) - a4(i,k,1)
       qmin = a4(i,k,1) - MIN(a4(i,k-1,1),a4(i,k,1),a4(i,k+1,1))
       dc(i,k) = SIGN(MIN(ABS(tmp),qmax,qmin), tmp)
       df2(i,k) = tmp
    ENDDO

    DO k=3,km-1
       DO i=1,im
          delq(i,k) =   a4(i,k+1,1) - a4(i,k,1)
          d4(i,k+1  ) = delp(i,k) + delp(i,k+1)
          d3(i,k  )   = d4(i,k)+delp(i,k+1)
          d4i(i,k+1 ) = 1.0_wp/d4(i,k+1)
          d3i(i,k )   = 1.0_wp/d3(i,k)
       END DO

     DO i=1,im
        c1  = (delp(i,k-1)+0.5_wp*delp(i,k))*d4i(i,k+1)
        c2  = (delp(i,k+1)+0.5_wp*delp(i,k))*d4i(i,k)
        tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) * d3i(i,k)
        qmax = MAX(a4(i,k-1,1),a4(i,k,1),a4(i,k+1,1)) - a4(i,k,1)
        qmin = a4(i,k,1) - MIN(a4(i,k-1,1),a4(i,k,1),a4(i,k+1,1))
        dc(i,k) = SIGN(MIN(ABS(tmp),qmax,qmin), tmp)
        df2(i,k) = tmp
     END DO

    !------------------------------------------------------------
    ! 4th order interpolation of the provisional cell edge value
    !------------------------------------------------------------

     DO i=1,im
        c1 = delq(i,k-1)*delp(i,k-1) * d4i(i,k)
        a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
        a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
        a4(i,k,2) = a4(i,k-1,1) + c1 + 2.0_wp/(d4(i,k-1)+d4(i,k+1)) *      &
             ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -            &
             delp(i,k-1)*a1*dc(i,k  ) )
      ENDDO
      IF (kord <=3) a4(1:im,k-1,3) = a4(1:im,k,2)
    ENDDO

    IF(kord > 3) CALL steepz(kbdim, im, km, a4, df2, dc, delq, delp, d4, d4i, d3i)

    ! Area preserving cubic with 2nd deriv. = 0 at the boundaries
    ! Top
    DO i=1,im
      d1 = delp(i,1)
      d2 = delp(i,2)
      qm = (d2*a4(i,1,1)+d1*a4(i,2,1)) / (d1+d2)
      dq = 2.0_wp*(a4(i,2,1)-a4(i,1,1)) / (d1+d2)
      c1 = (a4(i,3,2)-qm-d2*dq) / ( d2*(2.0_wp*d2*d2+d1*(d2+3.0_wp*d1)) )
      c3 = dq - 2.0_wp*c1*(d2*(5.0_wp*d1+d2)-3.0_wp*d1**2)
      a4(i,2,2) = qm - c1*d1*d2*(d2+3.0_wp*d1)
      a4(i,1,2) = d1*(8.0_wp*c1*d1**2-c3) + a4(i,2,2)
      dc(i,1) =  a4(i,1,1) - a4(i,1,2)
      ! No over- and undershoot condition
      cmax = MAX(a4(i,1,1), a4(i,2,1))
      cmin = MIN(a4(i,1,1), a4(i,2,1))
      a4(i,2,2) = MIN(cmax,MAX(cmin,a4(i,2,2)))
   END DO

   IF( iv == 0 ) THEN

      a4(1:im,1,2) = MAX(0.0_wp,a4(1:im,1,2))

   ELSEIF ( iv ==  1 ) THEN
      
      ! Monotone tracers:
      dc(1:im,1)   = 0.0_wp
      a4(1:im,1,2) = a4(1:im,1,1)
      a4(1:im,2,2) = a4(1:im,1,1)

  ELSEIF ( iv == -1 ) THEN

     DO i=1,im
        ! Winds:
        IF( a4(i,1,1)*a4(i,1,2) <=  0.0_wp ) THEN
           a4(i,1,2) = 0.0_wp
        ELSE
           a4(i,1,2) = SIGN(MIN(ABS(a4(i,1,1)),      &
                ABS(a4(i,1,2))), a4(i,1,1)  )
        ENDIF
     END DO
  ENDIF

  a4(1:im,1,3) = a4(1:im,2,2)
  IF (kord > 3) a4(1:im,2,3) = a4(1:im,3,2)

  DO i=1,im
    ! Bottom
    ! Area preserving cubic with 2nd deriv. = 0 at the surface
      d1 = delp(i,km)
      d2 = delp(i,km-1)
      qm = (d2*a4(i,km,1)+d1*a4(i,km-1,1)) / (d1+d2)
      dq = 2.0_wp*(a4(i,km-1,1)-a4(i,km,1)) / (d1+d2)
      c1 = (a4(i,km-1,2)-qm-d2*dq) / (d2*(2.0_wp*d2*d2+d1*(d2+3.0_wp*d1)))
      c3 = dq - 2.0_wp*c1*(d2*(5.0_wp*d1+d2)-3.0_wp*d1**2)
      a4(i,km,2) = qm - c1*d1*d2*(d2+3.0_wp*d1)
      a4(i,km,3) = d1*(8.0_wp*c1*d1**2-c3) + a4(i,km,2)
      dc(i,km) = a4(i,km,3) -  a4(i,km,1)
      ! No over- and under-shoot condition
      cmax = MAX(a4(i,km,1), a4(i,km-1,1))
      cmin = MIN(a4(i,km,1), a4(i,km-1,1))
      a4(i,km,2) = MAX(cmin,a4(i,km,2))
      a4(i,km,2) = MIN(cmax,a4(i,km,2))
   END DO
    ! Enforce constraint at the surface

   IF ( iv == 0 ) THEN

      ! Positive definite scalars:
      a4(1:im,km,3) = MAX(0.0_wp, a4(1:im,km,3))

   ELSEIF ( iv ==  1 ) THEN

      ! Monotone tracers:
      dc(1:im,km) = 0.0_wp
      a4(1:im,km,2) = a4(1:im,km,1)
      a4(1:im,km,3) = a4(1:im,km,1)

   ELSEIF ( iv == -1 ) THEN

      DO i=1,im
         ! Winds:
         IF( a4(i,km,1)*a4(i,km,3) <=  0.0_wp ) THEN
            a4(i,km,3) = 0.0_wp
         ELSE
            a4(i,km,3) = SIGN( MIN(ABS(a4(i,km,1)),      &
                 ABS(a4(i,km,3))), a4(i,km,1)  )
         ENDIF
      END DO
   ENDIF

   IF (kord > 3) a4(1:im,km-2,3) = a4(1:im,km-1,2)
   a4(1:im,km-1,3) = a4(1:im,km,2)

    ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )

    ! Top 2 and bottom 2 layers always use monotonic mapping
!CDIR UNROLL=2
    DO k=1,2
      DO i=1,im
        a4(i,k,4) = 3.0_wp*(2.0_wp*a4(i,k,1) - (a4(i,k,2)+a4(i,k,3)))
        ! Standard PPM constraint
        IF(dc(i,k) == 0.0_wp) THEN
          a4(i,k,2) = a4(i,k,1)
          a4(i,k,3) = a4(i,k,1)
          a4(i,k,4) = 0.0_wp
        ELSE
          da1  = a4(i,k,3) - a4(i,k,2)
          da2  = da1**2
          a6da = a4(i,k,4)*da1
          IF(a6da < -da2) THEN
            a4(i,k,4) = 3.0_wp*(a4(i,k,2)-a4(i,k,1))
            a4(i,k,3) = a4(i,k,2) - a4(i,k,4)
          ELSEIF(a6da > da2) THEN
            a4(i,k,4) = 3.0_wp*(a4(i,k,3)-a4(i,k,1))
            a4(i,k,2) = a4(i,k,3) - a4(i,k,4)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    IF(kord >= 7) THEN

      !----------------------------------------
      ! Huynh's 2nd constraint
      !----------------------------------------

!      DO k=1,km-2
!        DO i=1,im
! Method#1
!           h2(i,k+1) = delq(i,k+1) - delq(i,k)
! Method#2
!           h2(i,k+1) = 2.*(dc(i,k+2)/delp(i,k+2) - dc(i,k)/delp(i,k)) &
!                   / ( delp(i,k+1)+0.5*(delp(i,k)+delp(i,k+2)) )        &
!                   * delp(i,k+1)**2
! Method#3
!           h2(i,k+1) = dc(i,k+2) - dc(i,k)
!        ENDDO
!      ENDDO

       IF( kord == 7 ) THEN
          fac = 1.5_wp           ! original quasi-monotone
       ELSE
          fac = 0.125_wp         ! full monotone
       ENDIF
       
       ! mpuetz: 50% of the time in ppm2m spent here
       DO k=3, km-2
          DO i=1,im

             ! Right edges
             !        qmp   = a4(i,k,1) + 2.0*delq(i,k-1)
             !        lac   = a4(i,k,1) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
             !
             pmp   = 2.0_wp*dc(i,k)
             qmp   = a4(i,k,1) + pmp
             !          lac   = a4(i,k,1) + fac*h2(i,k-1) + dc(i,k)
             lac   = a4(i,k,1) + fac*(dc(i,k) - dc(i,k-2)) + dc(i,k)
             qmin  = MIN(a4(i,k,1), qmp, lac)
             qmax  = MAX(a4(i,k,1), qmp, lac)
             a4(i,k,3) = MIN(MAX(a4(i,k,3), qmin), qmax)
             ! Left  edges
             !        qmp   = a4(i,k,1) - 2.0*delq(i,k)
             !        lac   = a4(i,k,1) + fac*h2(i,k+1) - 0.5*delq(i,k)
             !
             qmp   = a4(i,k,1) - pmp
             !          lac   = a4(i,k,1) + fac*h2(i,k+1) - dc(i,k)
             lac   = a4(i,k,1) + fac*(dc(i,k+2) - dc(i,k)) - dc(i,k)
             qmin  = MIN(a4(i,k,1), qmp, lac)
             qmax  = MAX(a4(i,k,1), qmp, lac)
             a4(i,k,2) = MIN(MAX(a4(i,k,2), qmin), qmax)
             ! Recompute A6
             a4(i,k,4) = 3.0_wp*(2.0_wp*a4(i,k,1) - (a4(i,k,2)+a4(i,k,3)))
          END DO

          ! Additional constraint to prevent negatives when kord=7

          IF (iv /= -1 .AND. kord == 7) THEN
             DO i=1,im
                ! Positive definite constraint
                ! mpuetz: depending on the likelyhood of this condition
                ! this loop could be expensive
                IF( ABS(a4(i,k,3)-a4(i,k,2)) < -a4(i,k,4) ) THEN
                   fmin = a4(i,k,1)+0.25_wp*(a4(i,k,3)-a4(i,k,2))**2/a4(i,k,4)+a4(i,k,4)*r12
                   IF( fmin < 0.0_wp ) THEN
                      IF(a4(i,k,1) < a4(i,k,3) .AND. a4(i,k,1) < a4(i,k,2)) THEN
                         a4(i,k,3) = a4(i,k,1)
                         a4(i,k,2) = a4(i,k,1)
                         a4(i,k,4) = 0.0_wp
                      ELSEIF(a4(i,k,3) > a4(i,k,2)) THEN
                         a4(i,k,4) = 3.0_wp*(a4(i,k,2)-a4(i,k,1))
                         a4(i,k,3) = a4(i,k,2) - a4(i,k,4)
                      ELSE
                         a4(i,k,4) = 3.0_wp*(a4(i,k,3)-a4(i,k,1))
                         a4(i,k,2) = a4(i,k,3) - a4(i,k,4)
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
       ENDDO

    ELSE

      lmt = kord - 3
      lmt = MAX(0, lmt)
      IF (iv == 0) lmt = MIN(2, lmt)

      IF( kord /= 4) THEN
        DO k=3, km-2
          DO i=1,im
            a4(i,k,4) = 3.0_wp*(2.0_wp*a4(i,k,1) - (a4(i,k,2)+a4(i,k,3)))
          ENDDO
        ENDDO
      ENDIF

      IF(lmt == 0) THEN

        DO k=3, km-2
         ! Standard PPM constraint
          DO i=1,im
            IF(dc(i,k) == 0.0_wp) THEN
              a4(i,k,2) = a4(i,k,1)
              a4(i,k,3) = a4(i,k,1)
              a4(i,k,4) = 0.0_wp
            ELSE
              da1  = a4(i,k,3) - a4(i,k,2)
              da2  = da1**2
              a6da = a4(i,k,4)*da1
              IF(a6da < -da2) THEN
                a4(i,k,4) = 3.0_wp*(a4(i,k,2)-a4(i,k,1))
                a4(i,k,3) = a4(i,k,2) - a4(i,k,4)
              ELSEIF(a6da > da2) THEN
                a4(i,k,4) = 3.0_wp*(a4(i,k,3)-a4(i,k,1))
                a4(i,k,2) = a4(i,k,3) - a4(i,k,4)
              ENDIF
            ENDIF
          ENDDO
        ENDDO

      ELSEIF (lmt == 1) THEN

        DO k=3, km-2
          ! Improved full monotonicity constraint (Lin)
          ! Note: no need to provide first guess of A6 <-- a4(i,k,4)
          DO i=1,im
            qmp = 2*dc(i,k)
            a4(i,k,2) = a4(i,k,1)-SIGN(MIN(ABS(qmp),ABS(a4(i,k,2)-a4(i,k,1))), qmp)
            a4(i,k,3) = a4(i,k,1)+SIGN(MIN(ABS(qmp),ABS(a4(i,k,3)-a4(i,k,1))), qmp)
            a4(i,k,4) = 3.0_wp*( 2.0_wp*a4(i,k,1) - (a4(i,k,2)+a4(i,k,3)) )
          ENDDO
        ENDDO

      ELSEIF (lmt == 2) THEN

        DO k=3, km-2
          ! Positive definite constraint
          DO i=1,im
            IF( ABS(a4(i,k,3)-a4(i,k,2)) < -a4(i,k,4) ) THEN
              fmin = a4(i,k,1)+0.25_wp*(a4(i,k,3)-a4(i,k,2))**2/a4(i,k,4)+a4(i,k,4)*r12
              IF( fmin < 0.0_wp ) THEN
                IF(a4(i,k,1) < a4(i,k,3) .AND. a4(i,k,1) < a4(i,k,2)) THEN
                  a4(i,k,3) = a4(i,k,1)
                  a4(i,k,2) = a4(i,k,1)
                  a4(i,k,4) = 0.0_wp
                ELSEIF(a4(i,k,3) > a4(i,k,2)) THEN
                  a4(i,k,4) = 3.0_wp*(a4(i,k,2)-a4(i,k,1))
                  a4(i,k,3) = a4(i,k,2) - a4(i,k,4)
                ELSE
                  a4(i,k,4) = 3.0_wp*(a4(i,k,3)-a4(i,k,1))
                  a4(i,k,2) = a4(i,k,3) - a4(i,k,4)
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO

      ENDIF

    ENDIF

!CDIR UNROLL=2
    DO k=km-1,km
      DO i=1,im
        a4(i,k,4) = 3.0_wp*(2.0_wp*a4(i,k,1) - (a4(i,k,2)+a4(i,k,3)))
        ! Standard PPM constraint
        IF(dc(i,k) == 0.0_wp) THEN
          a4(i,k,2) = a4(i,k,1)
          a4(i,k,3) = a4(i,k,1)
          a4(i,k,4) = 0.0_wp
        ELSE
          da1  = a4(i,k,3) - a4(i,k,2)
          da2  = da1**2
          a6da = a4(i,k,4)*da1
          IF(a6da < -da2) THEN
            a4(i,k,4) = 3.0_wp*(a4(i,k,2)-a4(i,k,1))
            a4(i,k,3) = a4(i,k,2) - a4(i,k,4)
          ELSEIF(a6da > da2) THEN
            a4(i,k,4) = 3.0_wp*(a4(i,k,3)-a4(i,k,1))
            a4(i,k,2) = a4(i,k,3) - a4(i,k,4)
          ENDIF
        ENDIF
      ENDDO
    ENDDO

#ifdef _PROFILE
    CALL trace_stop ('ppm2m', 5)
#endif

  END SUBROUTINE ppm2m

  SUBROUTINE steepz(kbdim, im, km, a4, df2, dm, dq, dp, d4, d4i, d3i)

    INTEGER, INTENT(in) :: km                        ! Total levels
    INTEGER, INTENT(in) :: kbdim                    ! real leading dimension
    INTEGER, INTENT(in) :: im                       ! actual leading dimension
    REAL(wp), INTENT(in) :: dp(kbdim,km)            ! grid size
    REAL(wp), INTENT(in) :: dq(kbdim,km)            ! backward diff of q
    REAL(wp), INTENT(in) :: d4(kbdim,km)            ! backward sum:  dp(k)+ dp(k-1)
    REAL(wp), INTENT(in) :: d4i(kbdim,km)           ! inverse backward sum:  dp(k)+ dp(k-1)
    REAL(wp), INTENT(in) :: d3i(kbdim,km)           ! inverse backward sum:  dp(k+1) + dp(k) + dp(k-1)
    REAL(wp), INTENT(in) :: df2(kbdim,km)            ! first guess mismatch
    REAL(wp), INTENT(in) :: dm(kbdim,km)            ! monotonic mismatch

    REAL(wp), INTENT(inout) :: a4(kbdim,km,4)          ! first guess/steepened

    ! Local variables:
    INTEGER :: i, k
    REAL(wp) :: alfa(kbdim,km)
    REAL(wp) :: f(kbdim,km)
    REAL(wp) :: dg2

    ! Compute ratio of dq/dp
    ! Compute F
!CDIR EXPAND
    DO k=1,3
       DO i=1,im
          f(i,k+1) = (dq(i,k+1) * d4i(i,k+2) - dq(i,k) * d4i(i,k+1)) * d3i(i,k+1)
       END DO
    ENDDO
    k=3
    DO i=1,im
      IF(f(i,k+1)*f(i,k-1) < 0.0_wp .AND. df2(i,k) /= 0.0_wp) THEN
        dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2      &
             + d4(i,k)*d4(i,k+1) )
        alfa(i,k) = MAX(0.0_wp, MIN(0.5_wp, -0.1875_wp*dg2/df2(i,k)))
      ELSE
        alfa(i,k) = 0.0_wp
      ENDIF
    ENDDO
!CDIR OUTERUNROLL=8
    DO k=4,km-2
      DO i=1,im
        f(i,k+1) = (dq(i,k+1) * d4i(i,k+2) - dq(i,k) * d4i(i,k+1)) * d3i(i,k+1)
        IF(f(i,k+1)*f(i,k-1) < 0.0_wp .AND. df2(i,k) /= 0.0_wp) THEN
          dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2      &
               + d4(i,k)*d4(i,k+1) )
          alfa(i,k) = MAX(0.0_wp, MIN(0.5_wp, -0.1875_wp*dg2/df2(i,k)))
        ELSE
          alfa(i,k) = 0.0_wp
        ENDIF
        a4(i,k,2) = (1.0_wp-alfa(i,k-1)-alfa(i,k)) * a4(i,k,2) +      &
             alfa(i,k-1)*(a4(i,k,1)-dm(i,k))    +          &
             alfa(i,k)*(a4(i,k-1,1)+dm(i,k-1))
        a4(i,k-1,3) = a4(i,k,2)
      ENDDO
    ENDDO

  END SUBROUTINE steepz

#ifdef VECTOR

  SUBROUTINE fct_x(q, im, jm, jfirst, jlast, km, nq, ipx, qtmp)

    INTEGER, INTENT(in) :: im                  ! Longitudes
    INTEGER, INTENT(in) :: jm                  ! Total latitudes
    INTEGER, INTENT(in) :: km                  ! Number of levels
    INTEGER, INTENT(in) :: nq                  ! Number tracers

    INTEGER, INTENT(in) :: jfirst              ! Starting latitude
    INTEGER, INTENT(in) :: jlast               ! Finishing latitude

    REAL(wp), INTENT(inout) :: q(im,jfirst:jlast,km,nq) ! Field to adjust

    INTEGER, INTENT(out) :: ipx   ! Flag:  0 if Q not change, 1 if changed

    ! +1 in first index for explicit loop bank conflict prevention ...
    REAL(wp), INTENT(inout) :: qtmp((((jlast-jfirst+1)*km*nq)/2)*2+1,im)

    ! Local variables:
    REAL(wp), PARAMETER :: tiny = 1.0e-40_wp  ! A small number to pump up value
    REAL(wp) :: d0, d1, d2

    INTEGER :: i, j, k, iq, jm1, ip2, jmt
    INTEGER :: j1, j2, icnt
    INTEGER :: n, nl, nthreads, mythread, js, je, ks, ke
   
#ifdef _OPENMP
    nthreads = omp_get_num_threads()
    mythread = omp_get_thread_num()
#else
    nthreads = 1
    mythread = 0
#endif

    j1 = MAX( jfirst,   2 )
    j2 = MIN( jlast, jm-1 )
    jm1 = jm-1
    jmt = (jlast-jfirst+1)
    ipx = 0

    n=jmt*km*nq
    nl=(n+nthreads-1)/nthreads
    js=mythread*nl+1
    je=MIN(mythread*nl+nl,n)

    n=km
    nl=(n+nthreads-1)/nthreads
    ks=mythread*nl+1
    ke=MIN(mythread*nl+nl,n)

    ! Copy & swap direction for vectorization.

!CDIR NODEP
    DO iq = 1, nq
!CDIR NODEP
      DO k = ks, ke
        DO i = 1, im
          qtmp((((iq-1)*km+(k-1))*jmt)+1,i) = 0.0_wp
          qtmp((((iq-1)*km+(k-1))*jmt)+jlast-jfirst+1,i) = 0.0_wp
        ENDDO
        DO j = j1, j2
          DO i = 1, im
            qtmp((((iq-1)*km+(k-1))*jmt)+(j-jfirst+1),i) = q(i,j,k,iq)  
          ENDDO
        ENDDO
      ENDDO
    ENDDO
!$OMP BARRIER

    icnt=0
    loop1: DO i=1,im
      DO j=js,je
        IF(qtmp(j,i) < 0.0_wp) THEN
          icnt=1
          exit loop1
        ENDIF
      ENDDO
    ENDDO loop1
    IF (icnt == 1) THEN
      DO i=2,im-1
        DO j=js,je
          IF(qtmp(j,i) < 0.0_wp) THEN
            ipx =  1
            ! west
            d0 = MAX(0.0_wp,qtmp(j,i-1))
            d1 = MIN(-qtmp(j,i),d0)
            qtmp(j,i-1) = qtmp(j,i-1) - d1
            qtmp(j,i) = qtmp(j,i) + d1
            ! east
            d0 = MAX(0.0_wp,qtmp(j,i+1))
            d2 = MIN(-qtmp(j,i),d0)
            qtmp(j,i+1) = qtmp(j,i+1) - d2
            qtmp(j,i) = qtmp(j,i) + d2 + tiny
          ENDIF
        ENDDO
      ENDDO
    ENDIF
!$OMP BARRIER
    
    i=1
    DO j=js,je
      IF(qtmp(j,i) < 0.0_wp) THEN
        ipx =  1
        ! west
        d0 = MAX(0.0_wp,qtmp(j,im))
        d1 = MIN(-qtmp(j,i),d0)
        qtmp(j,im) = qtmp(j,im) - d1
        qtmp(j,i) = qtmp(j,i) + d1
        ! east
        d0 = MAX(0.0_wp,qtmp(j,i+1))
        d2 = MIN(-qtmp(j,i),d0)
        qtmp(j,i+1) = qtmp(j,i+1) - d2
        qtmp(j,i) = qtmp(j,i) + d2 + tiny
      ENDIF
    ENDDO
!$OMP BARRIER
    
    i=im
    DO j=js,je
      IF(qtmp(j,i) < 0.0_wp) THEN
        ipx =  1
        ! west
        d0 = MAX(0.0_wp,qtmp(j,i-1))
        d1 = MIN(-qtmp(j,i),d0)
        qtmp(j,i-1) = qtmp(j,i-1) - d1
        qtmp(j,i) = qtmp(j,i) + d1
        ! east
        d0 = MAX(0.0_wp,qtmp(j,1))
        d2 = MIN(-qtmp(j,i),d0)
        qtmp(j,1) = qtmp(j,1) - d2
        
        qtmp(j,i) = qtmp(j,i) + d2 + tiny
      ENDIF
    ENDDO
!$OMP BARRIER
    
    IF(ipx /= 0) THEN
      !-----------
      ! Final pass
      !-----------

    icnt=0
    loop2: DO i=1,im
      DO j=js,je
        IF(qtmp(j,i) < 0.0_wp) THEN
          icnt=1
          exit loop2
        ENDIF
      ENDDO
    ENDDO loop2
    IF (icnt == 1) THEN
      DO i=1,im-1
        DO j=js,je
          IF(qtmp(j,i) < 0.0_wp) THEN
            ! Take mass from east (essentially adjusting fx(i+1,j))
            qtmp(j,i+1) = qtmp(j,i+1) + qtmp(j,i)
            qtmp(j,i) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO
    ENDIF
!$OMP BARRIER

      ! Final sweep
    icnt=0
    loop3: DO i=1,im
      DO j=js,je
        IF(qtmp(j,i) < 0.0_wp) THEN
          icnt=1
          exit loop3
        ENDIF
      ENDDO
    ENDDO loop3
    IF (icnt == 1) THEN
      DO i=im,2,-1
        DO j=js,je
          IF(qtmp(j,i) < 0.0_wp) THEN
            ! Take mass from west (essentially adjusting fx(i,j))
            qtmp(j,i-1) = qtmp(j,i-1) + qtmp(j,i)
            qtmp(j,i) = 0.0_wp
          ENDIF
        ENDDO
        ! Note: qtmp(j,1) could still be negative
      ENDDO
    ENDIF
!$OMP BARRIER

    ENDIF

    ! transpose back ...

    DO iq = 1, nq
      DO k = ks, ke
        DO j = j1, j2
          DO i = 1, im
            q(i,j,k,iq) = qtmp((((iq-1)*km+(k-1))*jmt)+(j-jfirst+1),i)  
          ENDDO
        ENDDO
      ENDDO
    ENDDO
!$OMP BARRIER

    Multi_Tracer: DO iq=1,nq

      Vertical:  DO k=ks,ke

        ! Check Poles.
        IF ( jfirst == 1 ) THEN
          ip2 = 0
          ! SP
          IF(q(1,1,k,iq) < 0.0_wp) THEN
            CALL pfix(q(1,2,k,iq),q(1,1,k,iq),im,ipx)
          ELSE
            ! Check j=2
            DO i=1,im
              IF(q(i,2,k,iq) < 0.0_wp) THEN
                ip2 = 1
                EXIT
              ENDIF
            ENDDO
            IF(ip2 /= 0) CALL pfix(q(1,2,k,iq),q(1,1,k,iq),im,ipx)
          ENDIF
        ENDIF

        IF ( jlast == jm ) THEN
          ip2 = 0
          ! NP
          IF(q(1,jm,k,iq) < 0.0_wp) THEN
            CALL pfix(q(1,jm1,k,iq),q(1,jm,k,iq),im,ipx)
          ELSE
            ! Check j=jm1
            DO i=1,im
              IF(q(i,jm1,k,iq) < 0.0_wp) THEN
                ip2 = 1
                EXIT
              ENDIF
            ENDDO
            IF(ip2 /= 0) CALL pfix(q(1,jm1,k,iq),q(1,jm,k,iq),im,ipx)
          ENDIF
        ENDIF

      END DO Vertical

    END DO Multi_Tracer
!$OMP BARRIER

  END SUBROUTINE fct_x

#else

  SUBROUTINE fct_x(q, im, jm, jfirst, jlast, jlast_odd, ng, ipx)

    INTEGER, INTENT(in) :: im                  ! Longitudes
    INTEGER, INTENT(in) :: jm                  ! Total latitudes
    INTEGER, INTENT(in) :: jfirst              ! Starting latitude
    INTEGER, INTENT(in) :: jlast               ! Finishing latitude
    INTEGER, INTENT(in) :: jlast_odd
    INTEGER, INTENT(in) :: ng

    REAL(wp), INTENT(inout) :: q(im,jfirst-ng:jlast+ng) ! Field to adjust

    INTEGER, INTENT(out) :: ipx   ! Flag:  0 if Q not change, 1 if changed

    ! Local variables:
    REAL(wp), PARAMETER :: tiny = 1.0e-40_wp  ! A small number to pump up value
    REAL(wp) :: d0, d1, d2
    REAL(wp) :: qtmp(jfirst:jlast_odd,im)

    INTEGER :: i, j, jm1, ip2
    INTEGER :: j1, j2, icnt

    j1 = MAX( jfirst,   2 )
    j2 = MIN( jlast, jm-1 )
    jm1 = jm-1
    ipx = 0

    qtmp(:,:)=0.0_wp

    ! Copy & swap direction for vectorization.

    DO i=1,im
      DO j=j1,j2
        qtmp(j,i) = q(i,j)
      ENDDO
    ENDDO

    icnt=0
    IF (ANY(qtmp(:,:) < 0.0_wp)) icnt=1

    IF (icnt > 0) THEN
      DO i=2,im-1
        DO j=j1,j2
          IF(qtmp(j,i) < 0.0_wp) THEN
            ipx =  1
            ! west
            d0 = MAX(0.0_wp,qtmp(j,i-1))
            d1 = MIN(-qtmp(j,i),d0)
            qtmp(j,i-1) = qtmp(j,i-1) - d1
            qtmp(j,i) = qtmp(j,i) + d1
            ! east
            d0 = MAX(0.0_wp,qtmp(j,i+1))
            d2 = MIN(-qtmp(j,i),d0)
            qtmp(j,i+1) = qtmp(j,i+1) - d2
            qtmp(j,i) = qtmp(j,i) + d2 + tiny
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    i=1
    DO j=j1,j2
      IF(qtmp(j,i) < 0.0_wp) THEN
        ipx =  1
        ! west
        d0 = MAX(0.0_wp,qtmp(j,im))
        d1 = MIN(-qtmp(j,i),d0)
        qtmp(j,im) = qtmp(j,im) - d1
        qtmp(j,i) = qtmp(j,i) + d1
        ! east
        d0 = MAX(0.0_wp,qtmp(j,i+1))
        d2 = MIN(-qtmp(j,i),d0)
        qtmp(j,i+1) = qtmp(j,i+1) - d2
        qtmp(j,i) = qtmp(j,i) + d2 + tiny
      ENDIF
    ENDDO

    i=im
    DO j=j1,j2
      IF(qtmp(j,i) < 0.0_wp) THEN
        ipx =  1
        ! west
        d0 = MAX(0.0_wp,qtmp(j,i-1))
        d1 = MIN(-qtmp(j,i),d0)
        qtmp(j,i-1) = qtmp(j,i-1) - d1
        qtmp(j,i) = qtmp(j,i) + d1
        ! east
        d0 = MAX(0.0_wp,qtmp(j,1))
        d2 = MIN(-qtmp(j,i),d0)
        qtmp(j,1) = qtmp(j,1) - d2

        qtmp(j,i) = qtmp(j,i) + d2 + tiny
      ENDIF
    ENDDO


    IF(ipx /= 0) THEN
      !-----------
      ! Final pass
      !-----------

      icnt=0
      IF (ANY(qtmp(:,:) < 0.0_wp)) icnt=1
      
      IF (icnt > 0) THEN
        DO i=1,im-1
          DO j=j1,j2
            IF (qtmp(j,i) < 0.0_wp ) THEN
              ! Take mass from east (essentially adjusting fx(i+1,j))
              qtmp(j,i+1) = qtmp(j,i+1) + qtmp(j,i)
              qtmp(j,i) = 0.0_wp
            ENDIF
          ENDDO
        ENDDO
        
        icnt=0
        IF (ANY(qtmp(:,:) < 0.0_wp)) icnt=1
        
        ! Final sweep
        IF (icnt > 0) THEN
          DO i=im,2,-1
            DO j=j1,j2
              IF (qtmp(j,i) < 0.0_wp ) THEN
                ! Take mass from west (essentially adjusting fx(i,j))
                qtmp(j,i-1) = qtmp(j,i-1) + qtmp(j,i)
                qtmp(j,i) = 0.0_wp
              ENDIF
            ENDDO
            ! Note: qtmp(j,1) could still be negative
          ENDDO
        ENDIF
      ENDIF

      DO j=j1,j2
        DO i=1,im
          q(i,j) = qtmp(j,i)
          !         q(i,j) = max(0., qtmp(j,i))
        ENDDO
      ENDDO

    ENDIF

    ! Check Poles.
    IF ( jfirst == 1 ) THEN
      ip2 = 0
      ! SP
      IF(q(1,1) < 0.0_wp) THEN
        CALL pfix(q(1,2),q(1,1),im,ipx)
      ELSE
        ! Check j=2
        DO i=1,im
          IF(q(i,2) < 0.0_wp) THEN
            ip2 = 1
            EXIT
          ENDIF
        ENDDO
        IF(ip2 /= 0) CALL pfix(q(1,2),q(1,1),im,ipx)
      ENDIF
    ENDIF

    IF ( jlast == jm ) THEN
      ip2 = 0
      ! NP
      IF(q(1,jm) < 0.0_wp) THEN
        CALL pfix(q(1,jm1),q(1,jm),im,ipx)
      ELSE

        ! Check j=jm1
        DO i=1,im
          IF(q(i,jm1) < 0.0_wp) THEN
            ip2 = 1
            EXIT
          ENDIF
        ENDDO
        IF(ip2 /= 0) CALL pfix(q(1,jm1),q(1,jm),im,ipx)
      ENDIF
    ENDIF
  END SUBROUTINE fct_x

#endif

  SUBROUTINE fillz(im, i1, i2, km, nq, q, dp)

    INTEGER, INTENT(in) :: im                ! No. of longitudes
    INTEGER, INTENT(in) :: km                ! No. of levels
    INTEGER, INTENT(in) :: i1                ! Starting longitude
    INTEGER, INTENT(in) :: i2                ! Finishing longitude
    INTEGER, INTENT(in) :: nq                ! Total number of tracers
    REAL(wp), INTENT(in) ::  dp(im,km)       ! pressure thickness

    REAL(wp), INTENT(inout) :: q(im,km,nq)   ! tracer mixing ratio

    ! Local variables:
    INTEGER :: i, k, ic
    REAL(wp) :: qup, qly, dup

    DO ic=1,nq
      ! Top layer
      DO i=i1,i2
        IF( q(i,1,ic) < 0.0_wp) THEN
          q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
          q(i,1,ic) = 0.0_wp
        ENDIF
      ENDDO

      ! Interior
      DO k=2,km-1
        DO i=i1,i2
          IF( q(i,k,ic) < 0.0_wp ) THEN
            ! Borrow from above
            qup =  q(i,k-1,ic)*dp(i,k-1)
            qly = -q(i,k  ,ic)*dp(i,k  )
            dup =  MIN( 0.5_wp*qly, qup )        !borrow no more than 50%
            q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1)
            ! Borrow from below: q(i,k,ic) is still negative at this stage
            q(i,k+1,ic) = q(i,k+1,ic) + (dup-qly)/dp(i,k+1)
            q(i,k  ,ic) = 0.0_wp
          ENDIF
        ENDDO
      ENDDO

      ! Bottom layer
      k = km
      DO i=i1,i2
        IF( q(i,k,ic) < 0.0_wp) THEN
          ! Borrow from above
          qup =  q(i,k-1,ic)*dp(i,k-1)
          qly = -q(i,k  ,ic)*dp(i,k  )
          dup =  MIN( qly, qup )
          q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1)
          q(i,k,ic) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE fillz

  SUBROUTINE pfix(q, qp, im, ipx)

    INTEGER, INTENT(in) :: im        ! Longitudes

    REAL(wp), INTENT(inout) :: q(im)    ! Latitude-level field to adjust
    REAL(wp), INTENT(inout) :: qp(im)   ! Second latitude-level field to
                                        ! adjust (usually pole)

    INTEGER, INTENT(out) :: ipx   ! Flag:  0 if Q not change, 1 if changed

    ! Local variables:
    INTEGER :: i
    REAL(wp) :: summ, sump, pmean

    summ = 0.0_wp
    sump = 0.0_wp
    DO i=1,im
      summ = summ + q(i)
      sump = sump + qp(i)
    ENDDO

    pmean = (sump*gw(1) + summ*gw(2)) / (im*(gw(1)+gw(2)))

    DO i=1,im
      q(i) = pmean
      qp(i) = pmean
    ENDDO

    IF( qp(1) < 0.0_wp )  ipx = 1

  END SUBROUTINE pfix

  SUBROUTINE adj_fx(im, jm, km, jfirst, jlast, ak, bk, ffsl,  &
                     ps0, ps2, delp, fx3, cx, fy3, ng,     &
                     mg, tiny, n_adj)

    INTEGER, INTENT(in):: im
    INTEGER, INTENT(in):: jm
    INTEGER, INTENT(in):: km
    INTEGER, INTENT(in):: ng, mg
    INTEGER, INTENT(in):: jfirst, jlast
    INTEGER, INTENT(in):: n_adj
    REAL(wp), INTENT(in):: tiny
    REAL(wp), INTENT(in)::  ak(km+1)
    REAL(wp), INTENT(in)::  bk(km+1)
    REAL(wp), INTENT(in)::  ps2(im,jfirst-mg:jlast+mg)
    REAL(wp), INTENT(in)::   cx(im,jfirst-ng:jlast+ng,klev)

    LOGICAL, INTENT(in):: ffsl(jfirst-ng:jlast+ng,klev)

    REAL(wp), INTENT(inout):: ps0(im,jfirst:jlast)
    REAL(wp), INTENT(inout):: fx3(im,jfirst:jlast,klev)
    REAL(wp), INTENT(inout):: fy3(im,jfirst:jlast+mg,klev)
    REAL(wp), INTENT(inout):: delp(im,jfirst:jlast,klev)

    ! Local
    REAL(wp), PARAMETER :: fac = 0.25_wp
    REAL(wp) :: ps(im,jfirst-mg:jlast+mg)
    REAL(wp) :: fy(im,jfirst:jlast+mg)
    REAL(wp) :: fx(im+1)
    REAL(wp) :: dps(0:im)
    REAL(wp) :: dpy(im,jfirst-mg:jlast+mg)
    REAL(wp) :: er0, er1, er2
    INTEGER :: i,j,k,k2, it
    REAL(wp) :: tmpf, dh, dh1, dh2, tmp
    REAL(wp) :: dbk
    INTEGER :: js2g0, jn2g0

    js2g0  = MAX(2,jfirst)        ! No ghosting
    jn2g0  = MIN(jm-1,jlast)      ! No ghosting

!$OMP DO
    DO j=jfirst,jlast
      DO i=1,im
        ps(i,j) = ps0(i,j)
      ENDDO
    ENDDO
!$OMP END DO

    fx_iteration: DO it=1,n_adj

#ifndef NOMPI
!$OMP MASTER
      CALL ghost_update(ps, im, jm, 1, 1, jfirst, jlast, mg)
!$OMP END MASTER
!$OMP BARRIER
#endif

      !--- adjust fx ----
!$OMP DO
      DO k=1,klev
        k2=(k-1)*kstep+kstart
        IF(k<3) CYCLE
        dbk = bk(k2+1) - bk(k2)
        IF( dbk > 0.001_wp ) THEN
!CDIR NODEP
          DO j=js2g0,jn2g0
            DO i=1,im
              dps(i) = (ps(i,j) - ps2(i,j))*dbk
            ENDDO
            dps(0) = dps(im)
            DO i=1,im
              fx(i) = fac*(dps(i-1)-dps(i))
              tmpf = fx3(i,j,k) + fx(i)
              IF ( tmpf*fx3(i,j,k) > 0.0_wp ) THEN
                fx3(i,j,k) = tmpf
              ELSE
                fx(i)  =  fx3(i,j,k)
                fx3(i,j,k) = SIGN(MIN(ABS(tmpf), ABS(fx3(i,j,k))), fx3(i,j,k))
                fx(i) = fx3(i,j,k) - fx(i)
              ENDIF
            ENDDO
            fx(im+1) = fx(1)

            ! update delp
            DO i=1,im
              delp(i,j,k) = delp(i,j,k) + fx(i) - fx(i+1)
            ENDDO
          ENDDO     ! j-loop

          !--- adjust fy ----

          DO j=MAX(jfirst-1,1) ,MIN(jm,jlast+1)     ! Need ps at jlast+1
            DO i=1,im
              dpy(i,j) = (ps(i,j) - ps2(i,j))*dbk*gw(j)
            ENDDO
          ENDDO

          DO j=js2g0,MIN(jm,jlast+1)
            DO i=1,im
              fy(i,j) = fac*(dpy(i,j-1)-dpy(i,j))
              tmpf = fy3(i,j,k) + fy(i,j)
              IF ( tmpf*fy3(i,j,k) > 0.0_wp ) THEN
                fy3(i,j,k) = tmpf
              ELSE
                fy(i,j)  =  fy3(i,j,k)
                fy3(i,j,k) = SIGN(MIN(ABS(tmpf), ABS(fy3(i,j,k))), fy3(i,j,k))
                fy(i,j) = fy3(i,j,k) - fy(i,j)
              ENDIF
            ENDDO
          ENDDO

          ! update delp
          DO j=js2g0,jn2g0
            DO i=1,im
              delp(i,j,k) = delp(i,j,k) + (fy(i,j) - fy(i,j+1)) * rgw(j)
            ENDDO
          ENDDO

          ! Poles:
          IF ( jfirst == 1 ) THEN
            DO i=1,im
              delp(i,1,k) = delp(i,1,k) - fy(i,2)*rgw(1)
            ENDDO
            CALL xpavg(delp(1,1,k), im)
          ENDIF
          IF ( jlast == jm ) THEN
            DO i=1,im
              delp(i,jm,k) = delp(i,jm,k) + fy(i,jm)*rgw(jm)
            ENDDO
            CALL xpavg(delp(1,jm,k), im)
          ENDIF

        ENDIF
      ENDDO            ! k-loop
!$OMP END DO

      ! Update ps

#ifndef NOMPI
!$OMP MASTER
      CALL bcast_delp(delp)
!$OMP END MASTER
!$OMP BARRIER
#endif

!$OMP DO
      DO j=jfirst,jlast ! RJ was: js2g0,jn2g0
        DO i=1,im
          ps(i,j) = ak(1)
        ENDDO

        DO k=1,klev
          DO i=1,im
            ps(i,j) = ps(i,j) + delp(i,j,k)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO

    ENDDO fx_iteration

!$OMP DO
    DO k=1,klev
      k2=(k-1)*kstep+kstart
      DO j=js2g0,jn2g0
        IF(K>=3) THEN
          dbk = bk(k2+1) - bk(k2)
          IF( dbk > 0.001_wp ) THEN
            DO i=1,im
              dps(i) = (ps(i,j) - ps2(i,j))*dbk
            ENDDO
            dps(0) = dps(im)
            !
            i=1
            er0 =  dps(i-1)
            er1 =  dps(i)
            IF( er0*er1 < 0.0_wp ) THEN
              IF( er1 > 0.0_wp ) THEN
                dh = MIN(-er0, er1)
                fx3(i,j,k) = fx3(i,j,k) - dh
                delp(im,j,k) = delp(im,j,k) + dh
                delp(i,j,k) = delp(i,j,k) - dh
              ELSE
                dh = MIN(er0, -er1)
                fx3(i,j,k) = fx3(i,j,k) + dh
                delp(im,j,k) = delp(im,j,k) - dh
                delp(i,j,k) = delp(i,j,k) + dh
              ENDIF
            ENDIF

            er1 =  dps(1)
            er2 =  dps(2)
            IF( er1*er2 < 0.0_wp ) THEN
              IF( er2 > 0.0_wp ) THEN
                dh2 = MIN(-er1, er2)
                delp(i,j,k) = delp(i,j,k) + dh2
              ELSE
                dh2 = MIN(er1, -er2)
                delp(i,j,k) = delp(i,j,k) - dh2
              ENDIF
            ENDIF
            DO i=2,im-1
              er0 =  dps(i-1)
              er1 =  dps(i)
              er2 =  dps(i+1)
              tmp = 0.0_wp
              IF( er1*er2 < 0.0_wp ) THEN
                IF( er2 > 0.0_wp ) THEN
                  dh2 = MIN(-er1, er2)
                  tmp = tmp + dh2
                ELSE
                  dh2 = MIN(er1, -er2)
                  tmp = tmp - dh2
                ENDIF
              ENDIF
              IF( er0*er1 < 0.0_wp ) THEN
                IF( er1 > 0.0_wp ) THEN
                  dh1 = MIN(-er0, er1)
                  fx3(i,j,k) = fx3(i,j,k) - dh1
                  tmp = tmp - dh1
                ELSE
                  dh1 = MIN(er0, -er1)
                  fx3(i,j,k) = fx3(i,j,k) + dh1
                  tmp = tmp + dh1
                ENDIF
              ENDIF
              delp(i,j,k) = delp(i,j,k) + tmp
            ENDDO
            er0 =  dps(im-1)
            er1 =  dps(im)
            IF( er0*er1 < 0.0_wp ) THEN
              IF( er1 > 0.0_wp ) THEN
                dh1 = MIN(-er0, er1)
                fx3(im,j,k) = fx3(im,j,k) - dh1
                delp(im,j,k) = delp(im,j,k) - dh1
              ELSE
                dh1 = MIN(er0, -er1)
                fx3(im,j,k) = fx3(im,j,k) + dh1
                delp(im,j,k) = delp(im,j,k) + dh1
              ENDIF
            ENDIF
          ENDIF
        ENDIF ! K>=3
        IF( ffsl(j,k) ) THEN
          DO i=1,im
            fx3(i,j,k) = fx3(i,j,k)/SIGN(MAX(ABS(cx(i,j,k)),tiny),cx(i,j,k))
          ENDDO
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO

  END SUBROUTINE adj_fx

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

#ifndef NOMPI

  !-------------------------------------------------------------------------
  SUBROUTINE bcast_delp(delp)

    USE mo_mpi, ONLY: p_irecv, p_send, p_wait

    REAL(wp) :: delp(:,:,:)
    INTEGER :: src, k, n, idx

    IF (ldc%nprocb==1) RETURN ! Nothing to do

    DO k=1,plev
      IF(MOD(k-1,kstep)+1 == ldc%set_b) CYCLE ! Do not send to myself
      src = gdc(ldc%mapmesh(MOD(k-1,kstep)+1,ldc%set_a))%pe
      CALL p_irecv(delp(1,1,k),src,1000+k,ldc%nlon*ldc%ffsl%nlat)
    ENDDO

    DO k=kstart,plev,kstep
      DO n=1,ldc%nprocb
        IF(n  == ldc%set_b) CYCLE ! Do not send to myself
        idx = ldc%mapmesh(n,ldc%set_a)
        CALL p_send(delp(:,:,k),gdc(idx)%pe,1000+k)
      ENDDO
    ENDDO

    CALL p_wait

  END SUBROUTINE bcast_delp

  !-------------------------------------------------------------
  !-- RJ: The original ghost_update used a bad communication pattern,
  !       mpi_sendrecv is replaced by p_isend/p_recv in the following version
  SUBROUTINE ghost_update(q, im, jm, km, nq, jfirst, jlast, ng)

    USE mo_mpi, ONLY: p_isend, p_irecv, p_wait
    
    ! No buffer needed.
    
    INTEGER, INTENT(in):: im, jm, km, nq, jfirst, jlast
    INTEGER, INTENT(in):: ng         ! zones to be ghosted
    REAL(wp) :: q(im,jfirst-ng:jlast+ng,km,nq)
    INTEGER :: qsize

    REAL(wp) :: q_send1(im,ng,km,nq)
    REAL(wp) :: q_send2(im,ng,km,nq) 
    REAL(wp) :: q_recv1(im,ng,km,nq) 
    REAL(wp) :: q_recv2(im,ng,km,nq)

    qsize = im*ng*km*nq

!!$    IF ( jfirst > 1 ) CALL  p_irecv(q_recv1(1,1,1,1),ldc%ffsl%pe_s,9876,qsize) 
!!$    IF ( jlast < jm ) CALL  p_irecv(q_recv2(1,1,1,1),ldc%ffsl%pe_n,1234,qsize)
!!$
!!$    IF ( jfirst > 1 ) THEN
!!$      q_send1(:,1:ng,:,:) = q(:,jfirst:jfirst+ng-1,:,:)
!!$      CALL  p_isend(q_send1(1,1,1,1),ldc%ffsl%pe_s,1234,qsize)
!!$    END IF
!!$
!!$    IF ( jlast < jm ) THEN
!!$      q_send2(:,1:ng,:,:) = q(:,jlast-ng+1:jlast,:,:)
!!$      CALL  p_isend(q_send2(1,1,1,1),ldc%ffsl%pe_n,9876,qsize)
!!$    END IF

    IF ( jfirst > 1 ) CALL  p_irecv(q_recv1,ldc%ffsl%pe_s,9876,qsize) 
    IF ( jlast < jm ) CALL  p_irecv(q_recv2,ldc%ffsl%pe_n,1234,qsize)

    IF ( jfirst > 1 ) THEN
      q_send1(:,1:ng,:,:) = q(:,jfirst:jfirst+ng-1,:,:)
      CALL  p_isend(q_send1,ldc%ffsl%pe_s,1234,qsize)
    END IF

    IF ( jlast < jm ) THEN
      q_send2(:,1:ng,:,:) = q(:,jlast-ng+1:jlast,:,:)
      CALL  p_isend(q_send2,ldc%ffsl%pe_n,9876,qsize)
    END IF

    CALL p_wait

    IF ( jfirst > 1 )  q(:,jfirst-ng:jfirst-1,:,:) = q_recv1(:,1:ng,:,:)
    IF ( jlast < jm )  q(:,jlast+1:jlast+ng,:,:)   = q_recv2(:,1:ng,:,:)

  END SUBROUTINE ghost_update
  
#endif

END MODULE mo_tpcore



