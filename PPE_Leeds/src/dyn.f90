#ifdef __xlC__
@PROCESS HOT
@PROCESS ALIAS(NOARYOVRLP,NOPTEOVRLP)
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE dyn(petadot)

  ! Description:
  !
  ! Computes adiabatic tendencies and auxilliary hybrid variables.
  !
  ! Method:
  !
  ! The primary purpose is to compute adiabatic tendencies.
  ! Auxilliary variables connected with the vertical difference
  ! scheme are saved as they are subsequently required to compute
  ! input to physical parameterizations.
  !
  ! *dyn* is called from *gpc*. 
  ! Input is from long-term storage and modules *mo_hyb* and *mo_gaussgrid*
  ! Output is to long-term storage.
  !
  ! Externals:
  ! *pres* and *auxhyb* are called to calculate auxilliary variables.
  ! *geopot* is called to calculate full-level geopotentials.
  ! *locate*, *alloc* and *unloc* are called to manage storage.
  !
  ! *External documentation of the model equations, of
  ! the organization of the vertical calculation, and of the
  ! organization of the spectral horizontal calculation.
  !
  ! Authors:
  !
  ! A. J. Simmons, ECMWF, January 1982, original source
  ! U. Schlese, DKRZ, June 1991, changes for new advection scheme
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, August 1998, tendency diagnostics
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! L. Kornblueh, MPI, Februrary 2012, remove spitfire
  ! S. K. Cheedela, Jan 2010, options for single column model
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind,          ONLY: dp
  USE mo_scan_buffer,   ONLY: alnpr, alpha, alps, alpste, d, qte, rh,  &
                              t, tte, u, v, vo, vervel, vol, vom,      &
                              xite, xlte, dalpsl, dalpsm,              &
                              dtl, dtm, xtte
  USE mo_memory_gl,     ONLY: q, xl, xi
  USE mo_memory_g3a,    ONLY: apsm, geospm
  USE mo_memory_g3b,    ONLY: aps
  USE mo_control,       ONLY: nlev, nlevp1, nvclev, vct, ltimer,       &
                              lcolumn
  USE mo_geoloc,        ONLY: coriol_2d, cst_2d, rcsth_2d
  USE mo_hyb,           ONLY: alpham, alrrdic, ardprc, cetah, cpg,     &
                              delb, delpr, nlevm1, nlmsgl, nlmslp,     &
                              nplev, nplvp1, ralpha, rddelb, rdlnp0i,  &
                              rdt0ral, rlnpr, t0icao
  USE mo_physical_constants,     ONLY: rd, vtmpc1, cpd, rcpd, vtmpc2
  USE mo_decomposition, ONLY: ldc => local_decomposition
  USE mo_memory_gl,     ONLY: xt
  USE mo_column,        ONLY: nfor_div,nfor_omega,                     &
                              get_col_div,  &  ! prescribe divergence
                              get_col_omega    ! prescribe omega
  USE mo_advection
  USE mo_timer,         ONLY: timer_start, timer_stop, timer_dyn

  IMPLICIT NONE

  !  Array arguments (now indexed by continuous latitude index)

  REAL(dp) ,INTENT(OUT) :: petadot(ldc% nproma ,nlevp1, ldc% ngpblks)

  !  Local array bounds
  INTEGER :: nproma, ngpblks, nbdim

  !  Local scalars: 
  REAL(dp) :: dpde, sdoth, zalpha, zalphm, zardrc, zcorio, zcpdrrd,    &
              zcpg, zcst, zdelb, zdlp, zdt, zetdpde, zln, zlnpr,       &
              zrcsth, zrddlb, zrdwrp, zrgr, zrrd, zrtgr, ztq, zvcikp
  INTEGER :: ikp, jrow, jk, jl

  !  Local arrays: 
  REAL(dp) :: zgeos  (ldc%nproma)         !
  REAL(dp) :: ztbar  (ldc%nproma, nlev)   !
  REAL(dp) :: ztv    (ldc%nproma, nlev)   !
  REAL(dp) :: aph    (ldc%nproma ,nlevp1) ! pressure
  REAL(dp) :: zdpsl  (ldc%nproma        ) ! surface pressure gradient
  REAL(dp) :: zdpsm  (ldc%nproma        ) ! surface pressure gradient
  REAL(dp) :: delp   (ldc%nproma ,nlev  ) ! pressure difference across layers
  REAL(dp) :: zrdelp (ldc%nproma ,nlev  ) ! reciprocal of *pdelp.*
  REAL(dp) :: zvgrps (ldc%nproma ,nlev  ) ! v * grad(ps)
  REAL(dp) :: zsdiv  (ldc%nproma ,nlevp1) ! surface pressure tendency

  !  External subroutines 
  EXTERNAL :: auxhyb, geopot, pres
#ifdef _AIX
  EXTERNAL :: util_mem_zero_double
#endif

  !  Executable statements 

  !-- 1. Preliminary calculations

  !-- 1.1 Set local values

  zcpdrrd = cpd/rd
  zrrd    = 1.0_dp/rd

  !  Local array bounds

  ngpblks = ldc% ngpblks ! number of rows
  nbdim   = ldc% nproma

!CSD$ PARALLEL DO PRIVATE(nproma, dpde, sdoth, zalpha, zalphm, zardrc,     &
!CSD$&  zcorio, zcpg, zcst, zdelb, zdlp, zdt, zetdpde, zln, zlnpr, zrcsth, &
!CSD$&  zrddlb, zrdwrp, zrgr, zrtgr, ztq, zvcikp, ikp, jrow, jk, jl,       &
!CSD$&  zgeos, ztbar, ztv, aph, zdpsl, zdpsm, delp, zrdelp, zvgrps, zsdiv)

!$OMP PARALLEL

  IF (ltimer) CALL timer_start(timer_dyn)

!$OMP DO PRIVATE(nproma, dpde, sdoth, zalpha, zalphm, zardrc,             &
!$OMP  zcorio, zcpg, zcst, zdelb, zdlp, zdt, zetdpde, zln, zlnpr, zrcsth, &
!$OMP  zrddlb, zrdwrp, zrgr, zrtgr, ztq, zvcikp, ikp, jrow, jk, jl,       &
!$OMP  zgeos, ztbar, ztv, aph, zdpsl, zdpsm, delp, zrdelp, zvgrps, zsdiv)
  DO jrow = 1, ngpblks

    IF ( jrow == ngpblks ) THEN
      nproma = ldc% npromz
    ELSE
      nproma = ldc% nproma
    END IF

!-- 1.2 Compute surface pressure and its gradient

    aph(1:nproma,nlevp1) = EXP(alps(1:nproma,jrow))
    zdpsl(1:nproma) = aph(1:nproma,nlevp1)*dalpsl(1:nproma,jrow)
    zdpsm(1:nproma) = aph(1:nproma,nlevp1)*dalpsm(1:nproma,jrow)

!-- 1.3 Compute half-level pressures and auxilliary variables.

    CALL pres(aph,nbdim,aph(1,nlevp1),nproma)
    CALL auxhyb(delp,zrdelp,alnpr(:,:,jrow),alpha(:,:,jrow),   &
              aph,nbdim,nproma)

!-- 1.4 Compute v.grad(ps)

    DO jk = nplvp1, nlev
      DO jl = 1, nproma
        zvgrps(jl,jk) = u(jl,jk,jrow)*zdpsl(jl)                    &
                                        +(v(jl,jk,jrow)*zdpsm(jl))
      END DO
    END DO

!-- 2. Sum divergence and compute surface pressure tendency

   !>Prescribe divergence in SCM
   IF (lcolumn .AND. nfor_div(1)>0) CALL get_col_div(d)

!-- 2.1 Compute pressure-level sums

#ifdef _AIX
    CALL util_mem_zero_double(zsdiv(1,1), nproma)
#else
    zsdiv (1:nproma,1) = 0.0_dp
#endif

    DO jk = 1, nplev
      ikp = jk + 1
      zdlp = delpr(jk)

      DO jl = 1, nproma
        zsdiv(jl,ikp) = d(jl,jk,jrow)*zdlp + zsdiv(jl,jk)
      END DO
    END DO

!-- 2.2 Compute hybrid-level sums

    DO jk = nplvp1, nlev
      ikp = jk + 1
      zdelb = delb(jk)

      DO jl = 1, nproma
        zsdiv(jl,ikp) = d(jl,jk,jrow)*delp(jl,jk) +                &
                                  zdelb*zvgrps(jl,jk) + zsdiv(jl,jk)
      END DO
    END DO

!-- prescribe vertical velocity in SCM run

    IF (lcolumn .and. nfor_omega(1)>0) CALL get_col_omega(zsdiv,d,delp)

!-- 2.3 Tendency of logarithm of surface pressure

    DO jl = 1, nproma
      alpste(jl,jrow) = -zsdiv(jl,nlevp1)/aph(jl,nlevp1)
    END DO

!-- 3. Compute reference temperature and deviation
!      of virtual temprature

    DO jl = 1, nproma
      zgeos(jl) = rd*alps(jl,jrow)
    END DO

    DO jk = nlev, 1, -1

      DO jl = 1, nproma
        ztbar(jl,jk) = t0icao*EXP(alrrdic*(zgeos(jl)-                  &
                                       alpha(jl,jk,jrow)-rdlnp0i))
        zgeos(jl) = zgeos(jl) - alnpr(jl,jk,jrow)
        ztv(jl,jk) = t(jl,jk,jrow)*                                &
          (1._dp+vtmpc1*q(jl,jk,jrow)-(xl(jl,jk,jrow)+xi(jl,jk,jrow))) &
           - ztbar(jl,jk)
      END DO
    END DO

!-- 4. Compute vertical advection
!      and do zonal mean and box diagnostics.

!-- 4.1 Compute vertical advection

#ifdef _AIX
    CALL util_mem_zero_double(vom   (1,1,jrow),nproma)
    CALL util_mem_zero_double(vol   (1,1,jrow),nproma)
    CALL util_mem_zero_double(tte   (1,1,jrow),nproma)
    CALL util_mem_zero_double(qte   (1,1,jrow),nproma)
    CALL util_mem_zero_double(xlte  (1,1,jrow),nproma)
    CALL util_mem_zero_double(xite  (1,1,jrow),nproma)
    CALL util_mem_zero_double(vervel(1,1,jrow),nproma)
#else
    vom    (1:nproma,1,jrow) = 0.0_dp
    vol    (1:nproma,1,jrow) = 0.0_dp
    tte    (1:nproma,1,jrow) = 0.0_dp
    qte    (1:nproma,1,jrow) = 0.0_dp
    xlte   (1:nproma,1,jrow) = 0.0_dp
    xite   (1:nproma,1,jrow) = 0.0_dp
    vervel (1:nproma,1,jrow) = 0.0_dp
#endif

    DO jk = 1, nlevm1
      ikp = jk + 1
      zvcikp = vct(nvclev+ikp)

      DO jl = 1, nproma
        sdoth = 0.5_dp*(zvcikp*zsdiv(jl,nlevp1)-zsdiv(jl,ikp))
        vom(jl,jk ,jrow) = (u(jl,jk,jrow)-u(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,jk)) + vom(jl,jk,jrow)
        vom(jl,ikp,jrow) = (u(jl,jk,jrow)-u(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,ikp))
        vol(jl,jk ,jrow) = (v(jl,jk,jrow)-v(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,jk)) + vol(jl,jk,jrow)
        vol(jl,ikp,jrow) = (v(jl,jk,jrow)-v(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,ikp))
        tte (jl,jk,jrow) = (t(jl,jk,jrow)-t(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,jk)) + tte(jl,jk,jrow)
        tte (jl,ikp,jrow) =(t(jl,jk,jrow)-t(jl,ikp,jrow))* &
          (sdoth*zrdelp(jl,ikp))

        ! compute eta-dot for semi Lagrangian scheme

        zetdpde = 2._dp*sdoth
        dpde = (aph(jl,jk+2)-aph(jl,jk))/(cetah(jk+2)-cetah(jk))
        IF (iadvec == semi_lagrangian .AND. .NOT. lcolumn)         &
          petadot(jl,jk+1,jrow) = zetdpde/dpde
      END DO
    END DO

    IF (lcolumn) THEN
      DO jk = 1, nlevm1
        ikp = jk + 1
        zvcikp = vct(nvclev+ikp)

        DO jl = 1, nproma
          sdoth = 0.5_dp*(zvcikp*zsdiv(jl,nlevp1)-zsdiv(jl,ikp))

          ! vertical transport in single column model
          ! this is a poor man's replacement for the transport

            qte(jl,jk ,jrow)     = (q(jl,jk,jrow)-q(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,jk )) + qte(jl,jk,jrow)
            qte(jl,ikp,jrow)     = (q(jl,jk,jrow)-q(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,ikp))
            xlte(jl,jk ,jrow)  = (xl(jl,jk,jrow)-xl(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,jk )) + xlte(jl,jk,jrow)
            xlte(jl,ikp,jrow)  = (xl(jl,jk,jrow)-xl(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,ikp))
            xite(jl,jk ,jrow)  = (xi(jl,jk,jrow)-xi(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,jk )) + xite(jl,jk,jrow)
            xite(jl,ikp,jrow)  = (xi(jl,jk,jrow)-xi(jl,ikp,jrow))* &
              (sdoth*zrdelp(jl,ikp))
            xtte(jl,jk ,:,jrow)=                                   &
                                 (xt(jl,jk,:,jrow)-xt(jl,ikp,:,jrow))* &
                       (sdoth*zrdelp(jl,jk )) + xtte(jl,jk,:,jrow)
            xtte(jl,ikp,:,jrow)=                                   &
                                 (xt(jl,jk,:,jrow)-xt(jl,ikp,:,jrow))* &
                        (sdoth*zrdelp(jl,ikp))
        END DO
      END DO
    ENDIF
 
#ifdef _AIX
    CALL util_mem_zero_double(petadot(1,1     ,jrow),nproma)
    CALL util_mem_zero_double(petadot(1,nlevp1,jrow),nproma)
#else
    petadot(1:nproma,     1,jrow) = 0.0_dp
    petadot(1:nproma,nlevp1,jrow) = 0.0_dp
#endif

!-- 5. Compute energy conversion term for pressure levels
!      and do zonal mean and box diagnostics.

!-- 5.1 Compute energy conversion term for pressure levels

    DO jk = 1, nplev
      zardrc = ardprc(jk)
      zalphm = alpham(jk)

!OCL NOALIAS
      DO jl = 1, nproma
        ztq = (ztbar(jl,jk)+ztv(jl,jk))/(1._dp+vtmpc2*q(jl,jk,jrow))
        zdt = -ztq*(zsdiv(jl,jk)*zardrc+zalphm*d(jl,jk,jrow))
        tte(jl,jk,jrow) = tte(jl,jk,jrow) + zdt
      END DO
    END DO

!-- 6. Compute pressure-gradient terms,complete calculation
!      of energy conversion term

!-- 6.1 Hybrid levels

    DO jk = nplvp1, nlmsgl
      zdelb = delb(jk)
      zrddlb = rddelb(jk)
      zcpg = cpg(jk)

!OCL NOVREC,NOALIAS
      DO jl = 1, nproma
        zrgr = (zrddlb+zcpg*alnpr(jl,jk,jrow)*zrdelp(jl,jk))*      &
                                                          zrdelp(jl,jk)
        zrtgr = zrgr*ztv(jl,jk)
        zcst = cst_2d(jl,jrow)
        vom(jl,jk,jrow) = vom(jl,jk,jrow) - zcst*zrtgr*zdpsl(jl)
        vol(jl,jk,jrow) = vol(jl,jk,jrow) - zcst*zrtgr*zdpsm(jl)
        zrdwrp = (zrgr*zvgrps(jl,jk)-(zrdelp(jl,jk)*(zsdiv(jl,jk)*     &
                 alnpr(jl,jk,jrow)+alpha(jl,jk,jrow)*zdelb*    &
                 zvgrps(jl,jk))+                                       &
                 alpha(jl,jk,jrow)*d(jl,jk,jrow)))*rcpd
        ztq = (ztbar(jl,jk)+ztv(jl,jk))/(1._dp+vtmpc2*q(jl,jk,jrow))
        zdt = ztq*zrdwrp
        tte(jl,jk,jrow) = tte(jl,jk,jrow) + zdt

      END DO
    END DO

!-- 6.2 Sigma levels

    DO jk = nlmslp, nlev
      zalphm = alpham(jk)
      zlnpr = rlnpr(jk)
      zalpha = ralpha(jk)

      zrddlb = rddelb(jk)
!OCL NOALIAS
      DO jl = 1, nproma
        zrgr = zrddlb*zrdelp(jl,jk)
        zrtgr = zrgr*ztv(jl,jk)
        zcst = cst_2d(jl,jrow)
        vom(jl,jk,jrow) = vom(jl,jk,jrow) - zcst*zrtgr*zdpsl(jl)
        vol(jl,jk,jrow) = vol(jl,jk,jrow) - zcst*zrtgr*zdpsm(jl)
        zrdwrp = (zrgr*zvgrps(jl,jk)-(zrdelp(jl,jk)*(zsdiv(jl,jk)*     &
                 zlnpr+zalphm*                                         &
                 zvgrps(jl,jk))+zalpha*d(jl,jk,jrow)))*rcpd
        ztq = (ztbar(jl,jk)+ztv(jl,jk))/(1._dp+vtmpc2*q(jl,jk,jrow))
        zdt = ztq*zrdwrp
        tte(jl,jk,jrow) = tte(jl,jk,jrow) + zdt

      END DO
    END DO
 
!-- 6.3 Compute vertical velocity for mass-flux scheme

    DO jk = 1, nplev
      zardrc = ardprc(jk)
      zalphm = alpham(jk)
      DO jl = 1, nproma
        vervel(jl,jk,jrow) = &
          -(zsdiv(jl,jk)*zardrc+zalphm*d(jl,jk,jrow))*zcpdrrd
      END DO
    END DO

    DO jk = nplvp1, nlmsgl
      zdelb = delb(jk)
      zrddlb = rddelb(jk)
      zcpg = cpg(jk)
      DO jl = 1, nproma
        zrgr = (zrddlb+zcpg*alnpr(jl,jk,jrow)*zrdelp(jl,jk))*      &
                                                         zrdelp(jl,jk)
        vervel(jl,jk,jrow) = (zrgr*zvgrps(jl,jk)-zrdelp(jl,jk)*    &
                (zsdiv(jl,jk)*alnpr(jl,jk,jrow)+                   &
                alpha(jl,jk,jrow)*zdelb*zvgrps(jl,jk))-            &
                alpha(jl,jk,jrow)*d(jl,jk,jrow))*zrrd
      END DO
    END DO

    DO jk = nlmslp, nlev
      zalphm = alpham(jk)
      zlnpr = rlnpr(jk)
      zalpha = ralpha(jk)
      zrddlb = rddelb(jk)
      DO jl = 1, nproma
        zrgr = zrddlb*zrdelp(jl,jk)
        vervel(jl,jk,jrow) = (zrgr*zvgrps(jl,jk)-zrdelp(jl,jk)*    &
                (zsdiv(jl,jk)*zlnpr+zalphm*zvgrps(jl,jk))-             &
                 zalpha*d(jl,jk,jrow))*zrrd
      END DO
    END DO

    DO jk = 1, nlev
      DO jl = 1, nproma
        vervel(jl,jk,jrow) = vervel(jl,jk,jrow)                &
                                    * 0.5_dp*(aph(jl,jk)+aph(jl,jk+1))
      END DO
    END DO

!-- 7. Compute geopotential

!-- 7.1 Compute deviation of geopotential height at surface

    DO jl = 1, nproma
      zln = rd*alps(jl,jrow) - rdlnp0i
      zgeos(jl) = geospm(jl,jrow) + rdt0ral*EXP(alrrdic*zln)
    END DO

!-- 7.2 Compute deviation of geopotential height

    CALL geopot(rh(:,:,jrow),ztv,alnpr(:,:,jrow),              &
                 alpha(:,:,jrow), zgeos, nbdim, nproma)

!-- 8. Compute horizontal advection terms

    DO jk = 1, nlev
!OCL NOVREC,NOALIAS
      DO jl = 1, nproma
        zcorio = coriol_2d(jl,jrow)
        zrcsth = rcsth_2d(jl,jrow)
        rh(jl,jk,jrow)=zrcsth*(u(jl,jk,jrow)*u(jl,jk,jrow)+&
          v(jl,jk,jrow)*v(jl,jk,jrow)) + rh(jl,jk,jrow)
        IF(.NOT.lcolumn) THEN
          vom(jl,jk,jrow)=(vo(jl,jk,jrow)+zcorio)*             &
            v(jl,jk,jrow)+vom(jl,jk,jrow)
          vol(jl,jk,jrow) = -(vo(jl,jk,jrow)+zcorio)*          &
            u(jl,jk,jrow)+vol(jl,jk,jrow)
        END IF
        zdt = -u(jl,jk,jrow)*dtl(jl,jk,jrow)                   &
              - v(jl,jk,jrow)*dtm(jl,jk,jrow)
        tte (jl,jk,jrow) = tte(jl,jk,jrow) + zdt

      END DO
    END DO
 
!-- 10. Duplicate ps

    aps(1:nproma,jrow) = aph(1:nproma,nlevp1)
    apsm(1:nproma,jrow) = aph(1:nproma,nlevp1)

  END DO
!$OMP END DO
!CSD$ END PARALLEL DO

  IF (ltimer) CALL timer_stop(timer_dyn)

!$OMP END PARALLEL
 
END SUBROUTINE dyn
