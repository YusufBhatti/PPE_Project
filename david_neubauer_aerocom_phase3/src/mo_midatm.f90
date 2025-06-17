#ifdef __xlC__
@PROCESS HOT
#else
#define FSEL(a,b,c) MERGE(b,c,(a) >= 0._wp)
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_midatm

  USE mo_kind, ONLY: wp
#ifdef _PROFILE
  USE mo_profile,        ONLY : trace_start, trace_stop
#endif

  IMPLICIT NONE

CONTAINS

  SUBROUTINE gwspectrum ( krow, kproma, kbdim, klev  &
       ,paphm1,      papm1                       &
       ,ptm1,        pum1,        pvm1           &
       ,paprflux                                 &
       ,ptte,        pvol,        pvom   )


    !
    ! Description:
    ! 
    !   Hines parameterization from ccc/mam (Hines, 1997a,b):       
    !   physical tendencies of the prognostic variables u,v 
    !   due to vertical transports by a broad band spectrum
    !   of gravity waves. 
    !
    !   Note that diffusion coefficient and  heating rate 
    !   only calculated if iheatcal = 1.
    !
    !        *gwspectrum* is called from *physc*.
    !
    !  Authors:
    !    
    !   c. mclandress   ists   august 1995
    !   n. mcfarlane    cccma  may 1995
    !   m. charron      mpi    2000-2001
    !   e. manzini      mpi    february 2002 (re-write, based on cccgwd)                 
    !   h. schmidt      mpi    march 2003 
    !   h. schmidt      mpi    april 2010 (latitude dependent source function
    !                                      added)
    !
    !
    USE mo_physical_constants,ONLY: grav, rd, cpd, rhoh2o
    USE mo_gwspectrum,        ONLY: lextro, lfront, lozpr, iheatcal, &
                                    emiss_lev, rmscon, kstar, naz,  & 
                                    pcrit, pcons, lrmscon_lat, rmscon_lat
    USE mo_geoloc,            ONLY: ilat

    IMPLICIT NONE

    ! scalar argument with intent(IN)
    INTEGER ,INTENT(in) ::  krow, kproma, kbdim, klev

    !  Array arguments with intent(IN):
    ! Input 1D
    REAL(wp) ,INTENT(IN) :: paprflux(kbdim)     ! precipitation flux
    ! Input 2D
    REAL(wp) ,INTENT(IN) :: paphm1(kbdim,klev+1) ! half level pressure (t-dt)
    REAL(wp) ,INTENT(IN) :: papm1(kbdim,klev)    ! full level pressure (t-dt)
    REAL(wp) ,INTENT(IN) :: ptm1(kbdim,klev)     ! temperature (t-dt)
    REAL(wp) ,INTENT(IN) :: pum1(kbdim,klev)     ! zonal wind (t-dt)
    REAL(wp) ,INTENT(IN) :: pvm1(kbdim,klev)     ! meridional wind (t-dt)

    !  Array arguments with intent(InOut):
    ! - input/output 2d
    REAL(wp) ,INTENT(INOUT) :: ptte(kbdim,klev)  ! tendency of temperature
    REAL(wp) ,INTENT(INOUT) :: pvol(kbdim,klev)  ! tendency of meridional wind
    REAL(wp) ,INTENT(INOUT) :: pvom(kbdim,klev)  ! tendency of zonal wind

    !  Local arrays for ccc/mam hines gwd scheme:

    ! Important local parameter (passed to all subroutines):
    INTEGER, PARAMETER :: nazmth = 8  ! max azimuth array dimension size 

    REAL(wp) :: pressg(kproma)  ! Surface pressure (pascal)  
    REAL(wp) :: zpr(kproma)     ! precipitation (check dims: echam5 change)

    ! * Vertical positioning arrays and work arrays:                       
    REAL(wp) :: sgj(kproma,klev), shj(kproma,klev), shxkj(kproma,klev)
    REAL(wp) ::  dsgj(kproma,klev), dttdsf(kproma)

    REAL(wp) :: utendgw(kproma,klev) ! zonal tend, gravity wave spectrum (m/s^2)
    REAL(wp) :: vtendgw(kproma,klev) ! merid tend, gravity wave spectrum (m/s^2)
    REAL(wp) :: ttendgw(kproma,klev) ! temperature tend, gravity wave spectrum (K/s)
    REAL(wp) :: diffco(kproma,klev)  ! diffusion coefficient (m^2/s) 

    REAL(wp) ::  flux_u(kproma,klev) ! zonal momentum flux (pascals)
    REAL(wp) ::  flux_v(kproma,klev) ! meridional momentum flux (pascals) 

    REAL(wp) :: uhs(kproma,klev)      ! zonal wind (m/s), input for hines param
    REAL(wp) :: vhs(kproma,klev)      ! merid wind (m/s), input for hines param
    REAL(wp) :: bvfreq(kproma,klev)   ! background brunt vassala frequency (rad/s)
    REAL(wp) :: density(kproma,klev)  ! background density (kg/m^3)
    REAL(wp) :: visc_mol(kproma,klev) ! molecular viscosity (m^2/s) 
    REAL(wp) :: alt(kproma,klev)      ! background altitude (m)

    REAL(wp) :: rmswind(kproma)        ! rms gravity wave  wind, lowest level (m/s) 
    REAL(wp) :: anis(kproma,nazmth)    ! anisotropy factor (sum over azimuths = 1) 
    REAL(wp) :: k_alpha(kproma,nazmth) ! horizontal wavenumber of each azimuth (1/m)
    LOGICAL :: lorms(kproma)       ! .true. for rmswind /=0 at launching level 

    REAL(wp) :: m_alpha(kproma,klev,nazmth) ! cutoff vertical wavenumber (1/m)
    REAL(wp) :: mmin_alpha(kproma,nazmth)   ! minumum value of m_alpha
    REAL(wp) :: sigma_t(kproma,klev)        ! total rms gw wind (m/s)

    ! gw variances from orographic sources (for coupling to a orogwd)
    REAL(wp) :: sigsqmcw(kproma,klev,nazmth), sigmatm(kproma,klev)

    !
    ! Local scalars:
    INTEGER  :: jk, jl
    INTEGER  :: levbot     ! gravity wave spectrum lowest level
    REAL(wp)     :: rgocp, zpcons, hscal, ratio, paphm1_inv

#ifdef _PROFILE
  CALL trace_start ('gwspectrum', 20)
#endif

    IF (lextro) THEN

       !
       !--  Initialize the ccc/mam hines gwd scheme
       !

       utendgw(:,:) = 0.0_wp
       vtendgw(:,:) = 0.0_wp
       ttendgw(:,:) = 0.0_wp

       diffco(:,:) = 0.0_wp

       flux_u(:,:) = 0.0_wp
       flux_v(:,:) = 0.0_wp 

       uhs(:,:) = 0.0_wp
       vhs(:,:) = 0.0_wp

       ! Wind variances form orographic gravity waves
       ! Note: the code is NOT fully implemeted for this case!

       sigsqmcw(:,:,:) = 0.0_wp
       sigmatm(:,:)    = 0.0_wp

       ! precipitation (check the units!):
       zpcons = (1000.0_wp*86400.0_wp)/rhoh2o
       zpr(1:kproma)=zpcons*paprflux(1:kproma)

       rgocp=rd/cpd

       ! Vertical positioning arrays: 
       DO jk=1,klev
!IBM* novector
          DO jl=1,kproma
             paphm1_inv = 1._wp/paphm1(jl,klev+1)
             shj(jl,jk)=papm1(jl,jk)*paphm1_inv
             sgj(jl,jk)=papm1(jl,jk)*paphm1_inv
             dsgj(jl,jk)=(paphm1(jl,jk+1)-paphm1(jl,jk))*paphm1_inv
          END DO
          shxkj(1:kproma,jk) = sgj(1:kproma,jk)**rgocp
       END DO

       ! Surface pressure: 
       DO jl=1,kproma
          pressg(jl)=paphm1(jl,klev+1)
       END DO

       !
       !     * calculate b v frequency at all points
       !     * and smooth bvfreq.

       DO jk=2,klev
!IBM* novector
         DO jl=1,kproma
           dttdsf(jl)=(ptm1(jl,jk)/shxkj(jl,jk)-ptm1(jl,jk-1)/shxkj(jl,jk-1)) &
                     /(shj(jl,jk)-shj(jl,jk-1))
           dttdsf(jl)=MIN(dttdsf(jl), -5.0_wp/sgj(jl,jk))
           dttdsf(jl)=dttdsf(jl)*shxkj(jl,jk)

           bvfreq(jl,jk)=-dttdsf(jl)*sgj(jl,jk)/rd
         END DO
         bvfreq(1:kproma,jk) = SQRT(bvfreq(1:kproma,jk))

!IBM* novector
         DO jl=1,kproma
            bvfreq(jl,jk) = bvfreq(jl,jk)*grav/ptm1(jl,jk)
         END DO
       END DO

       bvfreq(:,1) = bvfreq(:,2)

       DO jk=2,klev
          DO jl=1,kproma
             ratio=5.0_wp*LOG(sgj(jl,jk)/sgj(jl,jk-1))
             bvfreq(jl,jk) = (bvfreq(jl,jk-1) + ratio*bvfreq(jl,jk))/(1.0_wp+ratio)
          END DO
       END DO

       !     * altitude and density at bottom.

       alt(:,klev) = 0.0_wp
!IBM* novector
       DO jl=1,kproma
          hscal = rd * ptm1(jl,klev) / grav
          density(jl,klev) = sgj(jl,klev) * pressg(jl) / (grav*hscal)
       END DO

       !     * altitude and density at remaining levels.

       DO jk=klev-1,1,-1
!IBM* novector
          DO jl=1,kproma
             hscal = rd * ptm1(jl,jk) / grav
             alt(jl,jk) = alt(jl,jk+1) + hscal * dsgj(jl,jk) / sgj(jl,jk)
             density(jl,jk) = sgj(jl,jk) * pressg(jl) / (grav*hscal)
          END DO
       END DO

       !
       !     * set molecular viscosity to a very small value.
       !     * if the model top is greater than 100 km then the actual
       !     * viscosity coefficient could be specified here.


       visc_mol(:,:) = 1.5e-5_wp


       ! use single value for azimuthal-dependent horizontal wavenumber:
       ! kstar = (old latitudinal dependence, introduce here if necessary)

       k_alpha(:,:) = kstar

       !     * defile bottom launch level (emission level of gws)

       levbot = klev-emiss_lev   

       !     * initialize switch for column calculation

       lorms(:) = .FALSE.

       !     * background wind minus value at bottom launch level.

       DO jk=1,levbot
         DO jl=1,kproma 
           uhs(jl,jk) = pum1(jl,jk) - pum1(jl,levbot)
           vhs(jl,jk) = pvm1(jl,jk) - pvm1(jl,levbot)
         END DO
       END DO

       !     * specify root mean square wind at bottom launch level.

       DO jl=1,kproma 
          anis(jl,:)   = 1.0_wp/REAL(naz,wp)
       END DO

       IF (lrmscon_lat) THEN
          !     * latitude dependent GW source
          DO jl=1,kproma
            rmswind(jl)=rmscon_lat(ilat(jl,krow))
          END DO
       ELSE
          DO jl=1,kproma
            rmswind(jl) = rmscon
          END DO
       ENDIF

       !     * gravity waves from fronts:
       IF (lfront) THEN
          CALL gw_fronts(krow, kproma, nazmth, rmswind, anis)
       ENDIF

       !     * modulation by precipitation:
       IF (lozpr) THEN
          DO jl=1,kproma 
             IF (zpr(jl) .GT. pcrit) THEN
                rmswind(jl) = rmscon + ( (zpr(jl)-pcrit)/zpr(jl) )*pcons
             ENDIF
          END DO
       ENDIF

       DO jl=1,kproma 
          IF (rmswind(jl) .GT. 0.0_wp) THEN
             lorms(jl) = .TRUE.
          ENDIF
       END DO


       !
       !     * calculate gw tendencies (note that diffusion coefficient and
       !     * heating rate only calculated if iheatcal = 1).
       !
       CALL hines_extro ( kproma, klev, nazmth,                          & 
                          utendgw, vtendgw, ttendgw, diffco,             &
                          flux_u, flux_v,                                & 
                          uhs, vhs, bvfreq, density, visc_mol, alt,      & 
                          rmswind, anis, k_alpha, sigsqmcw,              &
                          m_alpha,  mmin_alpha ,sigma_t, sigmatm,        & 
                          levbot, lorms)


       !   update tendencies: 
       !
       DO jk=1, klev
          DO jl=1,kproma
             pvom(jl,jk) = pvom(jl,jk) + utendgw(jl,jk)
             pvol(jl,jk) = pvol(jl,jk) + vtendgw(jl,jk)
          END DO
       END DO
       !
       IF ( iheatcal == 1 ) THEN  
          DO jk=1, klev
             DO jl=1,kproma   
                ptte(jl,jk) = ptte(jl,jk) + ttendgw(jl,jk)
             END DO
          END DO
       ENDIF
       !

    ENDIF
    !     * end of hines calculations.
#ifdef _PROFILE
     CALL trace_stop ('gwspectrum', 20)
#endif

    !-----------------------------------------------------------------------
  END SUBROUTINE gwspectrum

  SUBROUTINE hines_extro ( nlons, nlevs, nazmth,                          &
                           drag_u, drag_v, heat, diffco, flux_u, flux_v,  &
                           vel_u, vel_v, bvfreq, density, visc_mol, alt,  &
                           rmswind, anis, k_alpha, sigsqmcw,               &
                           m_alpha,  mmin_alpha, sigma_t, sigmatm,        &
                           lev2,  lorms)
    !
    !  main routine for hines' "extrowave" gravity wave parameterization based
    !  on hines' doppler spread theory. this routine calculates zonal
    !  and meridional components of gravity wave drag, heating rates
    !  and diffusion coefficient on a longitude by altitude grid.
    !  no "mythical" lower boundary region calculation is made. 
    !
    !  aug. 13/95 - c. mclandress
    !  sept. /95  - n. mcfarlane
    !  1995- 2002 - e. manzini
    !
    !  modifications:
    !
    !  output arguements:
    !
    !     * drag_u = zonal component of gravity wave drag (m/s^2).
    !     * drag_v = meridional component of gravity wave drag (m/s^2).
    !     * heat   = gravity wave heating (k/sec).
    !     * diffco = diffusion coefficient (m^2/sec)
    !     * flux_u = zonal component of vertical momentum flux (pascals)
    !     * flux_v = meridional component of vertical momentum flux (pascals)
    !
    !  input arguements:
    !
    !     * vel_u      = background zonal wind component (m/s).
    !     * vel_v      = background meridional wind component (m/s).
    !     * bvfreq     = background brunt vassala frequency (radians/sec).
    !     * density    = background density (kg/m^3) 
    !     * visc_mol   = molecular viscosity (m^2/s)
    !     * alt        = altitude of momentum, density, buoyancy levels (m)
    !     *              (note: levels ordered so that alt(i,1) > alt(i,2), etc.)
    !     * rmswind   = root mean square gravity wave wind at lowest level (m/s).
    !     * anis      = anisotropy factor (sum over azimuths is one)
    !     * lorms     = .true. for drag computation (column selector)
    !     * k_alpha    = horizontal wavenumber of each azimuth (1/m).
    !     * lev2       = index of last level (eg bottom) for drag calculation 
    !     *              (i.e., lev1 < lev2 <= nlevs).
    !     * nlons      = number of longitudes.
    !     * nlevs      = number of vertical levels.
    !     * nazmth     = azimuthal array dimension (nazmth >= naz).
    ! 
    !  ouput diagnostics:
    !     * m_alpha      = cutoff vertical wavenumber (1/m).
    !     * mmin_alpha   = minimum value of cutoff wavenumber.
    !     * sigma_t      = total rms horizontal wind (m/s).

    !  work arrays.
    !

    !     * v_alpha      = wind component at each azimuth (m/s) and if iheatcal=1
    !     *                holds vertical derivative of cutoff wavenumber.
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !     * ak_alpha     = spectral amplitude factor at each azimuth 
    !     *                (i.e.,{ajkj}) in m^4/s^2.
    !     * densb        = background density at bottom level.
    !     * bvfb         = buoyancy frequency at bottom level and
    !     *                work array for icutoff = 1.
    !
    !     * losigma_t    =  .true. for total sigma not zero

    USE mo_gwspectrum,     ONLY: iheatcal,                        &
                                 kstar, m_min,                    &
                                 naz, slope, f1, f2, f3, f5, f6,  &
                                 icutoff, alt_cutoff, smco, nsmax

    IMPLICIT NONE

    INTEGER :: nlons, nlevs, nazmth, lev2

    REAL(wp):: drag_u(nlons,nlevs),   drag_v(nlons,nlevs) 
    REAL(wp):: heat(nlons,nlevs),     diffco(nlons,nlevs)
    REAL(wp):: flux_u(nlons,nlevs),   flux_v(nlons,nlevs)
    REAL(wp):: flux(nlons,nlevs,nazmth)
    REAL(wp):: vel_u(nlons,nlevs),    vel_v(nlons,nlevs)
    REAL(wp):: bvfreq(nlons,nlevs),   density(nlons,nlevs)
    REAL(wp):: visc_mol(nlons,nlevs), alt(nlons,nlevs)
    REAL(wp):: rmswind(nlons),      bvfb(nlons),   densb(nlons)
    REAL(wp):: anis(nlons,nazmth)
    REAL(wp):: sigma_t(nlons,nlevs), sigsqmcw(nlons,nlevs,nazmth)
    REAL(wp):: sigma_alpha(nlons,nlevs,nazmth), sigmatm(nlons,nlevs)

    REAL(wp):: m_alpha(nlons,nlevs,nazmth), v_alpha(nlons,nlevs,nazmth)
    REAL(wp):: ak_alpha(nlons,nazmth),      k_alpha(nlons,nazmth)
    REAL(wp):: mmin_alpha(nlons,nazmth)    
    REAL(wp):: smoothr1(nlons,nlevs), smoothr2(nlons,nlevs)

    LOGICAL :: lorms(nlons), losigma_t(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER :: i, n, l, lev1, il1, il2, iprint

#ifdef _PROFILE
  CALL trace_start ('hines_extro', 21)
#endif
    !----------------------------------------------------------------------- 
    !

    ! range of longitude index:
    il1 = 1      
    il2 = nlons     

    lev1=1              ! top level index

    iprint = 0       !     * iprint     = 1 to print out various arrays.

    !
    !  buoyancy and density at bottom level.
    !
    DO i = il1,il2
       bvfb(i)  = bvfreq(i,lev2)
       densb(i) = density(i,lev2)
    END DO
    !
    !  initialize some variables
    !
    DO n = 1,naz
       DO l=lev1,lev2
          DO i=il1,il2
             m_alpha(i,l,n) =  m_min
          END DO
       END DO
    END DO
    !
    !  compute azimuthal wind components from zonal and meridional winds.
    !
    CALL hines_wind ( v_alpha,   & 
                      vel_u, vel_v, naz,   &
                      il1, il2, lev1, lev2, nlons, nlevs, nazmth )
    !
    !  calculate cutoff vertical wavenumber and velocity variances.
    !
    CALL hines_wavnum ( m_alpha, sigma_t, sigma_alpha, ak_alpha,   &
                        mmin_alpha, losigma_t,                     &
                        v_alpha, visc_mol, density, densb,         &
                        bvfreq, bvfb, rmswind, anis, lorms,        &
                        sigsqmcw, sigmatm,                         &
                        il1, il2, lev1, lev2, nlons, nlevs, nazmth)
    !
    !  smooth cutoff wavenumbers and total rms velocity in the vertical 
    !  direction nsmax times, using flux_u as temporary work array.
    !   
    IF (nsmax.GT.0)  THEN
       DO n = 1,naz
          DO l=lev1,lev2
             DO i=il1,il2
                smoothr1(i,l) = m_alpha(i,l,n)
             END DO
          END DO
          CALL vert_smooth (smoothr1,smoothr2, smco, nsmax,  &
                            il1, il2, lev1, lev2, nlons, nlevs )
          DO l=lev1,lev2
             DO i=il1,il2
                m_alpha(i,l,n) = smoothr1(i,l)
             END DO
          END DO
       END DO
       CALL vert_smooth ( sigma_t, smoothr2, smco, nsmax, &
                          il1, il2, lev1, lev2, nlons, nlevs )
    END IF
    !
    !  calculate zonal and meridional components of the
    !  momentum flux and drag.
    !
    CALL hines_flux ( flux_u, flux_v, flux, drag_u, drag_v,       &
                      alt, density, densb,                        &
                      m_alpha,  ak_alpha, k_alpha,                &
                      m_min, slope, naz,                          &
                      il1, il2, lev1, lev2, nlons, nlevs, nazmth, &
                      lorms )
    !
    !  cutoff drag above alt_cutoff, using bvfb as temporary work array.
    !
    IF (icutoff.EQ.1)  THEN     
       CALL hines_exp ( drag_u, bvfb, alt, alt_cutoff,  &
                        il1, il2, lev1, lev2, nlons, nlevs )
       CALL hines_exp ( drag_v, bvfb, alt, alt_cutoff,  &
                        il1, il2, lev1, lev2, nlons, nlevs )
    END IF
    !
    !  print out various arrays for diagnostic purposes.
    !
    IF (iprint.EQ.1)  THEN
       CALL hines_print ( flux_u, flux_v, drag_u, drag_v, alt,     &
                          sigma_t, sigma_alpha, v_alpha, m_alpha,  &
                          1, 1, il1, il2, lev1, lev2,              &
                          naz, nlons, nlevs, nazmth)
    END IF
    !
    !  if not calculating heating rate and diffusion coefficient then finished.
    !
    !  heating rate and diffusion coefficient.

    IF (iheatcal.EQ.1) THEN
       !
       !
       CALL hines_heat ( heat, diffco,                                 &
            alt, bvfreq, density, sigma_t, sigma_alpha,   &
            flux, visc_mol, kstar, f1, f2, f3, f5, f6,    &
            naz, il1, il2, lev1, lev2, nlons, nlevs,      &
            nazmth, losigma_t )
    END IF

#ifdef _PROFILE
     CALL trace_stop ('hines_extro', 21)
#endif

    !  finished.
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_extro

  SUBROUTINE hines_wavnum ( m_alpha, sigma_t, sigma_alpha, ak_alpha,     &
                            mmin_alpha, losigma_t,                       &
                            v_alpha, visc_mol, density, densb,           &
                            bvfreq, bvfb, rms_wind, anis, lorms,         &
                            sigsqmcw, sigmatm,                           &
                            il1, il2, levtop, levbot, nlons, nlevs, nazmth)
    !
    !  this routine calculates the cutoff vertical wavenumber and velocity
    !  variances on a longitude by altitude grid for the hines' doppler 
    !  spread gravity wave drag parameterization scheme.
    !  note: (1) only values of four or eight can be used for # azimuths (naz).
    !        (2) only values of 1.0, 1.5 or 2.0 can be used for slope (slope). 
    !        (3) if m_min not zero, only slope=1. can be used. 
    !
    !  aug. 10/95 - c. mclandress
    !  2000-2001  - m. charron
    !  2002       - e. manzini
    !
    !  output arguements:
    !
    !     * m_alpha      = cutoff wavenumber at each azimuth (1/m).
    !     * sigma_t      = total rms horizontal wind (m/s).
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !     * ak_alpha     = spectral amplitude factor at each azimuth 
    !     *                (i.e.,{ajkj}) in m^4/s^2.
    !     * losigma_t    =  .true. for total sigma not zero
    !     * mmin_alpha   = minimum value of cutoff wavenumber.


    !  input arguements:
    !
    !     * v_alpha  = wind component at each azimuth (m/s). 
    !     * visc_mol = molecular viscosity (m^2/s)
    !     * density  = background density (kg/m^3).
    !     * densb    = background density at model bottom (kg/m^3).
    !     * bvfreq   = background brunt vassala frequency (radians/sec).
    !     * bvfb     = background brunt vassala frequency at model bottom.
    !     * rms_wind = root mean square gravity wave wind at lowest level (m/s).
    !     * anis      = anisotropy factor (sum over azimuths is one)
    !     * lorms       = .true. for drag computation at lowest level 

    !     * levbot   = index of lowest vertical level.
    !     * levtop   = index of highest vertical level 
    !     *            (note: if levtop < levbot then level index 
    !     *             increases from top down).
    !     * il1      = first longitudinal index to use (il1 >= 1).
    !     * il2      = last longitudinal index to use (il1 <= il2 <= nlons).

    !     * nlons    = number of longitudes.
    !     * nlevs    = number of vertical levels.
    !     * nazmth   = azimuthal array dimension (nazmth >= naz).
    !

    !
    !  work arrays:
    !
    !     * i_alpha    = hines' integral at a single level.
    !
    !     * do_alpha = .true. for the azimuths and longitudes for
    !                   which to continue to compute the drag above  
    !                   the lowest level


    USE mo_io_units,       ONLY: nerr
    USE mo_exception,    ONLY: finish 
    USE mo_gwspectrum,   ONLY: kstar, m_min, slope, f1, f2, f3, naz

    IMPLICIT NONE

    INTEGER :: il1, il2, levtop, levbot, nlons, nlevs, nazmth
    REAL(wp):: m_alpha(nlons,nlevs,nazmth)
    REAL(wp):: sigma_alpha(nlons,nlevs,nazmth)
    REAL(wp):: sigalpmc(nlons,nlevs,nazmth)
    REAL(wp):: sigsqh_alpha(nlons,nlevs,nazmth)
    REAL(wp):: sigsqmcw(nlons,nlevs,nazmth)
    REAL(wp):: sigma_t(nlons,nlevs)
    REAL(wp):: sigmatm(nlons,nlevs)
    REAL(wp):: ak_alpha(nlons,nazmth)
    REAL(wp):: v_alpha(nlons,nlevs,nazmth)
    REAL(wp):: visc_mol(nlons,nlevs)
    REAL(wp):: f2mod(nlons,nlevs)
    REAL(wp):: density(nlons,nlevs),  densb(nlons)
    REAL(wp):: bvfreq(nlons,nlevs),   bvfb(nlons),  rms_wind(nlons)
    REAL(wp):: anis(nlons,nazmth) 
    REAL(wp):: i_alpha(nlons,nazmth), mmin_alpha(nlons,nazmth)

    LOGICAL :: lorms(nlons), losigma_t(nlons,nlevs), do_alpha(nlons,nazmth)
    !
    ! internal variables.
    !
    INTEGER :: ilorms(nlons),ialpha(nlons)
    INTEGER :: i, j, l, n, istart, lend, lincr, lbelow, nlorms, nalpha

    REAL(wp):: m_sub_m_turb, m_sub_m_mol, m_trial, mmsq
    REAL(wp):: visc, visc_min, sp1, f2mfac

    REAL(wp):: n_over_m(nlons), sigfac(nlons), vtmp1(nlons), vtmp2(nlons)
    !-----------------------------------------------------------------------     
    !
#ifdef _PROFILE
     CALL trace_start ('hines_wavnum', 22)
#endif

    visc_min = 1.e-10_wp

    sp1 = slope + 1.0_wp
    mmsq = m_min**2

    !
    !  indices of levels to process.
    !
    IF (levbot > levtop)  THEN
       istart = levbot - 1     
       lend   = levtop         
       lincr  = -1
    ELSE
       WRITE (nerr,*) ' Error: level index not increasing downward '
       CALL finish('hines_wavnum','Run terminated')
    END IF


    !   initialize logical flags and arrays
    DO l=1,nlevs
       losigma_t(:,l) = lorms(:)
    ENDDO
    DO n=1,nazmth
       do_alpha(:,n) = lorms(:)
    ENDDO

    nlorms = 0
    DO i = il1,il2
       IF (lorms(i)) THEN
          nlorms = nlorms + 1
          ilorms(nlorms) = i
       END IF
    END DO

    sigsqh_alpha(:,:,:) = 0.0_wp ! mpuetz: is this really necessary -> 15% of wavnum
    i_alpha(:,:) = 0.0_wp

    !
    ! calculate azimuthal variances at bottom level using anisotropy factor
    !
    DO n = 1,naz
       DO i = il1,il2
          sigsqh_alpha(i,levbot,n) = anis(i,n)* rms_wind(i)**2
       END DO
    END DO
    !
    !  velocity variances at bottom level.
    !
    CALL hines_sigma ( sigma_t, sigma_alpha,     &
                       sigsqh_alpha, naz, levbot,     &
                       il1, il2, nlons, nlevs, nazmth)

    CALL hines_sigma ( sigmatm, sigalpmc,     &
                       sigsqmcw, naz, levbot,     &
                       il1, il2, nlons, nlevs, nazmth)
    !
    !  calculate cutoff wavenumber and spectral amplitude factor 
    !  at bottom level where it is assumed that background winds vanish
    !  and also initialize minimum value of cutoff wavnumber.
    !
    IF ( ABS(slope-1.0_wp) < EPSILON(1.0_wp) ) THEN
       DO n = 1,naz
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
                m_alpha(i,levbot,n) =  bvfb(i) /    &
                     ( f1 * sigma_alpha(i,levbot,n)    &
                     + f2 * sigma_t(i,levbot) )
                ak_alpha(i,n)   = 2.0_wp * sigsqh_alpha(i,levbot,n)    &
                     / ( m_alpha(i,levbot,n)**2 - mmsq )
                mmin_alpha(i,n) = m_alpha(i,levbot,n)
          END DO
       END DO
    ELSE
       DO n = 1,naz
!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
             vtmp1(j) = m_alpha(i,levbot,n)
          END DO

          vtmp1(1:nlorms) = vtmp1(1:nlorms)**(-sp1)

!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
                m_alpha(i,levbot,n) =  bvfb(i) /    &
                     ( f1 * sigma_alpha(i,levbot,n)    &
                     + f2 * sigma_t(i,levbot) )
                ak_alpha(i,n)   = sigsqh_alpha(i,levbot,n)    &
                     * vtmp1(j) * sp1
                mmin_alpha(i,n) = m_alpha(i,levbot,n)
          END DO
       END DO
    ENDIF
    !
    !  calculate quantities from the bottom upwards, 
    !  starting one level above bottom.
    !
    DO l = istart,lend,lincr
       !
       !  level beneath present level.
       !
       lbelow = l - lincr 
       !
       !  calculate n/m_m where m_m is maximum permissible value of the vertical
       !  wavenumber (i.e., m > m_m are obliterated) and n is buoyancy frequency.
       !  m_m is taken as the smaller of the instability-induced 
       !  wavenumber (m_sub_m_turb) and that imposed by molecular viscosity
       !  (m_sub_m_mol). since variance at this level is not yet known
       !  use value at level below.
       !
       nlorms = 0
       DO i = il1,il2
          IF (losigma_t(i,lbelow)) THEN
             nlorms = nlorms + 1
             ilorms(nlorms) = i
          END IF
       END DO

!IBM* ASSERT(NODEPS)
       DO j = 1,nlorms
          i = ilorms(j)

          f2mfac=sigmatm(i,lbelow)**2
          f2mod(i,lbelow) =1.0_wp+ 2.0_wp*f2mfac / (f2mfac+sigma_t(i,lbelow)**2 )

          visc = MAX ( visc_mol(i,l), visc_min )
          n_over_m(i) = f2 *f2mod(i,lbelow)*sigma_t(i,lbelow)
          vtmp1(j) = bvfreq(i,l) / n_over_m(i)
          vtmp2(j) = bvfreq(i,l)*kstar/visc
       END DO

#ifdef  __xlC__
       call vcbrt(vtmp2,vtmp2,nlorms)
#else
       vtmp2(1:nlorms) = vtmp2(1:nlorms)**0.33333333_wp
#endif
       DO j = 1,nlorms
          i = ilorms(j)
          m_sub_m_turb = vtmp1(j)
          m_sub_m_mol  = vtmp2(j)/f3
          IF (m_sub_m_turb.GE.m_sub_m_mol) n_over_m(i) = bvfreq(i,l) / m_sub_m_mol
       END DO
       !
       !  calculate cutoff wavenumber at this level.
       !
       DO n = 1,naz
          ! set the default
          m_alpha(il1:il2,l,n) = m_min

          nalpha = 0
          DO j = 1,nlorms
             i = ilorms(j)
             IF (do_alpha(i,n)) THEN
                nalpha = nalpha + 1
                ialpha(nalpha) = i
             END IF
          END DO

!IBM* ASSERT(NODEPS)
          DO j = 1,nalpha
             i = ialpha(j)
             !
             !  calculate trial value (variance at this level is not yet known:
             !  use value at level below). if trial value negative or larger
             !  minimum value (not permitted) then set it to minimum value.
             !
!IBM* mem_delay(v_alpha(i,l,n),20)
!IBM* mem_delay(sigma_alpha(i,lbelow,n),20)
!IBM* mem_delay(sigalpmc(i,lbelow,n),20)
             m_trial = bvfb(i) / ( f1 * ( sigma_alpha(i,lbelow,n) + sigalpmc(i,lbelow,n)) + n_over_m(i) + v_alpha(i,l,n) )

             m_trial = MIN(m_trial,mmin_alpha(i,n))
             IF (m_trial <= 0.0_wp) m_trial = mmin_alpha(i,n)

             m_alpha(i,l,n) = m_trial

             !  do not permit cutoff wavenumber to be less than minimum  value.

             m_alpha(i,l,n) = MAX(m_alpha(i,l,n),m_min)
             !
             !  reset minimum value of cutoff wavenumber if necessary.
             !
             mmin_alpha(i,n) = MIN(mmin_alpha(i,n),m_alpha(i,l,n))
          END DO
       END DO
       !
       !  calculate the hines integral at this level.
       !
       CALL hines_intgrl ( i_alpha,                                     &
                           v_alpha, m_alpha, bvfb, m_min, slope, naz,   &
                           l, il1, il2, nlons, nlevs, nazmth,           &
                           lorms, do_alpha )

       !
       !  calculate the velocity variances at this level.
       !
!IBM* novector
       DO i = il1,il2
          sigfac(i) = densb(i) / density(i,l) * bvfreq(i,l) / bvfb(i) 
       END DO
       DO n = 1,naz
          DO i = il1,il2
             sigsqh_alpha(i,l,n) = sigfac(i) * ak_alpha(i,n) * i_alpha(i,n)
          END DO
       END DO
       CALL hines_sigma ( sigma_t, sigma_alpha, sigsqh_alpha, naz, l, &
                          il1, il2, nlons, nlevs, nazmth )

       CALL hines_sigma ( sigmatm, sigalpmc, sigsqmcw, naz, l,   &
                          il1, il2, nlons, nlevs, nazmth )

       !
       !  if total rms wind zero (no more drag) then set drag to false
       !
       DO i=il1,il2
          IF ( sigma_t(i,l) < EPSILON(1.0_wp)) THEN
             losigma_t(i,l) = .FALSE.
          ENDIF
       ENDDO
       !
       !  end of level loop.
       !
    END DO
#ifdef _PROFILE
     CALL trace_stop ('hines_wavnum', 22)
#endif
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_wavnum

  SUBROUTINE hines_wind (v_alpha,vel_u,vel_v,  &
                         naz,il1,il2,lev1,lev2,nlons,nlevs,nazmth)
    !
    !  this routine calculates the azimuthal horizontal background wind components 
    !  on a longitude by altitude grid for the case of 4 or 8 azimuths for
    !  the hines' doppler spread gwd parameterization scheme.
    !
    !  aug. 7/95 - c. mclandress
    !
    !  output arguement:
    !
    !     * v_alpha   = background wind component at each azimuth (m/s). 
    !     *             (note: first azimuth is in eastward direction
    !     *              and rotate in counterclockwise direction.)
    !
    !  input arguements:
    !
    !     * vel_u     = background zonal wind component (m/s).
    !     * vel_v     = background meridional wind component (m/s).
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1      = first altitude level to use (lev1 >=1). 
    !     * lev2      = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !
    !  constants in data statements.
    !
    !     * cos45 = cosine of 45 degrees.               
    !     * umin  = minimum allowable value for zonal or meridional 
    !     *         wind component (m/s).
    !
    !  subroutine arguements.

    INTEGER  :: naz, il1, il2, lev1, lev2
    INTEGER  :: nlons, nlevs, nazmth
    REAL(wp) :: v_alpha(nlons,nlevs,nazmth)
    REAL(wp) ::  vel_u(nlons,nlevs), vel_v(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER  :: i, l
    REAL(wp) :: u, v, vmu, vpu, cos45, umin
    !----------------------------------------------------------------------- 

    cos45 = 0.7071068_wp 
    umin  = 0.001_wp 

    SELECT CASE (naz) 
       !
    CASE(4)  !  case with 4 azimuths.

       DO l = lev1,lev2
!CDIR NODEP
          DO i = il1,il2
             !
             u = vel_u(i,l)
             IF (ABS(u) < umin)  u = umin
             v = vel_v(i,l)
             IF (ABS(v) < umin)  v = umin 
             !
             v_alpha(i,l,1) = u     ! east
             v_alpha(i,l,2) = v     ! north
             !
             v_alpha(i,l,3) = - u   ! west
             v_alpha(i,l,4) = - v   ! south
             !
          END DO
       END DO
       !
    CASE (8)   !  case with 8 azimuths.
       !
       DO l = lev1,lev2
!CDIR NODEP
          DO i = il1,il2
             !
             u = vel_u(i,l)
             IF (ABS(u) < umin)  u = SIGN(umin,u) 
             v = vel_v(i,l)
             IF (ABS(v) < umin)  v = SIGN(umin,v)
             vpu = v + u
             IF (ABS(vpu) < umin)  vpu = SIGN(umin,vpu)
             vmu = v - u
             IF (ABS(vmu) < umin)  vmu = SIGN(umin,vmu)
             !
             v_alpha(i,l,1) = u                  ! east
             v_alpha(i,l,2) = cos45 * vpu        ! north east
             v_alpha(i,l,3) = v                  ! north
             v_alpha(i,l,4) = cos45 * vmu        ! north west
             !
             v_alpha(i,l,5) = - v_alpha(i,l,1)   ! west
             v_alpha(i,l,6) = - v_alpha(i,l,2)   ! south west
             v_alpha(i,l,7) = - v_alpha(i,l,3)   ! south
             v_alpha(i,l,8) = - v_alpha(i,l,4)   ! south east
             !
          END DO
       END DO
       !
    END SELECT
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_wind

  SUBROUTINE hines_flux ( flux_u, flux_v, flux, drag_u, drag_v,        & 
                          alt, density, densb,                         &
                          m_alpha, ak_alpha, k_alpha,                  &
                          m_min, slope, naz,                           &
                          il1, il2, lev1, lev2, nlons, nlevs, nazmth,  &
                          lorms )
    !
    !  calculate zonal and meridional components of the vertical flux 
    !  of horizontal momentum and corresponding wave drag (force per unit mass)
    !  on a longitude by altitude grid for the hines' doppler spread 
    !  gwd parameterization scheme.
    !  note: only 4 or 8 azimuths can be used.
    !
    !  aug. 6/95 - c. mclandress
    !       2001 - m. charron
    !
    !  output arguements:
    !
    !     * flux_u = zonal component of vertical momentum flux (pascals)
    !     * flux_v = meridional component of vertical momentum flux (pascals)
    !     * drag_u = zonal component of drag (m/s^2).
    !     * drag_v = meridional component of drag (m/s^2).
    !
    !  input arguements:
    !
    !     * alt       = altitudes (m).
    !     * density   = background density (kg/m^3).
    !     * densb     = background density at bottom level (kg/m^3).
    !     * m_alpha   = cutoff vertical wavenumber (1/m).
    !     * ak_alpha  = spectral amplitude factor (i.e., {ajkj} in m^4/s^2).
    !     * k_alpha   = horizontal wavenumber (1/m).
    !     * slope     = slope of incident vertical wavenumber spectrum.
    !     * m_min     = minimum allowable cutoff wavenumber (1/m)
    !     *             for spectral slope of one.
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1      = first altitude level to use (lev1 >=1). 
    !     * lev2      = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !     * lorms     = .true. for drag computation (column selector)
    !
    !  constant in data statement.
    !
    !     * cos45 = cosine of 45 degrees.               
    !
    !  subroutine arguements.
    !

    INTEGER :: naz, il1, il2, lev1, lev2, lev2p
    INTEGER :: nlons, nlevs, nazmth
    REAL(wp) ::  slope, m_min
    REAL(wp) ::  flux_u(nlons,nlevs), flux_v(nlons,nlevs)
    REAL(wp) ::  flux(nlons,nlevs,nazmth)
    REAL(wp) ::  drag_u(nlons,nlevs), drag_v(nlons,nlevs)
    REAL(wp) ::  alt(nlons,nlevs),    density(nlons,nlevs), densb(nlons)
    REAL(wp) ::  m_alpha(nlons,nlevs,nazmth)
    REAL(wp) ::  ak_alpha(nlons,nazmth), k_alpha(nlons,nazmth)

    LOGICAL :: lorms(nlons)
    !
    !  internal variables.
    !
    INTEGER :: i, l, lev1p, lev2m, k
    REAL(wp) ::  cos45, dendz, dendz2
    !-----------------------------------------------------------------------
    cos45 = 0.7071068_wp   
    !
    lev1p = lev1 + 1
    lev2m = lev2 - 1
    lev2p = lev2 + 1
    !
    !  sum over azimuths for case where slope = 1.
    !
    IF ( ABS(slope-1.0_wp) .LT. EPSILON(1.0_wp) )  THEN
       !
       !  case with 4 azimuths.
       !
       IF (naz.EQ.4)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                flux(i,l,:) = ak_alpha(i,:)*k_alpha(i,:)*(m_alpha(i,l,:)-m_min)
                flux_u(i,l) = flux(i,l,1) - flux(i,l,3)
                flux_v(i,l) = flux(i,l,2) - flux(i,l,4)
             END DO
          END DO
       END IF
       !
       !  case with 8 azimuths.
       !
       IF (naz.EQ.8)  THEN
          DO l = lev1,lev2
             DO k = 1, nazmth
                DO i = il1,il2
                  flux(i,l,k) = ak_alpha(i,k)*k_alpha(i,k)*(m_alpha(i,l,k)-m_min)
                END DO
             END DO
             DO i = il1,il2
                flux_u(i,l) = flux(i,l,1) - flux(i,l,5) + cos45 *     &
                     ( flux(i,l,2) - flux(i,l,4) - flux(i,l,6) + flux(i,l,8) )
                flux_v(i,l) = flux(i,l,3) - flux(i,l,7) + cos45 *     &
                     ( flux(i,l,2) + flux(i,l,4) - flux(i,l,6) - flux(i,l,8) )
             END DO
          END DO
       END IF

    END IF
    !
    !  sum over azimuths for case where slope not equal to 1.
    !
    IF ( ABS(slope-1.0_wp) .GT. EPSILON(1.0_wp) )  THEN
       !
       !  case with 4 azimuths.
       !
       IF (naz.EQ.4)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                flux(i,l,:) = ak_alpha(i,:)*k_alpha(i,:)*m_alpha(i,l,:)**slope
                flux_u(i,l) = flux(i,l,1) - flux(i,l,3)
                flux_v(i,l) = flux(i,l,2) - flux(i,l,4)
             END DO
          END DO
       END IF
       !
       !  case with 8 azimuths.
       !
       IF (naz.EQ.8)  THEN
          DO l = lev1,lev2
             DO k = 1, nazmth
                DO i = il1,il2
                  flux(i,l,k) = ak_alpha(i,k)*k_alpha(i,k)*m_alpha(i,l,k)**slope
                END DO
             END DO
             DO i = il1,il2
                flux_u(i,l) = flux(i,l,1) - flux(i,l,5) + cos45 *     &
                     ( flux(i,l,2) - flux(i,l,4) - flux(i,l,6) + flux(i,l,8) )
                flux_v(i,l) = flux(i,l,3) - flux(i,l,7) + cos45 *     &
                     ( flux(i,l,2) + flux(i,l,4) - flux(i,l,6) - flux(i,l,8) )
             END DO
          END DO
       END IF

    END IF
    !
    !  calculate flux from sum.
    !
    DO l = lev1,lev2
       DO i = il1,il2
          flux_u(i,l) = flux_u(i,l) * densb(i) / slope
          flux_v(i,l) = flux_v(i,l) * densb(i) / slope
       END DO
       DO k = 1, nazmth
          DO i = il1,il2
            flux(i,l,k) = flux(i,l,k) * densb(i) / slope
          END DO
       END DO
    END DO
    !
    !  calculate drag at intermediate levels
    !      
    DO l = lev1p,lev2m
!IBM* NOVECTOR
       DO i = il1,il2
          IF (lorms(i)) THEN
             dendz2 = density(i,l) * ( alt(i,l-1) - alt(i,l) ) 
             drag_u(i,l) = - ( flux_u(i,l-1) - flux_u(i,l) ) / dendz2
             drag_v(i,l) = - ( flux_v(i,l-1) - flux_v(i,l) ) / dendz2       
          ENDIF
       END DO
    END DO
    !
    !  calculate drag at intermediate levels using centered differences (not used) 
    !ccc       dendz2 = density(i,l) * ( alt(i,l+1) - alt(i,l-1) )
    !ccc       drag_u(i,l) = - ( flux_u(i,l+1) - flux_u(i,l-1) ) / dendz2
    !ccc       drag_v(i,l) = - ( flux_v(i,l+1) - flux_v(i,l-1) ) / dendz2


    !  drag at first and last levels using one-side differences.
    ! 
!IBM* NOVECTOR
    DO i = il1,il2
       IF (lorms(i)) THEN
          dendz = density(i,lev1) * ( alt(i,lev1) - alt(i,lev1p) ) 
          drag_u(i,lev1) =  flux_u(i,lev1)  / dendz
          drag_v(i,lev1) =  flux_v(i,lev1)  / dendz
       ENDIF
    END DO
!IBM* NOVECTOR
    DO i = il1,il2
       IF (lorms(i)) THEN
          dendz = density(i,lev2) * ( alt(i,lev2m) - alt(i,lev2) )
          drag_u(i,lev2) = - ( flux_u(i,lev2m) - flux_u(i,lev2) ) / dendz
          drag_v(i,lev2) = - ( flux_v(i,lev2m) - flux_v(i,lev2) ) / dendz
       ENDIF
    END DO
    IF (nlevs .GT. lev2) THEN
!IBM* NOVECTOR
       DO i = il1,il2
          IF (lorms(i)) THEN
             dendz = density(i,lev2p) * ( alt(i,lev2) - alt(i,lev2p) )
             drag_u(i,lev2p) = -  flux_u(i,lev2)  / dendz
             drag_v(i,lev2p) = - flux_v(i,lev2)  / dendz
          ENDIF
       END DO
    ENDIF
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_flux

  SUBROUTINE hines_heat ( heat, diffco,                                 &
                          alt, bvfreq, density, sigma_t, sigma_alpha,   &
                          flux, visc_mol, kstar, f1, f2, f3, f5, f6,    &
                          naz, il1, il2, lev1, lev2, nlons, nlevs,      & 
                          nazmth, losigma_t )
    !
    !  this routine calculates the gravity wave induced heating and 
    !  diffusion coefficient on a longitude by altitude grid for  
    !  the hines' doppler spread gravity wave drag parameterization scheme.
    !
    !  This routine can be used for nonzero minimum cutoff wavenumber (m_min)
    !  only in the case of spectral slope=1, in which case m_min is not needed
    !  since its vertical derivative is zero.
    !
    !  aug. 6/95 - c. mclandress
    !  2001      - m. charron  
    !
    !  output arguements:
    !
    !     * heat   = gravity wave heating (k/sec).
    !     * diffco = diffusion coefficient (m^2/sec)
    !
    !  input arguements:
    !
    !
    !     * bvfreq      = background brunt vassala frequency (rad/sec).
    !     * density     = background density (kg/m^3).
    !
    !     * sigma_t     = total rms horizontal wind (m/s).
    !
    !     * visc_mol    = molecular viscosity (m^2/s).
    !     * kstar       = typical gravity wave horizontal wavenumber (1/m).
    !     * slope       = slope of incident vertical wavenumber spectrum.
    !     * f1,f2,f3,f5,f6 = hines's fudge factors.
    !
    !     * il1         = first longitudinal index to use (il1 >= 1).
    !     * il2         = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1        = first altitude level to use (lev1 >=1). 
    !     * lev2        = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons       = number of longitudes.
    !     * nlevs       = number of vertical levels.
    !     * nazmth      = azimuthal array dimension (nazmth >= naz).
    !     * losigma_t   = .true. for total sigma not zero
    !
    USE mo_physical_constants,    ONLY: cpd

    IMPLICIT NONE

    INTEGER ::  naz, il1, il2, lev1, lev2, nlons, nlevs, nazmth
    REAL(wp)::  kstar, f1, f2, f3, f5, f6
    REAL(wp)::  heat(nlons,nlevs), diffco(nlons,nlevs)
    REAL(wp)::  alt(nlons,nlevs), bvfreq(nlons,nlevs), density(nlons,nlevs) 
    REAL(wp)::  sigma_t(nlons,nlevs),  sigma_alpha(nlons,nlevs,nazmth)
    REAL(wp)::  flux(nlons,nlevs,nazmth), visc_mol(nlons,nlevs)
    LOGICAL ::  losigma_t(nlons,nlevs)
    !
    ! internal variables.
    !
    INTEGER  :: i, l, n, lev1p, lev2m
    REAL(wp) :: m_sub_m_turb, m_sub_m_mol, m_sub_m, heatng, dendz2
    REAL(wp) :: visc, visc_min

    REAL(wp) :: dfdz(nlons,nlevs,nazmth)
    !-----------------------------------------------------------------------   

    visc_min = 1.e-10_wp 

    lev1p = lev1 + 1
    lev2m = lev2 - 1   

    DO l = lev1p,lev2m
       DO i = il1,il2
          IF (losigma_t(i,l)) THEN
             dendz2 = density(i,l) * ( alt(i,l-1) - alt(i,l) )
             visc    = MAX ( visc_mol(i,l), visc_min )
             m_sub_m_turb = bvfreq(i,l) / ( f2 * sigma_t(i,l) )
             m_sub_m_mol  = (bvfreq(i,l)*kstar/visc)**0.33333333_wp/f3
             m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
!CDIR UNROLL=8
             dfdz(i,l,:) = ( flux(i,l-1,:) - flux(i,l,:) ) / dendz2 &
                  & * ( f1*sigma_alpha(i,l,:) + bvfreq(i,l)/m_sub_m )
          ENDIF
       END DO
    END DO

    DO i = il1,il2
       IF (losigma_t(i,lev1)) THEN
          dendz2 = density(i,lev1) * ( alt(i,lev1) - alt(i,lev1p) )
          visc    = MAX ( visc_mol(i,lev1), visc_min )
          m_sub_m_turb = bvfreq(i,lev1) / ( f2 * sigma_t(i,lev1) )
          m_sub_m_mol  = (bvfreq(i,lev1)*kstar/visc)**0.33333333_wp/f3
          m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
!CDIR UNROLL=8
          dfdz(i,lev1,:) = -flux(i,lev1,:) / dendz2 &
               & * ( f1*sigma_alpha(i,lev1,:) + bvfreq(i,lev1)/m_sub_m )
       ENDIF
    END DO

    DO i = il1,il2
       IF (losigma_t(i,lev2)) THEN
          dendz2 = density(i,lev2) * ( alt(i,lev2m) - alt(i,lev2) )
          visc    = MAX ( visc_mol(i,lev2), visc_min )
          m_sub_m_turb = bvfreq(i,lev2) / ( f2 * sigma_t(i,lev2) )
          m_sub_m_mol  = (bvfreq(i,lev2)*kstar/visc)**0.33333333_wp/f3
          m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
!CDIR UNROLL=8
          dfdz(i,lev2,:) = ( flux(i,lev2m,:) - flux(i,lev2,:) ) / dendz2 &
               & * ( f1*sigma_alpha(i,lev2,:) + bvfreq(i,lev2)/m_sub_m )
       ENDIF
    END DO
    !
    !  heating and diffusion.

    !
    !  maximum permissible value of cutoff wavenumber is the smaller 
    !  of the instability-induced wavenumber (m_sub_m_turb) and 
    !  that imposed by molecular viscosity (m_sub_m_mol).
    !
    !
    DO l = lev1,lev2
       DO i = il1,il2
          IF (losigma_t(i,l)) THEN
             visc    = MAX ( visc_mol(i,l), visc_min )
             m_sub_m_turb = bvfreq(i,l) / ( f2 * sigma_t(i,l) )
             m_sub_m_mol  = (bvfreq(i,l)*kstar/visc)**0.33333333_wp/f3
             m_sub_m      = MIN ( m_sub_m_turb, m_sub_m_mol )
             heatng = 0.0_wp
             DO n=1,naz
                heatng = heatng - f5 * dfdz(i,l,n)
             ENDDO
             diffco(i,l) = f6 * heatng**0.33333333_wp / m_sub_m**1.33333333_wp
             heat(i,l)   = heatng / cpd
          ENDIF
       END DO
    END DO

    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_heat

  SUBROUTINE hines_sigma (sigma_t,sigma_alpha,sigsqh_alpha,  &
                          naz,lev,il1,il2,nlons,nlevs,nazmth)
    !
    !  this routine calculates the total rms and azimuthal rms horizontal 
    !  velocities at a given level on a longitude by altitude grid for 
    !  the hines' doppler spread gwd parameterization scheme.
    !  note: only four or eight azimuths can be used.
    !
    !  aug. 7/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * sigma_t      = total rms horizontal wind (m/s).
    !     * sigma_alpha  = total rms wind in each azimuth (m/s).
    !
    !  input arguements:
    !
    !     * sigsqh_alpha = portion of wind variance from waves having wave
    !     *                normals in the alpha azimuth (m/s).
    !     * naz       = actual number of horizontal azimuths used (must be 4 or 8).
    !     * lev       = altitude level to process.
    !     * il1       = first longitudinal index to use (il1 >= 1).
    !     * il2       = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * nlons     = number of longitudes.
    !     * nlevs     = number of vertical levels.
    !     * nazmth    = azimuthal array dimension (nazmth >= naz).
    !
    !  subroutine arguements.
    !

    INTEGER  :: lev, naz, il1, il2
    INTEGER  :: nlons, nlevs, nazmth
    REAL(wp) ::  sigma_t(nlons,nlevs)
    REAL(wp) :: sigma_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: sigsqh_alpha(nlons,nlevs,nazmth)
    !
    !  internal variables.
    !
    INTEGER  :: i, n
    REAL(wp) :: sum_even, sum_odd 
    !-----------------------------------------------------------------------     
    !
    !  calculate azimuthal rms velocity for the 4 azimuth case.
    !
    IF (naz.EQ.4)  THEN
!CDIR NODEP
       DO i = il1,il2
          sigma_alpha(i,lev,1) = sigsqh_alpha(i,lev,1)+sigsqh_alpha(i,lev,3)
          sigma_alpha(i,lev,2) = sigsqh_alpha(i,lev,2)+sigsqh_alpha(i,lev,4)
       END DO

       sigma_alpha(il1:il2,lev,1) =  SQRT(sigma_alpha(il1:il2,lev,1))
       sigma_alpha(il1:il2,lev,2) =  SQRT(sigma_alpha(il1:il2,lev,2))

       DO i = il1,il2
          sigma_alpha(i,lev,3) = sigma_alpha(i,lev,1)
          sigma_alpha(i,lev,4) = sigma_alpha(i,lev,2)
       END DO
    END IF
    !
    !  calculate azimuthal rms velocity for the 8 azimuth case.
    !
    IF (naz.EQ.8)  THEN
!CDIR NODEP
       DO i = il1,il2
          sum_odd  = ( sigsqh_alpha(i,lev,1) + sigsqh_alpha(i,lev,3)   &
               + sigsqh_alpha(i,lev,5) + sigsqh_alpha(i,lev,7) ) * 0.5_wp
          sum_even = ( sigsqh_alpha(i,lev,2) + sigsqh_alpha(i,lev,4)   &
               + sigsqh_alpha(i,lev,6) + sigsqh_alpha(i,lev,8) ) * 0.5_wp
          sigma_alpha(i,lev,1) = sigsqh_alpha(i,lev,1) + sigsqh_alpha(i,lev,5) + sum_even
          sigma_alpha(i,lev,2) = sigsqh_alpha(i,lev,2) + sigsqh_alpha(i,lev,6) + sum_odd
          sigma_alpha(i,lev,3) = sigsqh_alpha(i,lev,3) + sigsqh_alpha(i,lev,7) + sum_even
          sigma_alpha(i,lev,4) = sigsqh_alpha(i,lev,4) + sigsqh_alpha(i,lev,8) + sum_odd
       END DO

       sigma_alpha(il1:il2,lev,1) = SQRT(sigma_alpha(il1:il2,lev,1))
       sigma_alpha(il1:il2,lev,2) = SQRT(sigma_alpha(il1:il2,lev,2))
       sigma_alpha(il1:il2,lev,3) = SQRT(sigma_alpha(il1:il2,lev,3))
       sigma_alpha(il1:il2,lev,4) = SQRT(sigma_alpha(il1:il2,lev,4))

       DO i = il1,il2
          sigma_alpha(i,lev,5) = sigma_alpha(i,lev,1)
          sigma_alpha(i,lev,6) = sigma_alpha(i,lev,2)
          sigma_alpha(i,lev,7) = sigma_alpha(i,lev,3)
          sigma_alpha(i,lev,8) = sigma_alpha(i,lev,4)
       END DO
    END IF
    !
    !  calculate total rms velocity.
    !
    DO i = il1,il2
       sigma_t(i,lev) = 0.0_wp
    END DO
    DO n = 1,naz
       DO i = il1,il2
          sigma_t(i,lev) = sigma_t(i,lev) + sigsqh_alpha(i,lev,n)
       END DO
    END DO
    sigma_t(il1:il2,lev) = SQRT(sigma_t(il1:il2,lev))
    !
    !-----------------------------------------------------------------------     
  END SUBROUTINE hines_sigma

  SUBROUTINE hines_intgrl (i_alpha,                                     &
                           v_alpha, m_alpha, bvfb, m_min, slope, naz,   &
                           lev, il1, il2, nlons, nlevs, nazmth,         &
                           lorms, do_alpha)
    !
    !  this routine calculates the vertical wavenumber integral
    !  for a single vertical level at each azimuth on a longitude grid
    !  for the hines' doppler spread gwd parameterization scheme.
    !  note: (1) only spectral slopes of 1, 1.5 or 2 are permitted.
    !        (2) the integral is written in terms of the product qm
    !            which by construction is always less than 1. series
    !            solutions are used for small |qm| and analytical solutions
    !            for remaining values.
    !
    !  aug. 8/95 - c. mclandress
    !  2001      - m. charron
    !  2003      - l. kornblueh
    !
    !  output arguement:
    !
    !     * i_alpha = hines' integral.
    !
    !  input arguements:
    !
    !     * v_alpha = azimuthal wind component (m/s). 
    !     * m_alpha = azimuthal cutoff vertical wavenumber (1/m).
    !     * bvfb    = background brunt vassala frequency at model bottom.
    !     * m_min   = minimum allowable cutoff vertical wavenumber (1/m)
    !     *           for spectral slope of one.
    !     * slope   = slope of initial vertical wavenumber spectrum 
    !     *           (must use slope = 1., 1.5 or 2.)
    !     * naz     = actual number of horizontal azimuths used.
    !     * lev     = altitude level to process.
    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical levels.
    !     * nazmth  = azimuthal array dimension (nazmth >= naz).
    !     * lorms     = .true. for drag computation (column selector)

    !
    !  constants in data statements:
    !
    !     * qmin = minimum value of q_alpha (avoids indeterminant form of integral)
    !     * qm_min = minimum value of q_alpha * m_alpha (used to avoid numerical
    !     *          problems).
    !

    USE mo_io_units,     ONLY: nout, nerr
    USE mo_exception,  ONLY: finish 

    IMPLICIT NONE

    INTEGER  :: lev, naz, il1, il2, nlons, nlevs, nazmth
    REAL(wp) :: i_alpha(nlons,nazmth)
    REAL(wp) :: v_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: m_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: bvfb(nlons), rbvfb(nlons), slope, m_min

    LOGICAL :: lorms(nlons), do_alpha(nlons,nazmth)
    !
    !  internal variables.
    !
    INTEGER  :: i, n, j
    REAL(wp) :: q_alpha, qm, qmm, sqrtqm, q_min, qm_min, ztmp

    !  variables for sparse vector optimization
    INTEGER :: ixi(nlons*nazmth), ixnaz(nlons*nazmth)
    INTEGER :: itx(nlons*nazmth)
    INTEGER :: ilorms(nlons),icond(nlons)

    REAL(wp) :: vinp (nlons*nazmth), vout (nlons*nazmth)
    REAL(wp) :: winp (nlons*nazmth), wout (nlons*nazmth)
    REAL(wp) :: xinp (nlons*nazmth), xout (nlons*nazmth)
    INTEGER  :: vlen,nlorms,ic,ix,nerror

    !-----------------------------------------------------------------------
    !
    !  initialize local scalar and arrays
#ifdef _PROFILE
     CALL trace_start ('hines_intgrl', 23)
#endif

    q_min = 1.0_wp 
    qm_min = 0.01_wp 

!IBM *NOVECTOR
    DO i = il1,il2
      rbvfb(i)=1.0_wp/bvfb(i)
    ENDDO
    nlorms = 0
    DO i = il1,il2
       IF (lorms(i)) THEN
          nlorms = nlorms + 1
          ilorms(nlorms) = i
       END IF
    END DO
    !
    !  for integer value slope = 1.
    !
    IF ( ABS(slope-1.0_wp) < EPSILON(1.0_wp) )  THEN
       DO n = 1,naz
!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
             icond(j) = INT(FSEL(m_min - m_alpha(i,lev,n),0._wp,1._wp))
          END DO

          vlen = 0
          DO j = 1,nlorms
             i = ilorms(j)
             IF (icond(j).GT.0) THEN
                vlen = vlen + 1
                ixi(vlen) = i
             ELSE
                i_alpha(i,n) = 0.0_wp
                do_alpha(i,n) = .FALSE.
             END IF
          END DO

!IBM* ASSERT(NODEPS)
          DO j = 1,vlen
             i = ixi(j)

             q_alpha = v_alpha(i,lev,n) * rbvfb(i)
             vinp(j) = q_alpha**2
             qm      = q_alpha * m_alpha(i,lev,n)
             winp(j) = 1.0_wp-qm
             qmm     = q_alpha * m_min
             xinp(j) = 1.0_wp-qmm
!             IF (ABS(q_alpha) .LT. q_min .OR. ABS(qm).LT. qm_min)) ic = ic + 1
             ztmp     = FSEL(ABS(qm)-qm_min  ,0._wp,1._wp)
             ztmp     = FSEL(vinp(j)-q_min**2,ztmp ,1._wp)
             icond(j) = INT(ztmp)
          END DO

          ! if |qm| is small (rare event) we must do a Taylor expansion

          ic = 0
          IF (vlen.GT.0) THEN
             DO j = 1,vlen
                IF (icond(j).LE.0) CYCLE
                ic = ic + 1
                itx(ic) = ixi(j)
             END DO
          END IF

          vout(1:vlen) = 1._wp/vinp(1:vlen)
          wout(1:vlen) = LOG(winp(1:vlen))
          if (vlen > 0) xout(1:vlen) = LOG(xinp(1:vlen))   ! don't fuse vectorized LOGs

          DO j = 1,vlen
             i = ixi(j)

             i_alpha(i,n) = - (wout(j) - winp(j) - xout(j) + xinp(j) ) * vout(j)
             !
             !  If i_alpha negative due to round off error, set it to zero
             !
             i_alpha(i,n) = MAX( i_alpha(i,n) , 0.0_wp )
          END DO

          ! do Taylor expansion for small |qm|
!IBM* ASSERT(NODEPS,ITERCNT(1))
          DO j = 1,ic
             i = itx(j)
             q_alpha = v_alpha(i,lev,n) * rbvfb(i)
             qm      = q_alpha * m_alpha(i,lev,n)
             qmm     = q_alpha * m_min
             IF ( ABS(q_alpha) < EPSILON(1.0_wp) )  THEN
                i_alpha(i,n) = ( m_alpha(i,lev,n)**2  - m_min**2 ) * 0.5_wp
             ELSE
                i_alpha(i,n) = ( qm**2 * 0.5_wp + qm**3  / 3.0_wp  + qm**4 * 0.25_wp + qm**5  * 0.2_wp   &
                     - qmm**2 * 0.5_wp - qmm**3 / 3.0_wp - qmm**4 * 0.25_wp - qmm**5 * 0.2_wp ) &
                     / q_alpha**2
             END IF
             i_alpha(i,n) = MAX( i_alpha(i,n) , 0.0_wp )
          END DO

       END DO ! n
    END IF
    !
    !  for integer value slope = 2.
    !
    IF ( ABS(slope-2.0_wp) < EPSILON(1.0_wp) )  THEN
       DO n = 1,naz
!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
             q_alpha = v_alpha(i,lev,n) * rbvfb(i)
             vinp(j) = q_alpha**3
             qm = q_alpha * m_alpha(i,lev,n)
             winp(j) = 1.0_wp-qm
!            IF (ABS(q_alpha) .LT. q_min .OR. ABS(qm) .LT. qm_min) ic = ic + 1
             ztmp = FSEL(ABS(qm)-qm_min,0._wp,1._wp)
             ztmp = FSEL(ABS(q_alpha)-q_min,ztmp,1._wp)
             icond(j) = INT(ztmp)
          END DO

          ! if |qm| is small (rare event) we must do a Taylor expansion

          ic = 0
          IF (nlorms.GT.0) THEN
             DO j = 1,nlorms
                IF (icond(j).LE.0) CYCLE
                ic = ic + 1
                itx(ic) = ilorms(j)
             END DO
          END IF

          vout(1:nlorms) = 1._wp/vinp(1:nlorms)
          wout(1:nlorms) = LOG(winp(1:nlorms))

!IBM* ASSERT(NODEPS)
          DO j = 1,nlorms
             i = ilorms(j)
             i_alpha(i,n) = - ( wout(j) + 1.0_wp - winp(j) + winp(j)**2 * 0.5_wp  ) * vout(j)
          END DO

          ! do Taylor expansion for small |qm|
!IBM* ASSERT(NODEPS,ITERCNT(1))
          DO j = 1,ic
             i = itx(j)
             q_alpha = v_alpha(i,lev,n) * rbvfb(i)
             qm = q_alpha * m_alpha(i,lev,n)
             IF ( ABS(q_alpha) < EPSILON(1.0_wp) )  THEN
                i_alpha(i,n) = m_alpha(i,lev,n)**3 / 3.0_wp
             ELSE
                i_alpha(i,n) = ( qm**3/3._wp + qm**4/4._wp + qm**5/5._wp    &
                     + qm**6/6._wp ) / q_alpha**3
             END IF
          END DO
       END DO
    END IF
    !
    !  for real value slope = 1.5
    !
    IF ( ABS(slope-1.5_wp) < EPSILON(1.0_wp) )  THEN
       ic = 0 
       DO n = 1,naz
          DO j = 1,nlorms
             i = ilorms(j)
             !
             q_alpha = v_alpha(i,lev,n) * rbvfb(i)
             qm = q_alpha * m_alpha(i,lev,n)
             !
             !  if |qm| is small then use first 4 terms series of taylor
             !  series expansion of integral in order to avoid
             !  indeterminate form of integral,
             !  otherwise use analytical form of integral.
             !
             IF (ABS(q_alpha) .LT. q_min .OR. ABS(qm) .LT. qm_min)  THEN
                ! taylor series expansion is a very rare event.
                ! do sparse processing separately
                ic = ic+1
                ixi(ic) = i
                ixnaz(ic) = n
             ELSE
                qm     = ABS(qm)
                sqrtqm = SQRT(qm)
                IF (q_alpha .GE. 0.0_wp)  THEN
                   i_alpha(i,n) = ( LOG( (1.0_wp+sqrtqm)/(1.0_wp-sqrtqm) )  &
                        -2.0_wp*sqrtqm*(1.0_wp+qm/3.0_wp) ) / q_alpha**2.5_wp
                ELSE
                   i_alpha(i,n) = 2.0_wp * ( ATAN(sqrtqm) + sqrtqm*(qm/3.0_wp-1.0_wp) ) &
                        / ABS(q_alpha)**2.5_wp
                ENDIF
             ENDIF
             !
          END DO
       END DO
       ! taylor series expansion is a very rare event.
       ! do sparse processing here separately 
!CDIR NODEP
       DO ix = 1, ic
          n = ixnaz(ix)
          i = ixi(ix)
          q_alpha = v_alpha(i,lev,n) * rbvfb(i)
          qm = q_alpha * m_alpha(i,lev,n)   
          IF ( ABS(q_alpha) < EPSILON(1.0_wp) )  THEN
             i_alpha(i,n) = m_alpha(i,lev,n)**2.5_wp / 2.5_wp
          ELSE
             i_alpha(i,n) = ( qm/2.5_wp + qm**2/3.5_wp   &
                  + qm**3/4.5_wp + qm**4/5.5_wp )   &
                  * m_alpha(i,lev,n)**1.5_wp / q_alpha
          END IF
       ENDDO
    END IF
    !
    !  if integral is negative (which in principle should not happen) then
    !  print a message and some info since execution will abort when calculating
    !  the variances.
    !
    DO n = 1,naz
       DO i = il1, il2
          icond(i-il1+1) = INT(FSEL(i_alpha(i,n),0._wp,1._wp))
       END DO
       IF (il2 >= il1) THEN  ! mpuetz: don't fuse with previous loop
          nerror = 0
          DO i = il1, il2
             nerror = nerror + icond(i-il1+1)
          END DO
       END IF

       IF (nerror > 0) THEN

          ! find the first error and dump

          DO i = il1, il2
             IF (icond(i-il1+1) > 0) EXIT
          END DO
          WRITE (nout,*) 
          WRITE (nout,*) '******************************'
          WRITE (nout,*) 'hines integral i_alpha < 0 '
          WRITE (nout,*) '  longitude i=',i
          WRITE (nout,*) '  azimuth   n=',n
          WRITE (nout,*) '  level   lev=',lev
          WRITE (nout,*) '  i_alpha =',i_alpha(i,n)
          WRITE (nout,*) '  v_alpha =',v_alpha(i,lev,n)
          WRITE (nout,*) '  m_alpha =',m_alpha(i,lev,n)
          WRITE (nout,*) '  q_alpha =',v_alpha(i,lev,n)*rbvfb(i)
          WRITE (nout,*) '  qm      =',v_alpha(i,lev,n)*rbvfb(i)*m_alpha(i,lev,n)
          WRITE (nout,*) '******************************'
          WRITE (nerr,*) ' Error: Hines i_alpha integral is negative  '
          CALL finish(' hines_intgrl','Run terminated')
       END IF
  
    END DO
#ifdef _PROFILE
     CALL trace_stop ('hines_intgrl', 23)
#endif
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_intgrl

  SUBROUTINE hines_print (flux_u, flux_v, drag_u, drag_v, alt, sigma_t,    &
                          sigma_alpha, v_alpha, m_alpha,                   &
                          iu_print, iv_print,                              &
                          ilprt1, ilprt2, levprt1, levprt2, naz,           &
                          nlons, nlevs, nazmth)
    !
    !  print out altitude profiles of various quantities from
    !  hines' doppler spread gravity wave drag parameterization scheme.
    !  (note: only for naz = 4 or 8). 
    !
    !  aug. 8/95 - c. mclandress
    !
    !  input arguements:
    !
    !     * iu_print = 1 to print out values in east-west direction.
    !     * iv_print = 1 to print out values in north-south direction.

    !     * ilprt1   = first longitudinal index to print.
    !     * ilprt2   = last longitudinal index to print.
    !     * levprt1  = first altitude level to print.
    !     * levprt2  = last altitude level to print.
    !
    USE mo_io_units, ONLY: nout

    INTEGER  :: naz, ilprt1, ilprt2, levprt1, levprt2
    INTEGER  :: nlons, nlevs, nazmth
    INTEGER  :: iu_print, iv_print
    REAL(wp) :: flux_u(nlons,nlevs), flux_v(nlons,nlevs)
    REAL(wp) :: drag_u(nlons,nlevs), drag_v(nlons,nlevs)
    REAL(wp) :: alt(nlons,nlevs), sigma_t(nlons,nlevs)
    REAL(wp) :: sigma_alpha(nlons,nlevs,nazmth)
    REAL(wp) :: v_alpha(nlons,nlevs,nazmth), m_alpha(nlons,nlevs,nazmth)
    !
    !  internal variables.
    !
    INTEGER :: n_east, n_west, n_north, n_south
    INTEGER ::  i, l
    !-----------------------------------------------------------------------
    !
    !  azimuthal indices of cardinal directions.
    !
    n_east = 1
    IF (naz.EQ.4)  THEN
       n_west  = 3       
       n_north = 2
       n_south = 4       
    ELSE IF (naz.EQ.8)  THEN
       n_west  = 5       
       n_north = 3
       n_south = 7       
    END IF
    !
    !  print out values for range of longitudes.
    !
    DO i = ilprt1,ilprt2
       !
       !  print east-west wind, sigmas, cutoff wavenumbers, flux and drag.
       !
       IF (iu_print.EQ.1)  THEN
          WRITE (nout,*) 
          WRITE (nout,'(a,i3)') 'hines gw (east-west) at longitude i =',i
          WRITE (nout,6005) 
6005      FORMAT (15x,' u ',2x,'sig_e',2x,'sig_t',3x,'m_e',  &
               &            4x,'m_w',4x,'fluxu',5x,'gwdu')
          DO l = levprt1,levprt2
             WRITE (nout,6701) alt(i,l)/1.e3_wp, v_alpha(i,l,n_east),   &
                  &                          sigma_alpha(i,l,n_east), sigma_t(i,l),  &
                  &                          m_alpha(i,l,n_east)*1.e3_wp,   &
                  &                          m_alpha(i,l,n_west)*1.e3_wp,  &
                  &                          flux_u(i,l)*1.e5_wp, drag_u(i,l)*24.0_wp*3600.0_wp
          END DO
6701      FORMAT (' z=',f7.2,1x,3f7.1,2f7.3,f9.4,f9.3)
       END IF
       !
       !  print north-south winds, sigmas, cutoff wavenumbers, flux and drag.
       !
       IF (iv_print.EQ.1)  THEN
          WRITE(nout,*) 
          WRITE(nout,'(a,i3)') 'hines gw (north-south) at longitude i =',i
          WRITE(nout,6006) 
6006      FORMAT (15x,' v ',2x,'sig_n',2x,'sig_t',3x,'m_n',   &
               &            4x,'m_s',4x,'fluxv',5x,'gwdv')
          DO l = levprt1,levprt2
             WRITE (nout,6701) alt(i,l)/1.e3_wp, v_alpha(i,l,n_north),    &
                  &                          sigma_alpha(i,l,n_north), sigma_t(i,l),   &
                  &                          m_alpha(i,l,n_north)*1.e3_wp,    &
                  &                          m_alpha(i,l,n_south)*1.e3_wp,   &
                  &                          flux_v(i,l)*1.e5_wp, drag_v(i,l)*24.0_wp*3600.0_wp
          END DO
       END IF
       !
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_print

  SUBROUTINE hines_exp (darr, data_zmax, alt, alt_exp,     &
                        il1, il2, lev1, lev2, nlons, nlevs)
    !
    !  this routine exponentially damps a longitude by altitude array 
    !  of darr above a specified altitude.
    !
    !  aug. 13/95 - c. mclandress
    !
    !  output arguements:
    !
    !     * darr = modified data array.
    !
    !  input arguements:
    !
    !     * darr    = original data array.
    !     * alt     = altitudes.
    !     * alt_exp = altitude above which exponential decay applied.

    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1    = first altitude level to use (lev1 >=1). 
    !     * lev2    = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical
    !
    !  input work arrays:
    !
    !     * data_zmax = data values just above altitude alt_exp.
    !

    USE mo_io_units,       ONLY: nerr
    USE mo_exception,    ONLY: finish 

    INTEGER  :: il1, il2, lev1, lev2, nlons, nlevs
    REAL(wp) :: alt_exp
    REAL(wp) :: darr(nlons,nlevs), data_zmax(nlons), alt(nlons,nlevs)
    !
    ! internal variables.
    !
    REAL(wp) :: vtmp(nlons)
    INTEGER  :: ialt(nlons)
    INTEGER  :: levbot, levtop, lincr, i, l, nalt, j
    REAL(wp) :: hscale
    !-----------------------------------------------------------------------     
    
    hscale = 5.e3_wp 

    !  index of lowest altitude level (bottom of drag calculation).
    !
    levbot = lev2
    levtop = lev1
    lincr  = 1
    IF (levbot > levtop)  THEN
       levbot = lev1
       levtop = lev2
       lincr  = -1
    ELSE
       WRITE (nerr,*) ' Error: level index not increasing downward '
       CALL finish('hines_exp','Run terminated')
    END IF
    !
    !  data values at first level above alt_exp.
    !
    DO i = il1,il2
       DO l = levtop,levbot,lincr
          IF (alt(i,l) .GE. alt_exp)  THEN
             data_zmax(i) = darr(i,l) 
          END IF
       END DO
    END DO
    !
    !  exponentially damp field above alt_exp to model top at l=1.
    !
    DO l = 1,lev2 
       nalt = 0
       DO i = il1,il2
          IF (alt(i,l) .GE. alt_exp)  THEN
             nalt = nalt + 1
             ialt(nalt) = i
          END IF
       END DO
       DO j = 1,nalt
          i = ialt(j)
          vtmp(j) = (alt_exp-alt(i,l))/hscale
       END DO

       vtmp(1:nalt) = EXP(vtmp(1:nalt))

       DO j = 1,nalt
          i = ialt(j)
          darr(i,l) = data_zmax(i) * vtmp(j)
       END DO
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE hines_exp

  SUBROUTINE vert_smooth (darr, work, coeff, nsmooth,         &
                          il1, il2, lev1, lev2, nlons, nlevs)
    !
    !  smooth a longitude by altitude array in the vertical over a
    !  specified number of levels using a three point smoother. 
    !
    !  note: input array darr is modified on output!
    !
    !  aug. 3/95 - c. mclandress
    !
    !  output arguement:
    !
    !     * darr    = smoothed array (on output).
    !
    !  input arguements:
    !
    !     * darr    = unsmoothed array of data (on input).
    !     * work    = work array of same dimension as darr.
    !     * coeff   = smoothing coefficient for a 1:coeff:1 stencil.
    !     *           (e.g., coeff = 2 will result in a smoother which
    !     *           weights the level l gridpoint by two and the two 
    !     *           adjecent levels (l+1 and l-1) by one).
    !     * nsmooth = number of times to smooth in vertical.
    !     *           (e.g., nsmooth=1 means smoothed only once, 
    !     *           nsmooth=2 means smoothing repeated twice, etc.)
    !     * il1     = first longitudinal index to use (il1 >= 1).
    !     * il2     = last longitudinal index to use (il1 <= il2 <= nlons).
    !     * lev1    = first altitude level to use (lev1 >=1). 
    !     * lev2    = last altitude level to use (lev1 < lev2 <= nlevs).
    !     * nlons   = number of longitudes.
    !     * nlevs   = number of vertical levels.
    !
    !  subroutine arguments.
    !
    IMPLICIT NONE

    INTEGER  :: nsmooth, il1, il2, lev1, lev2, nlons, nlevs
    REAL(wp) :: coeff
    REAL(wp) :: darr(nlons,nlevs), work(nlons,nlevs)
    !
    !  internal variables.
    !
    INTEGER  :: i, l, ns, lev1p, lev2m
    REAL(wp) :: sum_wts
    !-----------------------------------------------------------------------     
    !
    !  calculate sum of weights.
    !
#ifdef _PROFILE
     CALL trace_start ('vert_smooth', 24)
#endif
    sum_wts = coeff + 2.0_wp
    !
    lev1p = lev1 + 1
    lev2m = lev2 - 1
    !
    !  smooth nsmooth times
    !
    DO ns = 1,nsmooth
       !
       !  copy darr into work array.
       !
       DO l = lev1,lev2
          DO i = il1,il2
             work(i,l) = darr(i,l)
          END DO
       END DO
       !
       !  smooth array work in vertical direction and put into darr.
       !
       DO l = lev1p,lev2m
          DO i = il1,il2
             darr(i,l) = (work(i,l+1)+coeff*work(i,l)+work(i,l-1) ) / sum_wts 
          END DO
       END DO
    END DO
#ifdef _PROFILE
     CALL trace_stop ('vert_smooth', 24)
#endif
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vert_smooth

  SUBROUTINE  gw_fronts(krow, kproma, nazmth, rmswind, ani)
    !
    !
    !  may 22/2000 - m. charron
    !
    !  output arguements:
    !
    !     * rmswind    = root mean square gravity wave wind at lowest level (m/s).
    !     * ani        = anisotropy factor (sum over azimuths is one)
    !
    !  input arguments:
    !
    !     * klon     = number of longitudes 
    !     * nazmth     = azimuthal array dimension (nazmth >= naz).

    USE mo_gwspectrum,     ONLY: rmscon, rms_front, front_thres, naz
    !
    !  subroutine arguements.
    !
    INTEGER,  INTENT(IN)  :: krow, kproma, nazmth
    REAL(wp), INTENT(OUT) :: rmswind(kproma)
    REAL(wp), INTENT(OUT) :: ani(kproma,nazmth)
    !
    !  internal variables.
    !
    INTEGER  :: jl, jdir
    INTEGER  :: opp(nazmth)
    REAL(wp) :: angle(kproma), gen(kproma)


    IF ( naz == 8 ) THEN
       opp(1)=5
       opp(2)=6
       opp(3)=7
       opp(4)=8
       opp(5)=1
       opp(6)=2
       opp(7)=3
       opp(8)=4
    ELSE IF ( naz == 4 ) THEN
       opp(1)=3
       opp(2)=4
       opp(3)=1
       opp(4)=2
    END IF

    CALL calculate_gen ( krow, kproma, gen, angle )

    DO jl=1,kproma
       IF ( gen(jl) >= front_thres*2.7777778E-14_wp ) THEN
          rmswind(jl) = rms_front
          jdir=INT( angle(jl)/360.0_wp*REAL(naz,wp) ) + 1
          ani(jl,    :    ) = 0.0_wp
          ani(jl,    jdir ) = 0.5_wp
          ani(jl,opp(jdir)) = 0.5_wp
       ELSE
          rmswind(jl) = rmscon
          ani(jl,:)   = 1.0_wp/REAL(naz,wp)
       END IF
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE gw_fronts

  SUBROUTINE calculate_gen(krow, il2, gen, angle)
    !
    ! determine the value of the frontogenesis function
    !
    !  may 22/2000   - m. charron - first version
    !  march 29/2001 - m. charron - reorganized as a module
    !
    !  output arguement:
    !
    !     * gen    = frontogenesis function in (K/m)^2/hour
    !     * angle  = orientation (in degrees) of the gradiant of temp
    !
    !  input arguements:
    !
    !     * il2     = last longitudinal index to use (1 <= il2 <= nlons).
    
    USE mo_gaussgrid,  ONLY: gl_twomu, gl_sqcst
    USE mo_control,    ONLY: nvclev, vct, nlev
    USE mo_memory_g1a, ONLY: dalpslm1, dalpsmm1, vom1, dm1, alpsm1, tm1
    USE mo_memory_g2a, ONLY: dtlm1, dtmm1, um1, vm1, dudlm1, dvdlm1
    USE mo_gwspectrum, ONLY: emiss_lev, naz
    USE mo_geoloc,     ONLY: ilat

    !
    !  subroutine arguements.
    !
    INTEGER,                   INTENT(IN)  :: krow, il2
    REAL(wp), DIMENSION(il2),  INTENT(OUT) :: gen, angle
    !
    !  internal variables.
    !
    INTEGER                  :: jl, jglat, jrow, iplev
    REAL(wp), DIMENSION(il2) :: dalpsdl, dalpsdm, dudl, dvdl, ps  , vort, div
    REAL(wp), DIMENSION(il2) :: dtdl   , dtdm   , uu  , vv  , uup1, uum1, vvp1
    REAL(wp), DIMENSION(il2) :: vvm1   , ttp1   , ttm1
    REAL(wp), DIMENSION(il2) :: mu, cstmu, zrcst
    REAL(wp)                 :: za, zb, zap1, zbp1
    REAL(wp)                 :: zaplus, zaminus, zbplus, zbminus
    REAL(wp)                 :: term_1, term_2, term_3l, term_3m, term_4l, term_4m
    REAL(wp)                 :: term_5, term_6, term_7 , term_8
    REAL(wp),PARAMETER       :: a=6371000.0_wp, pr=100000.0_wp, kappa=2.0_wp/7.0_wp, pi=3.141592654_wp

!------------------------------------------------------------------------------

    jrow=krow
    iplev=nlev-emiss_lev

    ps     (:)=EXP(alpsm1(:,jrow))
    za        =vct     (iplev)
    zap1      =vct     (iplev+1)
    zb        =vct     (iplev+nvclev)
    zbp1      =vct     (iplev+nvclev+1)
    zaplus    =zap1+za
    zaminus   =zap1-za
    zbplus    =zbp1+zb
    zbminus   =zbp1-zb
    DO jl=1,il2
      jglat     =ilat(jl,jrow)
      mu(jl)    =0.5_wp*gl_twomu(jglat)
      zrcst(jl) =1.0_wp/gl_sqcst(jglat)
    ENDDO
    cstmu  (:)=1.0_wp-mu(:)*mu(:)
    dtdl   (:)=a*cstmu(:)*dtlm1 (:,iplev,jrow)
    dtdm   (:)=a*dtmm1 (:,iplev,jrow)
    dalpsdl(:)=a*cstmu(:)*dalpslm1(:,jrow)
    dalpsdm(:)=a*dalpsmm1(:,jrow)

    uu  (:) = um1 (:,iplev  ,jrow)
    vv  (:) = vm1 (:,iplev  ,jrow)
    uup1(:) = um1 (:,iplev+1,jrow)
    uum1(:) = um1 (:,iplev-1,jrow)
    vvp1(:) = vm1 (:,iplev+1,jrow)
    vvm1(:) = vm1 (:,iplev-1,jrow)
    ttp1(:) = tm1 (:,iplev+1,jrow)
    ttm1(:) = tm1 (:,iplev-1,jrow)
    vort(:) = vom1(:,iplev  ,jrow)
    div (:) = dm1 (:,iplev  ,jrow)

    dudl(:) = dudlm1(:,iplev  ,jrow)*zrcst(:)
    dvdl(:) = dvdlm1(:,iplev  ,jrow)*zrcst(:)

    DO jl=1,il2
      term_1  = (2.0_wp*pr/(zaplus+ps(jl)*zbplus))**kappa
      term_2  = 1.0_wp/(zaminus+ps(jl)*zbminus)
      term_3l = ps(jl)*dalpsdl(jl)*zbplus*term_2
      term_3m = ps(jl)*dalpsdm(jl)*zbplus*term_2
      term_4l = ( dtdl(jl) - 0.25_wp*(ttp1(jl)-ttm1(jl)) * term_3l ) * term_1
      term_4m = ( dtdm(jl) - 0.25_wp*(ttp1(jl)-ttm1(jl)) * term_3m ) * term_1
      term_5  = 0.25_wp * (uup1(jl) - uum1(jl))
      term_6  = 0.25_wp * (vvp1(jl) - vvm1(jl))
      term_7  = 1.0_wp/(a*a*cstmu(jl))
      term_8  = a * term_7
      gen(jl) = -term_7*term_4l**2                                             &
           *(term_8*(dudl(jl)-term_5*term_3l)-mu(jl)*vv(jl)*term_8)            &
           -cstmu(jl)/(a*a) * term_4m**2                                       &
           *(1.0_wp/a*(a*div(jl)-dudl(jl)/cstmu(jl)-term_6*term_3m)            &
           +mu(jl)*vv(jl)*term_8)-term_7* term_4l*term_4m                      &
           *( 1.0_wp/a*(dvdl(jl)-term_6*term_3l)                               &
           +cstmu(jl)/a*(dvdl(jl)/cstmu(jl)-a*vort(jl)-term_5*term_3m)         &
           +2.0_wp*mu(jl)*uu(jl)/a)
      angle(jl) = ATAN2(cstmu(jl)*dtdm(jl),dtdl(jl))*180.0_wp/pi+180.0_wp/REAL(naz,wp)
      IF ( angle(jl) < 0.0_wp ) angle(jl) = angle(jl) + 360.0_wp
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE calculate_gen

END MODULE mo_midatm
