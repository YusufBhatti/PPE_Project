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
SUBROUTINE scan1

  ! Description:
  !
  ! First three scans over the latitude lines in controlling the computations
  ! in Fourier space and the adiabatic part of grid point
  ! calculations including different transport schemes
  !
  ! Method:
  !
  ! This subroutine reads symmetric and antisymmetric parts
  ! of *Fourier coefficients produced by the scans controlling the
  ! inverse *Legendre transforms. It then scans over the latitude
  ! lines three times to perform the 2nd part of computations in
  ! *Fourier space and computes adiabatic tendencies in grid-point
  ! space including transport.
  !
  ! *scan1* is called from *stepon*
  !
  ! Externals.:
  ! *ffti*      inverse *fourier transforms.
  ! *tf1*       time filter (part 1)
  ! *tf2*       time filter (part 2)
  ! *dyn*       adiabatic tendencies (except transported variables)
  !
  ! *gpc*       grid point computations.
  ! *maxwind*   compute maximum wind
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, February 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! T. Diehl, DKRZ, July 1999, parallel version
  ! A. Rhodin, MPI, July 1999, parallel version    (Transpositions)
  ! U. Schlese, DKRZ, October 1999 , preliminary SPITFIRE version 
  !                                  subroutine renamed
  ! I. Kirchner, MPI, July 2000, revision tendency diagnostics
  ! I. Kirchner, MPI, November 2000, date/time control
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! L. Kornblueh, MPI, July 2002, parallelization and blocking of HD model
  ! L. Kornblueh, MPI, August 2002, longitudinal wind velocity derivatives 
  !                                 from Fourier space   
  ! L. Kornblueh, MPI, April 2003, added port test
  ! L. Kornblueh, MPI, October 2003, reactivated AMIP2 global tests
  ! D. Klocke, MPI, Nov 2010, add lnwp for NWP restarts
  ! L. Kornblueh, MPI, Februrary 2012, remove spitfire
  ! S.K. Cheedela, MPI, March 2010, added/renewed scm options
  ! 
  ! for more details see file AUTHORS
  !
  USE mo_kind,          ONLY: dp
  USE mo_physical_constants,     ONLY: earth_radius
  USE mo_control,       ONLY: lcolumn, nlev, nlevp1, &
                              ldiagamip, ltimer, lnwp
  USE mo_port_test,     ONLY: lport
  USE mo_diag_dynamics, ONLY: init_diag_dynamics, diag_dynamics, &
                              print_diag_dynamics
  USE mo_diag_amip2,    ONLY: amip2_global_diag
  USE mo_geoloc,        ONLY: budw_2d
  USE mo_memory_gl,     ONLY: q, xl, xi, xt
  USE mo_memory_g1a,    ONLY: alpsm1, dm1, qm1, tm1, vom1, xlm1, xim1, &
                              xtm1, dalpslm1, dalpsmm1
  USE mo_memory_g1b,    ONLY: alpsf, df, qf, tf, vof, xlf, xif, xtf
  USE mo_memory_g2a,    ONLY: um1, vm1, dtlm1, dtmm1, dudlm1, dvdlm1
  USE mo_memory_g2b,    ONLY: uf, vf, dudlf, dvdlf
  USE mo_memory_g3b,    ONLY: aps
  USE mo_legendre,      ONLY: rnmd
  USE mo_hyb,           ONLY: apsurf
  USE mo_radiation_parameters,  ONLY: ighg
  USE mo_co2,           ONLY: lco2_scenario
  USE mo_scan_buffer,   ONLY: m_bufscan, &
                              vo, d, t, alps, u, v, qte, xlte, xite, xtte, &
                              dalpsl, &
                              dalpsm, dudl, dvdl, dtl, dtm
  USE mo_memory_ls,     ONLY: ld, ltp, lu0, lvo
  USE mo_tracdef,       ONLY: ntrac
  USE mo_mpi,           ONLY: p_parallel_io
  USE mo_decomposition, ONLY: ldc=>local_decomposition
  USE mo_call_trans,    ONLY: legendre_to_spectral, legendre_to_fourier, &
                              fourier_to_gridpoint, fourier_to_legendre, &
                              gridpoint_to_fourier,                      &
                              test_scan_buffer, test_memory_f
  USE mo_global_op,     ONLY: sum_latit_sl, sum_global
  USE mo_column,        ONLY: get_col_states, get_col_dyns, get_col_trans, &
                              update_column, column_diagnostics_ls
  USE mo_test_trans
  USE mo_time_control,  ONLY: lstart, time_step_len, l_diagdyn,  l_trigrad 
  USE mo_advection        
  USE mo_semi_lagrangian, ONLY: setup_semi_lagrangian, extend_semi_lagrangian,&
                                semi_lagrangian_transport, mass_fixer
  USE mo_tpcore,          ONLY: setup_tpcore, tpcore_transport, &
                                tpcore_tendencies
  USE mo_timer,           ONLY: timer_start, timer_stop, &
                                timer_tpcore,   &
                                timer_slt, timer_ffti, timer_fftd
  USE mo_exception,       ONLY: message, message_text
  USE mo_greenhouse_gases,ONLY: interpolate_ghg, interpolate_co2
  USE mo_station_diag,    ONLY : lostation, station_diag
  USE mo_cosp_offline,    ONLY: locospoffl, cosp_offline 
  USE mo_submodel,        ONLY: lanysubmodel
  USE mo_submodel_interface, ONLY: scan1_subm
  IMPLICIT NONE

  ! Local arrays 
  !
  REAL(dp):: etadot(ldc% nproma, nlevp1, ldc% ngpblks)    ! vertical velocity 
!             ulz   (      nlev,  ngl),    &! dyn, maxwind
!             vmaxz (      nlev,  ngl)      ! dyn, maxwind

  REAL(dp):: alpha(nacnst), &
             hw1(nacnst),hw1p(nacnst),hw2(nacnst),hw3(nacnst),                 &
             hw1lat(nacnst,nalat), hw2lat(nacnst,nalat), hw3lat(nacnst,nalat), &
             qfcst(nalond,nalev,nacnst,nalatd,2)
#ifdef SLDIAG
  REAL(dp):: hqfcst(nalond,nalev,nacnst,nalat)
#endif

  !  Local parameters: 
  INTEGER, PARAMETER :: itermn = 1, itermx = 4
 
  !  Local scalars: 
  REAL(dp):: hwps, hwpsm1, pscor, psm1cor, ra, zpref, ztodt

  INTEGER :: iter, jglat, m, jrow

  ! Bounds on this PE
  INTEGER :: ngpblks, nproma

  REAL(dp) :: znaps(SIZE(aps,1),SIZE(aps,2))
  REAL(dp) :: znapsm1(SIZE(aps,1),SIZE(aps,2))

  LOGICAL       :: ldo_advection ! perform transport of tracers

  ! Variables needed by the port test

  INTEGER, ALLOCATABLE :: iseed(:) 
  INTEGER :: isize
  INTEGER :: ipos(3)   
  REAL(dp):: zrn(2)
  LOGICAL, SAVE :: lstart_port = .TRUE.

  !  External subroutines 
  EXTERNAL :: dyn, ewd, sym2, fftd, ffti, gpc, ltd, si2,  sym1, tf1, tf2, &
       maxwind, zonal_winds_gp
           
  !  Executable statements 

  ! bounds on theis PE

  ngpblks = ldc% ngpblks ! number of rows

  ! flag to run transport routine

  ldo_advection = iadvec /= no_advection .AND. .NOT. lcolumn

  ztodt = time_step_len

  ra = 1._dp/earth_radius

  !-- 1.2 Preset spectral components

!$OMP PARALLEL SECTIONS
!$OMP SECTION
  lvo(:,:,:) = 0.0_dp
!$OMP SECTION
  ld(:,:,:)  = 0.0_dp
  lu0(:,:)   = 0.0_dp
!$OMP SECTION
  ltp(:,:,:) = 0.0_dp
  lu0(:,:)   = 0.0_dp
!$OMP END PARALLEL SECTIONS

  IF (l_diagdyn) CALL init_diag_dynamics

  ! Prepare scan buffer (see bufsc1)

  CALL m_bufscan

  CALL test_scan_buffer ('before sym2')

  ! calculations in legendre space and fft taken out of main loop:

  !-- 2.3 Second part of computations in *fourier space
  !       (symmetric/antisymmetric recombination,
  !        zonal derivatives)

  !-- 2.3.1 Computation of *Fourier components from their symmetric
  !         and antisymmetric parts and retrievial of u and v.

  CALL test_memory_f    ('before sym2')

  CALL sym2

  CALL test_scan_buffer ('after sym2')

  CALL legendre_to_fourier

  !-- 2.3.2 East-west derivatives

  CALL ewd

  !-- 2.4 Inverse fast *fourier transforms

  IF (ltimer) CALL timer_start(timer_ffti)
  CALL ffti
  IF (ltimer) CALL timer_stop(timer_ffti)

  CALL fourier_to_gridpoint

  CALL test_scan_buffer ('after fourier_to_gridpoint')

  !-- 2.5 First part of grid point computations

  !>Required initialization and relaxation in single column model
  IF (lcolumn) CALL get_col_states  !relaxation can be moved after dyn

  !-- 2.5.1 Explicit dynamics for eulerian variables

  ! Compute t-dt values at model start.

  IF (lstart .OR. lnwp) THEN
!$OMP PARALLEL SECTIONS
!$OMP SECTION
     vom1   (:,:,:) =   vo (:,:,:)
     dm1    (:,:,:) =    d (:,:,:)
     qm1    (:,:,:) =    q (:,:,:)
     xlm1   (:,:,:) =    xl(:,:,:)
!$OMP SECTION
     xim1   (:,:,:) =    xi(:,:,:)
     tm1    (:,:,:) =    t (:,:,:)
     dtlm1  (:,:,:) =  dtl (:,:,:)
     dtmm1  (:,:,:) =  dtm (:,:,:)
!$OMP SECTION
     um1    (:,:,:) =    u (:,:,:)
     vm1    (:,:,:) =    v (:,:,:)
     dudlm1 (:,:,:) = dudl (:,:,:)
     dvdlm1 (:,:,:) = dvdl (:,:,:)
!$OMP SECTION
     IF (ntrac > 0) THEN
        xtm1(:,:,:,:) = xt(:,:,:,:)
     END IF
     alpsm1  (:,:)  = alps  (:,:)
     dalpslm1(:,:)  = dalpsl(:,:)
     dalpsmm1(:,:)  = dalpsm(:,:)
!$OMP END PARALLEL SECTIONS
  END IF

  IF (lport .AND. lstart_port) THEN
    lstart_port = .FALSE.
    CALL RANDOM_SEED   (size=isize)
    ALLOCATE (iseed(isize))
    CALL SYSTEM_CLOCK  (iseed(1))
    iseed(:) = iseed(1)
    CALL RANDOM_SEED   (put=iseed)
    CALL RANDOM_NUMBER (zrn)
    DEALLOCATE (iseed)
    
    ipos(1) = INT(REAL((SIZE(tm1,1)-1),dp)*zrn(1))+1
    ipos(2) = nlev
    ipos(3) = INT(REAL((SIZE(tm1,3)-1),dp)*zrn(2))+1
    
    WRITE (message_text,'(a,3i4,e30.19,e30.19)')    &
         'Perturbation on T    ', ipos, &
         tm1(ipos(1),ipos(2),ipos(3)),   &
         2.0_dp*SPACING(tm1(ipos(1),ipos(2),ipos(3))) 
    CALL message('scan1', TRIM(message_text))
    
    tm1(ipos(1),ipos(2),ipos(3)) = tm1(ipos(1),ipos(2),ipos(3)) &
         + 2.0_dp*SPACING(tm1(ipos(1),ipos(2),ipos(3))) 
    
    WRITE (message_text,'(a,3i4,e30.19)') &
         'Perturbed value of T ', ipos, &
         tm1(ipos(1),ipos(2),ipos(3))
    CALL message('scan1', TRIM(message_text))
    
  END IF
  
  ! Blank tendencies

!$OMP PARALLEL SECTIONS
!$OMP SECTION
     qte (:,:,:)   = 0.0_dp
!$OMP SECTION
     xlte(:,:,:)   = 0.0_dp
!$OMP SECTION
     xite(:,:,:)   = 0.0_dp
!$OMP SECTION
     xtte(:,:,:,:) = 0.0_dp
!$OMP END PARALLEL SECTIONS

  ! Eulerian advection, energy conversion term
  ! and computation of vertical velocity (etatdot)

  CALL dyn (etadot)
  CALL zonal_winds_gp
  
  !>Handling of large scale forcing in Single Column
  IF (lcolumn) CALL get_col_dyns

  CALL test_scan_buffer ('after dyn')

  ! Second part of time filter

  CALL tf2

  CALL test_scan_buffer ('after tf2')

  IF (l_diagdyn .AND. .NOT. lcolumn) CALL diag_dynamics

  CALL test_scan_buffer ('after diagnostics of the dynamic')

  !-- 2.8 advection ...

  IF (ldo_advection) THEN

    SELECT CASE (iadvec)
    CASE (semi_lagrangian)
      CALL setup_semi_lagrangian
    CASE (tpcore)
      CALL setup_tpcore  
    END SELECT

    IF (iadvec == semi_lagrangian) THEN
      IF (ltimer) CALL timer_start(timer_slt)
      CALL extend_semi_lagrangian

      IF (lstart .OR. lnwp) THEN
        iter = itermx
      ELSE
        iter = itermn
      END IF

#ifdef SLDIAG
      CALL semi_lagrangian_transport(ztodt,ra,iter,etadot, &
           qfcst, hw1lat, hw2lat, hw3lat, hqfcst)
#else
      CALL semi_lagrangian_transport(ztodt,ra,iter,etadot, &
           qfcst, hw1lat, hw2lat, hw3lat)
#endif
      IF (ltimer) CALL timer_stop(timer_slt)
    END IF

    !--  2.5.2 Zonal mass of dry air

    znapsm1(:,:) = budw_2d(:,:)*EXP(alpsm1(:,:))
    znaps(:,:)   = budw_2d(:,:)*aps(:,:)

    !-- 2.5.3 Global mass of dry air

    hwpsm1 = sum_global(znapsm1) 
    hwps   = sum_global(znaps)
    
    IF (iadvec == semi_lagrangian) THEN

      ! hwXlat are return values of the semi Lagrangian
 
      hw1(:) = sum_latit_sl ( hw1lat(:,:))
      hw2(:) = sum_latit_sl ( hw2lat(:,:))
      hw3(:) = sum_latit_sl ( hw3lat(:,:))
      
      DO m = 1, nacnst
        IF (hw3(m) > 0.0_dp) THEN
          alpha(m) = (hw1(m)-hw2(m))/hw3(m)
        ELSE
          alpha(m) = 0.0_dp
        END IF
        hw1p(m) = hw1(m)/hwps
      END DO
    END IF

    ! Coefficients for air mass correction

    zpref = apsurf
    psm1cor = zpref/hwpsm1
    pscor = zpref/hwps

    IF (l_diagdyn .AND. p_parallel_io) THEN
      CALL message('','')
      CALL message('',' Check of air mass correction:' )
      IF (iadvec == semi_lagrangian) THEN
        WRITE (message_text,'(a,f25.15,9x,f25.15,9x,f25.15)') '   hw1     = ', hw1
        CALL message('',message_text)
        WRITE (message_text,'(a,f25.15,9x,f25.15,9x,f25.15)') '   hw2     = ', hw2
        CALL message('',message_text)
        WRITE (message_text,'(a,f25.15,9x,f25.15,9x,f25.15)') '   hw3     = ', hw3
        CALL message('',message_text)
        WRITE (message_text,'(a,f25.15,9x,f25.15,9x,f25.15)') '   alpha   = ', alpha
        CALL message('',message_text)
        WRITE (message_text,'(a,f25.15,9x,f25.15,9x,f25.15)') '   hw1/ps  = ', hw1p
        CALL message('',message_text)
        WRITE (message_text,'(a,f25.15,a,f25.15)') '   psm1    = ', hwpsm1, &
                                            ' ps    = ', hwps              
        CALL message('',message_text)
        WRITE (message_text,'(a,f25.15,a,f25.15)') '   psm1cor = ', psm1cor, &
                                            ' pscor = ', pscor
        CALL message('',message_text)
      ELSE IF (iadvec == tpcore) THEN
        WRITE (message_text,'(a,f25.15,a,f25.15)') '   psm1    = ', hwpsm1, &
                                            ' ps    = ', hwps              
        CALL message('',message_text)
        WRITE (message_text,'(a,f25.15,a,f25.15)') '   psm1cor = ', psm1cor, &
                                            ' pscor = ', pscor
        CALL message('',message_text)
      END IF
    END IF
    
    IF (iadvec == tpcore) THEN
      IF (ltimer) CALL timer_start(timer_tpcore)
      CALL tpcore_transport
      IF (ltimer) CALL timer_stop(timer_tpcore)
    END IF

  ENDIF

  CALL test_scan_buffer ('before calculating advection tendencies')

!-- 2.8.1 Tendencies of advection and mass correction of dry air

  IF (ldo_advection) then

    SELECT CASE (iadvec)
    CASE (semi_lagrangian)
#ifdef SLDIAG
      CALL mass_fixer(qfcst,alpha,psm1cor,pscor,ztodt,hqfcst)
#else
      CALL mass_fixer(qfcst,alpha,psm1cor,pscor,ztodt)
#endif
    CASE (tpcore)
      CALL tpcore_tendencies (psm1cor, pscor) 
    END SELECT
  ENDIF

  !>Handling of large scale forcing for tracers in Single Column
  IF (lcolumn) CALL get_col_trans  !not necessary to be here

!-- 3. Final subscan

!-- 3.2 Start of mainloop for final subscan

  CALL test_scan_buffer ('before start of final subscan')

  !-- 3. Final subscan

  CALL si1_extended_const

  ! First part of time filtering for next time step

  CALL tf1

  IF ((l_trigrad) .AND. ighg.NE.0) CALL interpolate_ghg
  IF (lco2_scenario) CALL interpolate_co2

  !-- Initialisation of surface and soil temperatures.

  IF ( lstart ) THEN
    DO jrow = 1, ngpblks 
      CALL initemp(jrow)
    ENDDO
  ENDIF

  !-- 3.2 Start of mainloop for final subscan

!CSD$ PARALLEL DO PRIVATE(nproma,jglat,jrow) SCHEDULE(STATIC,1)
!$OMP PARALLEL DO PRIVATE(nproma,jglat,jrow) SCHEDULE(STATIC,1)
  DO jrow = 1, ngpblks          ! local number of rows

    IF ( jrow == ldc% ngpblks ) THEN
      nproma = ldc% npromz
    ELSE
      nproma = ldc% nproma
    END IF

    IF ( ldc% lreg ) THEN
      jglat = ldc%glat(jrow)    ! global index north -> south
    ELSE
      jglat = 1
    END IF

    !-- 3.4 2nd part of grid point computations

    IF(lcolumn) CALL column_diagnostics_ls  ! single column diagnostics largescale

    !-- 3.4.2 Physics 

    CALL gpc(jrow)

    !> Update in single column model
    IF (lcolumn) CALL update_column

    !-- 4. Grid point contributions to the semi implicit

    CALL si1(jrow)

!DIR$ CONCURRENT
    q     (1:nproma,:,jrow)   =  qm1(1:nproma,:,jrow)
!DIR$ CONCURRENT
    xl    (1:nproma,:,jrow)   =  xlm1(1:nproma,:,jrow) 
!DIR$ CONCURRENT
    xi    (1:nproma,:,jrow)   =  xim1(1:nproma,:,jrow) 
!DIR$ CONCURRENT
    xt    (1:nproma,:,:,jrow) =  xtm1(1:nproma,:,:,jrow)

    !-- 3.8 Compute contributions to symmetric and
    !       antisymmetric Fourier components.
    !-- 4.1 End of main loop for final subscan

  END DO
!$OMP END PARALLEL DO
!CSD$ END PARALLEL DO 

  CALL test_scan_buffer ('after end of final subscan')

!>>>mgs
  !-- Submodel global diagnostics
  IF (lanysubmodel) THEN
    CALL scan1_subm
  END IF
!<<<mgs

  !-- 3.5 Direct fast *fourier transforms

  CALL gridpoint_to_fourier 

  IF (ltimer) CALL timer_start(timer_fftd)
  CALL fftd
  IF (ltimer) CALL timer_stop(timer_fftd)

  !-- 3.6 Semi implicit adjustment (part 2)

  CALL si2

  ! second part of tendency diagnostics called after SI2

  CALL fourier_to_legendre

  !-- 3.8 Compute contributions to symmetric and
  !       antisymmetric *fourier components.

  CALL sym1
  !-- 4. Direct *Legendre transforms

  CALL ltd

  !-- Transposition Legendre -> spectral

  CALL legendre_to_spectral  

  !-- 5. Complete the scans

  !-- 5.1 Release space

  rnmd(:) = 0.0_dp ! Legendre fields

  !-- 5.2 Complete statistics

  CALL maxwind

!Store variables at station locations of previous time step
!The real winds are calculated inside station_diag for the station
  IF (lostation) THEN
     CALL station_diag
  END IF
  IF (locospoffl) THEN
     CALL cosp_offline
  END IF

  IF (.NOT. lcolumn) THEN
    IF (l_diagdyn) CALL print_diag_dynamics
    IF (ldiagamip) CALL amip2_global_diag
  END IF

  IF (lanysubmodel) CALL scan1_subm

!$OMP PARALLEL SECTIONS
!$OMP SECTION
  vom1  (:,:,:) = vof  (:,:,:)
  dm1   (:,:,:) = df   (:,:,:)
  tm1   (:,:,:) = tf   (:,:,:)
  dtlm1 (:,:,:) = dtl  (:,:,:)
!$OMP SECTION
  dtmm1 (:,:,:) = dtm  (:,:,:)
  qm1   (:,:,:) = qf   (:,:,:)
  xlm1  (:,:,:) = xlf  (:,:,:)
  xim1  (:,:,:) = xif  (:,:,:)
!$OMP SECTION
  um1   (:,:,:) = uf   (:,:,:)
  vm1   (:,:,:) = vf   (:,:,:)
  dudlm1(:,:,:) = dudlf(:,:,:)
  dvdlm1(:,:,:) = dvdlf(:,:,:)
!$OMP SECTION
  xtm1  (:,:,:,:) = xtf  (:,:,:,:)
  alpsm1(:,:)   = alpsf(:,:)
  dalpslm1(:,:) = dalpsl(:,:)
  dalpsmm1(:,:) = dalpsm(:,:)
!$OMP END PARALLEL SECTIONS

END SUBROUTINE scan1
