!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_m7.f90
!!
!! \brief
!! Module to provide interface to the M7 aerosol microphysics scheme
!!
!! \author Martin G. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin G. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# The original code is from J. Feichter, J. Wilson and E. Vignatti, JRC Ispra 
!! and was adapted for ECHAM by P. Stier, Oxford. Other contributions include 
!! D. O'Donnell, K. Zhang and others
!!   -# M.G. Schultz (FZ Juelich) - new code structure for integration into echam6-hammoz (2009-09-24)
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!!   -# H. Kokkola (FMI) - moved and renamed m7_interface into a general ham_subm_interface routine
!!
!! \limitations
!! None
!!
!! \details
!! This module contains the m7_interface routine and all individual routines
!! which make up M7. Parameter lists and flags are defined in mo_ham_m7ctl.
!! This module contains the following subroutines which used to be individual files.
!!       m7_cumulative_normal      (renamed from m7_cumnor)
!!       m7
!!       m7_averageproperties
!!       m7_kappa
!!       m7_equiz
!!       m7_equimix
!!       m7_equil
!!       m7_h2so4_cs
!!       m7_prod_cond
!!       m7_nuck
!!       m7_dnum
!!       m7_dconc
!!       m7_dconc_soa
!!       m7_coaset
!!       m7_concoag
!!       m7_delcoa 
!!
!! \bibliographic_references
!!   - Vignati E., Wilson J. and Stier P., M7: a size resolved aerosol mixture module 
!!     for the use in global aerosol models, JGR 109, D22 202, doi:10.1029/2003JD004 485, 2004.
!!   - Stier P. et al, The aerosol-climate model ECHAM5-HAM, ACP, 5, 1125–1156, 2005
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ham_m7

IMPLICIT NONE

PRIVATE

PUBLIC      :: m7, m7_cumulative_normal

CONTAINS


! ---------------------------------------------------------------------------

SUBROUTINE m7_cumulative_normal ( arg, presult, ccum )
  !
  !*******************************************************************************
  !
  !! CUMNOR computes the cumulative normal distribution.
  !
  !
  !     the integral from -infinity to x of
  !          (1/sqrt(2*pi)) exp(-u*u/2) du
  !
  !  Author:
  !  -------
  !  Original source:
  !
  !    W. J. Cody    Mathematics and Computer Science Division
  !                  Argonne National Laboratory
  !                  Argonne, IL 60439
  !
  !    DCDFLIB is attributed to Barry Brown, James Lovato, and Kathy Russell
  !            bwb@odin.mda.uth.tmc.edu.
  !
  !    Adopted to ECHAM/M7:
  !
  !    Philip Stier  (MPI-MET)                    2001
  !
  !
  !  Reference:
  !  ----------
  !
  !    W D Cody, 
  !    "ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special 
  !    Function Routines and Test Drivers"
  !    ACM Transactions on Mathematical Software,
  !    Volume 19, 1993, pages 22-32.
  !
  !  Parameters:
  !
  !     ARG --> Upper limit of integration.
  !                                        X is double precision
  !
  !     RESULT <-- Cumulative normal distribution.
  !                                        RESULT is double precision
  !
  !     CCUM <-- Complement of Cumulative normal distribution.
  !                                        CCUM is double precision
  !
  !
  ! Original Comments:
  !
  !
  ! This function evaluates the normal distribution function:
  !
  !                              / x
  !                     1       |       -t*t/2
  !          P(x) = ----------- |      e       dt
  !                 sqrt(2 pi)  |
  !                             /-oo
  !
  !   The main computation evaluates near-minimax approximations
  !   derived from those in "Rational Chebyshev approximations for
  !   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
  !   This transportable program uses rational functions that
  !   theoretically approximate the normal distribution function to
  !   at least 18 significant decimal digits.  The accuracy achieved
  !   depends on the arithmetic system, the compiler, the intrinsic
  !   functions, and proper selection of the machine-dependent
  !   constants.
  !
  !  Explanation of machine-dependent constants.
  !
  !   MIN   = smallest machine representable number.
  !
  !   EPS   = argument below which anorm(x) may be represented by
  !           0.5  and above which  x*x  will not underflow.
  !           A conservative value is the largest machine number X
  !           such that   1.0 + X = 1.0   to machine precision.
  !
  !  Error returns
  !
  !  The program returns  ANORM = 0     for  ARG .LE. XLOW.
  !
  !  Author: 
  !
  !    W. J. Cody
  !    Mathematics and Computer Science Division
  !    Argonne National Laboratory
  !    Argonne, IL 60439
  !
  !  Latest modification: March 15, 1992
  !
  USE mo_kind, ONLY: dp
  !
  IMPLICIT NONE
  !
  REAL(dp), PARAMETER, DIMENSION ( 5 ) :: a = (/ &
       2.2352520354606839287e00_dp, &
       1.6102823106855587881e02_dp, &
       1.0676894854603709582e03_dp, &
       1.8154981253343561249e04_dp, &
       6.5682337918207449113e-2_dp /)
  REAL(dp) :: arg
  REAL(dp), PARAMETER, DIMENSION ( 4 ) :: b = (/ &
       4.7202581904688241870e01_dp, &
       9.7609855173777669322e02_dp, &
       1.0260932208618978205e04_dp, &
       4.5507789335026729956e04_dp /)
  REAL(dp), PARAMETER, DIMENSION ( 9 ) :: c = (/ &
       3.9894151208813466764e-1_dp, &
       8.8831497943883759412e00_dp, &
       9.3506656132177855979e01_dp, &
       5.9727027639480026226e02_dp, &
       2.4945375852903726711e03_dp, &
       6.8481904505362823326e03_dp, &
       1.1602651437647350124e04_dp, &
       9.8427148383839780218e03_dp, &
       1.0765576773720192317e-8_dp /)
  REAL(dp) :: ccum
  REAL(dp), PARAMETER, DIMENSION ( 8 ) :: d = (/ &
       2.2266688044328115691e01_dp, &
       2.3538790178262499861e02_dp, &
       1.5193775994075548050e03_dp, &
       6.4855582982667607550e03_dp, &
       1.8615571640885098091e04_dp, &
       3.4900952721145977266e04_dp, &
       3.8912003286093271411e04_dp, &
       1.9685429676859990727e04_dp /)
  REAL(dp) :: del
!@@@ REAL(dp) :: dpmpar
  REAL(dp) :: eps
  INTEGER :: i
  REAL(dp) :: zmin
  REAL(dp), PARAMETER, DIMENSION ( 6 ) :: p = (/ &
       2.1589853405795699e-1_dp, &
       1.274011611602473639e-1_dp, &
       2.2235277870649807e-2_dp, &
       1.421619193227893466e-3_dp, &
       2.9112874951168792e-5_dp, &
       2.307344176494017303e-2_dp /)
  REAL(dp), PARAMETER, DIMENSION ( 5 ) :: q = (/ &
       1.28426009614491121e00_dp, &
       4.68238212480865118e-1_dp, &
       6.59881378689285515e-2_dp, &
       3.78239633202758244e-3_dp, &
       7.29751555083966205e-5_dp /)
  REAL(dp) :: presult
  REAL(dp), PARAMETER :: root32 = 5.656854248_dp
  REAL(dp), PARAMETER :: sixten = 16.0_dp
  REAL(dp) :: temp
  REAL(dp), PARAMETER :: sqrpi = 3.9894228040143267794e-1_dp
  REAL(dp), PARAMETER :: thrsh = 0.66291_dp
  REAL(dp) :: x
  REAL(dp) :: xden
  REAL(dp) :: xnum
  REAL(dp) :: y
  REAL(dp) :: xsq
  !
  !  Machine dependent constants
  !
  eps = EPSILON ( 1.0_dp ) * 0.5_dp
  !
  !@@@ Simplified calculation of the smallest machine representable number
  !    (Higher accuracy than needed!)
  !
  !@@@ min = dpmpar(2)

  zmin = EPSILON ( 1.0_dp )

  x = arg
  y = ABS ( x )

  IF ( y <= thrsh ) THEN
     !
     !  Evaluate  anorm  for  |X| <= 0.66291
     !
     IF ( y > eps ) THEN
        xsq = x * x
     ELSE
        xsq = 0.0_dp
     END IF

     xnum = a(5) * xsq
     xden = xsq
     DO i = 1, 3
        xnum = ( xnum + a(i) ) * xsq
        xden = ( xden + b(i) ) * xsq
     END DO
     presult = x * ( xnum + a(4) ) / ( xden + b(4) )
     temp = presult
     presult = 0.5_dp + temp
     ccum = 0.5_dp - temp
     !
     !  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
     !
  ELSE IF ( y <= root32 ) THEN

     xnum = c(9) * y
     xden = y
!DIR$ UNROLL
!CDIR UNROLL=7
     DO i = 1, 7
        xnum = ( xnum + c(i) ) * y
        xden = ( xden + d(i) ) * y
     END DO
     presult = ( xnum + c(8) ) / ( xden + d(8) )
     xsq = AINT ( y * sixten ) / sixten
     del = ( y - xsq ) * ( y + xsq )
     presult = EXP(-xsq*xsq*0.5_dp) * EXP(-del*0.5_dp) * presult
     ccum = 1.0_dp - presult

     IF ( x > 0.0_dp ) THEN
        temp = presult
        presult = ccum
        ccum = temp
     END IF
     !
     !  Evaluate  anorm  for |X| > sqrt(32).
     !
  ELSE

     presult = 0.0_dp
     xsq = 1.0_dp / ( x * x )
     xnum = p(6) * xsq
     xden = xsq
     DO i = 1, 4
        xnum = ( xnum + p(i) ) * xsq
        xden = ( xden + q(i) ) * xsq
     END DO

     presult = xsq * ( xnum + p(5) ) / ( xden + q(5) )
     presult = ( sqrpi - presult ) / y
     xsq = AINT ( x * sixten ) / sixten
     del = ( x - xsq ) * ( x + xsq )
     presult = EXP ( - xsq * xsq * 0.5_dp ) * EXP ( - del * 0.5_dp ) * presult
     ccum = 1.0_dp - presult  

     IF ( x > 0.0_dp ) THEN
        temp = presult
        presult = ccum
        ccum = temp
     END IF

  END IF

  IF ( presult < zmin ) THEN
     presult = 0.0_dp
  END IF

  IF ( ccum < zmin ) THEN
     ccum = 0.0_dp
  END IF

END SUBROUTINE m7_cumulative_normal

! ---------------------------------------------------------------------------

SUBROUTINE m7(kproma, kbdim,   klev, krow, &
              papp1,  ptp1,    prhp1,      &
              !>>dod soa
              pgas,                        &
              !<<dod
              pso4g,  pdso4g,              &
              paerml, paernl,              &
              pm6rp,  pm6dry,  prhop, pww, &
              pipr,   paclc,               &
              pdz,    pgrvol, ppbl)
  !
  ! Authors:
  ! ---------
  !
  ! 2000 E. Vignati (original code)
  ! 2001 P. Stier (f90-version, changes, comments)
  ! 2008 J. Kazil (changes, comments)
  ! 2009 D.O'Donnell, MPI-Met, Hamburg: introduction of secondary organics
  !                            and miscellaneous other changes (commented below)
  !
  ! Interface:
  ! ----------
  !
  ! *m7* is called from *ham_subm_interface*
  !
  ! Purpose
  ! ---------
  !
  ! *m7* Aerosol model for the system so4,bc,oc,ss,dust in 7 modes.
  !
  ! Externals and method
  ! --------------------
  !
  ! *m7_averageproperties* 
  !     calculates the average mass for all modes and the particle 
  !     dry radius and density for the insoluble modes.
  !  
  ! *m7_equiz*   
  !     calculates the ambient radius of sulfate particles
  !
  ! *m7_equimix* 
  !     calculates the ambient radius of so4, bc, oc, du, particles
  !
  ! *m7_equil* 
  !     calculates the ambient radius of so4, ss particles 
  !
  ! *m7_h2so4_cs*
  !     calculates the aerosol H2SO4 condensation sink
  !
  ! *m7_prod_cond*
  !     calculates new gas phase and aerosol H2SO4 after H2SO4 production/
  !     condensation
  !
  ! *m7_nuck*
  !     calculates new gas phase H2SO4 and newly produced aerosol after H2SO4
  !     nucleation
  !
  ! *m7_dnum*    
  !     calculates new aerosol numbers and masses after nucleation and
  !     coagulation
  !
  ! *m7_dconc*   
  !     repartitions aerosol number and mass between the the modes to account
  !     for condensational growth and the formation of an accumulation mode from
  !     the upper tail of the aitken mode and of a coarse mode from the upper
  !     tail of the accumulation mode
  
  USE mo_kind,           ONLY: dp
  USE mo_ham_m7ctl,      ONLY: nsnucl, nonucl, lnucl_stat, nwater
  USE mo_ham,            ONLY: subm_ngasspec, nsol
  USE mo_ham_subm_species, ONLY: isubm_so4g
  USE mo_ham_m7_nucl_diag, ONLY: sample_nucl_statistics
  !>>dod soa
  USE mo_ham,            ONLY: nsoa, naerocomp, nclass, nsol, lscond, lscoag
  !<<dod
  !>>dod timers
  USE mo_control,        ONLY: ltimer
  USE mo_hammoz_timer,   ONLY: timer_start,          timer_stop,            &
                               timer_ham_m7_avg,     timer_ham_m7_hygro,    &
                               timer_ham_m7_cs,      timer_ham_m7_cond,     &
                               timer_ham_m7_nucl,                           &
                               timer_ham_m7_dnum,    timer_ham_m7_redistr
  !<<dod

  IMPLICIT NONE
  
  !
  ! Input/output:
  !
  
  INTEGER :: kproma,                     & ! ECHAM5 geographic block number of locations
             kbdim,                      & !   "    geographic block maximum number of locations
             klev,                       & !   "    numer of levels
             krow                          !   "    geographic block number
  
  REAL(dp):: papp1(kbdim,klev),          & ! pressure [Pa] at t+1
             ptp1(kbdim,klev),           & ! temperature [K] at t+1
             prhp1(kbdim,klev)             ! relative humidity at t+1 [0,1]
  !>>dod soa
  REAL(dp) :: pgas(kbdim,klev,subm_ngasspec)    ! gas phase tracers 
  !<<dod

  REAL(dp):: pso4g(kbdim,klev)             ! on input : Gas phase H2SO4 concentration at time t-1 [molec. cm-3]
                                           ! on output: Gas phase H2SO4 concentration at time t+1 [molec. cm-3]
  
  REAL(dp):: pdso4g(kbdim,klev)            ! d[H2SO4(g)]/dt [cm-3 s-1] for processes preceding m7
  
  REAL(dp):: paerml(kbdim,klev,naerocomp), & ! aerosol mass for individual compounds [molec. cm-3 for sulfate and ug m-3 for bc, oc, ss, and dust]
             paernl(kbdim,klev,nclass),    & ! aerosol number for each mode [cm-3]
             pm6rp(kbdim,klev,nclass),     & ! mean mode actual radius (wet for soluble and dry for insoluble modes) [cm]
             pm6dry(kbdim,klev,nsol),    & ! dry radius for soluble modes [cm]
             prhop(kbdim,klev,nclass),     & ! mean mode particle density [g cm-3]
             pww(kbdim,klev,nclass)          ! aerosol water content for each mode [kg(water) m-3(air)]
  
  REAL(dp):: pipr(kbdim,klev),           & ! ion pair production rate [cm-3 s-1]
             paclc(kbdim,klev),          & ! cloud cover [0,1]
             pdz(kbdim,klev),            & ! layer thickness for diagnostics [m]
             pgrvol(kbdim,klev),         & ! grid box volume for diagnostics [m-3]
             ppbl(kbdim)                   ! Planetary boundary layer top level
  !
  ! Local variables:
  !
  
  REAL(dp):: zso4_5(kbdim,klev)          ! Number of H2SO4 molecules condensed on insoluble
  REAL(dp):: zso4_6(kbdim,klev)          ! mode x [molec. cm-3] (calculated in m7_prod_cond,
  REAL(dp):: zso4_7(kbdim,klev)          ! used in m7_concoag)
  
  REAL(dp):: zttn(kbdim,klev,naerocomp)    ! average mass for single compound in each mode [in molec. for sulfate and in ug for bc, oc, ss, and dust]
  REAL(dp):: zcs(kbdim,klev)             ! H2SO4 condensation sink of the aerosol population [s-1]
  REAL(dp):: zcsi(kbdim,klev,nclass)       ! H2SO4 condensation sink of the individual aerosol modes [s-1]
  REAL(dp):: zanew(kbdim,klev)           ! Number of nucleated particles [cm-3] over one timstep
  REAL(dp):: za4delt(kbdim,klev,naerocomp) ! Change in H2SO4 contents [cm-3] of the aerosol modes over one timstep
  
  !--- 0) Initialisations ------------------------------------------------------
  
  pm6dry(:,:,:)  = 0.0_dp 
  pm6rp(:,:,:)   = 0.0_dp
  zttn(:,:,:)    = 0.0_dp
  prhop(:,:,:)   = 0.0_dp
  pww(:,:,:)     = 0.0_dp 
  zso4_5(:,:)    = 0.0_dp
  zso4_6(:,:)    = 0.0_dp
  zso4_7(:,:)    = 0.0_dp
  za4delt(:,:,:) = 0.0_dp
  zanew(:,:)     = 0.0_dp
  
  !--- 1) Calculation of particle properties under ambient conditions ----------
  !
  !--- 1.1) Calculate mean particle mass for all modes 
  !         and dry radius and density for the insoluble modes.

  !>>dod timers
  IF (ltimer) CALL timer_start(timer_ham_m7_avg)
  !<<dod

  CALL m7_averageproperties(kproma, kbdim, klev, krow, paernl, paerml, zttn, pm6rp, prhop)
  
  !>>dod timers
  IF (ltimer) THEN 
     CALL timer_stop(timer_ham_m7_avg)
     CALL timer_start(timer_ham_m7_hygro)
  END IF
  !<<dod

  !--- 1.2) Calculate ambient count median radii and density 
  !         for lognormal distribution of particles.
  !
  
  IF (nwater == 1) THEN
  
     CALL m7_kappa(kproma, kbdim, klev, krow, prhp1, paernl, zttn, ptp1, &
                   pm6dry, pm6rp, prhop, pww)
  
  ELSE
    
     !         Sulfate particles:
    
     !>>dod soa 
     IF (nsoa == 0) &
     !<<dod
          CALL m7_equiz(kproma, kbdim, klev, krow,  &
                        papp1,  zttn,  ptp1,   &
                        prhp1,  pm6rp, pm6dry, &
                        prhop,  pww,   paernl  )
             
     !         Mixed particles with sulfate, black and organic carbon, and dust: 
    
     CALL m7_equimix(kproma, kbdim, klev, krow,  &
                     papp1,  zttn,  ptp1,   &
                     prhp1,  pm6rp, pm6dry, &
                     prhop,  pww,   paernl  )
    
     !         Accumulation and coarse mode particles in presence of
     !         sea salt particles:
     !         dod: deleted zhplus and added krow
     CALL m7_equil(kproma, kbdim,  klev,  krow,  prhp1, paerml, paernl, &
                   pm6rp,  pm6dry, pww,   prhop,  ptp1    )
    
  ENDIF

  !>>dod timers
  IF (ltimer) THEN
     CALL timer_stop(timer_ham_m7_hygro)
     CALL timer_start(timer_ham_m7_cs)
  END IF
  !<<dod
  
  !         Calculate the H2SO4 condensation sink of the entire aerosol
  !         population and of the individual modes:
  
  IF (lscond) CALL m7_h2so4_cs(kproma,kbdim,klev,ptp1,papp1,paernl,pm6rp,zcs,zcsi)

  !>>dod timers
  IF (ltimer) CALL timer_stop(timer_ham_m7_cs)
  !<<dod

  !--- 1.4) Sample the cloud-free volume as function of T, RH, [H2SO4(g)],'
  !         H2SO4 condensation sink, and ionization rate:

  IF (lnucl_stat) CALL sample_nucl_statistics(kproma,kbdim,klev, &
                        ptp1,prhp1,pso4g,pipr,zcs,paclc,pgrvol)
  
  !--- 2) Production/condensation ----------------------------------------------
  !
  !       Calculate H2SO4(g) production/condensation over a time step, and
  !       distribute the condensing gas phase H2SO4 onto the individual aerosol
  !       modes:
  
  !>>dod timers
  IF (ltimer) CALL timer_start(timer_ham_m7_cond)
  !<<dod

  IF (lscond) CALL m7_prod_cond(kproma, kbdim,  klev,   krow, &
                                pso4g,  pdso4g, paerml,       &
                                zso4_5, zso4_6, zso4_7,       &
                                zcs,    zcsi,   paclc,  pdz)

  !>>dod timers
  IF (ltimer) THEN
     CALL timer_stop(timer_ham_m7_cond)
     CALL timer_start(timer_ham_m7_nucl)
  END IF
  !<<dod

  !--- 3) Nucleation -----------------------------------------------------------
  !
  !       Calculate the nucleation rate, the number of nucleated particles, and
  !       change in gas phase H2SO4:

  IF (nsnucl+nonucl.gt.0) CALL m7_nuck(kproma,  kbdim,  klev,   krow,          &
                                       papp1,   ptp1,   prhp1,  paclc,   pipr, &
                                       pso4g,   zcs,    zanew,  za4delt, pdz, ppbl)

  !>>dod timers
  IF (ltimer) THEN
     CALL timer_stop(timer_ham_m7_nucl)
     CALL timer_start(timer_ham_m7_dnum)
  END IF
  !<<dod

  !--- 4) Coagulation ----------------------------------------------------------
  !
  !       Calculate change in particle number concentrations due coagulation:

  IF ((nsnucl+nonucl.gt.0).or.lscoag) CALL m7_dnum(kproma, kbdim,  klev,  krow, &
                                                   paerml, paernl, ptp1,  &
                                                   papp1,  pm6rp,  prhop, &
                                                   !>>dod soa: added zttn
                                                   zso4_5, zso4_6, zso4_7, zanew, za4delt, zttn)
                                                   !<<dod  
  !>>dod timers
  IF (ltimer) THEN
     CALL timer_stop(timer_ham_m7_dnum)
     CALL timer_start(timer_ham_m7_avg)
  END IF
  !<<dod

  !--- 5) Recalculation of particle properties under ambient conditions --------
  !
  !--- 5.1) Recalculate mean masses for all modes 
  !         and dry radius and density for the insoluble modes.
 
  CALL m7_averageproperties(kproma, kbdim, klev, krow, paernl, paerml, zttn, pm6rp, prhop) 

  !--- 5.2) Calculate ambient count median radii and density 
  !         for lognormal distribution of particles.
  
  !>>dod timers
  IF (ltimer) THEN
     CALL timer_stop(timer_ham_m7_avg)
     CALL timer_start(timer_ham_m7_hygro)
  END IF
  !<<dod

  IF (nwater == 1) THEN

     CALL m7_kappa(kproma, kbdim, klev, krow, prhp1, paernl, zttn, ptp1, &
                   pm6dry, pm6rp, prhop, pww)

  ELSE
    
     !>>dod soa 
     IF (nsoa == 0) &
     !<<dod
          !         Sulfate particles:
    
          CALL m7_equiz(kproma, kbdim, klev, krow,  &
                        papp1,  zttn,  ptp1,   &
                        prhp1,  pm6rp, pm6dry, &
                        prhop,  pww,   paernl  )
    
     !         Mixed particles with sulfate, black and organic carbon, and dust:
    
     CALL m7_equimix(kproma, kbdim, klev, krow,  &
                     papp1,  zttn,  ptp1,   &
                     prhp1,  pm6rp, pm6dry, &
                     prhop,  pww,   paernl  )
    
     !         Accumulation and coarse mode particles in presence of
     !         sea salt particles:  
     !         dod: deleted zhplus and added krow     
     CALL m7_equil(kproma, kbdim,  klev,   krow,  prhp1, paerml, paernl,  &
                   pm6rp,  pm6dry, pww,   prhop,  ptp1     )
    
  ENDIF
  
  !>>dod timers
  IF (ltimer) THEN
     CALL timer_stop(timer_ham_m7_hygro)
     CALL timer_start(timer_ham_m7_redistr)
  END IF
  !<<dod

  !--- 6) Repartitition particles among the modes ------------------------------
  
  !>>dod merged HAM R2 versions of m7_dconc
  IF (lscond .OR. lscoag .OR. nsoa == 1) THEN
          
     CALL m7_dconc(kproma, kbdim, klev, krow, paerml, paernl, pm6dry)
     
  END IF
  !<<dod
  !>>dod timers
  IF (ltimer) THEN
     CALL timer_stop(timer_ham_m7_redistr)
     CALL timer_start(timer_ham_m7_avg)
  END IF
  !<<dod
  
  !--- 7) Recalculation of particle properties under ambient conditions --------
  !
  !--- 7.1) Calculate mean particle mass for all modes 
  !         and dry radius and density for the insoluble modes:

  CALL m7_averageproperties(kproma, kbdim, klev, krow, paernl, paerml, zttn, pm6rp, prhop)
  
  !--- 7.2) Calculate ambient count median radii and density 
  !         for lognormal distribution of particles.
  !
  
  !>>dod timers
  IF (ltimer) THEN
     CALL timer_stop(timer_ham_m7_avg)
     CALL timer_start(timer_ham_m7_hygro)
  END IF
  !<<dod

  IF (nwater == 1) THEN

     CALL m7_kappa(kproma, kbdim, klev, krow, prhp1, paernl, zttn, ptp1, &
                   pm6dry, pm6rp, prhop, pww)

  ELSE
    
     !         Sulfate particles:
    
     !>>dod soa 
     IF (nsoa == 0) &
     !<<dod
          CALL m7_equiz(kproma, kbdim, klev,  krow, &
                        papp1,  zttn,  ptp1,   &
                        prhp1,  pm6rp, pm6dry, &
                        prhop,  pww,   paernl  )
    
     !         Mixed particles with sulfate, black and organic carbon, and dust: 
    
     CALL m7_equimix(kproma, kbdim, klev, krow,  &
                     papp1,  zttn,  ptp1,   &
                     prhp1,  pm6rp, pm6dry, &
                     prhop,  pww,   paernl  ) 
    
     !         Accumulation and coarse mode particles in presence of
     !         sea salt particles:
    
     CALL m7_equil(kproma, kbdim,  klev,   krow,  prhp1, paerml, paernl, &
                   pm6rp,  pm6dry, pww,   prhop,  ptp1    )
    
  ENDIF
  
  !>>dod timers
  IF (ltimer) CALL timer_stop(timer_ham_m7_hygro)
  !<<dod
     
  !>>dod soa
  pgas(1:kproma,:,isubm_so4g) = pso4g(1:kproma,:)
  !<<dod

END SUBROUTINE m7


SUBROUTINE m7_averageproperties(kproma, kbdim, klev, krow, paernl, paerml, pttn, pm6rp, prhop)
  !
  !  Author:
  !  --------
  !  E. Vignati, JRC/EI (original source)                10/2000
  !  P. Stier, MPI      (f90-version, changes, comments)    2001
  !  D. O'Donnell, MPI-M (rewritten without hardcoding of modes and species) 2008
  !
  !  Purpose:                                                           
  !  ---------                                                           
  !  Calculation of the mean particle mass (pttn).
  !     [molecules cm-3] for the sulphate mass
  !     [ug m-3] for the other compounds
  !
  !  Calculation of the (dry) radius and the density 
  !  of the particles of the insoluble modes.
  !
  !  Interface:
  !  ----------
  !  m7_averageproperties is called from m7

  USE mo_kind,           ONLY: dp
  USE mo_math_constants, ONLY: pi
  USE mo_ham_m7ctl,      ONLY: iso4ns,                                     &
                               critn,  ram2cmr,                            &
                               cdconv, cmin_aernl, cmin_aerml, cminvol
  USE mo_ham,            ONLY: subm_aerounitconv, subm_aero_idx, immr2ug
  USE mo_ham,            ONLY: naerocomp, aerocomp, sizeclass, nclass

  IMPLICIT NONE

  !---subroutine interface
  INTEGER,  INTENT(IN)  :: kproma, kbdim, klev, krow      ! grid parameters
  REAL(dp), INTENT(IN)  :: paerml(kbdim,klev,naerocomp)     ! aerosol mass per species and mode
  REAL(dp), INTENT(IN)  :: paernl(kbdim,klev,nclass)        ! aerosol number per mode
  REAL(dp), INTENT(INOUT) :: pttn(kbdim,klev,naerocomp)     ! mean mass per particle for each compound in each mode 
  !                                                       ! [in molec. for sulphate and in g for others]
  REAL(dp), INTENT(INOUT) :: pm6rp(kbdim,klev,nclass)       ! mode count median radius (wet cmr for soluble
                                                          ! modes, dry cmr for insoluble modes)
  REAL(dp), INTENT(INOUT) :: prhop(kbdim,klev,nclass)       ! mode mean particle density [g cm-3]
  !--- Local data:
  !    parameters
  !>>dod 
  REAL(dp), PARAMETER :: fac_ug = 1.E-12_dp   ! conversion of ug m-3 -> g cm-3
  REAL(dp), PARAMETER :: fac_molec = 1._dp    ! sulphate: no unit conversion
  REAL(dp), PARAMETER :: z4piover3 = 4._dp*pi/3._dp
  !<<dod

  INTEGER :: jclass, jn
  INTEGER :: jspec, jl
  
  REAL(dp) :: zunitfac
  REAL(dp) :: zdens
  REAL(dp) :: zaltern_val !SF #458

  REAL(dp) :: zinsvol(kbdim,klev,nclass)        ! mean volume for single particle per mode
  REAL(dp) :: zinsmas(kbdim,klev,nclass)        ! mean mass for single particle per mode

  REAL(dp) :: ztmp1(kbdim,klev) !SF #458 temporary variable

  LOGICAL :: ll1(kbdim,klev) !SF #458 temporary logical array

  !---executable procedure

  zinsmas(:,:,:) = 0._dp
  zinsvol(:,:,:) = 0._dp

  !--- 1) Calculate mean particle masses [g]:
  !
  !       To be able to compute a intra-modal coagulation coefficient for the nucleation
  !       mode for the case of no pre-existing particles but coagulation of freshly formed
  !       particles during the timestep, pttn is set to the mass of the critical cluster
  !       for this case. This allows to calculate an ambient radius of the 
  !       freshly formed particles and subsequently the calculation of the coagulation 
  !       coefficient. This mass is "virtual" as it is not added to the mode but used 
  !       only for the described computation of the coagulation coefficient. 
  !       !@@@ Check whether this is always fulfilled. 

  DO jn=1,naerocomp
     jspec = aerocomp(jn)%spid
     jl    = subm_aero_idx(jspec)
     IF (subm_aerounitconv(jl) == immr2ug) THEN
        zunitfac = fac_ug
     ELSE
        zunitfac = fac_molec
     END IF
     
     jclass = aerocomp(jn)%iclass

     !>>SF #458 (replacing WHERE statements)
     zaltern_val = MERGE(critn, 0._dp, (jn == iso4ns))

     ll1(1:kproma,:) = (paernl(1:kproma,:,jclass) > cmin_aernl) &
                 .AND. (paerml(1:kproma,:,jn)     > cmin_aerml)

     ztmp1(1:kproma,:) = MERGE(paernl(1:kproma,:,jclass), 1._dp, ll1(1:kproma,:)) !SF 1. is a dummy val.

     pttn(1:kproma,:,jn) = MERGE( &
                                zunitfac*paerml(1:kproma,:,jn)/ztmp1(1:kproma,:), &
                                zaltern_val, &
                                ll1(1:kproma,:))
     !<<SF #458 (replacing WHERE statements)

     !---in case of fp underflow
     pttn(1:kproma,:,jn) = MAX(pttn(1:kproma,:,jn),0._dp)

     !---for calculation of count median radii and density for insoluble modes:
     !   calculation of mode mass and volume
     IF (.NOT. sizeclass(jclass)%lsoluble) THEN
        zdens = cdconv*aerocomp(jn)%species%density
        zinsmas(1:kproma,:,jclass) = zinsmas(1:kproma,:,jclass) + pttn(1:kproma,:,jn)
        zinsvol(1:kproma,:,jclass) = zinsvol(1:kproma,:,jclass) + pttn(1:kproma,:,jn)/zdens
     END IF

  END DO

  !---calculate count median radii and density for insoluble modes:
  
  DO jclass=1,nclass
     IF (.NOT. sizeclass(jclass)%lsoluble) THEN

        !>>SF #458 (replacing WHERE statements)
        ll1(1:kproma,:) = (paernl(1:kproma,:,jclass)  > cmin_aernl) &
                    .AND. (zinsvol(1:kproma,:,jclass) > cminvol)

        ztmp1(1:kproma,:) = MERGE(zinsvol(1:kproma,:,jclass), 1._dp, ll1(1:kproma,:)) !SF 1. is a dummy val.

        prhop(1:kproma,:,jclass) = MERGE( &
                                        zinsmas(1:kproma,:,jclass) / ztmp1(1:kproma,:), &
                                        prhop(1:kproma,:,jclass), &
                                        ll1(1:kproma,:))

        ztmp1(1:kproma,:) = ram2cmr(jclass)*((ztmp1(1:kproma,:)/z4piover3)**(1._dp/3._dp))

        pm6rp(1:kproma,:,jclass) = MERGE(ztmp1(1:kproma,:), pm6rp(1:kproma,:,jclass), ll1(1:kproma,:))
        !<<SF #458 (replacing WHERE statements)

     END IF
  END DO

END SUBROUTINE m7_averageproperties


SUBROUTINE m7_kappa(kproma, kbdim, klev, krow, prelhum, paernl, pttn, ptp1, &
                    pm6dry, pm6rp, prhop, pww)

  ! m7_kappa implements the parameterisation of the hygroscopic
  ! growth of aerosols, as measured by the growth factor (gf), as a
  ! function of  the hygroscopicity parameter kappa and 
  ! the ambient temperature and relative humidity. The relation is:
  !            [    A   ]      (gf**3 -1)
  !      RH*exp[- ----- ] = ----------------
  !            [  Rd*gf ]   gf**3 - (1-kappa)
  !
  ! source:
  ! A single parameter representation of hygroscopic growth and cloud
  ! condensation nucleus activity
  ! M.D. Petters and S.M. Kreidenweis
  ! ACP 7, 1961-1971, 2007
  !
  ! Equation (11) from the Petters & Kreidenweis paper, given above, expresses
  ! the hygroscopic growth factor (gf) as a function of aerosol dry radius (Rd),
  ! temperature (T), relative humidity (rh) and a substance property denoted kappa
  ! that encapsulates hygroscopic properties of that substance. This transcendental
  ! equation is solved offline for various Rd, T, rh and kappa of atmospheric
  ! relevance and the results stored in a lookup table. 
  ! m7_kappa uses the ambient T and rh, uses the mode count median 
  ! dry radius and the volume-weighted average kappa to find the entry point into
  ! the GF lookup table. Using the growth factor the wet median radius and water
  ! uptake of the mode are calculated.
  !   
  ! Author:
  ! -------
  ! Declan O'Donnell, MPI-M, 2008

  ! Changes for release HAM 2.3: 
  ! ---------------------------
  ! 1. Bugfix (redmine issue #357) 
  ! 2. Changes to lookup table to improve RH lookup and extend to higher RH (redmine issue #260)
  !    SEE mo_ham_kappa for more details (better still, read the redmine issue)

  USE mo_kind,            ONLY: dp
  USE mo_physical_constants, ONLY: avo
  USE mo_math_constants,     ONLY: pi
  USE mo_species,         ONLY: speclist
  USE mo_ham_m7ctl,       ONLY: cmr2ram, ram2cmr,   &
                                cdconv, cminvol, cminrad
  USE mo_ham,             ONLY: subm_aerounitconv, subm_aero_idx, immr2molec, &
                                naerocomp, aerocomp, sizeclass, nclass, nsol

  USE mo_ham_kappa  

  IMPLICIT NONE
  
  !--- subroutine interface
  INTEGER, INTENT(IN) :: kproma, kbdim, klev, krow       ! grid parameters

  REAL(dp), INTENT(IN) :: prelhum(kbdim,klev)            ! relative humidity
  REAL(dp), INTENT(IN) :: ptp1(kbdim,klev)               ! temperature at t+dt
  REAL(dp), INTENT(IN) :: paernl(kbdim,klev,nclass)        ! aerosol number density per mode [ cm-3 ]
  REAL(dp), INTENT(IN) :: pttn(kbdim,klev,naerocomp)     ! aerosol mean particle mass per species and mode
  REAL(dp), INTENT(INOUT) :: pm6dry(kbdim,klev,nsol)     ! dry aerosol count median radius [cm]
  REAL(dp), INTENT(INOUT) :: pm6rp(kbdim,klev,nclass)      ! aerosol count median radius [cm]
                                                         ! (wet radius for soluble modes)
  REAL(dp), INTENT(INOUT) :: prhop(kbdim,klev,nclass)      ! aerosol density [g cm-3]
  REAL(dp), INTENT(OUT) :: pww(kbdim,klev,nclass)          ! aerosol water [kg m-3]

  !--- local data
  REAL(dp), PARAMETER :: z4piover3 = 4._dp*pi/3._dp
  REAL(dp), PARAMETER :: cepssec   = 1.E-25_dp           ! used to avoid divide by zero

  !--- local variables
  REAL(dp) :: zkappa(kbdim,klev,nclass)                    ! hygroscopicity factor
  REAL(dp) :: zdrymass(kbdim,klev,nclass)                  ! volume of dry  material
  REAL(dp) :: zdryvol(kbdim,klev,nclass)                   ! volume of total dry material
  REAL(dp) :: zm6dry_safe(kbdim,klev)                      ! tmp duplicate of pm6rdry where pm6rdry > cminrad
  REAL(dp) :: zww_tmp(kbdim,klev)                          ! tmp duplicate of pww
  REAL(dp) :: zm6rp_tmp(kbdim,klev)                        ! tmp duplicate of pm6rp
  REAL(dp) :: zrhop_tmp(kbdim,klev)                        ! tmp duplicate of prhop
  REAL(dp) :: zgf                                          ! growth factor
  REAL(dp) :: zunitfac, zmass, zdensity, zvolume
  REAL(dp) :: zwetvol, zdvol, zwatervol, zwatermass
  REAL(dp) :: zTrange, zRHrange, zKrange, zRDrange
  REAL(dp) :: ztsteps, zrhsteps, zksteps, zrdsteps
  REAL(dp) :: zt1, zrh1, zk1, zr1, zr2, zr3

  LOGICAL  :: ll1(kbdim,klev) 

  INTEGER  :: ix_rd, ix_k                                ! index to radius and kappa in lookup table
  INTEGER  :: ix_t(kbdim,klev)                           ! index to temperature in lookup table
  INTEGER  :: ix_rh(kbdim,klev)                          ! index to RH in lookup table
  INTEGER  :: jk, jl, jn, jclass, jspec                  ! loop counters

  !--- executable procedure

  !---initialisations
  zdrymass(:,:,:) = 0._dp
  zdryvol(:,:,:) = 0._dp
  zkappa(:,:,:) = 0._dp

  !---range of possible values of temperature, RH, kappa and radius in the lookup table:
  zTrange = T_max-T_min
  zRHrange = rh_max-rh_min
  zKrange = kappa_max-kappa_min
  zRDrange = ln_Rd_max - ln_Rd_min

  ztsteps = REAL((N_t-1),dp)                   ! number of steps on the temp axis 
  zrhsteps = REAL((N_rh-1),dp)                 ! number of steps on the RH axis 
  zksteps = REAL((N_kappa-1),dp)               ! number of steps on the kappa axis  
  zrdsteps = REAL((N_Rd-1),dp)                 ! number of steps on the radius axis

  !---since temperature and RH are independent of mode and species we can
  !   calculate the coordinate value for the lookup table just once
  !   rather than per mode, so we do this for the first loop iteration only
  DO jk=1,klev
     DO jl=1,kproma
        
        zt1 = MAX(ptp1(jl,jk), T_min)               ! limit to at least T_min
        zt1 = MIN(zt1, T_max)                       ! limit to at most T_max
        ix_t(jl,jk) = 1 + NINT(ztsteps*(zt1-T_min)/zTrange) ! linear interpolation from actual
                                                           ! temperature to lookup table entry
  !>>dod redmine #260     
     END DO
  END DO

  IF (lut_kappa_version == 1) THEN
     DO jk=1,klev
        DO jl=1,kproma
           zrh1 = MAX(prelhum(jl,jk), rh_min)       ! limit to at least RH_min
           zrh1 = MIN(zrh1, rh_max)                 ! limit to at most T_max
           ix_rh(jl,jk) = 1 + NINT(zrhsteps*(zrh1-rh_min)/zRHrange) ! linear interpolation from actual
                                                                    ! RH to lookup table entry
        END DO
     END DO
  ELSE                                              ! lookup table version 2
     DO jk=1,klev
        DO jl=1,kproma
           ix_rh(jl,jk) = CEILING( &                ! see redmine #260
                          (LOG(1._dp-(1._dp-rh_step_compress)*(prelhum(jl,jk)-rh_min)*inv_rh_init_step)) &
                          / log_rh_step_compress )
           ix_rh(jl,jk) = MAX(ix_rh(jl,jk),1)
           ix_rh(jl,jk) = MIN(ix_rh(jl,jk),N_rh)
        END DO
     END DO
  END IF
  !<<dod
  

  DO jn=1,naerocomp

     jclass  = aerocomp(jn)%iclass
     jspec = aerocomp(jn)%spid
     jl    = subm_aero_idx(jspec)

     IF (sizeclass(jclass)%lsoluble) THEN
       
        IF (subm_aerounitconv(jl) == immr2molec) THEN
           zunitfac = speclist(jspec)%moleweight / avo
        ELSE
           zunitfac = 1._dp
        END IF

        zdensity = cdconv * speclist(jspec)%density 

        DO jk=1,klev
           DO jl=1,kproma

              zmass = zunitfac * pttn(jl,jk,jn)
              zvolume = zmass / zdensity
              zdrymass(jl,jk,jclass) = zdrymass(jl,jk,jclass) + zmass

              zdryvol(jl,jk,jclass) = zdryvol(jl,jk,jclass) + zvolume
       
           END DO
        END DO

     END IF

  END DO

  DO jclass = 1,nclass
     IF (sizeclass(jclass)%lsoluble) THEN
        pm6dry(1:kproma,:,jclass) = ram2cmr(jclass)*(zdryvol(1:kproma,:,jclass)/z4piover3)**(1._dp/3._dp)
     END IF
  END DO

  DO jn=1,naerocomp

     jclass = aerocomp(jn)%iclass
     jspec  = aerocomp(jn)%spid
     jl     = subm_aero_idx(jspec)
     
     IF (sizeclass(jclass)%lsoluble) THEN
       
        IF (subm_aerounitconv(jl) == immr2molec) THEN
           zunitfac = speclist(jspec)%moleweight / avo
        ELSE
           zunitfac = 1._dp
        END IF

        zdensity = cdconv * speclist(jspec)%density 

        DO jk=1,klev
           DO jl=1,kproma
              zkappa(jl,jk,jclass) = zkappa(jl,jk,jclass) + &
                                     zunitfac * speclist(jspec)%kappa * &
                                     (pttn(jl,jk,jn)/zdensity) / (zdryvol(jl,jk,jclass) + cepssec) 
           END DO
        END DO
     END IF

  END DO

  DO jclass = 1,nclass
     IF (sizeclass(jclass)%lsoluble) THEN

        ll1(1:kproma,:) = ( pm6dry(1:kproma,:,jclass) > cminrad )
        zm6dry_safe(1:kproma,:) = MERGE(pm6dry(1:kproma,:,jclass), cminrad, ll1(1:kproma,:))
        !SFNote: cminrad is dummy in the above line. This is merely a protection against LOG(0.) later on

        DO jk = 1,klev
           DO jl = 1,kproma

              zk1 = MAX (zkappa(jl,jk,jclass), kappa_min)       ! limit to at least kappa_min
                 
              zk1 = MIN(zk1, kappa_max)        ! limit to at most kappa_max
              ix_k = 1 + NINT(zksteps*(zk1-kappa_min)/zKrange) ! linear interpolation from actual
                                                                 ! kappa to lookup table entry

              zr1 = 1.E-2*zm6dry_safe(jl,jk)                  ! LOG of radius in m.
              zr2 = LOG(zr1)         
               
              zr3 = MAX(zr2, ln_Rd_min)                       ! limit to at least Rd_min
              zr3 = MIN(zr3, ln_Rd_max)                       ! limit to at most Rd_max 
              ix_rd = 1 + NINT(zrdsteps*(zr3-ln_Rd_min) / zRDrange )   ! to lookup table entry
               
              !---get the growth factor
              zgf = gf(ix_rd, ix_t(jl,jk), ix_rh(jl,jk), ix_k)

              !---wet radius = dryradius * growth factor
              zm6rp_tmp(jl,jk) = zm6dry_safe(jl,jk) * zgf
               
              !---calculate water uptake and wet count median radius
              !   We assume that the volume of dissolved material is the same as the volume of 
              !   the dry material.
             
              zdvol = z4piover3*(cmr2ram(jclass)*zm6dry_safe(jl,jk))**3
              zwetvol = z4piover3*(cmr2ram(jclass)*zm6rp_tmp(jl,jk))**3

              zwatervol = zwetvol - zdvol

              zwatermass = 1._dp * zwatervol ! water in g per particle

              !---aerosol water in kg m-3
              zww_tmp(jl,jk) = 1.E3_dp*paernl(jl,jk,jclass)* zwatermass 

              zrhop_tmp(jl,jk) = (zdrymass(jl,jk,jclass) + zwatermass) / zwetvol

           END DO
        END DO

        pww(1:kproma,:,jclass)   = MERGE(zww_tmp(1:kproma,:),   0._dp,                     ll1(1:kproma,:))
        pm6rp(1:kproma,:,jclass) = MERGE(zm6rp_tmp(1:kproma,:), pm6dry(1:kproma,:,jclass), ll1(1:kproma,:))
        prhop(1:kproma,:,jclass) = MERGE(zrhop_tmp(1:kproma,:), 1._dp,                     ll1(1:kproma,:))

     END IF
  END DO

END SUBROUTINE m7_kappa


SUBROUTINE m7_equiz(kproma,  kbdim, klev, krow,  &
                    papp1,   pttn,  ptp1,   &
                    prelhum, pm6rp, pm6dry, &
                    prhop,   pww,   paernl  )
  !
  !   *m7_equiz*   calculates the ambient radii of the sulphate particles
  !
  !    Authors:
  !    --------
  !    J. Wilson, E. Vignati, JRC/EI (original source)                05/2000
  !    P. Stier, MPI                 (f90-version, changes, comments)    2001
  !
  !    Purpose:
  !    --------
  !    This routine calculates the ambient radii for sulfate particles
  !    with mass of ttn molecules, converts them to count mean radii and
  !    stores them in the array with address pm6rp.
  !    It additionally calculates the ambient particle density.
  !
  !    Method:
  !    -------
  !    The calculations of the ambient particle properties are based on 
  !    parameterisations of the mass of sulfate and density derived 
  !    by Julian Wilson from a regression analysis of results of solving 
  !    the generalised Kelvin equation using (F. J. Zeleznik, J. Phys. Chem.
  !    Ref. Data 20, 1157, 1991), for an H2SO4-H2O mixture, in the 
  !    following parameter ranges: 

  !       1e2   < pttn    < 1E11  [molecules]
  !       0.2   < prelhum < 0.9   [1]
  !       240   < ptp1    < 330   [K]
  !       10000 < papp1   < 100000  [Pa]
  !
  !    Due to the limitations of the parametrisation, the ambient temperature
  !    is restricted to a minimum of 240 K within this subroutine. 
  !      
  !    Interface:
  !    ----------
  !    *m7_equiz* is called from *m7* and *m7_dconc*
  !
  !    Externals:
  !    ----------
  !    none
  !
  !
  USE mo_kind,          ONLY: dp
  USE mo_physical_constants, ONLY: avo
  USE mo_math_constants,     ONLY: pi
  USE mo_ham_m7ctl,     ONLY: dh2so4, dh2o,                       &
                              wvb, gmb, wh2so4, ram2cmr, cmin_aerml
  USE mo_species,       ONLY: speclist
  USE mo_ham_species,   ONLY: id_so4
  USE mo_ham,           ONLY: naerocomp, aerocomp, sizeclass, nclass, nsol
  !
  ! pttn      = average mass for single compound in each mode 
  !             [in molec. for sulphate and in ug for bc, oc, ss, and dust]
  ! pm6rp     = count mean radius under ambient conditions [cm]
  ! pm6dry    = count mean radius under dry conditions [cm]
  ! paernl    = aerosol number for each mode [cm-3]
  ! pww       = aerosol water content for each mode [kg(water) m-3(air)]
  ! zwso4     = percentage by mass of sulfate in a H2O-H2SO4 particle 
  !             containing pttn molecules of sulfate under ambient conditions
  ! zvso4     = volume of pttn molecules of sulfate [cm3]
  ! zmso4     = mass of pttn molecules of sulfate [g]
  ! zdso4h2o  = density of sulfate-h2o fraction of a particle with average 
  !             mass [g.cm-3]
  ! zmso4h2o  = mass of sulfate-h2o fraction of a particle with average mass [g]
  ! zvso4h2o  = volume of sulfate-h2o fraction of a particle with average 
  !             mass [cm3]

  IMPLICIT NONE

  INTEGER :: kproma, kbdim, klev, krow

  REAL(dp):: papp1(kbdim,klev),        ptp1(kbdim,klev),       &
             prelhum(kbdim,klev)

  REAL(dp):: pttn(kbdim,klev,naerocomp), prhop(kbdim,klev,nclass), &
             pm6dry(kbdim,klev,nsol),  pm6rp(kbdim,klev,nclass), &
             pww(kbdim,klev,nclass),     paernl(kbdim,klev,nclass)

  !--- Local variables:

  INTEGER :: jk, jl, jclass, jn

  REAL(dp) :: zaerelse(kbdim,klev,nsol)
                                         
  REAL(dp):: zwso4,       zvso4,       zmso4,                  &
             zvso4h2o,    zmso4h2o,    zdso4h2o,               &
             zapp1,       ztk,         zrh

  REAL(dp):: ztk2,        zln3,        zln32,                  &
             zlnm,        zss2,        zlnm2

  INTEGER :: ixso4
  !------------------------------------------------------------------------------

  zaerelse(:,:,:) = 0._dp

  !--- 1) Determine mass of non sulfate compounds in a particle: ---------

  DO jn = 1,naerocomp                               ! for each mode+species
     jclass = aerocomp(jn)%iclass                    ! check if the species is NOT sulphate 
                                                  ! if so, add to the total mode non-S mass
     IF (aerocomp(jn)%spid /= id_so4 .AND. &
         sizeclass(jclass)%lsoluble) THEN        
        zaerelse(1:kproma,:,jclass) = zaerelse(1:kproma,:,jclass) + pttn(1:kproma,:,jn)
     END IF

  END DO

  DO 100 jclass=1,nsol
     ixso4 = speclist(id_so4)%iaerocomp(jclass)    !!mgs!! im7table(jclass, id_so4)
     DO 90 jk=1,klev
        DO 80 jl=1,kproma

           !       other compounds than sulfate:

           IF (pttn(jl,jk,ixso4) > 0.0_dp .AND. zaerelse(jl,jk,jclass) < cmin_aerml) THEN
              !
              !--- 2.1) Calculation of the ambient particle properties: -----------
              !      
              !--- Constrain ambient temperature to conditions for which the 
              !    parametrisation of the liquid water content works:
              
              ! Temperature:
              ztk = ptp1(jl,jk)
              ztk = MAX(ztk , 240._dp)

              ! Relative Humidity:
              zrh = prelhum(jl,jk)
              zrh = MAX(zrh , 0.05_dp)
              zrh = MIN(zrh , 0.90_dp)

              !--- Assign auxiliary variables:

              zapp1=papp1(jl,jk)
              zlnm = LOG(pttn(jl,jk,ixso4))
              zlnm2 = zlnm*zlnm
              zss2 = zrh**2
              ztk2 = ztk*ztk
              zln3 = zlnm/3.0_dp
              zln32 = zln3*zln3
              !
              !--- Percentage by weight of sulfate in the particle [%]:
              !    (Here we ignore any insoluble mass.)
              !
              zwso4 = wvb(1) + wvb(2)*zlnm + wvb(3)*zrh*zlnm + wvb(4)*ztk*zlnm + &
                      wvb(5)*zrh/ztk + wvb(6)*zlnm2*zrh + wvb(7)*zlnm2*ztk +     &
                      wvb(8)*zlnm*zss2 + wvb(9)*zlnm*ztk2 + wvb(10)*zlnm2*zss2 + &
                      wvb(11)*zlnm2*ztk2 + wvb(12)*zss2/ztk2 + wvb(13)*zlnm2 +   &
                      wvb(14)*zlnm2*zlnm + wvb(15)*zlnm2*zlnm2 +                 &
                      wvb(16)*zss2*zrh/(ztk2*ztk) + wvb(17)*LOG(zrh*ztk/zapp1)

              !--- Dry mass of sulfate in an average particle [g]:

              zmso4 = pttn(jl,jk,ixso4)*wh2so4/avo

              !--- Dry volume of sulfate in an average particle [cm3]:
              !    Any temperature or pressure dependency of the 
              !    sulfate density is ingored.

              zvso4 = zmso4/dh2so4

              !--- Mass of sulfate + water in an average particle [g]:

              zmso4h2o = zmso4/(zwso4/100.0_dp)

              !--- Density of the sulfate-water fraction of an average particle [g cm-3]:
              !@@@ Check: changed zwvso4 into zwso4 (now the mass!)

              zdso4h2o = gmb(1) + gmb(2)*zwso4 + gmb(3)*zln3 + gmb(4)*zrh +        &
                         gmb(5)*ztk + gmb(6)*zln32 + gmb(7)*zln3/zrh +             &
                         gmb(8)*zln3/ztk + gmb(9)*ztk2
         
              !--- Limits for zdso4h2o: H2O(0.99) and pure H2SO4 (1.841):
              !
              zdso4h2o=MAX(zdso4h2o,dh2o)
              zdso4h2o=MIN(zdso4h2o,dh2so4)
              
              !--- Volume of sulfate-water fraction of an average particle [cm3]:

              zvso4h2o = zmso4h2o/zdso4h2o

              !--- 2.2) Calculatiion of the particle radii: ----------------------------

              !--- 2.2.1) Dry count mean radius [cm]:
 
              pm6dry(jl,jk,jclass)=((zvso4)*0.75_dp/pi)**(1._dp/3._dp)*ram2cmr(jclass)

              !--- 2.2.2) Equilibrium wet count mean radius [cm]:
             
              pm6rp(jl,jk,jclass) =((zvso4h2o)*0.75_dp/pi)**(1._dp/3._dp)*ram2cmr(jclass)
              
              !--- 2.3) Assignment of the particle density [g cm-3]: -------------------
 
              prhop(jl,jk,jclass)=zdso4h2o

              !--- 2.4) Store aerosol water for each mode [kg(water) m-3(air)]:

              pww(jl,jk,jclass)=(zmso4h2o-zmso4)*paernl(jl,jk,jclass)*1.E3_dp

           END IF
80      END DO
90   END DO
100 END DO

END SUBROUTINE m7_equiz


SUBROUTINE m7_equimix(kproma,  kbdim, klev, krow,  &
                      papp1,   pttn,  ptp1,   &
                      prelhum, pm6rp, pm6dry, &
                      prhop,   pww,   paernl  )
  !
  !    *m7_equimix*   calculates the ambient radii of the particles with
  !                   sulphate, b/o carbon and dust.
  !
  !    Authors:
  !    --------
  !    J. Wilson and E. Vignati, JRC (original source)               05/2000
  !    P. Stier, MPI                 (f90-version, changes, comments)   2001
  !
  !    Purpose:
  !    --------
  !    This routine calculates the ambient radii for mixed particles without
  !    sea salt with mass of ttn molecules, converts them to count mean radii and
  !    stores them in the array with address pm6rp.
  !    It additionally calculates the ambient particle density.
  !
  !    Method:
  !    -------
  !    The calculations of the ambient particle properties are based on 
  !    parameterisations of the mass of sulfate and density derived 
  !    by Julian Wilson from a regression analysis of results of solving 
  !    the generalised Kelvin equation using (F. J. Zeleznik, J. Phys. Chem.
  !    Ref. Data 20, 1157, 1991), for an H2SO4-H2O mixture, in the 
  !    following parameter ranges: 

  !       1e2   < pttn    < 1E11  [molecules]
  !       0.2   < prelhum < 0.9   [1]
  !       240   < ptp1    < 330   [K]
  !       10000 < papp1   < 100000  [Pa]
  !
  !    Due to the limitations of the parametrisation, the ambient temperature
  !    is restricted to a minimum of 240 K within this subroutine. 
  !
  !    For this application to mixed aerosols with an insoluble core we
  !    assume the H2O uptake by the particle to be that of a pure 
  !    H2SO4 / H2O particle containing the H2SO4 mass of the mixed aerosol.
  !
  !    Interface:
  !    ----------
  !    *m7_equimix* is called from *m7*
  !
  !    Externals:
  !    ----------
  !    none
  !
  USE mo_kind,         ONLY: dp
  USE mo_physical_constants, ONLY: avo
  USE mo_math_constants,     ONLY: pi
  USE mo_ham_m7ctl,    ONLY: dh2so4, dh2o, wvb, gmb, wh2so4, &
                             ram2cmr, cmin_aerml, cdconv 
  USE mo_species,      ONLY: speclist
  USE mo_ham_species,  ONLY: id_ss, id_so4
  USE mo_ham,          ONLY: naerocomp, aerocomp, sizeclass, nclass, nsol

  IMPLICIT NONE
  !
  ! pttn      = average mass for single compound in each mode 
  !             [in molec. for sulphate and in ug for bc, oc, ss, and dust]
  ! pm6rp     = count mean radius under ambient conditions
  ! pm6dry    = count mean radius under dry conditions
  ! paernl    = aerosol number for each mode [cm-3]
  ! pww       = aerosol water content for each mode [kg(water) m-3(air)]
  ! zwso4     = percentage by mass of sulfate in a H2O-H2SO4 particle 
  !             containing pttn molecules of sulfate under ambient conditions
  ! zvso4     = volume of pttn molecules of sulfate [cm3]
  ! zmso4     = mass of pttn molecules of sulfate [g]
  ! zdso4h2o  = density of sulfate-h2o fraction of a particle with average 
  !             mass [g.cm-3]
  ! zmso4h2o  = mass of sulfate-h2o fraction of a particle with average mass [g]
  ! zvso4h2o  = volume of sulfate-h2o fraction of a particle with average 
  !             mass [cm3]
  ! zinsvol   = total volume of insoluble compounds in a single particle of 
  !             average mass [cm3]
  ! zinsmass  = total mass of insoluble compounds in a single 
  !             particle of average mass [cm3]

  INTEGER :: kproma, kbdim, klev, krow 

  REAL(dp):: papp1(kbdim,klev),        ptp1(kbdim,klev),       &
             prelhum(kbdim,klev)

  REAL(dp):: pttn(kbdim,klev,naerocomp), prhop(kbdim,klev,nclass), &
             pm6dry(kbdim,klev,nsol),  pm6rp(kbdim,klev,nclass), &
             pww(kbdim,klev,nclass),   paernl(kbdim,klev,nclass)
  !
  ! Local variables:
  !

  INTEGER :: jk, jl, jclass
  
  REAL(dp):: zseasalt,                                         &
             zwso4,       zvso4,       zmso4,                  &
             zvso4h2o,    zmso4h2o,    zdso4h2o,               &
             zapp1,       zrh,         ztk
                    
  REAL(dp):: zlnm2,       zln3,        zln32,      ztk2,       &
             zlnm,        zss2       

  REAL(dp) :: zinsmas(kbdim,klev,nsol)
  REAL(dp) :: zinsvol(kbdim,klev,nsol)
  REAL(dp) :: zdens
  INTEGER :: ixso4, ixss, jspec, jn
  !------------------------------------------------------------------------------
  !
  zinsmas(:,:,:) = 0._dp
  zinsvol(:,:,:) = 0._dp

  !--- 1) Determine mass of non sulfate compounds in a particle: ---------

  DO jn = 1,naerocomp                      ! for each mode+species
     jclass = aerocomp(jn)%iclass           ! check if the species is neither sulphate nor seasalt
     jspec = aerocomp(jn)%spid           ! if so, add to the total mode 'others' mass
     IF (sizeclass(jclass)%lsoluble .AND. (.NOT. speclist(jspec)%lelectrolyte)) THEN        
        zdens = cdconv*speclist(jspec)%density
        zinsmas(1:kproma,:,jclass) = zinsmas(1:kproma,:,jclass) + pttn(1:kproma,:,jn)
        zinsvol(1:kproma,:,jclass) = zinsvol(1:kproma,:,jclass) + pttn(1:kproma,:,jn)/zdens
        zinsmas(1:kproma,:,jclass) = MAX(zinsmas(1:kproma,:,jclass),0._dp)
        zinsvol(1:kproma,:,jclass) = MAX(zinsvol(1:kproma,:,jclass),0._dp)
     END IF

  END DO

!CDIR unroll=5
  DO 100 jclass=1,nsol
     ixso4 = speclist(id_so4)%iaerocomp(jclass)     !!mgs!! im7table(jclass,id_so4)
     ixss  = speclist(id_ss)%iaerocomp(jclass)      !!mgs!! im7table(jclass, id_ss)

     DO 90 jk=1,klev
        DO 80 jl=1,kproma

           !--- 1) Split particle quantities into soluble (sea salt) and --------------
           !       non soluble (organic carbon + black carbon + dust) parts:
           !       (N.B. densities are assumed independent of temperature & pressure)

           IF (ixss > 0) THEN 
              zseasalt = pttn(jl,jk,ixss)
           ELSE
              zseasalt = 0._dp
           END IF

           !--- 2) Calculation of the particle properties in the absense of sea salt: --
           !  
             IF (pttn(jl,jk,ixso4) > 0.0_dp .AND. zseasalt < cmin_aerml ) THEN

              !--- 2.1) Calculation of the ambient particle properties: ----------------
              !         
              !--- Constrain ambient temperature and relative humidity to 
              !    conditions for which the parametrisation of the liquid 
              !    water content works:
              !
              ! Temperature:
              ztk = ptp1(jl,jk)
              ztk = MAX(ztk , 240._dp)

              ! Relative Humidity:
              zrh = prelhum(jl,jk)
              zrh = MAX(zrh , 0.05_dp)
              zrh = MIN(zrh , 0.90_dp)
              
              !--- Assign auxiliary variables:

              zapp1=papp1(jl,jk)
              zlnm = LOG(pttn(jl,jk,ixso4))
              zlnm2 = zlnm*zlnm
              zss2 = zrh**2
              ztk2 = ztk*ztk
              zln3 = zlnm/3.0_dp
              zln32 = zln3*zln3

              !--- Percentage by weight of sulfate in the particle [%]:
              !    (Here we ignore any insoluble mass.)

              zwso4 =  wvb(1) + wvb(2)*zlnm + wvb(3)*zrh*zlnm + wvb(4)*ztk*zlnm +  &
                       wvb(5)*zrh/ztk + wvb(6)*zlnm2*zrh + wvb(7)*zlnm2*ztk +      &
                       wvb(8)*zlnm*zss2 + wvb(9)*zlnm*ztk2 + wvb(10)*zlnm2*zss2 +  &
                       wvb(11)*zlnm2*ztk2 + wvb(12)*zss2/ztk2 + wvb(13)*zlnm2 +    &
                       wvb(14)*zlnm2*zlnm + wvb(15)*zlnm2*zlnm2 +                  &
                       wvb(16)*zss2*zrh/(ztk2*ztk) + wvb(17)*LOG(zrh*ztk/zapp1)

              !--- Dry mass of sulfate in an average particle [g]:

              zmso4 = pttn(jl,jk,ixso4)*wh2so4/avo

              !--- Dry volume of sulfate in an average particle [cm3]:
              !    Any temperature or pressure dependency of the 
              !    sulfate density is ingored.

              zvso4 = zmso4/dh2so4

              !--- Mass of sulfate + water in an average particle [g]:

              zmso4h2o = zmso4/(zwso4/100.0_dp)

              !--- Density of the sulfate-water fraction of an average particle [g cm-3]:
              !@@@ Check: changed zwvso4 into zwso4 (now the mass!)

              zdso4h2o = gmb(1) + gmb(2)*zwso4 + gmb(3)*zln3 + gmb(4)*zrh +        &
                         gmb(5)*ztk + gmb(6)*zln32 + gmb(7)*zln3/zrh +             &
                         gmb(8)*zln3/ztk + gmb(9)*ztk2
         
              !--- Limits for zdso4h2o: H2O(0.99) and pure H2SO4 (1.841):
              !
              zdso4h2o=MAX(zdso4h2o,dh2o)
              zdso4h2o=MIN(zdso4h2o,dh2so4)
              
              !--- Volume of sulfate-water fraction of an average particle [cm3]:

              zvso4h2o = zmso4h2o/zdso4h2o

              !--- 2.2) Calculatiion of the particle radii: ----------------------------

              !--- 2.2.1) Dry count mean radius [cm]:
 
              pm6dry(jl,jk,jclass)=((zvso4+zinsvol(jl,jk,jclass))*0.75_dp/pi)**(1._dp/3._dp)*ram2cmr(jclass)

              !--- 2.2.2) Equilibrium wet count mean radius [cm]:
             
              pm6rp(jl,jk,jclass) =((zvso4h2o+zinsvol(jl,jk,jclass))*0.75_dp/pi)**(1._dp/3._dp)*ram2cmr(jclass)

              !--- 2.3) Calculation of the particle density [g cm-3]:-------------------
 
              prhop(jl,jk,jclass)=(zmso4h2o+zinsmas(jl,jk,jclass))/                  &
                                (zvso4h2o+zinsvol(jl,jk,jclass))                 

              !--- 2.4) Store aerosol water for each mode [kg(water) m-3(air)]:

              pww(jl,jk,jclass)=(zmso4h2o-zmso4)*paernl(jl,jk,jclass)*1.E3_dp

           END IF

80      END DO
90   END DO
100 END DO

END SUBROUTINE m7_equimix


SUBROUTINE m7_equil (kproma, kbdim,  klev,   krow, prelhum, paerml, paernl, &
                     pm6rp,  pm6dry, pww,    prhop,  ptp1    )
  !
  !**** *m7_equil* calculates the ambient radii of accumulation 
  !                and coarse mode particles in presence of sea salt
  !     
  !
  !  Authors:
  !  --------
  !  E. Vignati, JRC/EI (original source)                    01/2000
  !  P. Stier, MPI      (f90-version, changes, comments)        2001
  !  D. O'Donnell, MPI-M, see Revision information below        2007-2008

  !  Purpose:
  !  --------
  !  This routine calculates the equilibrium radius of sea salt
  !  particles, and of sea salt particles mixed with sulphate for 
  !  accumulation and coarse modes. 
  !
  !  Interface:
  !  ----------
  !  *m7_equil* is called from *m7*
  !
  !  References:
  !  -----------
  !  Jacobson, M.Z., Tabazadeh, A., Turco, R.P., (1996). Simulating
  !     equilibrium within aerosols and nonequilibrium between gases
  !     and aerosols. 
  !  Tang, I.N., (1997). Thermodynamic and optical properties of 
  !     mixed-salt aerosols of atmospheric importance. 
  
  !  Revision information
  !  --------------------
  !  Based on: m7_equil from HAM_me_v1 (svn rev:266) 
  !  Major changes: changed to become more generic with respect to non-electrolytic species
  !                 i.e. should not require updating if more such species are added to the model
  !  Minor changes: 1. removed computation of [H+]: this is not used anywhere (this also affects 
  !                 the subroutine call in m7)
  !                 2. cleaned up mode and species indexing - fixed relations between mode numbers
  !                    and species identities no longer required
 
  !--- Local variables: to be completed

  USE mo_kind,            ONLY: dp
  USE mo_physical_constants, ONLY: avo
  USE mo_math_constants,     ONLY: pi
  USE mo_ham,              ONLY: naerocomp, nclass, nsol, subm_naerospec, subm_aerospec
  USE mo_ham_m7ctl,       ONLY: iaccs, icoas,                         &
                                wnacl,   wna2so4, wnahso4, wh2so4,    &
                                dna2so4, dnahso4, dh2so4,  dh2o,      &
                                crh, ram2cmr, cmin_aerml, cmin_aernl, cdconv
  USE mo_ham_subm_species, ONLY: isubm_wat

  USE mo_species,         ONLY: speclist
  USE mo_ham_species,     ONLY: id_ss, id_so4    !!mgs!! ham_aerospec, ham_naerospec 

  IMPLICIT NONE
  !
  !--- Parameter list:
  !
  !
  !--- Local variables:
  !
  ! !@@@ To be completed!
  !
  ! Indices for the arrays zmm and zmmr:
  !  -----------------------------------
  ! |     ions               ion pair   |
  ! |  (mole m-3)            (mole m-3) |
  ! | i   zmm(i)             zmmr(i)    |
  ! | 1   Na+                NaCl       |
  ! | 2   Cl-                NaHSO4     |
  ! | 3   SO4--              Na2SO4     |
  ! | 4   HSO4-              H2-SO4     |
  ! | 5   H+                            |
  !  -----------------------------------

  !--- Parameters:
  
  INTEGER :: kproma, kbdim, klev, krow  

  REAL(dp):: prelhum(kbdim,klev),        ptp1(kbdim,klev)
  
  REAL(dp):: paerml(kbdim,klev,naerocomp), paernl(kbdim,klev,nclass), &
             pm6rp(kbdim,klev,nclass),     pm6dry(kbdim,klev,nsol), &
             prhop(kbdim,klev,nclass),     pww(kbdim,klev,nclass)

  !--- Local variables:

  INTEGER :: jl,            jk,           jclass, jn, jm

  REAL(dp):: zaw,         zvolw,         zdryvol,      zdrymass,  &
             zdryvol_mean,  zambvol_mean

  !>>dod soa removed zmol
  REAL(dp):: zmm(4),      zmmr(4),       zmo(4),     zmmt(4)
  !<<dod

  REAL(dp) :: zothermass(kbdim,klev)
  REAL(dp) :: zothervol(kbdim,klev)
  REAL(dp) :: zdss 
  INTEGER :: ixss, ixso4 
  !-------------------------------------------------------------------------------
  !   
  !
  !---initialisations
  zothermass(:,:) = 0._dp
  zothervol(:,:) = 0._dp
  zdss = cdconv*speclist(id_ss)%density

  !---calculation of wet count median radii and water mass for soluble modes
  DO 50 jclass=iaccs, icoas
     
     ixso4 = speclist(id_so4)%iaerocomp(jclass)
     ixss  = speclist(id_ss)%iaerocomp(jclass)

     !---calculate mass [g m-3] and volume [cm3 m-3] of 
     !   non-electrolyte species (bc, oc, du, soa) in the mode 
     DO jn=1,subm_naerospec     !!mgs!!     ham_naerospec
        IF ((.NOT. speclist(subm_aerospec(jn))%lelectrolyte) .AND. (jn /= isubm_wat)) THEN
           jm = speclist(subm_aerospec(jn))%iaerocomp(jclass)
           DO jk=1,klev
              DO jl=1,kproma
                 !--- mass in g m-3
                 zothermass(jl,jk) = zothermass(jl,jk) + 1.e-6_dp*paerml(jl,jk,jm)  
                 !--- vmr in cm3 m-3
                 zothervol(jl,jk) = zothervol(jl,jk) + &
                                   1.e-6_dp*paerml(jl,jk,jm)/(cdconv*speclist(subm_aerospec(jn))%density)
              END DO
           END DO
        END IF
     END DO

     DO 60 jk=1,klev
        DO 70 jl=1,kproma

           IF ((paerml(jl,jk,ixss) > cmin_aerml) .AND. (paernl(jl,jk,jclass)>cmin_aernl))  THEN
              !
              !--- 1) Dry Calculations: --------------------------------------------------
              !
              !--- 1.1) Calculate initial concentrations of the compounds in mole/m+3:
              !
              !--- Na, Cl:
              !
              zmm(1)=paerml(jl,jk,ixss)*1.E-6_dp / wnacl       ! n(Na)    [mole m-3]
              zmm(2)=paerml(jl,jk,ixss)*1.E-6_dp / wnacl       ! n(Cl)    [mole m-3]
              !      [         g m-3           ] / [g mole-1]           = [mole m-3]
              !
              !--- SO4--:
              !
              zmm(3)=paerml(jl,jk,ixso4) *1.E+6_dp / avo         ! n(H2SO4) [mole m-3]
              !      [      m-3                ] / [mole-1]             = [mole m-3]
              !
              !--- HSO4-:
              !
              zmm(4)=0._dp                                        ! n(HSO4)  [mole m-3]
              !
              !--- 1.2) Calculation of the concentration of the different species:
              !         The ions are supposed to be arranged such that
              !         sodium is associated to sulphate first in
              !         Na2SO4, the remaining sodium is associated to Cl
              !         in form of NaCl.
              !
              zmmt(1)=zmm(1)                        ! n(Na)  
              zmmt(3)=zmm(3)                        ! n(SO4)    
              !
              zmmr(3)=MIN(zmmt(1)/2._dp , zmmt(3))     ! n(Na2SO4) 
              zmmt(1)=zmmt(1)-2._dp*zmmr(3)             ! Remaining n(Na) after association
                                                    ! with Na2SO4: n(Na)=n(Na)-2*n(Na2SO4)
              !
              zmmr(1)=MIN(zmm(2),zmmt(1))           ! n(NaCl) 
              !
              zmm(2)=zmmr(1)                        ! n(Cl) bound in NaCl, the rest is
                                                    ! assumed to evaporate in presence of SO4
              !
              zmmr(2)=0._dp                            ! n(NaHSO4)
              !
              zmmr(4)=zmm(3)-zmmr(2)-zmmr(3)        ! n(H2-SO4)(t)=n(H2SO4)(t0)-n(NaHSO4)-n(Na2SO4)
              !                                     ! as n(H2SO4)(t0)=n(SO4--)(t0)
              !
              !--- 1.3) Total aerosol dry volume [cm3/m3]:
              !
              zdryvol= zmmr(1)*wnacl/zdss + &
                       zmmr(2)*wnahso4/dnahso4           + &
                       zmmr(3)*wna2so4/dna2so4           + &
                       zmmr(4)*wh2so4/dh2so4             + &
                       zothervol(jl,jk)

              !
              !--- 1.4) Mean aerosol dry volume [cm+3]:
              zdryvol_mean = zdryvol    / (paernl(jl,jk,jclass)*1.E6_dp)
              ! [cm+3]     = [cm+3/m+3] /  [         m-3            ]
              !
              !--- 1.5) Dry radius [cm]:
              !
              pm6dry(jl,jk,jclass)=((3._dp/(4._dp*pi))*zdryvol_mean)**(1._dp/3._dp) * &
                                  ram2cmr(jclass)
              !
              !--- 1.6) Total aerosol dry mass [gr/m3]:
              !
              zdrymass= zmmr(1)*wnacl                    + &
                        zmmr(2)*wnahso4                  + &
                        zmmr(3)*wna2so4                  + &
                        zmmr(4)*wh2so4                   + &
                        zothermass(jl,jk)
              !
              !
              !--- 2) Wet calculations: --------------------------------------------------
              !
              !--- Set threshold for relative humidity:
              !    If RH is smaller than the Critical Relative Humidity 
              !    (currently crh=0.45) the equilibrium radius is set to the dry radius:

              IF (prelhum(jl,jk) < crh) THEN
                 !
                 pww(jl,jk,jclass)     = 0._dp
                 !
                 pm6rp(jl,jk,jclass) = pm6dry(jl,jk,jclass)
                 !
                 prhop(jl,jk,jclass) = zdrymass/zdryvol
                 !
              ELSE
                 !
                 !--- 2.1) Calculate thermodynamic properties under ambient conditions
                 !
                 !--- 2.1.1) Water activity:
                 !
                 zaw=prelhum(jl,jk)
                 !
                 !--- 2.1.2) Molality as function of the water activity:
                 !           Currently sulfate is assumed to be fully dissociated,
                 !           i.e. zmmr(2)=0. and zmo(2) is not calculated.
                 !
                 !           Changed reference to Jacobson et al. (1996):
                 !
                 !--- NaCl:
                 
                 zmo(1)=(-1.918004E2_dp+2.001540E3_dp*zaw-8.557205E3_dp*zaw**2           &
                         +1.987670E4_dp*zaw**3-2.717192E4_dp*zaw**4+2.187103E4_dp*zaw**5 &
                         -9.591577E3_dp*zaw**6+1.763672E3_dp*zaw**7                    )**2 

                 !--- NaHSO4:

                 zmo(2)=(+4.662777E0_dp-1.128472E1_dp*zaw+7.049464E1_dp*zaw**2           &
                         -2.788050E2_dp*zaw**3+6.103105E2_dp*zaw**4-7.409417E2_dp*zaw**5 & 
                         +4.614577E2_dp*zaw**6-1.150735E2_dp*zaw**7                    )**2 

                 !--- Na2SO4:

                 zmo(3)=(-3.295311E3_dp+3.188349E4_dp*zaw-1.305168E5_dp*zaw**2           &
                         +2.935608E5_dp*zaw**3-3.920423E5_dp*zaw**4+3.109519E5_dp*zaw**5 &
                         -1.356439E5_dp*zaw**6+2.510249E4_dp*zaw**7                    )**2

                 !--- H2-SO4:
                 
                 zmo(4)=(+5.611895_dp-1.387446E1_dp*zaw+1.750682E1_dp*zaw**2             &
                         +7.138146E1_dp*zaw**3-3.109173E2_dp*zaw**4+4.662288E2_dp*zaw**5 &
                         -3.128612E2_dp*zaw**6+7.76097E1_dp*zaw**7                     )**2
                 !
                 !
                 !--- 2.2) Calculation of the water content in kg water/m3 air:
                 !         (zmmr[mole/m3(air)]/zmo[mole/kg(water)] =zww[kg(water)/m3(air)]
                 !
                 pww(jl,jk,jclass)=zmmr(1)/zmo(1)+zmmr(2)/zmo(2)+zmmr(3)/zmo(3)+zmmr(4)/zmo(4)
                 !
                 !--- 2.3) Calculate the molality of the ions
                 !
                 !--- 2.3.1) For Na+, Cl-, SO4--, HSO4- :
                 !
                 !>>dod soa removed these calculations - not used anywhere
                 !<<dod
                 
                 !--- 2.4) Calculation of the wet radius
                 !
                 !--- 2.4.1) Total water volume  [cm3/m3]:
                 !
                 zvolw=pww(jl,jk,jclass)/(dh2o*1.E-3_dp)  ![cm3/m3]=[kg/m3]/([g/cm3]*[1.E-3 kg/g])
                 !
                 !
                 !--- 2.4.2) Mean aerosol ambient volume:
                 zambvol_mean = (zdryvol+zvolw) / (paernl(jl,jk,jclass)*1.E6_dp)
                 ! [cm+3]     = [  cm+3/m+3   ] / [         m-3             ]                 
                 
                 !--- 2.4.3) Equilibrium wet count mean radius [cm]:
                 !
                 pm6rp(jl,jk,jclass)=((3._dp/(4._dp*pi))*zambvol_mean)**(1._dp/3._dp) * &
                                     ram2cmr(jclass)
                 !
                 !--- 2.4.4) Calculation of the particle density (g cm-3):
                 !
                 prhop(jl,jk,jclass)=(zdrymass+zvolw*dh2o)/(zvolw+zdryvol)
                 !

              END IF !(prelhum(jl,jk) < crh)

           END IF !((paerml(jl,jk,ixss) > cmin_aerml) .AND. (paernl(jl,jk,jclass) > cmin_aernl)) 

70      END DO
60   END DO
50 END DO

END SUBROUTINE m7_equil




SUBROUTINE m7_h2so4_cs(kproma,kbdim,klev,ptp1,papp1,paernl,pm6rp,pcs,pcsi)

!
! *m7_h2so4_cs* calculates the H2SO4 condensation sink (s-1) of the given
! aerosol size distribution and of the individual modes.
!
! Method:
! -----------------
! The transfer of sulfate to the particles is based on Fuchs (1959). Soluble
! and insoluble particles are distinguished by a different accomodation
! coefficient "caccso4" defined in mo_ham_m7ctl. (Currently 1.0 for soluble and
! 0.3 for insoluble modes). 
!
! References:
! -----------
! Fuchs, N.A. (1959). Evaporation and droplet growth in gaseous media; Pergamon, New York, pp72.
!
! Authors:
!
! 2001 Philip Stier <philip.stier@atm.ox.ac.uk> AOPP, University of Oxford
! 2000 J. Wilson, E. Vignati, JRC/EI (original source)

USE mo_kind,         ONLY: dp
USE mo_time_control, ONLY: time_step_len
USE mo_ham_m7ctl,    ONLY: wh2so4,rerg,caccso4
USE mo_math_constants, ONLY: pi
USE mo_ham,          ONLY: nclass

IMPLICIT NONE

!
! Input
!

INTEGER :: kproma, kbdim, klev

REAL(dp):: ptp1(kbdim,klev) ! Temperature (K)
REAL(dp):: papp1(kbdim,klev) ! Pressure (Pa)
REAL(dp):: paernl(kbdim,klev,nclass) ! Aerosol number for each mode (cm-3)
REAL(dp):: pm6rp(kbdim,klev,nclass) ! Mean mode actual radius (wet radius for soluble modes and dry radius for insoluble modes) (cm)

!
! Output
!

REAL(dp):: pcs(kbdim,klev)       ! H2SO4 condensation sink of the aerosol population (s-1)
REAL(dp):: pcsi(kbdim,klev,nclass) ! H2SO4 condensation sink of the individual aerosol modes (s-1)

!
! Local variables:
!

INTEGER :: jclass,jk,jl

REAL(dp):: zpbyone,zde2,zvelb,zxibc,zm6rp,zf1,zqtmst

!--- 0) Initizalizations: ------------------------------------------------

pcs(:,:)=0.0_dp

pcsi(:,:,:) = 0.0_dp  

zqtmst=1._dp/time_step_len

!--- 1) Calculate condensation rate for cm diameter sulphate aerosols: ---

DO jclass=1,nclass
   DO jk=1,klev
      DO jl=1,kproma
         IF (pm6rp(jl,jk,jclass).GT.0._dp) THEN                          ! try re-code avgprops so that this is never the case. (call 1 of avgprops)
            
            !--- Diffusion coefficient (Reference???):
            
            zpbyone=1000.0_dp / (papp1(jl,jk)/100.0_dp)
            
            zde2=0.073_dp * zpbyone * (ptp1(jl,jk) / 298.15_dp)**1.5_dp  ! merge this and the previous statement a la coaset. In 
                                                                         ! fact consider re-using it in coaset. NB multiply by 4. 
                                                                         ! in subsequent usage, zde2 is always multiplied by 4.
            
            !--- Mean molecule velocity (Moore, 1962 (S+P equ. 8.2)):
            
            zvelb=SQRT(8.0_dp * rerg * ptp1(jl,jk) / pi / wh2so4)       ! merge the constants. isolate sqrt(t) since it usually doesn't vectorise.
            
            !--- ???Fuchs???
            
            zxibc=8.0_dp * zde2 / pi / zvelb                            ! again merge constants. consider separate calc of 3*zxibc 
            
            ! Use count median radius:
            
            zm6rp=pm6rp(jl,jk,jclass)
            
            !--- Distance from particle up to which the kinetic regime applies:
            
            zf1=( (zm6rp + zxibc)**3 - (zm6rp**2 + zxibc**2)**1.5_dp ) / &
                (3.0_dp * zm6rp * zxibc) - zm6rp
            
            !--- Diffusive flux to single particle surface:
            !    (Elisabetta's thesis: fraction in equ. 2.26)
            
            pcsi(jl,jk,jclass)=(4.0_dp * pi * zde2 * zm6rp ) /                      &
                             ((4.0_dp * zde2) / (zvelb * zm6rp * caccso4(jclass)) +   &
                             (zm6rp/(zm6rp+zf1))                              )
            
            !--- Total diffusive flux to all particles in the respective mode:
            !    (per concentration of gas phase sulfate)
            
            pcsi(jl,jk,jclass)=pcsi(jl,jk,jclass) * paernl(jl,jk,jclass)
            
            !--- Total diffusive flux to all particles in all modes:
            !    (per concentration of gas phase sulfate)
            
            pcs(jl,jk)=pcs(jl,jk)+ pcsi(jl,jk,jclass)  
            
         END IF
      END DO
   END DO
END DO

END SUBROUTINE m7_h2so4_cs


  
  
SUBROUTINE m7_prod_cond(kproma, kbdim,  klev,  krow, &
                        pso4g,  dpso4g, paerml,      &
                        pso4_5, pso4_6, pso4_7,      &
                        pcs,    pcsi,   paclc, pdz)
  !
  ! Authors:
  ! --------
  !
  ! 2000 J. Wilson and E. Vignati(original code)
  ! 2001 P. Stier (f90-version, changes, comments)
  ! 2008 J. Kazil
  ! - cloud correction
  ! - new time integration scheme
  ! - comments
  !
  ! Interface:
  ! ----------
  !
  ! *m7_prod_cond* is called from *m7* before *m7_nuck*
  !
  ! Purpose:
  ! --------
  !
  ! *prod_cond* performs the production/condensation part of the time
  ! integration of the differential equation for sulfate production/condensation
  ! /nucleation over a time step, removes the condensing H2SO4 from the gas
  ! phase and distributes it onto the individual aerosol modes, given their
  ! H2SO4 condensation sinks. Inside clouds, all H2SO4 is removed from the gas
  ! phase.
  
  USE mo_kind,          ONLY: dp
  USE mo_time_control,  ONLY: time_step_len,delta_time
  USE mo_physical_constants, ONLY: avo
  USE mo_ham_m7ctl,     ONLY: iso4ns,iso4ks,iso4as,iso4cs,cmin_aernl
  USE mo_ham,           ONLY: naerocomp, mw_so4, nclass
  USE mo_ham_streams,   ONLY: d_cond_so4
  
  IMPLICIT NONE
  
  !
  ! Input/output:
  !
  
  INTEGER :: kproma, kbdim, klev, krow  ! ECHAM indices
  
  REAL(dp):: pso4g(kbdim,klev)          ! [H2SO4(g)] [molec. cm-3]
  REAL(dp):: dpso4g(kbdim,klev)         ! [H2SO4(g)] production rate [molec. cm-3]
  REAL(dp):: paerml(kbdim,klev,naerocomp) ! aerosol mass for individual compounds [molec. cm-3 for sulfate and ug m-3 for bc, oc, ss, and dust]
  REAL(dp):: pso4_5(kbdim,klev)         ! H2SO4 condensing on insoluble mode 5 [molec. cm-3]
  REAL(dp):: pso4_6(kbdim,klev)         ! H2SO4 condensing on insoluble mode 6 [molec. cm-3]
  REAL(dp):: pso4_7(kbdim,klev)         ! H2SO4 condensing on insoluble mode 7 [molec. cm-3]
  REAL(dp):: pcs(kbdim,klev)            ! H2SO4 condensation sink of the aerosol population [s-1]
  REAL(dp):: pcsi(kbdim,klev,nclass)      ! H2SO4 condensation sink of the individual aerosol modes [s-1]
  REAL(dp):: paclc(kbdim,klev)          ! Cloud cover [0,1]
  REAL(dp):: pdz(kbdim,klev)            ! layer thickness dz for diagnostics [m]
  
  !
  ! Local variables:
  !
  
  INTEGER :: jl,jk
  
  REAL(dp):: zqtmst
  REAL(dp):: zh2so4_cf0 ! Start [H2SO4(g)] [molec. cm-3] in the cloud-free part of the grid box
  REAL(dp):: zh2so4_cf1 ! End [H2SO4(g)] [molec. cm-3] in the cloud-free part of the grid box
  REAL(dp):: zh2so4_cy0 ! Start [H2SO4(g)] [molec. cm-3] in the cloudy part of the grid box
  REAL(dp):: zh2so4_cy1 ! End [H2SO4(g)] [molec. cm-3] in the cloudy part of the grid box
             
  REAL(dp):: zfcond_cf  ! [H2SO4(g)] condensing during the time step in the cloud-free part of the grid box
  REAL(dp):: zfcond_cy  ! [H2SO4(g)] condensing during the time step in the cloudy part of the grid box
  
  REAL(dp):: cc         ! Corrected cloud cover [0,1]
  
  ! Initialisations:
  
  zqtmst=1.0_dp/time_step_len
  
  ! Calculation of the new sulfate aerosol masses and of the mass of sulfate
  ! condensing on the respective modes:
  
  DO jk=1,klev
    DO jl=1,kproma
      
      IF (pcs(jl,jk) > 1.0E-10_dp) THEN  ! Regular H2SO4 condensation sink of the aerosol
        
        ! Safety check on the cloud fraction:
        cc = min(paclc(jl,jk),1.0_dp)
        cc = max(cc,0.0_dp)
        
        !
        ! Cloud-free part of the grid box:
        !
        
        ! Start [H2SO4(g)]:
        zh2so4_cf0 = pso4g(jl,jk)
        
        ! End [H2SO4(g)]:
        zh2so4_cf1 = &
        (zh2so4_cf0 - dpso4g(jl,jk)/pcs(jl,jk))*exp(-pcs(jl,jk)*time_step_len) &
        + dpso4g(jl,jk)/pcs(jl,jk)
        
        ! Safety check:
        zh2so4_cf1 = max(zh2so4_cf1,0.0_dp)
        
        !
        ! Cloudy part of the grid box:
        !
        
        ! Start [H2SO4(g)]:
        zh2so4_cy0 = pso4g(jl,jk)
        
        ! End [H2SO4(g)]:
        zh2so4_cy1 = 0.0_dp
        
        ! 
        ! New grid box-averaged [H2SO4(g)]:
        !
        
        pso4g(jl,jk) = (1.0_dp-cc)*zh2so4_cf1 + cc*zh2so4_cy1
        
        ! [H2SO4(g)] that condensed onto aerosol particles during the time step
        ! in the cloud-free and cloudy part of the grid box:
        
        zfcond_cf = (1.0_dp-cc)*(zh2so4_cf0 - zh2so4_cf1 + dpso4g(jl,jk)*time_step_len)
        zfcond_cy = cc*(zh2so4_cy0 - zh2so4_cy1 + dpso4g(jl,jk)*time_step_len)
        
        ! Distribute H2SO4(g) condensing in the cloud-free part of the grid box
        ! on the soluble and insoluble aerosol modes according to their
        ! condensation sinks:
        
        paerml(jl,jk,iso4ns) = paerml(jl,jk,iso4ns) + pcsi(jl,jk,iso4ns)/pcs(jl,jk)*zfcond_cf
        paerml(jl,jk,iso4ks) = paerml(jl,jk,iso4ks) + pcsi(jl,jk,iso4ks)/pcs(jl,jk)*zfcond_cf
        paerml(jl,jk,iso4as) = paerml(jl,jk,iso4as) + pcsi(jl,jk,iso4as)/pcs(jl,jk)*zfcond_cf
        paerml(jl,jk,iso4cs) = paerml(jl,jk,iso4cs) + pcsi(jl,jk,iso4cs)/pcs(jl,jk)*zfcond_cf
        
        ! Number of H2SO4 molecules condensing on the insoluble modes:
        ! (Transfer from insoluble to soluble modes is calculated in m7_concoag)
        
        pso4_5(jl,jk) = pcsi(jl,jk,5)/pcs(jl,jk)*zfcond_cf
        pso4_6(jl,jk) = pcsi(jl,jk,6)/pcs(jl,jk)*zfcond_cf
        pso4_7(jl,jk) = pcsi(jl,jk,7)/pcs(jl,jk)*zfcond_cf
        
        ! Commit the H2SO4(g) condensing in the cloudy part of the grid box
        ! to the largest aerosol mode, assuming that that mode contributes
        ! aerosol particles that activated and became cloud droplets:
        
        IF (pcsi(jl,jk,iso4cs) > 0.0_dp) THEN
          paerml(jl,jk,iso4cs) = paerml(jl,jk,iso4cs) + zfcond_cy
        ELSEIF (pcsi(jl,jk,iso4as) > 0.0_dp) THEN
          paerml(jl,jk,iso4as) = paerml(jl,jk,iso4as) + zfcond_cy
        ELSEIF (pcsi(jl,jk,iso4ks) > 0.0_dp) THEN
          paerml(jl,jk,iso4ks) = paerml(jl,jk,iso4ks) + zfcond_cy
        ELSEIF (pcsi(jl,jk,iso4ns) > 0.0_dp) THEN
          paerml(jl,jk,iso4ns) = paerml(jl,jk,iso4ns) + zfcond_cy
        ELSEIF (pcsi(jl,jk,7) > 0.0_dp) THEN
          pso4_7(jl,jk) = pso4_7(jl,jk) + zfcond_cy
        ELSEIF (pcsi(jl,jk,6) > 0.0_dp) THEN
          pso4_6(jl,jk) = pso4_6(jl,jk) + zfcond_cy
        ELSEIF (pcsi(jl,jk,5) > 0.0_dp) THEN
          pso4_5(jl,jk) = pso4_5(jl,jk) + zfcond_cy
        ENDIF
        
        ! Vertically integrate mass of condensed sulfate for diagnostics,
        ! converting [molec. cm-3] to [kg(SO4) m-2]:
        
        d_cond_so4(jl,krow) = d_cond_so4(jl,krow) &
         + (((zfcond_cf + zfcond_cy)*mw_so4*1.E3_dp)/avo)*pdz(jl,jk)*zqtmst*delta_time
        
      ELSE ! H2SO4 condensation sink of the aerosol extremely small:
        
        ! New grid box-averaged [H2SO4(g)]:
        pso4g(jl,jk) = pso4g(jl,jk) + time_step_len*dpso4g(jl,jk)
        
        ! Number of H2SO4 molecules condensing on the insoluble modes:
        ! (Transfer from insoluble to soluble modes is calculated in m7_concoag)
        pso4_5(jl,jk) = 0.0_dp
        pso4_6(jl,jk) = 0.0_dp
        pso4_7(jl,jk) = 0.0_dp
        
      END IF
    END DO
  END DO
  
END SUBROUTINE m7_prod_cond


SUBROUTINE m7_nuck(kproma,  kbdim,  klev,   krow,          &
                   papp1,   ptp1,   prhp1,  paclc,   pipr, &
                   ph2so4,  pcs,    panew,  pa4delt, pdz,  &
                   ppbl)
  ! 
  ! Authors:
  ! --------
  !
  ! 2000 J. Wilson and E. Vignati (original code)
  ! 2001-2003 P. Stier
  ! - F90-version
  ! - Changes
  ! - Comments
  ! - Modularisation
  ! - Implemented H2SO4/H2O nucleation scheme (Vehkamaeki et al., JGR 2002)
  ! 2008 J. Kazil
  ! - Removed H2SO4/H2O nucleation (Kulmala et al., JGR 1998)
  ! - Implemented neutral and charged H2SO4/H2O nucleation (Kazil and Lovejoy, ACP 2007)
  ! - Implemented organic nucleation parameterizations:
  !   activation nucleation (Kulmala et al., ACP 2006)
  !   kinetic nucleation (Laakso et al., ACP 2004)
  ! - Implemented a cloud correction
  ! - Implemented new time integration scheme
  ! 
  ! Interface:
  ! ----------
  !
  ! *m7_nuck* is called from *m7* after *m7_prod_cond*
  ! 
  ! Purpose:
  ! --------
  !
  ! *m7_nuck* performes the nucleation part of the time integration of the
  ! differential equation for sulfate production/condensation/nucleation over a
  ! time step.
  !
  ! Method:
  ! -------
  !
  ! First, the grid box averaged H2SO4 gas phase concentration is converted to
  ! the H2SO4 gas phase concentration in the cloud-free portion of the grid box,
  ! assuming that [H2SO4(g)] = 0 in the cloudy part. Then the nucleation rate J
  ! is calculated for the cloud-free portion of the grid box and the resulting
  ! change in gas phase H2SO4 integrated over the time step. From this change,
  ! the number of newly formed aerosol particles is calculated. Finally, the
  ! concentrations of H2SO4(g) and newly formed aerosol particles are converted
  ! to grid box averaged values.
  !
  ! Note: The removal rate of H2SO4 from the gas phase due to nucleation and its
  ! conversion rate to freshly nucleated particles are approximated with n*J,
  ! where n is the H2SO4 content of the smallest newly formed particles. These
  ! are approximations which neglect the loss of freshly nucleated particles
  ! containing < n H2SO4 molecules onto pre-existing aerosol, and the
  ! contribution to new particle formation by coagulation.
  !
  ! References:
  ! -----------
  !
  ! Vignati, E., Modelling Interactions between Aerosols and 
  ! Gaseous Compounds in the Polluted Marine Atmosphere. PhD-Thesis,
  ! RISO National Laborartory Copenhagen, Riso-R-1163(EN), 1999
  
  USE mo_time_control, ONLY: time_step_len, delta_time
  USE mo_kind,         ONLY: dp
  USE mo_physical_constants, ONLY: avo
  USE mo_ham_m7ctl,    ONLY: nsnucl, nonucl
  USE mo_ham,          ONLY: naerocomp, mw_so4
  USE mo_ham_streams,  ONLY: d_nuc_so4
  USE mo_ham_m7_nucl,  ONLY: nucl_vehkamaeki, nucl_kazil_lovejoy, &
                             nucl_activation, nucl_kinetic
  USE mo_memory_g3b,   ONLY: forest
  
  IMPLICIT NONE
  
  !
  ! Input/output:
  !
  
  INTEGER :: kproma, kbdim, klev, krow   ! ECHAM indices
  
  REAL(dp):: papp1(kbdim,klev)           ! Atmospheric pressure at time t+1 [Pa]
  REAL(dp):: ptp1(kbdim,klev)            ! Atmospheric temperature at time t+1 [K]
  REAL(dp):: prhp1(kbdim,klev)           ! Relative humidity at time t+1 [0,1]
  REAL(dp):: paclc(kbdim,klev)           ! Cloud cover [0,1]
  REAL(dp):: pipr(kbdim,klev)            ! Ion pair production rate [cm-3 s-1]
  REAL(dp):: ph2so4(kbdim,klev)          ! Gas phase H2SO4 concentration [molec. cm-3]
  REAL(dp):: pcs(kbdim,klev)             ! H2SO4 condensation sink [s-1]
  REAL(dp):: panew(kbdim,klev)           ! Concentration of nucleated particles [cm-3] over a timestep
  REAL(dp):: pa4delt(kbdim,klev,naerocomp) ! pa4delt(:,:,1) = Number of H2SO4 per cm-3 added
                                         ! to the nucleation mode due to H2SO4 nucleation over a timestep
  REAL(dp):: pdz(kbdim,klev)             ! Layer thickness dz for diagnostics [m]
  REAL(dp):: ppbl(kbdim)                 ! Planetary boundary layer top level
  
  !
  ! Local variables:
  !
  
  INTEGER :: jk,jl
  
  REAL(dp):: ztmst,zqtmst,zeps
  
  REAL(dp):: zrh(kbdim,klev)             ! Relative humidity [%]
  REAL(dp):: zh2so4_cf0(kbdim,klev)      ! Initial gas phase H2SO4 concentration in the cloud-free portion of the grid box [molec. cm-3]
  REAL(dp):: zh2so4_cf1(kbdim,klev)      ! Final gas phase H2SO4 concentration in the cloud-free portion of the grid box [molec. cm-3]
  REAL(dp):: zsnrate(kbdim,klev)         ! H2SO4/H2O nucleation rate in the cloud-free portion of the grid box [cm-3 s-1]
  REAL(dp):: zsns(kbdim,klev)            ! Number of H2SO4 molecules new particles from H2SO4/H2O nucleation
  REAL(dp):: zonrate(kbdim,klev)         ! Organic aerosol nucleation rate in the cloud-free portion of the grid box [cm-3 s-1]
  REAL(dp):: zons(kbdim,klev)            ! Number of H2SO4 molecules new particles from organic nucleation
  
  REAL(dp):: znucrate                    ! Total nucleation (particle formation) rate
  REAL(dp):: zsnucloss                   ! Total H2SO4(g) loss due to nucleation
  
  ! Initialisations
  
  ztmst = time_step_len
  
  zqtmst = 1.0_dp/ztmst
  
  zeps = EPSILON(1.0_dp)
  
  ! Relative humidity [%]:
  
  zrh(1:kproma,:) = 100.0_dp*prhp1(1:kproma,:)
  
  ! Gas phase H2SO4 concentration in the cloud-free portion of the grid box:
  
  DO jk = 1, klev
    DO jl = 1, kproma
      
      IF (paclc(jl,jk) >= 1.0_dp) THEN
        zh2so4_cf0(jl,jk) = 0.0_dp
      ELSEIF (paclc(jl,jk) <= 0.0_dp) THEN
        zh2so4_cf0(jl,jk) = ph2so4(jl,jk)
      ELSE
        zh2so4_cf0(jl,jk) = ph2so4(jl,jk)/(1.0_dp-paclc(jl,jk))
      ENDIF
      
    ENDDO
  ENDDO
  
  ! H2SO4/H2O nucleation:
  
  IF (nsnucl == 1) THEN
    CALL nucl_vehkamaeki(kproma,kbdim,klev,ptp1,prhp1,zh2so4_cf0,zsnrate,zsns)
  ELSE IF (nsnucl == 2) THEN
    CALL nucl_kazil_lovejoy(kproma,kbdim,klev,ptp1,zrh,zh2so4_cf0,pcs,pipr,zsnrate,zsns)
  ELSE
    zsns(1:kproma,:)    = 0.0_dp
    zsnrate(1:kproma,:) = 0.0_dp
  END IF
  
  ! Organic nucleation:
  
  IF (nonucl == 1) THEN
    CALL nucl_activation(kproma,kbdim,klev,zh2so4_cf0,forest(:,krow),ppbl,zonrate,zons)
  ELSEIF (nonucl == 2) THEN
    CALL nucl_kinetic(kproma,kbdim,klev,zh2so4_cf0,forest(:,krow),ppbl,zonrate,zons)
  ELSE
    zons(1:kproma,:)    = 0.0_dp
    zonrate(1:kproma,:) = 0.0_dp
  END IF
  
  DO jk = 1, klev
    DO jl = 1, kproma
      
      ! Total nucleation (particle formation) rate:
      znucrate = zsnrate(jl,jk) + zonrate(jl,jk)
      
      IF (znucrate > zeps) THEN
        
        ! H2SO4(g) loss due to nucleation:
        zsnucloss = zsns(jl,jk)*zsnrate(jl,jk) + zons(jl,jk)*zonrate(jl,jk)
        
        ! Calculate the gas phase concentration after nucleation:
        ! N(t)   =   N(0) - (zsnucloss  *  dt)/( 1  +   pcs *  dt)
        ! [cm-3] = [cm-3] - ([cm-3 s-1] * [s])/([1] + [s-1] * [s]) 
        zh2so4_cf1(jl,jk) = zh2so4_cf0(jl,jk) - ztmst*zsnucloss/(1.0_dp+ztmst*pcs(jl,jk))
        
        ! Safety check:
        zh2so4_cf1(jl,jk) = max(zh2so4_cf1(jl,jk),0.0_dp)
        zh2so4_cf1(jl,jk) = min(zh2so4_cf1(jl,jk),zh2so4_cf0(jl,jk))
        
        ! Calculate concentration of nucleated H2SO4 from the net loss of gas
        ! phase H2SO4 in the cloud-free grid box portion, and average it over
        ! the grid box:
        pa4delt(jl,jk,1) = (zh2so4_cf0(jl,jk)-zh2so4_cf1(jl,jk))*(1.0_dp-paclc(jl,jk))
        
        ! Calculate the number of nucleated particles, making sure that total
        ! H2SO4 (gas + aerosol contributions) is conserved:
        panew(jl,jk) = pa4delt(jl,jk,1)*znucrate/zsnucloss
        
        ! Calculate the new grid box averaged gas phase H2SO4 concentration:
        ph2so4(jl,jk) = ph2so4(jl,jk)-pa4delt(jl,jk,1)
        
        ! Vertically integrate the mass of nucleated sulfate for diagnostics;
        ! convert [molec. cm-3] to [kg(S04) m-2]:
        d_nuc_so4(jl,krow) = d_nuc_so4(jl,krow) &
         + pa4delt(jl,jk,1)*mw_so4*1.E+03_dp/avo*pdz(jl,jk)*zqtmst*delta_time
        
      END IF
    
    END DO
  END DO
  
END SUBROUTINE m7_nuck


SUBROUTINE m7_dnum(kproma, kbdim,  klev,   krow,  &
                   paerml, paernl, ptp1,  &
                   papp1,  pm6rp,  prhop, &
                   !>>dod soa: added pttn
                   pso4_5, pso4_6, pso4_7, panew, pa4delt, pttn)
                   !<<dod
  !
  !  *m7_dnum*  changes gas phase sulfate, aerosol numbers and masses
  !             due to nucleation and coagulation
  !
  !  Authors:
  !  ---------
  !  J. Wilson  and E. Vignati, JRC/EI (original source)                 09/2000
  !  P. Stier, MPI                     (f90-version, changes, comments)     2001 
  !  D.O'Donnell MPI-M, optimisation
  !
  !  Version: 
  !  --------- 
  !  This version is equivalent to the version dnum2 of the m7 boxmodel. 
  !
  !  Purpose
  !  ---------
  !  This routine calculates new gas phase sulfate and aerosol
  !  numbers and masses after timestep ztmst.
  !
  !  Interface
  !  -----------
  !  *m7_dnum* is called from *m7*
  !
  !  Externals
  !  -----------
  !  *m7_coaset*   calculates the coagulation kernels
  !  *m7_delcoa*   integrates equations for the changes in aerosol numbers
  !                dn/dt=c -an2 -bn over one timestep and calculates the 
  !                corresponding changes in aerosol masses
  !  *m7_concoag*  calculates particle numbers and mass moved from the 
  !                insoluble to the mixed modes due to coagulation with
  !                smaller mixed particles and condensation of H2SO4.
  !
  !  Warning:
  !  --------
  !  For optimization purposes currently only "physically reasonable" elements of the
  !  coagulation kernel zcom are calculated in m7_concoag. These elements are specified
  !  in the matrix locoagmask in mo_ham_m7ctl.   Check carefully and adapt locoagmask 
  !  accordingly before changes in the code below.
  
  USE mo_kind,       ONLY: dp
  USE mo_ham_m7ctl,  ONLY: ncoag
  USE mo_ham,        ONLY: naerocomp, nclass, lscoag
  !>>dod timers
  USE mo_control,    ONLY: ltimer
  USE mo_hammoz_timer, ONLY: timer_start, timer_stop, timer_ham_m7_coaset, timer_ham_m7_delcoa
  !<<dod

  IMPLICIT NONE
  
  !--- Input/output:
  !
  !  paerml     = total aerosol mass for each compound 
  !               [molec. cm-3 for sulphate and ug m-3 for bc, oc, ss, and dust]
  !  paernl     = aerosol number for each mode [cm-3]
  !  ptp1       = atmospheric temperature at time t+1 [K]
  !  papp1      = atmospheric pressure at time t+1 [Pa]
  !  pm6rp      = mean mode actual radius (wet radius for soluble modes 
  !               and dry radius for insoluble modes) [cm]
  !  prhop      = mean mode particle density [g cm-3]
  !  pso4_x     = mass of sulphate condensed on insoluble mode x []
  !  pa4delt    = On input: change in H2SO4 mass of the respective mode over one
  !               timstep due to:
  !                - nucleation of H2SO4 (calculated in m7_nuck).
  !               On output: change in H2SO4 mass of the respective mode over one
  !               timstep due to:
  !                - nucleation of H2SO4 (calculated in m7_nuck)
  !                - coagulation (calculated in m7_concoag)
  !  panew      = On input: number of nucleated particles panew, calculated in m7_nuck
  !  pttn       = mean mass in one particle (in moleccules for SO4 and in g for others)
  !
  !--- Local variables:
  !
  ! zcom            = general coagulation coefficient []
  !                   (calculated in m7_coaset)
  ! za              = effectively used coagulation coefficient for 
  !                   unimodal coagulation []
  ! zb              = effectively used coagulation coefficient for 
  !                   inter-modal coagulation []
  ! zbfractx(:,:,y) = fraction of the total number of particles removed by 
  !                   coagulation from mode x (finally calculated in m7_delcoa)
  !                   that is moved to mode y / y+1 (modes 5,6,7 and mode 1,2 resp.) [1]
  !                   !@@@ CAUTION: Inconsistent usage!!!
  
  !--- Input/output:
  
  INTEGER :: kproma, kbdim, klev, krow
  
  REAL(dp):: ptp1(kbdim,klev), &
             papp1(kbdim,klev), &
             pso4_5(kbdim,klev), &
             pso4_6(kbdim,klev), &
             pso4_7(kbdim,klev)
  
  REAL(dp):: paerml(kbdim,klev,naerocomp), &
             paernl(kbdim,klev,nclass), &
             pm6rp(kbdim,klev,nclass), &
             prhop(kbdim,klev,nclass)
  
  REAL(dp):: panew(kbdim,klev), &
             pa4delt(kbdim,klev,naerocomp)
  
  !>>dod soa
  REAL(dp) :: pttn(kbdim,klev,naerocomp)
  !<<dod soa

  !--- Local Variables:
  
  INTEGER :: jl, jk
  
  REAL(dp):: za(kbdim,klev,nclass), zb(kbdim,klev,nclass)

  !>>dod opt deleted zbfract6 and zbfract7 (not used)
  REAL(dp):: zbfract1(kbdim,klev,nclass-1),zbfract2(kbdim,klev,nclass-1),      &
             zbfract5(kbdim,klev,3) 
  !<<dod
        
  !>>dod opt changed from 4-D (kbdim,klev,nclass,nclass) to 3-D (kbdim,klev,ncoag)
  REAL(dp):: zcom(kbdim,klev,ncoag)
  !<<dod

  !--- 0) Initialisations: ----------------------------------------------
  
  !--- 1) Calculate  coagulation coefficients: --------------------------

  !>>dod timers
  IF (ltimer) CALL timer_start(timer_ham_m7_coaset)
  !<<dod
  !>>dod opt added krow for debugging 
  IF (lscoag) CALL m7_coaset(kproma, kbdim, klev,  krow, paernl, ptp1, &
                             papp1,  pm6rp, prhop, zcom          )

  !>>dod timers
  IF (ltimer) CALL timer_stop(timer_ham_m7_coaset)
  !<<dod

  !--- 3) Assign coagulation coefficients (unimodal coagulation)---------
  !       and the normalised fraction of particle numbers moved 
  !       between the modes (inter-modal coagulation):
  !
  !       The general equation for dn/dt for each mode is:
  !
  !       dn/dt=-za*n^2 - zb*n + zc
  !
  !       where za=unimodal coagulation coefficient (zcom(mod))
  !             zb=inter-modal coagulation with higher modes
  !                (zcom(mod) * n(jclass+1))
  !             zc=particle formation rate 
  !                (=panew/ztmst if jclass=1, zero for higher modes) 
  !
  !             zb is zero when n(jclass+1)=zero, or jclass=naerocomp
  
  IF (lscoag) THEN
    
     DO jk=1,klev
        DO jl=1,kproma 
        
           !---  3.1) Unimodal coagulation coefficients:
           !@@@Coag:
           !>>dod opt changed indexing of zcom
           za(jl,jk,1)=zcom(jl,jk,1)/2._dp    ! Unimodal coagulation
           za(jl,jk,2)=zcom(jl,jk,8)/2._dp    ! only allowed for modes
           za(jl,jk,3)=zcom(jl,jk,14)/2._dp    ! 1,2,3 and 5.
           za(jl,jk,4)=0._dp
           za(jl,jk,5)=zcom(jl,jk,16)/2._dp
           za(jl,jk,6)=0._dp
           za(jl,jk,7)=0._dp
           !<<dod

           !---  3.2) Inter-modal coagulation - soluble modes
         
           !--- Sum all higher mode coagulation terms for 
           !    soluble modes 1,2,3,4:
        
           !--- 3.2.1) Mode 1:
        
           !--- Number of particles (zbfract1(:,:,x)) that are moved 
           !    from mode 1 to the mode x+1 :
           !    !@@@ Clumsy! Change to x - also in concoag!!!
           !>>dod opt changed indexing of zcom
           zbfract1(jl,jk,1)=zcom(jl,jk,2)*paernl(jl,jk,2)
           zbfract1(jl,jk,2)=zcom(jl,jk,3)*paernl(jl,jk,3)
           zbfract1(jl,jk,3)=zcom(jl,jk,4)*paernl(jl,jk,4)
           zbfract1(jl,jk,4)=zcom(jl,jk,5)*paernl(jl,jk,5)
           zbfract1(jl,jk,5)=zcom(jl,jk,6)*paernl(jl,jk,6)
           zbfract1(jl,jk,6)=zcom(jl,jk,7)*paernl(jl,jk,7)
           !<<dod

           !--- Sum of all particles that are moved from mode 1:
        
           zb(jl,jk,1)=zbfract1(jl,jk,1)+zbfract1(jl,jk,2)+            &
                       zbfract1(jl,jk,3)+zbfract1(jl,jk,4)+            &
                       zbfract1(jl,jk,5)+zbfract1(jl,jk,6)
        
           !>>dod moved IF to try to improve prefetch
           !<<dod

           !--- 3.2.2) Mode 2:
        
           !--- Number of particles (zbfract1(:,:,x)) that are moved 
           !    from mode 2 to the mode x+1 :
        
           !>>dod changed zcom indexing...also commented out pointless calculation
           zbfract2(jl,jk,2)=zcom(jl,jk,9)*paernl(jl,jk,3)
           zbfract2(jl,jk,3)=zcom(jl,jk,10)*paernl(jl,jk,4)
           ! zbfract2(jl,jk,4)=zcom(jl,jk,11)*paernl(jl,jk,5)
           zbfract2(jl,jk,4)=0._dp
           zbfract2(jl,jk,5)=zcom(jl,jk,12)*paernl(jl,jk,6)
           zbfract2(jl,jk,6)=zcom(jl,jk,13)*paernl(jl,jk,7)
           !<<dod

           !--- Sum of all particles that are moved from mode 2:
        
           zb(jl,jk,2)=zbfract2(jl,jk,2)+zbfract2(jl,jk,3)+            &
                       zbfract2(jl,jk,4)+zbfract2(jl,jk,5)+            &
                       zbfract2(jl,jk,6)
        
           !>>dod moved IF to try to improve prefetch
           !<<dod

           !--- 3.2.3) Mode 3 and Mode 4 - considered ineffective:
        
           zb(jl,jk,3)=0.0_dp
        
           zb(jl,jk,4)=0.0_dp
        
           !--- 3.3) Inter-modal coagulation - insoluble modes
           !
           !         For the insoluble modes coagulation with soluble modes
           !         is a sink as they are transfered to the corresponding
           !         mixed/solublemode. Therefore, terms with a lower mode 
           !         number or the same mode number are included. 
           !         (!@@@ There are still some switches for testing.)
        
           !--- 3.3.1) Mode 5:
        
           !--- Number of particles (zbfract5(:,:,x)) that are moved 
           !    from mode 5 to the mode x:
        
           zbfract5(jl,jk,1)=0._dp
           zbfract5(jl,jk,2)=zcom(jl,jk,11)*paernl(jl,jk,2)
           zbfract5(jl,jk,3)=zcom(jl,jk,15)*paernl(jl,jk,3)
        
           !--- Sum of all particles that are moved from mode 5:
        
           zb(jl,jk,5)=zbfract5(jl,jk,1)+zbfract5(jl,jk,2)+zbfract5(jl,jk,3)
        
           !--- Normalize number of particles by the total number of 
           !    particles moved from mode 1:
        
           IF (zb(jl,jk,1).GT.0.0_dp) THEN
              zbfract1(jl,jk,1)=zbfract1(jl,jk,1)/zb(jl,jk,1)
              zbfract1(jl,jk,2)=zbfract1(jl,jk,2)/zb(jl,jk,1)
              zbfract1(jl,jk,3)=zbfract1(jl,jk,3)/zb(jl,jk,1)
              zbfract1(jl,jk,4)=zbfract1(jl,jk,4)/zb(jl,jk,1)
              zbfract1(jl,jk,5)=zbfract1(jl,jk,5)/zb(jl,jk,1)
              zbfract1(jl,jk,6)=zbfract1(jl,jk,6)/zb(jl,jk,1)
           END IF
        
           !--- Normalize particle numbers by the total number of 
           !    particles moved from mode 2:
        
           IF (zb(jl,jk,2).GT.0.0_dp) THEN
              zbfract2(jl,jk,2)=zbfract2(jl,jk,2)/zb(jl,jk,2)
              zbfract2(jl,jk,3)=zbfract2(jl,jk,3)/zb(jl,jk,2)
              zbfract2(jl,jk,4)=zbfract2(jl,jk,4)/zb(jl,jk,2)
              zbfract2(jl,jk,5)=zbfract2(jl,jk,5)/zb(jl,jk,2)
              zbfract2(jl,jk,6)=zbfract2(jl,jk,6)/zb(jl,jk,2)
           END IF
        
           !--- Normalize number of particles by the total number of 
           !    particles moved from mode 5:
        
           IF (zb(jl,jk,5).GT.0.0_dp) THEN
              zbfract5(jl,jk,1)=zbfract5(jl,jk,1)/zb(jl,jk,5)
              zbfract5(jl,jk,2)=zbfract5(jl,jk,2)/zb(jl,jk,5)
              zbfract5(jl,jk,3)=zbfract5(jl,jk,3)/zb(jl,jk,5)
           END IF
        
           !>>dod removed calculation for modes 6 and 7, which were not used
        
        END DO
     END DO
    
  ELSE
    
    za(:,:,:)       = 0._dp
    zb(:,:,:)       = 0._dp
    zbfract1(:,:,:) = 0._dp
    zbfract2(:,:,:) = 0._dp
    zbfract5(:,:,:) = 0._dp
    !>>dod removed modes 6 and 7
    
 END IF !(lscoag)

 !>>dod timers
 IF (ltimer) CALL timer_start(timer_ham_m7_delcoa)
 !<<dod
 !>>dod soa added krow and pttn to call of m7_delcoa
 CALL m7_delcoa(kproma,   kbdim,  klev,     krow,      paerml,   &
                paernl,   pm6rp,  pttn,     pa4delt,   panew,    &
                za,       zb,     zbfract1, zbfract2, &
                zbfract5, pso4_5, pso4_6,   pso4_7    )
 !>>dod timers
 IF (ltimer) CALL timer_stop(timer_ham_m7_delcoa)
 !<<dod

END SUBROUTINE m7_dnum


SUBROUTINE m7_dconc(kproma, kbdim, klev, krow, paerml, paernl, pm6dry)
  !
  !    *m7_dconc*  changes aerosol numbers and masses to account for
  !                    condensational growth of the mode mean radii
  !  
  !    Authors:
  !    --------
  !    D. O'Donnell. MPI-M, rewrite to reduce CPU and introduce SOA   Jan 2007
  !    J. Wilson and E. Vignati, JRC (original source)            May 2000
  !    P. Stier, MPI-MET (f90 version, changes, comments)             2001
  !
  !    Purpose:
  !    --------
  !    This routine repartitions aerosol number and mass between the
  !    the modes to account for condensational growth and the formation
  !    of an accumulation mode from the upper tail of the aitken mode.
  !
  !    Interface:
  !    ----------
  !    *m7_dconc* is called from *m7*
  !
  !    Method:
  !    -------
  !    The routine calculates the cumulative number and mass distribution of the
  !    modes up to the respective mode boundary:
  !
  !                        / x                              _
  !                 N      |       1           1   ln(R)-ln(R)  2
  !    N(0,x) = ---------  |   --------  exp(- - ( ----------- )   ) d ln(R) 
  !             ln(sigma)  |   sqrt(2PI)       2    ln(sigma)
  !                        / 0 
  !                         
  !                         /tx                   2
  !                        |        1            t
  !           =     N      |     --------  exp(- - ) d t 
  !                        |     sqrt(2PI)       2 
  !                        /-inf
  ! 
  !    where:                   
  !
  !                        _
  !               ln(R)-ln(R)
  !    t      =   -----------
  !                ln(sigma)
  !
  !    and:
  !                        _
  !               ln(x)-ln(R)
  !    tx     =   -----------
  !                ln(sigma)
  !    _
  !    R is the Count Mean Radius or the Mass Mean Radius.
  !
  !    Practically, the routine m7_cumulative_normal calculates the fraction of the number and
  !    mass distribution for each mode lying below the respective upper mode boundary (1).
  !    In a next step the net fraction of each mode lying between the upper and lower
  !    mode boundaries are summed up (2) and the numbers and masses exceeding the mode
  !    boundaries are transfered to the neighbouring larger mode (3). 
  !    Finally, these quantities are stored in the respective arrays
  !    paernl and paerml (4).
  !    The repartititioning is currently only done for the soluble modes as it is
  !    assumed that insoluble modes are rather transfered to the soluble modes 
  !    and grow as soluble particles.
  !
  !    When secondary organic aerosols (SOA) are present, the volatile nature of SOA means that,
  !    unlike previous versions of M7, aerosol mass can evaporate into the vapour phase. This 
  !    means that the total aerosol mass can shrink (even before considering deposition processes) and
  !    that accordingly it must be possible to transfer aerosol mass and number not only from one
  !    mode into the neighbouring larger mode but also into the neighbouring smaller mode.
  !    Accordingly we also evaluate N(0,x) using the above integrals but with
  !                   _
  !                ln(R)-ln(Rl)
  !    t       =   ------------
  !                 ln(sigma)
  !
  !    where Rl = lower mode boundary, and 
  !                   _
  !                ln(R)-ln(x)
  !    tx      =   ------------
  !                 ln(sigma)
  !
  !    This is subject to limits in the model, however, since sea salt is considered to reside only in 
  !    the accumulation and coarse modes and hence cannot be transferred 'downwards' to the Aitken mode,
  !    and likewise for black carbon and nonvolatiles (oc) to the nucleation mode. Because the mixed modes
  !    are considered to be internal mixtures, this means that in the presence of appreciable amounts
  !    of sea salt, no transfer is done from accumulation to Aitken mode and similarly for transfer from
  !    Aitken to nucleation mode in the presence of bc or oc.
  !
  !    Changes compared to previous m7_dconc
  !    -------------------------------------
  !    Functional changes:
  !    SOA is added. 
  !    Other changes:
  !    This subroutine rewrites m7_dconc to a large extent.
  !    It has been heavily restructured in order to eliminate redundant
  !    gridpoint loops.
  !    The use of the notion of 'equivalent molecules sulphate' of a species, as used in m7_donc, 
  !    is dispensed with, since this is, to within a multiplicative constant, (dh2so4*avo/wh2so4),
  !    the same as the vmr of the species.
  !
  !    Externals:
  !    ----------
  !    None
  !
  !--- Parameter list:
  !
  !    paerml(kbdim,klev,naerocomp)= total aerosol mass for each compound 
  !                                [molec. cm-3 for sulfate and ug m-3 for others]
  !    paernl(kbdim,klev,nclass)   = aerosol number for each mode [cm-3]  !
  !    sigmaln(jclass)             = log of standard deviation of mode jclass 
  !    crdiv                     = threshold radii between the different modes [cm]
  !                                crdiv(jclass) is the lower bound and crdiv(jclass+1) is 
  !                                the upper bound of the respective mode
  !
  !--- Local Variables:
  !
  !    zfconn(:,:,jnum,jclass)     = absolute fraction of the number of particles in mode jclass, 
  !                                with CMD=2*pm6dry(jclass) and a geometric standard 
  !                                deviation zrcsig, that are smaller than crdiv(jnum+1).
  !                                i.e. with 0 < CMD < crdiv(jnum+1) 
  !    zfconm(:,:,jnum,jclass)     = absolute fraction of the mass in mode jclass, 
  !                                with CMD=2*pm6dry(jclass) and a geometric standard 
  !                                deviation zrcsig, that are smaller than crdiv(jnum+1).
  !                                i.e. with 0 < CMD < crdiv(jnum+1) 
  
  USE mo_kind,        ONLY: dp
  USE mo_ham_m7ctl,   ONLY: sigmaln, crdiv,                &
                            inucs, iaits, iaccs, icoas,    &
                            iso4ns,iso4ks, iso4as, iso4cs, &
                            ibcks, ibcas, ibccs,           &
                            iocks, iocas, ioccs,           &
                            issas, isscs,                  &
                            iduas, iducs,                  &
                            wh2so4,                        &
                            dh2so4,dbc, doc,               &
                            dnacl, ddust,                  & 
                            cmr2ram,                       &
                            cmin_aernl, cminrad
  USE mo_physical_constants, ONLY: avo
  USE mo_math_constants,     ONLY: pi
  USE mo_ham,         ONLY: naerocomp, aerocomp, nsoaspec, nsoa, sizeclass, nclass, nsol
  !>>dod replaced isoa_ix
  USE mo_species,     ONLY: speclist
  USE mo_ham_soa,     ONLY: soaprop
  !<<dod

  IMPLICIT NONE


  !------- subroutine interface ----------------------------------------------------------------
  INTEGER, INTENT(IN)  :: kproma, kbdim, klev,                    & ! Grid parameters  
                          krow                                      ! For debugging
               
  REAL(dp), INTENT(IN) :: pm6dry(kbdim,klev,nsol)                   ! dry CMR
  REAL(dp), INTENT(INOUT) :: paerml(kbdim,klev,naerocomp)             ! Aerosol mass 
  REAL(dp), INTENT(INOUT) :: paernl(kbdim,klev,nclass)                ! Aerosol number
  !---------------------------------------------------------------------------------------------

  !------- local data --------------------------------------------------------------------------
  !
  !------- local parameters
  !
  REAL(dp), PARAMETER :: z4piover3 = 4._dp*pi/3._dp                
  REAL(dp), PARAMETER :: zfacso4m = wh2so4/avo                      ! conversion factor sulphate mass
  REAL(dp), PARAMETER :: zfacso4v = wh2so4/avo/dh2so4               ! conversion factor sulphate VMR
  REAL(dp), PARAMETER :: zfac2cgs = 1.E-12_dp                  ! conversion factor ug m-3 to g cm-3
  
  !------- local variables 
  INTEGER  :: jclass, jl, jk, jm                           ! loop counters

  REAL(dp) :: ztotmass(kbdim,klev),                   &  ! total mass mixing ratio of the mode
              ztotvmr(kbdim,klev),                    &  ! total volume mixing ration of the mode
              zavdens,                                &  ! average density of the mode 
              zttnj,                                  &  ! average particle mass
              zavnj,                                  &  ! average mass of particle with radius zdm
              zmrj,                                   &  ! mass remaining in mode
              znt,                                    &  ! number transferred between modes
              zmt                                        ! mass transferred between modes

  !>>dod deleted arrays zsumm, sumn
  !<<dod
  !csld to resolve bug burden -> re-introduction of the variables zsumm/-n
  ! in the buggy version, zsumm 
  REAL(dp) :: zsumn(kbdim,klev,nclass),                 &  ! new number concentration
              zsumm(kbdim,klev,naerocomp)                  ! new mass concentration 
 
  REAL(dp) :: zfconn(kbdim,klev,nclass), zfconm(kbdim,klev,nclass)
                                                                 
  REAL(dp) :: zdummy,                                 &  ! Unused output from m7_cumulative_normal
              zr1,                                    &  ! Mode lower boundary (cm)
              zr2,                                    &  ! Mode upper boundary (cm) 
              zmr1,                                   &  ! Mass equivalent of zr1
              zmr2,                                   &  ! Mass equivalent of zr2
              zdm,                                    &  ! Radius of average mass of the mode
              zdm3,                                   &  ! zdm**3
              zdens                                      ! density

  REAL(dp) :: zt                                         ! see comments above

  REAL(dp) :: zeps

  REAL(dp) :: zcdf1, zr23, zmrad

  INTEGER  :: ispec, iaero, isoans, isoaks, isoaas, isoacs

  !---------------------------------------------------------------------------------------------


  !------- executable procedure ----------------------------------------------------------------

  !--- 0) Initialisations: ----------------------------------------------------------

  zfconm(:,:,:) = 1._dp
  zfconn(:,:,:) = 1._dp

  zeps = EPSILON(1._dp)

  !
  !--- 1) Identify how much the mode jclass has grown into the next higher mode 
  !       All compound masses are expressed in cgs units in this subroutine
  !       For sulphate, held in molecules per cm-3 in paerml, we1 have to convert units

  DO jclass=1,nsol-1

     !--- Total volume mixing ratio the mode - assumed to be the sum of the vmr's of the 
     !    individual compounds 

     ! Start with inorganics, bc and nonvolatiles (oc)

     SELECT CASE(jclass)
     CASE(inucs)
                                                            ! Sulphate mmr in molecules cm-3 -> g cm-3
        ztotmass(1:kproma,:) = paerml(1:kproma,:,iso4ns)*zfacso4m
                                                            ! Sulphate in molecules cm-3 -> vmr
        ztotvmr(1:kproma,:) = paerml(1:kproma,:,iso4ns)*zfacso4v

     CASE(iaits)
        ztotmass(1:kproma,:) = paerml(1:kproma,:,iso4ks)*zfacso4m + &
                               zfac2cgs*(paerml(1:kproma,:,ibcks) + paerml(1:kproma,:,iocks))

        ztotvmr(1:kproma,:) = paerml(1:kproma,:,iso4ks)*zfacso4v + &
                              zfac2cgs*(paerml(1:kproma,:,ibcks)/dbc+paerml(1:kproma,:,iocks)/doc)  

     CASE(iaccs)
        ztotmass(1:kproma,:) = paerml(1:kproma,:,iso4as)*zfacso4m + &
                               zfac2cgs*(paerml(1:kproma,:,ibcas) + paerml(1:kproma,:,iocas) + &
                                          paerml(1:kproma,:,issas) + paerml(1:kproma,:,iduas) )

        ztotvmr(1:kproma,:) = paerml(1:kproma,:,iso4as)*zfacso4v + &
                              zfac2cgs*(paerml(1:kproma,:,ibcas)/dbc+paerml(1:kproma,:,iocas)/doc+  &
                                         paerml(1:kproma,:,issas)/dnacl+paerml(1:kproma,:,iduas)/ddust) 

     END SELECT

     !--- Add SOA for all modes

     IF (nsoa == 1 .AND. sizeclass(jclass)%lsoainclass) THEN
        DO jm = 1,nsoaspec                                      ! SOA in ug m-3 -> vmr
           ispec = soaprop(jm)%spid_soa
           iaero = speclist(ispec)%iaerocomp(jclass)
           zdens = 1.E-3_dp*aerocomp(iaero)%species%density
           ztotmass(1:kproma,:) = ztotmass(1:kproma,:) + zfac2cgs*paerml(1:kproma,:,iaero)
           ztotvmr(1:kproma,:) = ztotvmr(1:kproma,:) + &
                                 zfac2cgs*paerml(1:kproma,:,iaero)/zdens
        END DO
     END IF

     !--- Set minimum radius and maximum radius:

     zr1 = crdiv(jclass)     
     zr2 = crdiv(jclass+1)

     !--- Radius of average mass for a lognormal distribution
     !    @@@ TO DO: Move this to mo_ham_m7ctl 
     
     zdm = EXP((LOG(zr1)+LOG(zr2))/2.0_dp)*cmr2ram(jclass)
     zdm3 = zdm**3

     !--- For calculating average mass of particles of radius zr1 and zr2 within the gridpoint loop
         
     zr23 = zr2**3                                ! Mode upper boundary
     zmr1 = z4piover3*(zr1**3)                    ! Mode volume lower boundary

     DO jk=1,klev
        DO jl=1,kproma

           zmrad = pm6dry(jl,jk,jclass)

           IF ((paernl(jl,jk,jclass) > cmin_aernl) .AND. (zmrad > cminrad) &
                .AND. ztotvmr(jl,jk) > 0._dp) THEN

              !--- Average mass of a particle in the mode 
              zttnj = ztotmass(jl,jk)/paernl(jl,jk,jclass)

              !--- Average particle density:
              zavdens = ztotmass(jl,jk)/ztotvmr(jl,jk)
 
              !--- mass of a particle with radius zdm
              zavnj = z4piover3*zdm3*zavdens 

              !--- mass of a particle with radius zr2
              zmr2 = z4piover3*zr23*zavdens 

              !--- If the average mass contained in the mode is larger than the average mass
              !    for a particle with radius zdm, the transfer of number and mass is done,
              !    from mode jclass to mode jclass+1. Not for the coarse mode, obviously, but
              !    this is taken care of by the fact that zfconn =1 for the coarse mode

              IF (zttnj > zavnj) THEN

                 IF (jclass /= icoas) THEN 

                       
                    !--- Calculate the cumulative of the log-normal number distribution:
                    !    Fraction of particles having radius < upper mode boundary
                    !    Input to m7_cumnor is t as described in the comments at the beginning of the
                    !    subroutine. t = [log(crdiv(jnum+1)-log(pm6dry(jclass))]/sigmaln(jclass)
                    !    t is the distance of the CMD of the mode from the upper threshold mode  
                    !    diameter in terms of geometric standard deviations:
                    !    Output from m7_cumnor is zfconn(:,:,jclass), see description under 'parameters'
                    !    in the comments at the beginning of the subroutine
                    
                    zt = LOG(zr2/zmrad)/sigmaln(jclass)
                    !DIR$ INLINE
                    CALL m7_cumulative_normal(zt, zcdf1, zdummy)
                       
                    ! if mode cmr >> upper mode boundary (this can happen due to SOA condensation) 
                    ! m7_cumulative_normal can return cumulative distribution functions < zeps at the mode
                    ! boundary. Then no transfer would take place even though the particles have grown
                    ! by more than an order of magnitude. So the following Q&D is applied until someone
                    ! thinks of something better.

                    IF (zmrad > (zr2*10._dp)) THEN
                       zfconn(jl,jk,jclass) = 0._dp
                    ELSE
                       IF (zcdf1 > zeps) zfconn(jl,jk,jclass) = zcdf1
                    END IF

                    !--- Numbers transferred 
                    znt=(1.0_dp-zfconn(jl,jk,jclass))*paernl(jl,jk,jclass)
                    
                    !---mass transferred
                    IF (ztotmass(jl,jk) > 0._dp .AND. paernl(jl,jk,jclass) > 0._dp .AND. znt > 0._dp) THEN

                       !--- Mass remaining in the mode
                       zmrj = zfconn(jl,jk,jclass)*paernl(jl,jk,jclass)*zavnj
                       zmrj = MAX(zmrj, 0._dp)

                       !--- Mass transferred
                       zmt = ztotmass(jl,jk)-zmrj

                       zfconm(jl,jk,jclass)=zmrj/ztotmass(jl,jk)
                          
                       !--- If the average mass of particle transferred is smaller than the average mass
                       !    mass of particles with radius zmr2 then reduce the particles transferred
                       !    so that average mass per particle transferred is (zmr2*zavdens)
                    
                       IF ((zmt/znt) < zmr2 ) THEN
                          zfconn(jl,jk,jclass) = 1.0_dp - &
                                                 (zmt/zmr2)/paernl(jl,jk,jclass)
                       END IF

                    ELSE
                       zfconn(jl,jk,jclass) = 1._dp
                       zfconm(jl,jk,jclass) = 1._dp
                    END IF

                 END IF
                                
              END IF                   !--- end IF avg mass =< mass for radius zdm
           END IF                   !--- end IF enough particles

        END DO
     END DO
  END DO

  ! ---- number concentrations -------------------------------------------------------------------------------

  ! Nucleation mode: particles that stay in nucl. mode
                       
  zsumn(1:kproma,:,inucs) = paernl(1:kproma,:,inucs)*zfconn(1:kproma,:,inucs)   

  ! Aitken mode: particles that stay in aitken mode (do not transfer into accum. mode) plus
  !              particles transferred from the nucl. mode (considering upwards transfer)

  zsumn(1:kproma,:,iaits) = paernl(1:kproma,:,iaits)*zfconn(1:kproma,:,iaits) +        &
                             paernl(1:kproma,:,inucs)*(1._dp-zfconn(1:kproma,:,inucs)) 

  ! Accumulation mode: particles that stay in accum. mode plus:
  !                    particles transferred from the aitken soluble mode 

  zsumn(1:kproma,:,iaccs) = paernl(1:kproma,:,iaccs)*zfconn(1:kproma,:,iaccs) +        &
                             paernl(1:kproma,:,iaits)*(1._dp-zfconn(1:kproma,:,iaits)) 

  ! Coarse mode: particles that stay in coarse mode plus particles transferred from accum. mode
                       
  zsumn(1:kproma,:,icoas) = paernl(1:kproma,:,icoas) +  &  
                             paernl(1:kproma,:,iaccs)*(1._dp-zfconn(1:kproma,:,iaccs))

!csld bug burden
  ! Store the final number concentration
  paernl(1:kproma,:,inucs) = zsumn(1:kproma,:,inucs)
  paernl(1:kproma,:,iaits) = zsumn(1:kproma,:,iaits)
  paernl(1:kproma,:,iaccs) = zsumn(1:kproma,:,iaccs)
  paernl(1:kproma,:,icoas) = zsumn(1:kproma,:,icoas)
  !end csld bug burden
  ! ---- end number concentrations ---------------------------------------------------------------------------

  ! ---- sulphate mass ---------------------------------------------------------------------------------------

  ! Nucleation mode: mass that stays in nucl. mode 

  zsumm(1:kproma,:,iso4ns) = paerml(1:kproma,:,iso4ns)*zfconm(1:kproma,:,inucs)   

  ! Aitken mode: mass that stays in aitken mode plus mass transferred from the nucl. mode 

  zsumm(1:kproma,:,iso4ks) = paerml(1:kproma,:,iso4ks)*zfconm(1:kproma,:,iaits) + &
                              paerml(1:kproma,:,iso4ns)*(1._dp-zfconm(1:kproma,:,inucs)) 

  ! Accumulation mode: mass that stays in accum. mode plus mass transferred from the aitken mode

  zsumm(1:kproma,:,iso4as) = paerml(1:kproma,:,iso4as)*zfconm(1:kproma,:,iaccs) + &
                              paerml(1:kproma,:,iso4ks)*(1._dp-zfconm(1:kproma,:,iaits)) 
                       
  ! Coarse mode: mass that stays in coarse mode plus mass transferred from accum. mode
                       
  zsumm(1:kproma,:,iso4cs) = paerml(1:kproma,:,iso4cs) +                            &  
                              paerml(1:kproma,:,iso4as)*(1._dp-zfconm(1:kproma,:,iaccs))

 !csld bug burden
  ! Store the final sulphate mass concentrations
  paerml(1:kproma,:,iso4ns) = zsumm(1:kproma,:,iso4ns)
  paerml(1:kproma,:,iso4ks) = zsumm(1:kproma,:,iso4ks)
  paerml(1:kproma,:,iso4as) = zsumm(1:kproma,:,iso4as)
  paerml(1:kproma,:,iso4cs) = zsumm(1:kproma,:,iso4cs)
  !end csld
  ! ---- end sulphate mass -----------------------------------------------------------------------------------

  ! ---- SOA mass --------------------------------------------------------------------------------------------

  IF (nsoa == 1) THEN
     DO jm = 1,nsoaspec        
        ispec = soaprop(jm)%spid_soa

        ! Nucleation mode: mass that stays in nucl. mode
        IF (sizeclass(inucs)%lsoainclass) THEN
           isoans = speclist(ispec)%iaerocomp(inucs)
           zsumm(1:kproma,:,isoans) = paerml(1:kproma,:,isoans)*zfconm(1:kproma,:,inucs)     
        END IF
    
        ! Aitken mode: mass that stays in aitken mode plus mass transferred from the nucl. mode
        isoaks = speclist(ispec)%iaerocomp(iaits)

        zsumm(1:kproma,:,isoaks) = paerml(1:kproma,:,isoaks)*zfconm(1:kproma,:,iaits) 

        IF (sizeclass(inucs)%lsoainclass) THEN
           zsumm(1:kproma,:,isoaks) = zsumm(1:kproma,:,isoaks) + &
                                       paerml(1:kproma,:,isoans)*(1._dp-zfconm(1:kproma,:,inucs)) 
        END IF

        ! Accumulation mode: mass that stays in accum. mode plus:
        !                    mass transferred from the aitken mode 
        isoaas = speclist(ispec)%iaerocomp(iaccs)

        zsumm(1:kproma,:,isoaas) = paerml(1:kproma,:,isoaas)*zfconm(1:kproma,:,iaccs) + &
                                    paerml(1:kproma,:,isoaks)*(1._dp-zfconm(1:kproma,:,iaits)) 
                       
        ! Coarse mode: mass that stays in coarse mode plus mass transferred from accum. mode
        isoacs = speclist(ispec)%iaerocomp (icoas)

        zsumm(1:kproma,:,isoacs) = paerml(1:kproma,:,isoacs) +  &  
                                    paerml(1:kproma,:,isoaas)*(1._dp-zfconm(1:kproma,:,iaccs))

!csld bug burden
        IF (sizeclass(inucs)%lsoainclass) THEN
           isoans = speclist(ispec)%iaerocomp(inucs)
           paerml(1:kproma,:,isoans)=zsumm(1:kproma,:,isoans)
        ENDIF
        paerml(1:kproma,:,isoaks)=zsumm(1:kproma,:,isoaks) 
        paerml(1:kproma,:,isoaas)=zsumm(1:kproma,:,isoaas)
        paerml(1:kproma,:,isoacs)=zsumm(1:kproma,:,isoacs)
!end csld
     END DO
  END IF

  ! ---- end SOA mass ----------------------------------------------------------------------------------------
  
  ! ---- Black carbon (bc) and nonvolatile organic carbon (oc) mass ------------------------------------------

  ! Aitken mode: mass that stays in aitken mode 
                       
  zsumm(1:kproma,:,ibcks) = paerml(1:kproma,:,ibcks)*zfconm(1:kproma,:,iaits)

  zsumm(1:kproma,:,iocks) = paerml(1:kproma,:,iocks)*zfconm(1:kproma,:,iaits)     

  ! Accumulation mode: mass that stays in accum. mode plus: 
  !                    mass transferred from the aitken mode 

  zsumm(1:kproma,:,ibcas) = paerml(1:kproma,:,ibcas)*zfconm(1:kproma,:,iaccs) + &
                             paerml(1:kproma,:,ibcks)*(1._dp-zfconm(1:kproma,:,iaits)) 

  zsumm(1:kproma,:,iocas) = paerml(1:kproma,:,iocas)*zfconm(1:kproma,:,iaccs) + &
                             paerml(1:kproma,:,iocks)*(1._dp-zfconm(1:kproma,:,iaits)) 
                                                                     

  ! Coarse mode: mass that stays in coarse mode plus mass transferred from accum. mode
                       
  zsumm(1:kproma,:,ibccs) = paerml(1:kproma,:,ibccs) +  &
                             paerml(1:kproma,:,ibcas)*(1._dp-zfconm(1:kproma,:,iaccs))

  zsumm(1:kproma,:,ioccs) = paerml(1:kproma,:,ioccs) + &
                             paerml(1:kproma,:,iocas)*(1._dp-zfconm(1:kproma,:,iaccs))
 
!csld bug burden
   ! Store the final mass concentrations
  paerml(1:kproma,:,ibcks) = zsumm(1:kproma,:,ibcks)
  paerml(1:kproma,:,iocks) = zsumm(1:kproma,:,iocks)  
  paerml(1:kproma,:,ibcas) = zsumm(1:kproma,:,ibcas)
  paerml(1:kproma,:,iocas) = zsumm(1:kproma,:,iocas)  
  paerml(1:kproma,:,ibccs) = zsumm(1:kproma,:,ibccs)
  paerml(1:kproma,:,ioccs) = zsumm(1:kproma,:,ioccs) 
  !end csld
  ! ---- end black carbon (bc) and nonvolatile organic carbon (oc) mass --------------------------------------


  ! ---- sea salt and dust mass ------------------------------------------------------------------------------

  ! Accumulation mode: mass that stays in accum. mode 
                       
  zsumm(1:kproma,:,issas) = paerml(1:kproma,:,issas)*zfconm(1:kproma,:,iaccs)     

  zsumm(1:kproma,:,iduas) = paerml(1:kproma,:,iduas)*zfconm(1:kproma,:,iaccs)     

  ! Coarse mode: mass that stays in coarse mode plus mass transferred from accum. mode
                       
  zsumm(1:kproma,:,isscs) = paerml(1:kproma,:,isscs) + &
                             paerml(1:kproma,:,issas)*(1._dp-zfconm(1:kproma,:,iaccs))

  zsumm(1:kproma,:,iducs) = paerml(1:kproma,:,iducs) + &
                             paerml(1:kproma,:,iduas)*(1._dp-zfconm(1:kproma,:,iaccs))

!csld bug burden
  ! ---- store the mass concentrations
  paerml(1:kproma,:,issas) = zsumm(1:kproma,:,issas)
  paerml(1:kproma,:,iduas) = zsumm(1:kproma,:,iduas)  
  paerml(1:kproma,:,isscs) = zsumm(1:kproma,:,isscs)
  paerml(1:kproma,:,iducs) = zsumm(1:kproma,:,iducs)  
  !end csld
  ! ---- end sea salt and dust mass --------------------------------------------------------------------------

END SUBROUTINE m7_dconc


SUBROUTINE m7_coaset(kproma, kbdim, klev,  krow, paernl, ptp1, &
                     papp1,  pm6rp, prhop, pcom          )

  !
  ! *m7_coaset*  calculates the coagulation kernels between the modes
  !
  ! Authors:
  ! ---------
  ! J. Wilson  and E. Vignati, JRC/EI (original source)                09/2000
  ! P. Stier, MPI                     (f90-version, changes, comments)    2001 
  !
  ! Modifications:
  ! --------------
  ! Philip Stier, MPI                             2001
  ! Declan O'Donnell, MPI-M, 2009: performance improvements
  !
  ! Purpose
  ! ---------
  ! This routine calculates the coaglation kernels between particles
  ! with the count median radii of the three modes.
  ! Coagulation allowed between the following modes:
  ! soluble modes:   1+1=1, 2+2=2, 1+2=2, 1+3=3, 1+4=4, 2+3=3, 2+4=4
  ! insoluble modes: 2i+2i=2i
  ! mixed modes:     1+2i=2, 1+4i=4, 2+2i=2, 2+4i=4, 3+2i=3, 4+2i=4
  !
  ! Interface:
  ! -----------
  ! *m7_coaset* is called from *m7_dnum*
  !
  ! Externals:
  ! -----------
  ! none
  !
  ! Reference:
  ! -----------
  ! The calculations are based on:
  ! Fuchs, N.A. (1964). The Mechanics of Aerosols. Pergamon Press. Oxford. 
  ! (Chapter VII, Section 49)
  !
  !  Warning:
  !  --------
  !  For optimization purposes currently only "physically reasonable" elements of the
  !  coagulation kernel pcom are calculated in m7_concoag. These elements are specified
  !  in the matrix locoagmask in mo_ham_m7ctl. Check carefully and adapt locoagmask 
  !  accordingly  before changes in the code below.
  !
  USE mo_ham,         ONLY: nclass
  USE mo_kind,        ONLY: dp
  USE mo_ham_m7ctl,   ONLY: bk, ncoag, coag_modes, &
                            cmin_aernl, cminrad, cminrho
  USE mo_math_constants, ONLY: pi, sqrt2
  !

  IMPLICIT NONE
  !
  !--- Parameter list:
  !
  !  paernl            = aerosol number for each mode [cm-3]
  !  ptp1              = atmospheric temperature at time t+1 [K]
  !  papp1             = atmospheric pressure at time t+1 [Pa]
  !  pm6rp             = mean mode actual radius (wet radius for soluble modes 
  !                      and dry radius for insoluble modes) [cm]
  !  prhop             = mean mode particle density [g cm-3]
  !  pcom(:,:,jm) = Coagulation coefficient 
  !
  !--- List of local variables:
  !
  ! zwlc              = mean free pathlength []
  ! zairvisc          = air viscosity []
  ! zrpav             = average radius of the two interacting modes
  ! zpvx              = volume of the xth interacting mode
  ! zpmx              = mass of the xth interacting mode
  !
  ! zrknudx           = knudsen number of the xth interacting mode
  ! zpd2x             = particle diffusion of the xth interacting mode
  !
  !---subroutine interface 
  INTEGER,  INTENT(IN) :: kproma, kbdim, klev, krow
  REAL(dp), INTENT(IN) :: ptp1(kbdim,klev)
  REAL(dp), INTENT(IN) :: papp1(kbdim,klev)
  REAL(dp), INTENT(IN) :: pm6rp(kbdim,klev,nclass)
  REAL(dp), INTENT(IN) :: prhop(kbdim,klev,nclass)
  REAL(dp), INTENT(IN) :: paernl(kbdim,klev,nclass)
  REAL(dp), INTENT(OUT) :: pcom(kbdim,klev,ncoag)

  !--- Local variables:
  !    Parameters:
  REAL(dp), PARAMETER :: zbkover6pi = bk / (6._dp*pi)
  REAL(dp), PARAMETER :: z4piover3 = 4._dp*pi/3._dp
  REAL(dp), PARAMETER :: z16pi = 16._dp*pi
  REAL(dp), PARAMETER :: z8overpi = 8._dp/pi
  REAL(dp), PARAMETER :: z8bkoverpi = 8._dp*bk/pi
  REAL(dp), PARAMETER :: qsqrt2 = 1._dp / sqrt2

  !    Local variables
  INTEGER :: jm, jm2, jm1, jl, jk, jclass
  !
  REAL(dp) :: zbtketc(kbdim,klev)
  REAL(dp) :: zwlc(kbdim,klev)
  REAL(dp) :: zpd(kbdim,klev,nclass)
  REAL(dp) :: zcv2(kbdim,klev,nclass)
  REAL(dp) :: zh2(kbdim,klev,nclass)

  LOGICAL  :: lsuffaero(kbdim,klev,nclass)

  REAL(dp):: zairvisc,    zeps,      &
             zrpvm,       z2rpvm,      zrpav,       zpv1,      &
             zpm1,       &
             zcv2av,      zrknud,      &
             ze,         &
             zpdav,       zxd,         &
             zh2av,       zcoc,        zhu1,      &
             zhu2,        zhu

  !---executable procedure
  zeps = EPSILON(1._dp)

  !---1) gridpoint properties
  DO jk=1,klev
     DO jl=1,kproma

        !--- Mean free pathlength ? (from Knudsen Number below):
        !    Parametrisation?
        !    Grid point calculation, can be done once rather than 16 times.
        ! old:
        !  zpbyone=1000.0_dp / (papp1(jl,jk)/100.0_dp)
        !  zwlc(jl,jk)=6.6e-6_dp * ptp1(jl,jk) / 293.15_dp * zpbyone
        ! new:
        zwlc(jl,jk) = 2.251E-3_dp * ptp1(jl,jk) / papp1(jl,jk)

        !--- Viscosity:        (291.15..?)
        zairvisc=1.827e-4_dp * (ptp1(jl,jk) / 291.15_dp)**0.74_dp

        !--- 
        zbtketc(jl,jk) = zbkover6pi * ptp1(jl,jk) / zairvisc

     END DO
  END DO

  !---2) Per-mode properties
  DO jclass=1,nclass
     DO jk=1,klev
        DO jl=1,kproma

           lsuffaero(jl,jk,jclass) = (paernl(jl,jk,jclass) > cmin_aernl .AND. &
                                      pm6rp(jl,jk,jclass)  > cminrad .AND. &
                                      prhop(jl,jk,jclass)  > cminrho)

           IF (lsuffaero(jl,jk,jclass)) THEN

              !---mode radius
              zrpvm = pm6rp(jl,jk,jclass)

              !--- Volume and mass of mode:
              zpv1=z4piover3 * zrpvm**3
              zpm1=zpv1 * prhop(jl,jk,jclass)

              !--- Squared mean particle velocity of mode:
              zcv2(jl,jk,jclass)=z8bkoverpi *  ptp1(jl,jk) /  zpm1

              !---Knudsen number of the mode:
              zrknud=0.5_dp*zwlc(jl,jk)/zrpvm

              !---Diffusivity of the mode:
              ze=EXP(-0.43_dp/zrknud)
              zpd(jl,jk,jclass)=zbtketc(jl,jk) * (1.0_dp + zrknud*2.492_dp + zrknud*0.84_dp*ze) / zrpvm

              !--- Average mean free path of particles in the mode:
              zxd=z8overpi * zpd(jl,jk,jclass) / SQRT(zcv2(jl,jk,jclass))

              !--- Mean distance from surface after mean free path (Eq. 49.13):
              z2rpvm = 2._dp*zrpvm
              zh2(jl,jk,jclass)=(((z2rpvm + zxd)**3 -                        &
                                SQRT((z2rpvm*z2rpvm + zxd*zxd)**3) ) /             &
                               (6.0_dp*zrpvm*zxd) - z2rpvm           ) * sqrt2

           ELSE
              zpd(jl,jk,jclass) = 0._dp
              zcv2(jl,jk,jclass) = 0._dp
              zh2(jl,jk,jclass) = 0._dp
           END IF

        END DO
     END DO
  END DO

  !--- 1) Calculation of the coagulation coefficient: ---------------------------
  !
  DO jm=1,ncoag

     jm1 = coag_modes(jm)%mode1
     jm2 = coag_modes(jm)%mode2

     DO jk=1,klev
        DO jl=1,kproma

           IF (lsuffaero(jl,jk,jm1) .AND. lsuffaero(jl,jk,jm2)) THEN
              !--- Average radius of the modes:
              zrpav=zeps+0.5_dp*(pm6rp(jl,jk,jm1)+pm6rp(jl,jk,jm2)) 

              !--- Fuchs: G_r (below Eq. 49.27):
              zcv2av=SQRT(zcv2(jl,jk,jm1) + zcv2(jl,jk,jm2))

              !--- Average diffusivity of the modes:
              zpdav=0.5_dp*(zpd(jl,jk,jm1) + zpd(jl,jk,jm2)) 
                    
              !--- Fuchs: delta_r !@@@ (why division by sqrt2?)
              !>>dod deleted the division by sqrt(2) (redmine #58)
              zh2av=SQRT(zh2(jl,jk,jm1)*zh2(jl,jk,jm1) + zh2(jl,jk,jm2)*zh2(jl,jk,jm2)) 
              !<<dod
                 
              !--- 1.2) Calculation of the coagulation coefficient pcom (Eq. 49.26):
              !         Factor 16 instead factor 8 as in Fuchs as his formulation
              !         applies for the inter-modal coagulation. This is taken into
              !         account in the assignment of the inter-modal coagulation
              !         coefficient.

              zcoc=z16pi * zpdav * zrpav

              !--- Calculation of beta=1/zhu (Eq. 49.27):
              zhu1=4.0_dp * zpdav / (zcv2av * zrpav)
              zhu2=zrpav / (zrpav + 0.5_dp*zh2av )
              zhu=zhu1 +  zhu2

              !--- Coagulation coefficient following (Eq.49.26):
              pcom(jl,jk,jm)=zcoc / zhu

           ELSE
              pcom(jl,jk,jm)  = 0._dp
           END IF
        END DO
     END DO


  END DO

END SUBROUTINE m7_coaset


SUBROUTINE m7_concoag (kproma,   kbdim,   klev, krow,               &
                       paerml,   paernl,  pm6rp,  pa4delt, panli,   &
!>>dod soa added pttn, removed bc, oc, dust averages
                       pa4av1,   pa4av2,  pttn,                     &
                       pso4_5,  pso4_6,   pso4_7,                   &
!<<dod
                       pbfract1, pbfract2                            )
  !
  !   *m7_concoag*
  !
  !   Author:
  !   ----------
  !   E. Vignati, JRC/EI     (original source)                09/2000
  !   P. Stier, MPI          (f90-version, changes, comments)    2001 

  !   Version:
  !   ----------
  !   This version is equivalent to the version concoa_n of the boxmodel. 
  !
  !   Purpose
  !   ----------
  !   m7_concoag transfers aerosol mass and numbers from the insoluble
  !   to the soluble modes.
  !
  !   Interface:
  !   ----------
  !   *m7_concoag* is called from *m7_delcoa*
  !
  !   Externals
  !   ----------
  !   none

  USE mo_kind,      ONLY: dp
  USE mo_ham_m7ctl, ONLY: ibcks,   ibcki,   iocks,  iocki,   &
                          iduas,   iducs,   iduai,  iduci,   &
                          iso4as,  iso4cs,                   &
                          iaiti,   iacci,   icoai,           &
                          iaits,   iaccs,   icoas, cmin_aernl
  !>>dod soa
  USE mo_ham,       ONLY: nsoa, nsoaspec, naerocomp, nclass
  USE mo_ham_soa,   ONLY: soaprop
  !>>dod deleted isoa_ix
  USE mo_species,   ONLY: speclist
  !<<dod
  !<<dod
  
  IMPLICIT NONE 

  !--- Parameters:
  !
  ! paerml          = total aerosol mass for each compound 
  !                   [molec. cm-3 for sulphate and ug m-3 for bc, oc, ss, and dust]
  ! paernl          = aerosol number for each mode [cm-3]
  ! pm6rp           = mean mode actual radius (wet radius for soluble modes 
  !                   and dry radius for insoluble modes) [cm]
  ! pa4delt(:,:,:)  = change in H2SO4 mass of the respective mode over one timstep 
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in m7_nuck)
  !                      - coagulation (calculated here in m7_concoag)
  ! pxxavy          = average mass of species xx in mode y []!@@@
  !                   where xx is ss, du, bc, oc, or a4 for sulfate
  ! panli(:,:,x)    = total number of particles moved by inter-modal 
  !                   coagulation from mode x [cm-3]
  ! pbfractx(:,:,y) = fraction of the total number of particles removed by 
  !                   coagulation from mode x that is moved to mode y+1 [1] 
  !                   !@@@ Clumsy notation! Should be moved to mode y !!!
  ! pso4_x          = mass of sulphate condensed on insoluble mode x [molec. cm-3]
  !
  !--- Local variables / Constants:
  !
  ! zso4x    = available mass of sulfate from mode 1 and 2 
  !            condensing and coagulating on mode x (x = insoluble modes 5,6,7).  
  !
  ! zcrtcst  = Critical constant, i.e. number of sulfate molecules to cover 
  !            an average particle of the mode with a layer of the thickness
  !            determined by cLayerThickness in mo_ham_m7ctl.   Calculated by
  !            m7_coat.
  !
  ! =>         zso4x/zcrtcst is the total number of particles that could be moved
  !            from insoluble mode x to soluble modes.
  !
  ! zcrit_x  = total available number of particles in mode x that are moved from 
  !            insoluble mode x to the corresponding soluble mode.
  !>>dod soa added krow
  INTEGER :: kproma, kbdim, klev, krow
  !<<dod

  REAL(dp):: pso4_5(kbdim,klev),          pso4_6(kbdim,klev),                   &
             pso4_7(kbdim,klev),                                                &
             !>>dod soa removed bc, oc and du averages
             pa4av1(kbdim,klev),          pa4av2(kbdim,klev)
             !<<dod

  REAL(dp):: paerml(kbdim,klev,naerocomp),  paernl(kbdim,klev,nclass),              &
             pbfract1(kbdim,klev,nclass-1), pbfract2(kbdim,klev,nclass-1),          &
             panli(kbdim,klev,nclass),      pa4delt(kbdim,klev,naerocomp),          &
             pm6rp(kbdim,klev,nclass)

  !>>dod soa 
  REAL(dp) :: pttn(kbdim,klev,naerocomp)
  !<<dod

  ! Local variables:

  INTEGER :: jl,          jk,            jclass   
  !>>dod soa
  INTEGER :: jm, jn, jn2, jspec
  !<<dod

  REAL(dp):: zcrit_5,     zcrit_6,       zcrit_7,                               &
             zso45,       zso46,         zso47,                                 &
             zeps

  REAL(dp):: zm6rp(nclass), zcrtcst(nclass)

  !--- 0) Initializations:

  zeps=EPSILON(1._dp)


  !--- 1) Redistribution of mass and numbers after nucleation, coagulation ----
  !       and coagulation calculated in the preceeding subroutines:

  DO jk=1,klev
     DO jl=1,kproma

        !--- 1.1) Sum mass of sulphate added to modes 5, 6, and 7 due to 
        !         coagulation with modes 1 and 2 (1st term) and the mass
        !         of sulfate condensed on the insoluble mode x (pso4_x):
        
        zso45=panli(jl,jk,1)*pbfract1(jl,jk,4)*pa4av1(jl,jk)+pso4_5(jl,jk)

        zso46=panli(jl,jk,1)*pbfract1(jl,jk,5)*pa4av1(jl,jk)+                &
              panli(jl,jk,2)*pbfract2(jl,jk,5)*pa4av2(jl,jk)+pso4_6(jl,jk)

        zso47=panli(jl,jk,1)*pbfract1(jl,jk,6)*pa4av1(jl,jk)+                &
              panli(jl,jk,2)*pbfract2(jl,jk,6)*pa4av2(jl,jk)+pso4_7(jl,jk)

        !--- 1.2) Determine number of particles that can be sufficiently coated
        !         by the available sulfate to be transfered to the soluble modes:

        !    Optimization of the call of m7_coat to allow for unroll and 
        !    subsequent vectorization.
!DIR$ UNROLL
!CDIR UNROLL=7
        DO jclass = 1, nclass
          zm6rp(jclass) = pm6rp(jl,jk,jclass)
        END DO

        CALL m7_coat(zm6rp,zcrtcst)

        !@@@ Changed security check to allow for inconsistent radii:

        IF(paernl(jl,jk,iaiti) >= cmin_aernl .AND. zcrtcst(5)>zeps) THEN
           zcrit_5=MIN(paernl(jl,jk,iaiti), zso45/zcrtcst(5))
        ELSE
           zcrit_5=0._dp
        END IF
        IF(paernl(jl,jk,iacci) >= cmin_aernl .AND. zcrtcst(6)>zeps) THEN
           zcrit_6=MIN(paernl(jl,jk,iacci), zso46/zcrtcst(6))
        ELSE
           zcrit_6=0._dp
        END IF
        IF(paernl(jl,jk,icoai) >= cmin_aernl .AND. zcrtcst(7)>zeps) THEN
           zcrit_7=MIN(paernl(jl,jk,icoai), zso47/zcrtcst(7))
        ELSE
           zcrit_7=0._dp
        END IF

        !--- 1.3) Number of particles moved from the mode 5 to 2 due to
        !         interaction with 1 and due to condensation:
        
        paernl(jl,jk,iaits)=paernl(jl,jk,iaits)+zcrit_5
        paernl(jl,jk,iaiti)=paernl(jl,jk,iaiti)-zcrit_5
        
        !--- 1.4) Mass moved from mode 5 to 2:
        
        pa4delt(jl,jk,2)=pa4delt(jl,jk,2)+pso4_5(jl,jk)
        pa4delt(jl,jk,ibcks)=pa4delt(jl,jk,ibcks)+zcrit_5 * 1.E12_dp * pttn(jl,jk,ibcki)
        pa4delt(jl,jk,iocks)=pa4delt(jl,jk,iocks)+zcrit_5 * 1.E12_dp * pttn(jl,jk,iocki)

        !>>dod soa
        IF (nsoa ==1) THEN
           DO jm=1,nsoaspec
              IF (.NOT. soaprop(jm)%lvolatile) THEN
                 !>>dod deleted isoa_ix
                 jspec = soaprop(jm)%spid_soa
                 jn = speclist(jspec)%iaerocomp(iaits)
                 jn2 = speclist(jspec)%iaerocomp(iaiti)
                 !<<dod
                 pa4delt(jl,jk,jn) = pa4delt(jl,jk,jn)+zcrit_5* 1.E12_dp * pttn(jl,jk,jn2)
              END IF
           END DO
        END IF
        !<<dod

        !--- 1.5) Mass remaining in mode 5:
        
        paerml(jl,jk,ibcki)=paerml(jl,jk,ibcki)-zcrit_5 * 1.E12_dp * pttn(jl,jk,ibcki)
        paerml(jl,jk,iocki)=paerml(jl,jk,iocki)-zcrit_5 * 1.E12_dp * pttn(jl,jk,iocki)

        !>>dod soa
        IF (nsoa == 1) THEN
           DO jm=1,nsoaspec
              IF (.NOT. soaprop(jm)%lvolatile) THEN
                 !>>dod deleted isoa_ix
                 jspec = soaprop(jm)%spid_soa                 
                 jn = speclist(jspec)%iaerocomp(iaiti)
                 !<<dod
                 paerml(jl,jk,jn) = paerml(jl,jk,jn)-zcrit_5 * 1.E12_dp * pttn(jl,jk,jn)
              END IF
           END DO
        END IF
        !<<dod

        !--- 1.6) Number of particles moved from the mode 6 to 3:
        
        paernl(jl,jk,iaccs)=paernl(jl,jk,iaccs)+zcrit_6
        paernl(jl,jk,iacci)=paernl(jl,jk,iacci)-zcrit_6
        
        !--- 1.7) Mass moved from mode 6 to 3:
        !>>dod soa
        pa4delt(jl,jk,iso4as)=pa4delt(jl,jk,iso4as)+pso4_6(jl,jk)
        pa4delt(jl,jk,iduas)=pa4delt(jl,jk,iduas)+zcrit_6 * 1.E12_dp * pttn(jl,jk,iduai)
        !<<dod

        !--- 1.8) Mass remaining in mode 6:
        !>>dod soa
        paerml(jl,jk,iduai)=paerml(jl,jk,iduai)-zcrit_6 * 1.E12_dp * pttn(jl,jk,iduai)
        !<<dod

        !--- 1.9) Number of particles moved from the mode 7 to 4:

        paernl(jl,jk,icoas)=paernl(jl,jk,icoas)+zcrit_7
        paernl(jl,jk,icoai)=paernl(jl,jk,icoai)-zcrit_7

        !--- 1.10) Mass moved from mode 7 to 4:
        !>>dod soa
        pa4delt(jl,jk,iso4cs)=pa4delt(jl,jk,iso4cs)+pso4_7(jl,jk)
        pa4delt(jl,jk,iducs)=pa4delt(jl,jk,iducs)+zcrit_7 * 1.E12_dp * pttn(jl,jk,iduci)

        !--- 1.11) Mass remaining in mode 7:

        paerml(jl,jk,iduci)=paerml(jl,jk,iduci)-zcrit_7 * 1.E12_dp * pttn(jl,jk,iduci)

     END DO
  END DO

END SUBROUTINE m7_concoag


  SUBROUTINE m7_coat(pm6rp_lon_lev, pcrtcst)
        
    ! Purpose:
    ! ---------
    ! *m7_coat* calculates the number of sulfate 
    !           molecules required to coat a particle
    !           with cLayerThickness of sulfate
    !
    ! Author:
    ! ---------
    ! Philip Stier, MPI                          2001
    !
    ! Interface:
    ! ---------
    ! *m7_coat* is called from *m7_concoag*
    !

    USE mo_ham,          ONLY: nclass
    USE mo_kind,         ONLY: dp
    USE mo_math_constants, ONLY: pi
    USE mo_ham_m7ctl,    ONLY: so4_coating_threshold,        &
                               cmr2ras

    IMPLICIT NONE

    INTEGER   :: jclass 

    REAL(dp)  :: pm6rp_lon_lev(nclass) ! Ambient radii for current
                                       ! longitude and level [cm]
    REAL(dp)  :: pcrtcst(nclass)       ! Critical constant, i.e. number of
                                       ! sulfate to cover an average particle
                                       ! of the mode with a layer of the 
                                       ! thickness determined by cLayerThickness.
    REAL(dp)  :: zras(nclass)          ! Radius of average surface
                                       ! for a single particle [cm]
    REAL(dp)  :: zas(nclass)           ! Average surface
                                       ! for single particle [cm+2]
        
    REAL(dp), PARAMETER :: csurf_molec = 2.39E-15_dp ! Average cross-section 
                                              ! of a single H2SO4 molecule [cm+2]

    !--- 1) Calculate the radii of average surface for modes 5-7:

    zras(5) = pm6rp_lon_lev(5) * cmr2ras(5)
    zras(6) = pm6rp_lon_lev(6) * cmr2ras(6)
    zras(7) = pm6rp_lon_lev(7) * cmr2ras(7)

    DO jclass=5, 7
           
       !--- 2) Calculate the average surface of an particle for modes 5-7:

       zas(jclass)    = 4._dp * zras(jclass)**2 * pi

       !--- 3) Determine the number of sulfate molecules needed to form
       !       a cLayerThickness thick layer of sulfate on the particles
       !       in modes 5-7:
           
       pcrtcst(jclass) = (zas(jclass) / csurf_molec) * so4_coating_threshold

    END DO

  END SUBROUTINE m7_coat
 

 
  SUBROUTINE m7_delcoa(kproma,   kbdim,  klev,     krow, paerml,   &
                       !>>dod soa (added pttn)
                       paernl,   pm6rp,  pttn,     pa4delt,  panew,    &
                       !<<dod
                       pa,       pb,     pbfract1, pbfract2, &
                       pbfract5, pso4_5, pso4_6,   pso4_7    )
    ! 
    !    Authors: 
    !    --------- 
    !    E. Vignati and J. Wilson, JRC/EI (original source)                09/2000
    !    P. Stier, MPI                    (f90-version, changes, comments)    2001 
    !    D. O'Donnell, MPI-M, 2007-2008 
    ! 
    !    Version/History: 
    !    ----------------
    !    equivalent to the version delco_n2 of the m7 boxmodel 
    !    + use of the analytical solution 
    ! 
    !    Adapted to calculate mass transfer of non-volatile SOA.
    !    Also deleted calculation of average compound mass per mode, since this is already
    !    done in m7_averageproperties.
    
    !    Purpose 
    !    --------- 
    !    This routine calculates changes in number concentration of 
    !    each aerosol mode over the time step, due to coagulation with 
    !    the current mode and all higher ones. 
    ! 
    !    Method: 
    !    ----------- 
    !    *delcoa*  integrates for each mode dn/dt=c -a*n^2 -b*n  over ztmst  
    ! 
    !    The resulting particles are assumed to reside in the 
    !    mode of highest mode of the pair of particles colliding. 
    !    1+1=>1 1+2=>2 1+3=>3, 2+2=>2 2+3=>3, 3+3=>3. 
    !    zc is now non zero for mode 1 only (nucleation). 
    !    All formation of higher mode particles is handled in dconc. 
    ! 
    !    For climatological studies, 5 day accumulation mode concs are  
    !    within a factor of 2 of the full model.  
    ! 
    !    Interface 
    !    ----------- 
    !    *m7_delcoa* is called from *m7_dnum* 
    ! 
    !    Externals 
    !    ----------- 
    !    none 
    ! 
 
    USE mo_kind,          ONLY: dp
    USE mo_time_control,  ONLY: time_step_len
    USE mo_ham_m7ctl,     ONLY: cmin_aernl, cmin_aerml,          & 
                                iaiti,   ibcki,   iocki,   ibcks,      & 
                                iocks,   ibcas,   iocas,   ibccs,      &
                                ioccs,   iacci,                        &
                                !>>dod soa                              
                                iso4ns,  iso4ks,  iso4as,  iso4cs,     &
                                icoai,   iaits,   inucs,   iaccs,      &
                                icoas
                                !<<dod
    USE mo_ham,           ONLY: nsoa, nsoaspec, naerocomp, sizeclass, nclass,nsol
    !>>dod deleted isoa_ix
    USE mo_ham_soa,       ONLY: soaprop
    !<<dod
    USE mo_species,       ONLY: speclist
    USE mo_ham_species,   ONLY: id_so4, id_bc, id_oc, id_du
    !<<dod
    !>>dod timers
    USE mo_control,       ONLY: ltimer
    USE mo_hammoz_timer,  ONLY: timer_start, timer_stop, timer_ham_m7_concoag

    IMPLICIT NONE 

  !--- Parameter list:
  !
  ! paerml          = total aerosol mass for each compound 
  !                   [molec. cm-3 for sulphate; ug m-3 for others]
  ! paernl          = aerosol number for each mode [cm-3]
  ! pm6rp           = mean mode actual radius (wet radius for soluble modes 
  !                   and dry radius for insoluble modes) [cm]
  ! pa4delt(:,:,:)  = change in H2SO4 mass of the respective mode over one timstep 
  !                   due to:
  !                      - nucleation of H2SO4 (calculated in m7_nuck)
  !                      - coagulation (calculated in m7_concoag)
  ! panew           = number of nucleated particles (during 1 timestep) [1] 
  ! pa              = unimodal coagulation coefficient (zcom(mod)) []
  ! pb              = inter-modal coagulation with higher modes
  !                   (zcom(mod) * n(jclass+1)) []
  ! pbfractx(:,:,y) = fraction of the total number of particles removed by 
  !                   coagulation from mode x that is moved to mode y+1 [1]
  ! pso4_x          = mass of sulphate condensed on insoluble mode x [molec. cm-3]
  ! 
  !--- Local variables:
  !
  ! zansum          = aerosol number in the respective mode [cm-3]
  ! zxxsum          = aerosol mass for compound xx in the respective
  !                   mode, e.g. xx = bc, oc, a4 (sulfate) 
  !                   [g cm-3 for bc,oc and molec. cm-3 for sulfate]
  ! zxxav           = average mass of a sulfate particle in the respective
  !                   mode [molecules]
  ! zxxavy          = average mass of species xx in mode y []
  !                   where xx is ss, du, bc, oc, or a4 for sulfate
  !                   [molecules for sulfate and ug for others]
  ! zanli(:,:,:)    = Number of particles moved by the inter-modal 
  !                   coagulation []
  ! zansq(:,:,:)    = Number of particles moved by the intra-modal 
  !                   coagulation []
  ! zaernt(:,:,:)   = New particle number n(t+dt) after the integration
  !                   of the aerosol dynamics equation [cm-3]

  !--- Parameters:

  INTEGER :: kproma, kbdim, klev, krow
 
  REAL(dp):: paerml(kbdim,klev,naerocomp),   paernl(kbdim,klev,nclass),         &
             pa4delt(kbdim,klev,naerocomp),  panew(kbdim,klev),               &
             pm6rp(kbdim,klev,nclass)
  !>>dod soa
  REAL(dp) :: pttn(kbdim,klev,naerocomp)
  !<<dod

  REAL(dp):: pbfract1(kbdim,klev,nclass-1), pbfract2(kbdim,klev,nclass-1),     &
             pbfract5(kbdim,klev,3)
 
  REAL(dp):: pa(kbdim,klev,nclass), pb(kbdim,klev,nclass)
 
  REAL(dp):: pso4_5(kbdim,klev),           pso4_6(kbdim,klev),              & 
             pso4_7(kbdim,klev)
 
  ! Local variables: 
 
  INTEGER :: jclass, jl, jk, jn, jm, kmod

  !>>dod soa: removed arrays of average particle mass (replaced by pttn) 
  REAL(dp):: za4av1(kbdim,klev),           za4av2(kbdim,klev)
  !<<dod

  !>>dod replaced arrays zansq and zaernt with scalars
  REAL(dp):: zanli(kbdim,klev,nclass), zansq, zaernt 
  !<<dod

  REAL(dp):: zm6rp(nclass), zcrtcst(nclass)

  ! Auxiliary variables: 
 
  REAL(dp):: zansum, za4sum, za4av,                                         & 
             ztop,   zbot,   zatot,  zanse,  zanle,                         & 
             ze1,    zf1,    zf2,    zf3,    zf4,    zr1,                   &
             ztmst,  zc

  REAL(dp):: zamloss5,      zamloss6,      zamloss7,      zanloss5,         &
             zanloss6,      zanloss7,      ztloss,        zbfofac,   zms,   &
             zbfnfac,       zbftot,        zaerni,        zanytnt,          &
             zanytni,       zanytns,       zanytnm,       ztotn,            &
             zaerns,        zeps

  !>>dod soa
  REAL(dp) :: zbcmass, zocmass, zsoamass
  INTEGER :: jspec, isoans, isoaks, isoaas, isoacs, jn2
  !<<dod 

    REAL(dp), PARAMETER :: zmin_coag = 1.e-15 !gf #240 homogenized minimum value for coagulation coeffs

  !--- 0) Initialisations: ------------------------------------------------ 
 
  ztmst  = time_step_len 
  zeps = EPSILON(1._dp)
 
  za4av       = 0._dp 
  za4av1(:,:) = 0._dp
  za4av2(:,:) = 0._dp
  !>>dod soa deleted array initialisations but inserted zanli
  zanli(:,:,:) = 0._dp
  !<<dod
     
  !--- 1) Insoluble modes 
  ! 
  DO jk=1,klev 
     DO jl=1,kproma 

        !--- Dust modes insoluble: 
        !>>dod soa : do nothing.
        !<<dod

        !--- Aitken mode insoluble:
        !    Only considered process: 
        !    Coagulation and transfer from the insoluble aitken mode 
        !    to the soluble aitken and accumulation modes.

        zansum = paernl(jl,jk,iaiti)
        zansq = 0.0_dp 

        !--- Calculations only in presence of sufficient particles:

        IF (zansum > cmin_aernl) THEN 
 
           !--- 1.1) Case of no coagulation:
           !
           !         (pa(jl,jk,iaiti) < 1.e-15 .AND. pb(jl,jk,iaiti) < 1.e-15) 
           !
           !         => Nothing to be done
              
           !--- 1.2) Case with coagulation:

           IF (pa(jl,jk,iaiti) >= zmin_coag .OR. pb(jl,jk,iaiti) >= zmin_coag) THEN
                 
              !--- 1.2.1) Case of no inter-modal coagulation:
              !           dn/dt = -a*n**2 => 
              !           n(t)  = n0/(1 + n0*a*(t-t0))

              IF (pb(jl,jk,iaiti) < zmin_coag) THEN 
                 zaernt = zansum/(1.0_dp+zansum*pa(jl,jk,iaiti)*ztmst) 
                 !---zanli unchanged = 0
                 zansq = zansum-zaernt

              !--- 1.2.2) Case with inter- and intra-modal coagulation:
              !           dn/dt = -a*n**2 - b*n => 
              !           n(t)  = (b*n0*exp(-b(t-t0)))/((n0*a)(1-exp(-b(t-t0)))+b)
              ELSE   
                 !--- Calculate n(t+dt):
                 ze1=EXP(-pb(jl,jk,iaiti)*ztmst) 
                 ztop=pb(jl,jk,iaiti)*zansum*ze1 
                 zbot=zansum*pa(jl,jk,iaiti)*(1.0_dp-ze1)+pb(jl,jk,iaiti) 
                 zaernt=ztop/zbot 
                 !--- Limit n(t+dt) to available particle in the mode:
                 zaernt=MIN(zaernt, zansum) 
                 !--- Total change in particle numbers of the mode due to coagulation:
                 zatot=zansum-zaernt
                 !--- Contribution of the intra-modal coagulation:
                 zanse=zansum*zansum*pa(jl,jk,iaiti)
                 !--- Contribution of the inter-modal coagulation:
                 zanle=zansum*pb(jl,jk,iaiti) 
                 !--- Number of particles moved by the inter-modal coagulation:
                 zanli(jl,jk,iaiti)=zatot*zanle/(zanse+zanle) 
                 !--- Number of particles moved by the intra-modal coagulation:
                 zansq=zatot*zanse/(zanse+zanle) 

              END IF 

              !--- 1.2.3) Change masses of the insoluble aitken mode due to 
              !           intra-modal coagulation and the coagulation with the
              !           nucleation mode (transfers to the soluble modes
              !           of the particles coagulating with the nucleation mode 
              !           are done in m7_concoag):
 
              paerml(jl,jk,ibcki)=( zaernt + zansq + pbfract5(jl,jk,1)*zanli(jl,jk,iaiti) ) * & 
                                    pttn(jl,jk,ibcki)*1.e12_dp 

              paerml(jl,jk,iocki)=( zaernt + zansq + pbfract5(jl,jk,1)*zanli(jl,jk,iaiti) ) * & 
                                    pttn(jl,jk,iocki)*1.e12_dp 
             
              !--- 1.2.4) Change the numbers of the insoluble aitken mode due to 
              !           intra-modal coagulation:
 
              paernl(jl,jk,iaiti)=zaernt + pbfract5(jl,jk,1)*zanli(jl,jk,iaiti)
 
              !--- 1.2.5) Store changes in masses of compounds in the insoluble 
              !           aitken mode due to inter-modal coagulation:
              !           (zanli(:,:,x)   = total number of particles moved from mode x
              !            pbfract5(:,:,x)= fraction of the total number of particles
              !                             moved from mode 5 that is moved to mode x   )
 
              pa4delt(jl,jk,ibcks)=pbfract5(jl,jk,2)*zanli(jl,jk,iaiti)*pttn(jl,jk,ibcki)*1.e12_dp
              pa4delt(jl,jk,iocks)=pbfract5(jl,jk,2)*zanli(jl,jk,iaiti)*pttn(jl,jk,iocki)*1.e12_dp
              pa4delt(jl,jk,ibcas)=pbfract5(jl,jk,3)*zanli(jl,jk,iaiti)*pttn(jl,jk,ibcki)*1.e12_dp
              pa4delt(jl,jk,iocas)=pbfract5(jl,jk,3)*zanli(jl,jk,iaiti)*pttn(jl,jk,iocki)*1.e12_dp
 
              IF (nsoa == 1) THEN
                 DO jm=1,nsoaspec
                    IF (.NOT. soaprop(jm)%lvolatile) THEN
                       !>>dod deleted isoa_ix
                       jspec = soaprop(jm)%spid_soa
                       jn = speclist(jspec)%iaerocomp(iaiti)
                       !<<dod
                       zsoamass = 1.E12_dp*pttn(jl,jk,jn)
                       paerml(jl,jk,jn) = ( zaernt + zansq +       &
                                             pbfract5(jl,jk,1)*zanli(jl,jk,iaiti) ) * &
                                             zsoamass

                       !>>dod deleted isoa_ix
                       jn = speclist(jspec)%iaerocomp(iaits)
                       jn2 = speclist(jspec)%iaerocomp(iaccs)
                       
                       pa4delt(jl,jk,jn) = pbfract5(jl,jk,2)*zanli(jl,jk,iaiti) * zsoamass
                       pa4delt(jl,jk,jn2) = pbfract5(jl,jk,3)*zanli(jl,jk,iaiti) * zsoamass
                       !<<dod
                    END IF
                 END DO
              END IF

           END IF
        END IF
     END DO
  END DO


  ! 
  !--- 2) Soluble modes: --------------------------------------------------
  ! 
!CDIR unroll=5
  mode : DO jclass=1,nsol 
     level : DO jk=1,klev 
        longitude : DO jl=1,kproma 


           !--- Nucleation mode:
           IF (jclass .EQ. 1) THEN 
              zansum=paernl(jl,jk,jclass)+panew(jl,jk) 
              za4sum=paerml(jl,jk,jclass)+pa4delt(jl,jk,1)

           !--- Others:
           ELSE 
              zansum=paernl(jl,jk,jclass) 
              za4sum=paerml(jl,jk,jclass)
           END IF 

           zaernt=zansum 
           zansq=0.0_dp 

           !--- Calculations only in presence of sufficient particles:

           IF (zansum > cmin_aernl) THEN 

              za4av=za4sum/zansum 
              
              IF (jclass.EQ.1) THEN 
                 za4av1(jl,jk)=za4av 
              ELSE IF (jclass.EQ.2) THEN 
                 za4av2(jl,jk)=za4av 
              END IF

              !--- 2.1) Case of no coagulation:
              !         
              IF (pa(jl,jk,jclass) < zmin_coag .AND. pb(jl,jk,jclass) < zmin_coag) THEN  
                 !--- Nucleation in mode 1 only. 
                 !    Nothing to be done for other modes.

                 IF(jclass.EQ.1) THEN 
                    paerml(jl,jk,iso4ns)=za4sum
                    paernl(jl,jk,inucs)=zansum

                 END IF 
 
              !--- 2.2) Case with coagulation:

              ELSE 

                 !--- 2.2.1) Case of no nucleation:

                 !--- Not Mode 1 or Nucleation rate below 1/s:

                 IF ( (jclass .NE. 1) .OR. (panew(jl,jk)/ztmst < 1.0_dp) ) THEN

                    paernl(jl,jk,jclass)=zansum 

                    !--- 2.2.1a) Case of no inter-modal coagulation:
                    !            dn/dt = -a*n**2 => 
                    !            n(t)  = n0/(1 + n0*a*(t-t0))

                    IF (pb(jl,jk,jclass) < zmin_coag) THEN 
                       zaernt=zansum/(1.0_dp+zansum*pa(jl,jk,jclass)*ztmst) 
                       !---zanli unchanged (=0)
                       zansq=zansum-zaernt

                    !--- 2.2.1b) Case with inter- and intra-modal coagulation:
                    !            dn/dt = -a*n**2 - b*n => 
                    !            n(t)  = (b*n0*exp(-b(t-t0)))/((n0*a)(1-exp(-b(t-t0)))+b)

                    ELSE            
                       !--- Calculate n(t+dt):
                       ze1=EXP(-pb(jl,jk,jclass)*ztmst) 
                       ztop=pb(jl,jk,jclass)*zansum*ze1 
                       zbot=zansum*pa(jl,jk,jclass)*(1.0_dp-ze1)+pb(jl,jk,jclass) 
                       zaernt=ztop/zbot 
                       !--- Limit n(t+dt) to available particle in the mode:
                       zaernt=MIN(zaernt, zansum) 
                       !--- Total change in particle numbers of the mode due to coagulation:
                       zatot=zansum-zaernt
                       !--- Contribution of the intra-modal coagulation:
                       zanse=zansum*zansum*pa(jl,jk,jclass)
                       !--- Contribution of the inter-modal coagulation:
                       zanle=zansum*pb(jl,jk,jclass) 
                       !--- Number of particles moved by the inter-modal coagulation:
                       zanli(jl,jk,jclass)=zatot*zanle/(zanse+zanle) 
                       !--- Number of particles moved by the intra-modal coagulation:
                       zansq=zatot*zanse/(zanse+zanle) 
                    END IF 
                 
                 !--- 2.2.2) Case with nucleation:

                 ELSE IF ( (jclass .EQ. 1) .AND. (panew(jl,jk)/ztmst >= 1.0_dp) ) THEN

                    !--- 2.2.2a) Nucleation, inter- and intra-modal coagulation:
                    !            dn/dt = -a*n**2 - b*n + c => 
                    !            n(t)  = -(b/(2a)) + 
                    !                    R/2a * [ ((1 - (-2ax0-b+R)/(+2ax0+b+R))exp(-Rt)) /
                    !                             ((1 + (-2ax0-b+R)/(+2ax0+b+R))exp(-Rt))  ]
                    !            where:  R=SQRT(b**2+4ac)
                    !
                    !            If b/=0 then always a/=0. The only case where a would be 0
                    !            and b unequal zero is the case of no pre-existing particles 
                    !            in the nucleation mode but pre-existing particles in other
                    !            modes. For this case a is calculated for an assumed radius
                    !            of a critical cluster in m7_coaset. 

                    IF (pb(jl,jk,jclass) >= zmin_coag) THEN
                       !--- Calculate n(t):
                       !--- c:
                       zc=panew(jl,jk)/ztmst
                       !--- R:
                       zf1=pb(jl,jk,jclass)*pb(jl,jk,jclass)+4.0_dp*pa(jl,jk,jclass)*zc 
                       zr1=SQRT(zf1) 
                       !--- exp(-Rt):
                       ze1=EXP(-zr1*ztmst) 
                       !--- 2ax0+b:
                       zf2=2.0_dp*pa(jl,jk,jclass)*paernl(jl,jk,jclass)+pb(jl,jk,jclass) 
                       !--- Term in squared bracket:
                       zf3=ze1*(zr1-zf2)/(zr1+zf2) 
                       zf4=(1.0_dp-zf3)/(1.0_dp+zf3) 
                       !--- n(t):
                       zaernt=(zr1*zf4-pb(jl,jk,jclass))/2.0_dp/pa(jl,jk,jclass) 
                       !--- Limit n(t+dt) to available particle in the mode:
                       zaernt=MIN(zaernt, zansum) 
                       !--- Total change in particle numbers of the mode due to coagulation:
                       zatot=zansum-zaernt
                       !--- Contribution of the intra-modal coagulation:
                       zanse=zansum*zansum*pa(jl,jk,jclass)
                       !--- Contribution of the inter-modal coagulation:
                       zanle=zansum*pb(jl,jk,jclass) 
                       !--- Number of particles moved by the inter-modal coagulation:
                       zanli(jl,jk,jclass)=zatot*zanle/(zanse+zanle) 
                       !--- Number of particles moved by the intra-modal coagulation:
                       zansq=zatot*zanse/(zanse+zanle) 

                       !--- 2.2.2b) Nucleation and intra-modal coagulation:
                       !            dn/dt = -a*n**2 - b*n + c with b=0 =>
                       !            dn/dt = -a*n**2 + c => 
                       !            n(t)  = R/2a * [ ((1 - (-2ax0+R)/(+2ax0+R))exp(-Rt)) /
                       !                             ((1 + (-2ax0+R)/(+2ax0+R))exp(-Rt))  ]
                       !            where:  R=SQRT(4ac)
                       !
                       !            Can be shown to be equivalent to:
                       !
                       !            n(t)  = R1*((x0+R1)/(x0-R1)+exp(-SQRT(-4ac)t)) / 
                       !                       ((x0+R1)/(x0-R1)-exp(-SQRT(-4ac)t))
                       !            where R1=SQRT(c/a)
                       
                    ELSE IF (pb(jl,jk,jclass) < zmin_coag) THEN 
                       !--- c:
                       zc=panew(jl,jk)/ztmst
                       !--- R1:
                       zr1=SQRT(zc/pa(jl,jk,jclass)) 
                       !--- exp(-Rt):
                       ze1=EXP(-zr1*2.0_dp*pa(jl,jk,jclass)*ztmst)
                       !--- n(t):
                       zf1=(paernl(jl,jk,jclass)+zr1)/(paernl(jl,jk,jclass)-zr1) 
                       ztop=zr1*(zf1+ze1) 
                       zbot=zf1-ze1 
                       IF (zbot < zmin_coag) THEN 
                          zaernt=zansum 
                       ELSE 
                          zaernt=ztop/zbot 
                       END IF 
                       !--- Limit n(t+dt) to available particle in the mode:
                       zaernt=MIN(zaernt, zansum) 
                       !--- Number of particles moved by the inter-modal coagulation:
                       zanli(jl,jk,jclass)=0.0_dp 
                       !--- Number of particles moved by the intra-modal coagulation:
                       zansq=zansum-zaernt
                    END IF 
                 END IF 
                 !---2.2.3 New bit for insoluble/souble coagulation
                 !--- sum total insoluble+soluble paticles in mode jclass JJNW
                 IF (jclass .EQ. 1 .AND. zanli(jl,jk,jclass)>0.0_dp) THEN
                    zaerni=paernl(jl,jk,iaiti)+paernl(jl,jk,iacci)+paernl(jl,jk,icoai)
                    zaerns=zansum+paernl(jl,jk,iaits)
                    !                 zaerns=zansum
                    ztotn=zaerns+zaerni
                    IF (zaerns .gt. zaerni .and. zaerni .gt. zeps) THEN !SF (#136): use zeps instead of 0._dp
                       ! calculate analytical solution no of mixed particles for coagulation
                       ! between paernl(jl,jk,jclass) soluble particles and zaerni insouble of
                       ! the same dimensions
                       IF (zaerni .gt. 1.0_dp) then
                          zanytni=4.0_dp*zaerni/((2.0_dp+pa(jl,jk,jclass)*ztmst*ztotn)*                     &
                               (2.0_dp+pa(jl,jk,jclass)*ztmst*(ztotn-zaerni)))
                       ELSE
                          zanytni = 0.0_dp
                       ENDIF
                       zanytnt=2.0_dp*ztotn/(2.0_dp+pa(jl,jk,jclass)*ztmst*ztotn)
                       zanytns=4.0_dp*zaerns/((2.0_dp+pa(jl,jk,jclass)*ztmst*ztotn)*            &
                               (2.0_dp+pa(jl,jk,jclass)*ztmst*(ztotn-zaerns)))
                       zanytnm=zanytnt-(zanytni+zanytns)
                       zanytnm=min(zanytnm,zaerni)
                       zanytni=zaerni-zanytnm
                       zanytns=zaernt
 
                       IF (za4av > 0._dp) THEN
!DIR$ UNROLL
!CDIR UNROLL=7
                          DO kmod=1,nclass
                             zm6rp(kmod)=pm6rp(jl,jk,kmod)
                          END DO
                          CALL m7_coat(zm6rp,zcrtcst)
                          zamloss5=paernl(jl,jk,5)/zaerni*zanytnm*zcrtcst(5)
                          zanloss5=zamloss5/za4av
                          zamloss6=paernl(jl,jk,6)/zaerni*zanytnm*zcrtcst(6)
                          zanloss6=zamloss6/za4av
                          zamloss7=paernl(jl,jk,7)/zaerni*zanytnm*zcrtcst(7)
                          zanloss7=zamloss7/za4av
                       ELSE
                          zamloss5 = 0._dp   
                          zamloss6 = 0._dp   
                          zamloss7 = 0._dp   
                          zanloss5 = 0._dp   
                          zanloss6 = 0._dp   
                          zanloss7 = 0._dp   
                       END IF
                       ztloss=zanloss5+zanloss6+zanloss7
                       zms=zansq*0.95_dp

                       ztloss=min(ztloss,zansq*0.95_dp)                               
                       zbfofac=zanli(jl,jk,jclass)/(zanli(jl,jk,jclass)+ztloss)
                       zbfnfac=ztloss/(zanli(jl,jk,jclass)+ztloss)
                       zanli(jl,jk,jclass)=zanli(jl,jk,jclass)+ztloss
                       zansq=zansq-ztloss
                       zbftot=0.0_dp
!DIR$ UNROLL
!CDIR UNROLL=7
                       DO kmod=1,nclass
                          IF(kmod>jclass) THEN
                             pbfract1(jl,jk,kmod-jclass)=pbfract1(jl,jk,kmod-jclass)*zbfofac


                             IF (kmod.GE.5) THEN
                                pbfract1(jl,jk,kmod-jclass)=pbfract1(jl,jk,kmod-jclass)+                  &
                                                          zbfnfac*paernl(jl,jk,kmod)/zaerni

                             END IF
                             zbftot=zbftot+pbfract1(jl,jk,kmod-jclass)
                          END IF
                       END DO
!DIR$ UNROLL
!CDIR UNROLL=7
                       DO kmod=1,nclass
                          IF (kmod>jclass) THEN
!kai 
                             IF (abs(zbftot).gt.zeps) then 
                                pbfract1(jl,jk,kmod-jclass)=pbfract1(jl,jk,kmod-jclass)/zbftot
                             ELSE
                                pbfract1(jl,jk,kmod-jclass)=0._dp 
                             END IF 

                          END IF
                       END DO
                    END IF
                 END IF
                 !---- End of new inslouble/soluble caogulation routine JJNW
                 !--- 2.3) Change masses and numbers of the respective modes to account-----------
                 !         for intra-modal coagulation (zansq) and coagulation with
                 !         higher modes (zaernt):
                 !
                 !--- 2.3.1) Change mass of the sulfur compounds:

                 paerml(jl,jk,jclass)=(zaernt+zansq)*za4av 
 
                 !--- 2.3.2) Change mass of the carbon compounds:

                 IF (jclass.EQ.2) THEN 
                    paerml(jl,jk,ibcks)=(zaernt+zansq)*pttn(jl,jk,ibcks)*1.E12_dp
                    paerml(jl,jk,iocks)=(zaernt+zansq)*pttn(jl,jk,iocks)*1.E12_dp
                 ELSE IF (jclass.EQ.3) THEN 
                    paerml(jl,jk,ibcas)=(zaernt+zansq)*pttn(jl,jk,ibcas)*1.E12_dp
                    paerml(jl,jk,iocas)=(zaernt+zansq)*pttn(jl,jk,iocas)*1.E12_dp
                 ELSE IF (jclass.EQ.4) THEN 
                    paerml(jl,jk,ibccs)=(zaernt+zansq)*pttn(jl,jk,ibccs)*1.E12_dp
                    paerml(jl,jk,ioccs)=(zaernt+zansq)*pttn(jl,jk,ioccs)*1.E12_dp
                 END IF 

                 IF (nsoa == 1) THEN
                    DO jm=1,nsoaspec
                       IF (sizeclass(jclass)%lsoainclass .AND. (.NOT. soaprop(jm)%lvolatile)) THEN
                          !>>dod deleted isoa_ix
                          jspec = soaprop(jm)%spid_soa
                          jn = speclist(jspec)%iaerocomp(jclass)
                          !<<dod
                          paerml(jl,jk,jn) = (zaernt+zansq) * pttn(jl,jk,jn)*1.E12_dp
                       END IF
                    END DO
                 END IF

                 !--- 2.3.3) Particle numbers:
 
                 paernl(jl,jk,jclass)=zaernt             !@@@+zansq(jl,jk,jclass)/2.0 
 
                 !--- 2.4) Calculate changes in particle masses due to inter-modal --------------
                 !         coagulation:

                 !--- 2.4.1) Transfer of mass from mode 1 to other modes:

                 IF (jclass .EQ. 1) THEN 
                     
                    ! Mass from 1 to 2: 
 
                    pa4delt(jl,jk,iso4ks)=pbfract1(jl,jk,1)*zanli(jl,jk,1)*za4av 

                    ! Mass from 1 to 2 due to coag. with 5:
 
                    pa4delt(jl,jk,iso4ks)=pa4delt(jl,jk,iso4ks)+pbfract1(jl,jk,4)*zanli(jl,jk,1)*za4av 
  
                    ! Mass from 1 to 3:
 
                    pa4delt(jl,jk,iso4as)=pbfract1(jl,jk,2)*zanli(jl,jk,1)*za4av 

                    ! Mass from 1 to 3 due to coag. with 6:
 
                    pa4delt(jl,jk,iso4as)=pa4delt(jl,jk,iso4as)+pbfract1(jl,jk,5)*zanli(jl,jk,1)*za4av 

                    ! Mass from 1 to 4: 
 
                    pa4delt(jl,jk,iso4cs)=pbfract1(jl,jk,3)*zanli(jl,jk,1)*za4av 

                    ! Mass from 1 to 4 due to coag. with 7:
 
                    pa4delt(jl,jk,iso4cs)=pa4delt(jl,jk,iso4cs)+pbfract1(jl,jk,6)*zanli(jl,jk,1)*za4av 

                    IF (nsoa == 1 .AND. sizeclass(inucs)%lsoainclass) THEN
                       DO jm=1,nsoaspec
                          IF (.NOT. soaprop(jm)%lvolatile) THEN
                             !>>dod deleted isoa_ix
                             jspec = soaprop(jm)%spid_soa
                             isoans = speclist(jspec)%iaerocomp(inucs)
                             isoaks = speclist(jspec)%iaerocomp(iaits)
                             isoaas = speclist(jspec)%iaerocomp(iaccs)
                             isoacs = speclist(jspec)%iaerocomp(icoas)

                             zsoamass = 1.E12_dp*pttn(jl,jk,isoans)*zanli(jl,jk,1)

                             !---mass from nucs to aits
                             pa4delt(jl,jk,isoaks) = zsoamass *  ( pbfract1(jl,jk,1) +  &
                             !---add mass from nucs to aits due to coag with aiti                  
                                                                  pbfract1(jl,jk,4) )

                             !---mass from nucs to accs
                             pa4delt(jl,jk,isoaas) = zsoamass * ( pbfract1(jl,jk,2) + &
                             !---add mass from nucs to accs due to coag with acci
                                                                  pbfract1(jl,jk,5) )

                             !---mass from nucs to coas
                             pa4delt(jl,jk,isoacs) = zsoamass * ( pbfract1(jl,jk,3) + &
                             !---add mass from nucs to coas due to coag with coai
                                                                  pbfract1(jl,jk,6) )
                             !<<dod
                          END IF
                       END DO
                    END IF
 
                 !---  2.4.2) Transfer of mass from mode 2 to other modes:
 
                 ELSE IF (jclass .EQ. 2) THEN 
 
                    zbcmass = 1.E12_dp * pttn(jl,jk,ibcks) * zanli(jl,jk,2)
                    zocmass = 1.E12_dp * pttn(jl,jk,iocks) * zanli(jl,jk,2)

                    ! Mass from 2 to 3: 
                     
                    pa4delt(jl,jk,iso4as)=pa4delt(jl,jk,iso4as)+pbfract2(jl,jk,2)*zanli(jl,jk,2)*za4av 
                    pa4delt(jl,jk,ibcas)=pa4delt(jl,jk,ibcas)+                                & 
                                         pbfract2(jl,jk,2)*zbcmass
                    pa4delt(jl,jk,iocas)=pa4delt(jl,jk,iocas)+                                & 
                                         pbfract2(jl,jk,2)*zocmass
  
                    ! Mass from 2 to 3 due to coag. with 6: 
 
                    pa4delt(jl,jk,iso4as)=pa4delt(jl,jk,iso4as)+pbfract2(jl,jk,5)*zanli(jl,jk,2)*za4av 
                    pa4delt(jl,jk,ibcas)=pa4delt(jl,jk,ibcas)+                                & 
                                         pbfract2(jl,jk,5)*zbcmass
                    pa4delt(jl,jk,iocas)=pa4delt(jl,jk,iocas)+                                & 
                                         pbfract2(jl,jk,5)*zocmass

                    ! Mass from 2 to 4: 
 
                    pa4delt(jl,jk,iso4cs)=pa4delt(jl,jk,iso4cs)+pbfract2(jl,jk,3)*zanli(jl,jk,2)*za4av 
                    pa4delt(jl,jk,ibccs)=pa4delt(jl,jk,ibccs)+                                & 
                                         pbfract2(jl,jk,3)*zbcmass
                    pa4delt(jl,jk,ioccs)=pa4delt(jl,jk,ioccs)+                                & 
                                         pbfract2(jl,jk,3)*zocmass

                    ! Mass from 2 to 4 due to coag. with 7: 
 
                    pa4delt(jl,jk,iso4cs)=pa4delt(jl,jk,iso4cs)+pbfract2(jl,jk,6)*zanli(jl,jk,2)*za4av 
                    pa4delt(jl,jk,ibccs)=pa4delt(jl,jk,ibccs)+                                & 
                                         pbfract2(jl,jk,6)*zbcmass
                    pa4delt(jl,jk,ioccs)=pa4delt(jl,jk,ioccs)+                                & 
                                         pbfract2(jl,jk,6)*zocmass

 
                    IF (nsoa == 1) THEN
                       DO jm=1,nsoaspec
                          IF (.NOT. soaprop(jm)%lvolatile) THEN
                             !>>dod deleted isoa_ix
                             jspec = soaprop(jm)%spid_soa
                             isoaks = speclist(jspec)%iaerocomp(iaits)
                             isoaas = speclist(jspec)%iaerocomp(iaccs)
                             isoacs = speclist(jspec)%iaerocomp(icoas)

                             zsoamass = 1.E12_dp*pttn(jl,jk,isoaks)*zanli(jl,jk,2)

                             !---mass from aits to accs
                             pa4delt(jl,jk,isoaas) = pa4delt(jl,jk,isoaas) + &
                                                     zsoamass * ( pbfract2(jl,jk,2) +  &
                             !---add mass from aits to accs due to coag with acci                  
                                                                pbfract2(jl,jk,5) )

                             !---mass from aits to coas
                             pa4delt(jl,jk,isoacs) = pa4delt(jl,jk,isoacs) + &
                                                     zsoamass *  ( pbfract2(jl,jk,3) + &
                             !---add mass from aits to coas due to coag with coai
                                                                  pbfract2(jl,jk,6) )
                             !<<dod
                          END IF
                       END DO
                    END IF
 
                    ! Mass from 2 due to coagulation of 2 with 5 remains in 2:
                    !
                    !@@@ No effect as pbfract2(:,:,4)=0.!
                    !@@@ (Not needed as it is assumed that 5 coagulates with 2 
                    !@@@ and therefor the masses in 2 remain unchanged!)

                    pa4delt(jl,jk,iso4ks)=pa4delt(jl,jk,iso4ks)+pbfract2(jl,jk,4)*zanli(jl,jk,2)*za4av
                    pa4delt(jl,jk,ibcks)=pa4delt(jl,jk,ibcks)+                                &
                                         pbfract2(jl,jk,4)*zbcmass
                    pa4delt(jl,jk,iocks)=pa4delt(jl,jk,iocks)+                                &
                                         pbfract2(jl,jk,4)*zocmass
 
                    !>>dod soa
                    IF (nsoa == 1) THEN
                       DO jm=1,nsoaspec
                          IF (.NOT. soaprop(jm)%lvolatile) THEN
                             !>>dod deleted isoa_ix
                             jspec = soaprop(jm)%spid_soa
                             isoaks = speclist(jspec)%iaerocomp(iaits)

                             zsoamass = 1.E12_dp*pttn(jl,jk,isoaks)*zanli(jl,jk,2)

                             pa4delt(jl,jk,isoaks) = pa4delt(jl,jk,isoaks) + zsoamass * pbfract2(jl,jk,4)
                             !<<dod
                          END IF
                       END DO
                    END IF
                    !<<dod
                 END IF             !---end IF jclass==1 / ELSIF jclass==2
              END IF             !---end IF pa < zmin_coag .AND. pb < zmin_coag / ELSE
           END IF            !---end IF sufficient particles
        END DO longitude 
     END DO level 
  END DO mode 


  !--- 3) Calculate transfer from the insoluble to the soluble modes: -------------------------
  !>>dod timers
  IF (ltimer) CALL timer_start(timer_ham_m7_concoag)
  !<<dod
  CALL m7_concoag (kproma,   kbdim,   klev,   krow,              &
                   paerml,   paernl,  pm6rp,  pa4delt, zanli,    & 
                   !>>dod soa: added pttn, removed bc,oc and du averages
                   za4av1,   za4av2, pttn,                       & 
                   pso4_5,  pso4_6, pso4_7,                      &
                   !<<dod
                   pbfract1, pbfract2                            )
  !>>dod timers
  IF (ltimer) CALL timer_stop(timer_ham_m7_concoag)
  !<<dod

  !--- 4) Final change of the aerosol masses due to nucleation, -------------------------------
  !       inter-modal coagulation and condensation on the insoluble modes:
  !       (Nucleation mode already done above.)


  !>>dod soa 
  DO jclass=1,nclass
     jm = speclist(id_so4)%iaerocomp(jclass)    !!mgs!!   im7table(jclass,id_so4)
     IF (jm > 0 .AND. jclass /= inucs) paerml(1:kproma,:,jm) = paerml(1:kproma,:,jm) + &
                                                            pa4delt(1:kproma,:,jm)
     jm = speclist(id_bc)%iaerocomp(jclass)    !!mgs!!   im7table(jclass,id_bc)
     IF (jm > 0) paerml(1:kproma,:,jm)=paerml(1:kproma,:,jm)+pa4delt(1:kproma,:,jm) 
     jm = speclist(id_oc)%iaerocomp(jclass)    !!mgs!!   im7table(jclass,id_oc)
     IF (jm > 0) paerml(1:kproma,:,jm)=paerml(1:kproma,:,jm)+pa4delt(1:kproma,:,jm) 
     jm = speclist(id_du)%iaerocomp(jclass)    !!mgs!!   im7table(jclass,id_du)
     IF (jm > 0) paerml(1:kproma,:,jm)=paerml(1:kproma,:,jm)+pa4delt(1:kproma,:,jm) 
     
     IF (nsoa == 1) THEN
        DO jm=1,nsoaspec
           IF (sizeclass(jclass)%lsoainclass .AND. (.NOT. soaprop(jm)%lvolatile)) THEN
              !>>dod deleted isoa_ix
              jspec = soaprop(jm)%spid_soa
              jn = speclist(jspec)%iaerocomp(jclass)
              paerml(1:kproma,:,jn) = paerml(1:kproma,:,jn) + pa4delt(1:kproma,:,jn)
              !<<dod
           END IF
        END DO
     END IF
  END DO 
  !<<dod   

END SUBROUTINE m7_delcoa 



END MODULE mo_ham_m7
