!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_tools.f90
!!
!! \brief
!! mo_ham_tools hold auxiliary routines for the 
!! HAM aerosol model
!!
!! \author Philip Stier (MPI-Met)
!!
!! \responsible_coder
!! Philip Stier, philip.stier@physics.ox.ac.uk
!!
!! \revision_history
!!   -# Philip Stier (MPI-Met) - original code (2002)
!!   -# Philip Stier (MPI-Met) - added ham_logtail (2004)
!!   -# Betty Croft (Dalhousie University) - added scavenging coefficient bilinear interpolation (2008)
!!   -# Sylvaine Ferrachat (ETH Zurich) - cleanup and security (2011)
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
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

MODULE mo_ham_tools

  ! *mo_ham_tools* hold auxiliary routines for the 
  !                 HAM aerosol model

  USE mo_kind,               ONLY: dp
  USE mo_exception,          ONLY: finish

  IMPLICIT NONE

CONTAINS

  SUBROUTINE geoindex (plon, plat, klon, klat)

    !   *geoindex* calculates the corresponding lat-lon (klon,klat) index
    !              for given lat-lon [-90,90][0,360] coordinates (plon, plat)
    !
    !   Authors: 
    !   --------- 
    !   Philip Stier, MPI-MET            2002 
    !
    !   Externals 
    !   ----------- 
    !   none 
    

    USE mo_exception, ONLY: finish
    USE mo_control,   ONLY: ngl, nlon
    USE mo_gaussgrid, ONLY: philat, philon
    IMPLICIT NONE

    REAL(dp),INTENT(IN)  :: plon, plat   ! Coordinates in degrees
    INTEGER, INTENT(OUT) :: klon, klat   ! Corresponding indices

    REAL(dp)            :: zxdif, zydif
    INTEGER             :: jlon, jlat
      

    IF((NINT(plon)<0).OR.(NINT(plon)>360).OR.(NINT(plat)<-90).OR.(NINT(plat)>90)) THEN
       CALL finish('geoindex:', 'Coordinates out of range')
    END IF

    klon=-999
    klat=-999

    zxdif=360._dp
    zydif=180._dp

    !--- 1) Search for closest longitude in [0,360}:

    DO jlon = 1, nlon
       IF(ABS(plon-philon(jlon)) < zxdif) THEN
          klon=jlon
          zxdif=ABS(plon-philon(jlon))
       END IF
    END DO

    !--- 2) Search for closest latitude [90,-90]:

    DO jlat = 1, ngl
       IF( ABS(plat-philat(jlat)) < zydif ) THEN
          klat=jlat
          zydif=ABS(plat-philat(jlat)) 
       END IF
    END DO


    IF(klon==-999 .OR. klat==-999) CALL finish('geoindex:', 'No index found')

  END SUBROUTINE geoindex

!---------------------------------------------------------------------------------------------

  SUBROUTINE calc_daylength

    ! *dayfac* calculates relative daylength
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-Met          14/11/2002
    !
    ! Method:
    ! -------
    ! Adapted from the routines "solang" and "prerad"

    USE mo_kind,           ONLY: dp
    USE mo_control,        ONLY: nlon, ngl
    USE mo_gaussgrid,      ONLY: coslon, sinlon, gl_twomu, gl_sqcst
    USE mo_time_control,   ONLY: get_clock, current_date
    USE mo_radiation_parameters, ONLY: decl_sun_cur
    USE mo_decomposition,  ONLY: gdc => global_decomposition
    USE mo_transpose,      ONLY: scatter_gp
    USE mo_test_trans,     ONLY: test_gridpoint
    USE mo_ham_streams,    ONLY: daylength 

    IMPLICIT NONE

    INTEGER :: jl, jlat

    LOGICAL :: lo

    REAL(dp):: zclock
    REAL(dp):: czen1,  czen2,  czen3, ztim1, ztim2, ztim3, zsum

    REAL(dp):: zmu0(nlon), zrdayl(nlon)

    REAL(dp), POINTER :: zdaylength(:,:)


    !--- 0) Initialization:

    ALLOCATE (zdaylength(nlon,ngl))


    !--- 1) Compute orbital parameters for present time step:

    zclock = get_clock(current_date)

!!mgs: removed orbital calculation - can now make use of variables from mo_radiation_parameters

    czen1 = SIN(decl_sun_cur)
    czen2 = COS(decl_sun_cur)*COS(zclock)
    czen3 = COS(decl_sun_cur)*SIN(zclock)


    !--- 2) Calculate relative daylength (from routine solang):

    DO jlat = 1, ngl

       ztim1 =  czen1*0.5_dp*gl_twomu(jlat)
       ztim2 = -czen2*gl_sqcst(jlat)
       ztim3 =  czen3*gl_sqcst(jlat)

       DO jl = 1, nlon
          zmu0(jl) = ztim1 + ztim2*coslon(jl) + ztim3*sinlon(jl)
          lo = zmu0(jl) >= 0._dp
          zrdayl(jl) = MERGE(1._dp,0._dp,lo)
       END DO

       zsum = SUM(zrdayl(1:nlon))

       IF (ABS(zsum) > 0._dp) THEN
          zdaylength(:,jlat) = zsum/REAL(nlon,dp)
       ELSE
          zdaylength(:,jlat) = 0._dp
       END IF

    END DO


    !--- 3) Scatter in stream element:

    CALL scatter_gp(zdaylength, daylength, gdc)

    CALL test_gridpoint(daylength, 'daylength')

    DEALLOCATE (zdaylength)

  END SUBROUTINE calc_daylength

!---------------------------------------------------------------------------------------------

  SUBROUTINE ham_m7_logtail(kproma, kbdim,  klev,   krow,  kmod, &
                            ld_numb, pcmr, pr, pfrac)

    ! *ham_m7_logtail* calculates mass- or number-fraction larger than
    !                  the radius pr for one given mode of a superposition 
    !                  of nclass log-normal aerosol distributions
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-MET                       2004
    !
    ! Revision:
    ! ---------
    ! Sylvaine Ferrachat, ETH Zurich,             2013
    !    --> this routine computes now one mode at a time, which allows to reduce the computational load
    !        for all irrelevant modes
    !
    ! Method:
    ! -------

    !
    ! The calculation of the activated number fraction and mass fraction
    ! from the radius of activation is done by a transformation of the 
    ! log-normal distribution to the error function which is then computed
    ! using the routine m7_cumulative_normal:
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


    USE mo_ham_m7ctl,   ONLY: sigmaln, cmedr2mmedr
    USE mo_ham_m7,      ONLY: m7_cumulative_normal

    IMPLICIT NONE

    !--- Subroutine parameters:

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow, kmod

    LOGICAL, INTENT(in) :: ld_numb   !number vs mass switch

    REAL(dp), INTENT(in) :: pcmr(kbdim,klev) !count mean radius
    REAL(dp), INTENT(in) :: pr(kbdim,klev) !lower bound radius

    REAL(dp), INTENT(out) :: pfrac(kbdim,klev)

    !--- Local variables:

    INTEGER :: jl, jk

    REAL(dp) :: zt, zdummy, zeps, zfact


    !--- 0) 

    zeps=EPSILON(1.0_dp)

    !--- 1) 

    !>>SF
    IF (ld_numb) THEN !number calculation
       zfact = 1._dp
    ELSE !mass calculation
       zfact = cmedr2mmedr(kmod)
    ENDIF
    !<<SF

    DO jk=1, klev
       DO jl=1, kproma

          IF (pr(jl,jk)>zeps .AND. pcmr(jl,jk)>zeps) THEN

             !--- Transform number distribution to error function:

             zt=(LOG(pr(jl,jk))-LOG(pcmr(jl,jk)*zfact))/sigmaln(kmod)

             !--- Calculate the cumulative of the log-normal number distribution:

             CALL m7_cumulative_normal(zt,zdummy,pfrac(jl,jk))

             !--- Calculate the cumulative of the log-normal mass distribution:

          ELSE IF (pr(jl,jk)<zeps .AND. pcmr(jl,jk)>zeps) THEN

             pfrac(jl,jk)=1.0_dp

          ELSE

             pfrac(jl,jk)=0.0_dp

          END IF

       END DO !kproma
    END DO !klev

  END SUBROUTINE ham_m7_logtail

!>>gf
! ------------------------------------------------------------------------------------------------------

  SUBROUTINE ham_m7_invertlogtail(kproma, kbdim,  klev,  krow, kmod, &
                               pcmr, pxie, pcritrad )

    ! *ham_m7_invertlogtail* calculates the critical                                    
                                                                                      
    ! radius for one mode with a log-normal aerosol distribution       
        
    ! that contains a given number in the logtail       
   
    ! The input to the subroutine is pxie, the x for erf^-1(x)     
                                                                       
    ! Author:        

    ! -------          
                       
    ! Betty Croft, Dalhousie University                       2007    

    ! Updates:
    ! Grazia Frontoso, C2SM-ETHZ, adjusted for e6-h2          2013    
    ! Sylvaine Ferrachat, ETHZ                                2013
    !     ---> cleanup, optimization, and change so that it computes one mode only (allows to reduce the
    !          computational load for all irrelevant modes)
    ! 
    ! Method:     
                  
    !   The cumulative number in a logtail is:       
    !                                                      
    !                                              ln (R/Rg)                                  
    !   Tail = N - Fn(r) = N/2 - N/2 * erf  [ --------------------]  
    !                                        sqrt (2) * ln (sigma)                                
    ! 
    !   This can be solved for R by inverting the error function (inverf)  
    !                                                                       
    !   R = Rg *exp( sqrt (2) * ln (sigma) * inverf(1 - (Tail*2/N)))                              
    !
    !  The approximation for the inverse error function used is (valid for x in (0,1)):
    !      
    !                     -2      ln(1-x^2)             2     ln(1-x^2)            ln(1-x^2)      
    !  erf^-1 (x) = sqrt( ---  -  --------   + sqrt((( ---  + -------  ))^2    -   --------- ))   
    !                     pi*a        2                pi*a       2                    a
    !
    !                                                                     
    !   where a = 8(pi-3)/(3*pi*(4-pi))                                                                      

    USE mo_ham_m7ctl,   ONLY: sigmaln
    USE mo_ham_streams, ONLY: rwet, rdry 
    USE mo_math_constants, ONLY: pi

    IMPLICIT NONE

    !--- Subroutine parameters: 

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow, kmod

    REAL(dp), INTENT(in) :: pcmr(kbdim,klev)
    REAL(dp), INTENT(in) :: pxie(kbdim,klev)

    REAL(dp), INTENT(out) :: pcritrad(kbdim,klev)

    !--- Local variables:  
  
    INTEGER  :: jl, jk

    LOGICAL :: ll1(kbdim, klev), ll2(kbdim,klev), ll3(kbdim,klev)

    REAL(dp) :: za_rcp, zb(kbdim, klev), zc(kbdim, klev), &
                zx2(kbdim, klev), zy(kbdim, klev), zcritrad(1:kbdim, klev), &
                zpre_fact, zfact(kbdim, klev)

    za_rcp = (3._dp * pi * (4._dp-pi)) / (8._dp * (pi-3._dp))

    zpre_fact = SQRT(2._dp)*sigmaln(kmod)

    pcritrad(1:kproma,:) = 0._dp !initialization

    ll1(1:kproma,:) = (ABS(pxie(1:kproma,:)) < 1._dp)
    ll2(1:kproma,:) = (pxie(1:kproma,:) >= 1._dp)
    ll3(1:kproma,:) = ll1(1:kproma,:) .AND. (pxie(1:kproma,:) < 0._dp)

    zx2(1:kproma,:) = pxie(1:kproma,:)**2
    zx2(1:kproma,:) = MERGE(zx2(1:kproma,:), 0._dp, ll1(1:kproma,:)) !to avoid illegal operations
    zb(1:kproma,:)  = 2.0_dp/pi*za_rcp + (0.5_dp * LOG(1.0_dp-zx2(1:kproma,:)))       
    zc(1:kproma,:)  = LOG(1.0_dp-zx2(1:kproma,:))*za_rcp
    zy(1:kproma,:)  = SQRT(-zb(1:kproma,:) + SQRT(zb(1:kproma,:)**2-zc(1:kproma,:))) 
    zy(1:kproma,:)  = MERGE(-zy(1:kproma,:), zy(1:kproma,:), ll3(1:kproma,:))

    zfact(1:kproma,:)    = zpre_fact * zy(1:kproma,:)
    zcritrad(1:kproma,:) = EXP(zfact(1:kproma,:)) * pcmr(1:kproma,:)

    pcritrad(1:kproma,:) = MERGE(zcritrad(1:kproma,:), pcritrad(1:kproma,:), ll1(1:kproma,:))

    ! Minimal scavenging of the mode by using an artificial large critical radius:
    pcritrad(1:kproma,:) = MERGE(500.e-6_dp, pcritrad(1:kproma,:), ll2(1:kproma,:))

  END SUBROUTINE ham_m7_invertlogtail
!<<gf

! ------------------------------------------------------------------------------------------------------

  SUBROUTINE scavcoef_bilinterp(kproma, kbdim,  klev,  krow, ktop, &
                        pfprecip, pmr, &
                        X1, X2, Y1, Y2, &
                        Q11, Q12, Q21, Q22, &
                        pscavcoef)

    ! *scavcoef_bilinterp* performs bilinear interpolation of below-cloud
    ! aerosol scavenging coefficients
    !
    !     Author:
    !     -------
    !     Betty Croft, Dalhousie University           2008
    !
    !     Contributors:
    !     -------------
    !     Sylvaine Ferrachat, ETHZ, 2011 (cleanup + reinforced security)
    !

    !--- Subroutine parameters:

    INTEGER, INTENT(in)  :: kproma, kbdim,  klev,  krow, ktop
    REAL(dp), INTENT(in) :: Q11(kbdim,klev)
    REAL(dp), INTENT(in) :: Q12(kbdim,klev)
    REAL(dp), INTENT(in) :: Q21(kbdim,klev)
    REAL(dp), INTENT(in) :: Q22(kbdim,klev)
    REAL(dp), INTENT(in) :: pfprecip(kbdim,klev) ! precip rate (rain or snow)
    REAL(dp), INTENT(in) :: pmr(kbdim,klev)      ! aerosol wet radius
    REAL(dp), INTENT(in) :: X1(kbdim,klev)
    REAL(dp), INTENT(in) :: X2(kbdim,klev)
    REAL(dp), INTENT(in) :: Y1(kbdim,klev)
    REAL(dp), INTENT(in) :: Y2(kbdim,klev)

    REAL(dp), INTENT(out) :: pscavcoef(kbdim,klev)

    ! --- Local variables

    REAL(dp)  :: ztmp1(kbdim,klev), ztmp2(kbdim,klev)
    REAL(dp)  :: ztmp3(kbdim,klev)
    LOGICAL :: lint1(kbdim,klev),  lint2(kbdim,klev)
    LOGICAL :: lint3(kbdim,klev),  lint4(kbdim,klev)

! ----------------------------------------------------------------------------------------

    pscavcoef(1:kproma,:) = 0._dp

    lint1(1:kproma,:) = (X1(1:kproma,:) /= X2(1:kproma,:)) .AND. (Y1(1:kproma,:) == Y2(1:kproma,:))
    lint2(1:kproma,:) = (X1(1:kproma,:) == X2(1:kproma,:)) .AND. (Y1(1:kproma,:) /= Y2(1:kproma,:))
    lint3(1:kproma,:) = (X1(1:kproma,:) == X2(1:kproma,:)) .AND. (Y1(1:kproma,:) == Y2(1:kproma,:))
    lint4(1:kproma,:) = (X1(1:kproma,:) /= X2(1:kproma,:)) .AND. (Y1(1:kproma,:) /= Y2(1:kproma,:))

! --- Interpolation in x direction only (rain rates) -------------------------------------

    ztmp1(1:kproma,:) = MERGE( &
                             X2(1:kproma,:)-X1(1:kproma,:), &
                             1._dp, & !SF dummy value
                             lint1(1:kproma,:))

    ztmp1(1:kproma,:) = &
        (((X2(1:kproma,:)-pfprecip(1:kproma,:))/ztmp1(1:kproma,:))*Q11(1:kproma,:)) &
       +(((pfprecip(1:kproma,:)-X1(1:kproma,:))/ztmp1(1:kproma,:))*Q21(1:kproma,:))

    pscavcoef(1:kproma,:) = MERGE(ztmp1(1:kproma,:), pscavcoef(1:kproma,:), lint1(1:kproma,:))

! --- Interpolation in y direction only (aerosol radii) ----------------------------------

    ztmp1(1:kproma,:) = MERGE( &
                             Y2(1:kproma,:)-Y1(1:kproma,:), &
                             1._dp, & !SF dummy value
                             lint2(1:kproma,:))

    ztmp1(1:kproma,:) = &
        (((Y2(1:kproma,:)-pmr(1:kproma,:))/ztmp1(1:kproma,:))*Q21(1:kproma,:))   &
       +(((pmr(1:kproma,:)-Y1(1:kproma,:))/ztmp1(1:kproma,:))*Q22(1:kproma,:))

    pscavcoef(1:kproma,:) = MERGE(ztmp1(1:kproma,:), pscavcoef(1:kproma,:), lint2(1:kproma,:))

! --- No interpolation of below-cloud scavenging coefficents -----------------------------

    ztmp1(1:kproma,:) = Q11(1:kproma,:)

    pscavcoef(1:kproma,:) = MERGE(ztmp1(1:kproma,:), pscavcoef(1:kproma,:), lint3(1:kproma,:))

! --- Bilinear interpolation of below-cloud scavenging coefficients ----------------------             

    ztmp1(1:kproma,:) = MERGE( &
                             X2(1:kproma,:)-X1(1:kproma,:), &
                             1._dp, & !SF dummy value
                             lint4(1:kproma,:))

    ztmp2(1:kproma,:) = MERGE( &
                             Y2(1:kproma,:)-Y1(1:kproma,:), &
                             1._dp, & !SF dummy value
                             lint4(1:kproma,:))

    ztmp3(1:kproma,:) =                                          &
      (((Y2(1:kproma,:)-pmr(1:kproma,:))/ztmp2(1:kproma,:))*                   &
       ((((X2(1:kproma,:)-pfprecip(1:kproma,:))/ztmp1(1:kproma,:))*Q11(1:kproma,:))     &
       +((pfprecip(1:kproma,:)-X1(1:kproma,:))/ztmp1(1:kproma,:))*Q21(1:kproma,:))) +   &
      (((pmr(1:kproma,:)-Y1(1:kproma,:))/ztmp2(1:kproma,:))*                   &
       ((((X2(1:kproma,:)-pfprecip(1:kproma,:))/ztmp1(1:kproma,:))*Q12(1:kproma,:))     &
       +((pfprecip(1:kproma,:)-X1(1:kproma,:))/ztmp1(1:kproma,:))*Q22(1:kproma,:)))

    pscavcoef(1:kproma,:) = MERGE(ztmp3(1:kproma,:), pscavcoef(1:kproma,:), lint4(1:kproma,:))

   END SUBROUTINE scavcoef_bilinterp

END MODULE mo_ham_tools
