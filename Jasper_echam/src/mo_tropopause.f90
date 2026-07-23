#if defined (__SX__) || defined (ES)
#define VECTOR 1
#endif
#ifdef __xlC__
@PROCESS HOT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_tropopause

  !
  ! - Description:
  !
  !   This module contains parameters and routines needed for the
  !   WMO defined tropopause height.
  !
  ! Hauke: Is this routine HAMMONIA save?

  USE mo_kind, ONLY: wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: WMO_tropopause

CONTAINS

  SUBROUTINE WMO_tropopause(kproma, kbdim, klev, &
       ptm1, papm1, &
       ptropo, &
       ktrpwmo, ktrpwmop1)

    !
    !**** WMO_tropopause - calculation of the tropopause height.
    !
    !   programmed by Dieter Nodorp                 V1.0  March 95
    !   modified by Thomas Reichler                 V1.2  May 95
    !   Rewritten for use in ECHAM  Christine Land  V2.0  Jan 96
    !   Rewritten for use in ECHAM5 Monika Esch           Sep 2003
    !   Debugged                    Luis Kornblueh        Nov 2003        
    !
    !     purpose
    !     -------
    !     WMO_tropopause computes the tropopause height following the
    !     definition of the height of the tropopause as postulated
    !     by the WMO (1957).
    !
    !     interface
    !     ---------
    !     call WMO_tropopause(kproma, kbdim, klev, krow, ptm1, papm1, ptropo)
    !
    !          kproma: 'strip' elements
    !          kbdim:  'strip' length of array arguments 
    !          klev:   number of vertical levels
    !          ptm1:   temperature at full levels (t-dt)
    !          papm1:  pressure at full levels (t-dt)
    !          ptropo: height of the tropopause in Pa
    !
    !     method
    !     ------
    !   - p**k interpolation for the temperature Tm and
    !     the temperature gradient Gamma of each model layer
    !   - p**kappa interpolation for p-WMO
    !
    !     references
    !     ----------
    !
    !     WMO (1992): International meteorological vocabulary, Genf, 784pp.
    !
    !   " 1. The first tropopause is defined as the lowest level at which
    !        the lapse rate decreases to 2 deg C per kilometer or less,
    !        provided also the average lapse rate between this level and
    !        all higher levels within 2 kilometers does not exceed 2 deg C
    !        per kilometer."
    !
    !     Reichler (1995): Eine globale Klimatologie der Tropopausenhoehe
    !                 auf der Basis von ECMWF-Analysen, Diplomarbeit
    !                 Universitaet Augsburg, 145pp.
    !
    !     Dameris, M., D. Nodorp and R.Sausen, 1995: Correlation between
    !                 Tropopause Height Pressure and TOMS-Data for the
    !                 EASOE-Winter 1991/1992, Beitr. Phys. Atmosph., 68,
    !                 227-232.
    !
    !
    !     called by
    !     ---------
    !               physc
    !
    !     no externals
    !     ------------
    !
    !
    USE mo_physical_constants,    ONLY: grav, rd, cpd
    USE mo_echam_cloud_params,        ONLY: ncctop, nccbot
    USE mo_time_control, ONLY: lresume

    ! scalar arguments
    INTEGER, INTENT(in) :: kproma, kbdim, klev
    
    ! array arguments
    REAL(wp), INTENT(in)  :: ptm1(kbdim,klev), papm1(kbdim,klev)
    REAL(wp), INTENT(inout) :: ptropo(kbdim)

    INTEGER, INTENT(out) :: ktrpwmo(kbdim), ktrpwmop1(kbdim)

    ! local arrays
#ifdef VECTOR
    REAL(wp) :: ztropov(kproma)
#else
    REAL(wp) :: ztropo(kproma)
#endif

    REAL(wp) :: zpmk(kproma,klev), zpm(kproma,klev),   za(kproma,klev)
    REAL(wp) :: zb(kproma,klev),   ztm(kproma,klev),   zdtdz(kproma,klev)
    REAL(wp) :: zplimb(kproma), zplimt(kproma)

    REAL(wp) :: zptph, zp2km, zkappa, zzkap, zfaktor, zgwmo, zdeltaz
    REAL(wp) :: zag, zbg, zdp

    INTEGER :: jl, jk, jj 

    INTEGER :: iplimb     ! bottom level to search for the tropopause
    INTEGER :: iplimt     !   top  level to search for the tropopause

    LOGICAL, SAVE :: linit = .FALSE. ! for initialization
!$OMP THREADPRIVATE (linit)

    LOGICAL :: lo1, lo2, lo3

#ifdef VECTOR
    INTEGER :: n_still, n_still_new, n_still_still, n_still_aux
    INTEGER :: jlv
    INTEGER :: ind_jl(kproma), ind_jl_aux(kproma), kcountv(kproma)
    INTEGER :: ind_jl_temp(kproma), kcount_temp(kproma)
    LOGICAL :: jl_todo(kproma)

    REAL(wp) :: zptphv(kproma), zp2kmv(kproma), zasumv(kproma)
    REAL(wp) :: zptph_temp(kproma), zp2km_temp(kproma), zasum_temp(kproma)
#else
    INTEGER :: kcount

    REAL(wp) :: zasum, zamean
#endif
    REAL(wp) :: zpapm1(kproma,klev)

    IF (.NOT. lresume) THEN
      IF (.NOT. linit) THEN
        ptropo (1:kproma) = 20000.0_wp
        linit = .TRUE.
      END IF
    END IF
    
    zkappa   = rd/cpd                   ! 0.286
    zzkap    = 1.0_wp/zkappa
    zfaktor  = (-1.0_wp)*grav/rd           ! -9.81/287.0
    zgwmo    = -0.002_wp
    zdeltaz  = 2000.0_wp

    iplimb=nccbot
    iplimt=ncctop+2

    ! Calculate the height of the tropopause

    DO jl = 1, kproma
#ifdef VECTOR
      ztropov(jl) = -999.0_wp
#else
      ztropo(jl) = -999.0_wp
#endif
      zplimb(jl) = papm1(jl,iplimb)
      zplimt(jl) = papm1(jl,iplimt)
    ENDDO

    ! compute dt/dz

!LK may generate problem on NEC !cdir collapse
    DO jk = iplimt-2, iplimb+1
      DO jl = 1, kproma
        zpapm1(jl,jk)=papm1(jl,jk)**zkappa
      ENDDO
    ENDDO
!LK may generate problem on NEC !cdir collapse
    DO jk = iplimt-1, iplimb+1
      DO jl = 1, kproma

        ! ztm   lineare Interpolation in p**kappa
        ! gamma         dt/dp = a * kappa + papm1(jl,jk)**(kappa-1.)

        zpmk(jl,jk) = 0.5_wp*(zpapm1(jl,jk-1)+zpapm1(jl,jk))
        zpm(jl,jk)  = zpmk(jl,jk)**zzkap                   ! p centre
        za(jl,jk)   = (ptm1(jl,jk-1)-ptm1(jl,jk)) &
             /(zpapm1(jl,jk-1)-zpapm1(jl,jk))
        zb(jl,jk)   = ptm1(jl,jk)-(za(jl,jk)*zpapm1(jl,jk))
        ztm(jl,jk)  = za(jl,jk)*zpmk(jl,jk)+zb(jl,jk)      ! T centre
        
        zdtdz(jl,jk)=zfaktor*zkappa*za(jl,jk)*zpmk(jl,jk)/ztm(jl,jk)
      ENDDO
    ENDDO

#ifdef VECTOR

    n_still = kproma

    DO jl = 1, n_still
      ind_jl(jl) = jl
    ENDDO

    jk = iplimb+1

    ! vertical_loop:
    DO WHILE (n_still > 0 .AND. jk >= iplimt-1)

      n_still_new   = 0
      n_still_still = 0
      DO jlv = 1, n_still
        jl = ind_jl(jlv)
        ! First test: valid dt/dz ?
        IF (zdtdz(jl,jk) >  zgwmo .AND.  &     ! dt/dz > -2K/km
             zpm(jl,jk) <= zplimb(jl)) THEN    ! zpm not too low
          
          ! dtdz is valid - calculatye p_WMO by interpolation between
          ! actual and layer above linear in p^kappa (Dieters new method)

          n_still_still = n_still_still + 1
          ind_jl_aux(n_still_still) = ind_jl(jlv)
        ELSE
          n_still_new = n_still_new + 1
          ind_jl_temp(n_still_new) = ind_jl(jlv)
        ENDIF
      END DO

      n_still = n_still_still
 
      DO jlv = 1, n_still
        ind_jl(jlv) = ind_jl_aux(jlv)
      END DO
      n_still_still = 0
      DO jlv = 1, n_still
        jl = ind_jl(jlv)
        zag = (zdtdz(jl,jk)-zdtdz(jl,jk+1))/(zpmk(jl,jk)-zpmk(jl,jk+1))
        zbg = zdtdz(jl,jk+1)-zag*zpmk(jl,jk+1)
        IF ((zgwmo-zbg)/zag >= 0.0_wp) THEN
          zptph = ABS((zgwmo-zbg)/zag)**zzkap
        ELSE
          zptph = 0.0_wp
        ENDIF
        IF (zdtdz(jl,jk+1) >= zgwmo) zptph = zpm(jl,jk)
        zp2km = zptph+zdeltaz*zfaktor*zpm(jl,jk)/ztm(jl,jk) ! p at ptph + 2 km
        IF (zptph >= zplimt(jl) )THEN
          n_still_still = n_still_still + 1
          ind_jl_aux(n_still_still) = jl
          zptphv(n_still_still) = zptph
          zp2kmv(n_still_still) = zp2km
        ENDIF
      END DO

      jj = jk
      DO jlv = 1, n_still_still
        jl_todo(jlv)= .TRUE.
        zasumv(jlv) = 0.0_wp
        kcountv(jlv)= 0.0_wp
      END DO
      
      ! vertical_sub_loop
      n_still_aux = 0
      DO WHILE (n_still_still > 0 .AND.  jj >= iplimt-1)
        DO jlv = 1, n_still_still
          jl    = ind_jl_aux(jlv)
          zptph = zptphv(jlv)
          zp2km = zp2kmv(jlv)
          IF (zpm(jl,jj) <= zptph .AND.zpm(jl,jj) >= zp2km )THEN
            zasumv(jlv) = zasumv(jlv) + zdtdz(jl,jj)
            kcountv(jlv) = kcountv(jlv) + 1
            IF( zasumv(jlv) <= zgwmo* REAL(kcountv(jlv),wp) ) THEN
              n_still_new = n_still_new + 1
              ind_jl_temp(n_still_new) = ind_jl_aux(jlv)
              jl_todo(jlv) = .FALSE.
            ENDIF
          ELSE IF (zpm(jl,jj) < zp2km )THEN
            ztropov(jl) = zptph
            jl_todo(jlv) = .FALSE.
          ENDIF
          
        END DO
        
        n_still_aux = 0
        DO jlv = 1, n_still_still
          IF(jl_todo(jlv)) THEN
            n_still_aux = n_still_aux + 1
            ind_jl_temp(n_still_new+n_still_aux) = ind_jl_aux(jlv)
            zptph_temp(n_still_aux)  = zptphv(jlv)
            zp2km_temp(n_still_aux)  = zp2kmv(jlv)
            zasum_temp(n_still_aux)  = zasumv(jlv)
            kcount_temp(n_still_aux) = kcountv(jlv)
          ENDIF
        END DO
        
        n_still_still = n_still_aux
        
        DO jlv = 1, n_still_still
          ind_jl_aux(jlv) = ind_jl_temp(n_still_new+jlv)
          jl_todo(jlv)=.TRUE.
          zptphv(jlv) = zptph_temp(jlv)
          zp2kmv(jlv) = zp2km_temp(jlv)
          zasumv(jlv) = zasum_temp(jlv)
          kcountv(jlv)= kcount_temp(jlv)
        ENDDO
        
        jj = jj - 1
        
      END DO ! vertical_sub_loop
      
      n_still_new = n_still_new + n_still_aux
      n_still = n_still_new

      if ( n_still > kproma ) n_still = kproma
      
      DO jlv = 1, n_still
        ind_jl(jlv) = ind_jl_temp(jlv)
      END DO
      
      jk = jk - 1
      
    END DO ! vertical_loop

#else

    nproma_loop: DO jl = 1, kproma
      vertical_loop: DO jk = iplimb+1, iplimt-1, -1
        ! First test: valid dt/dz ?
        IF (zdtdz(jl,jk) >  zgwmo .AND.  &     ! dt/dz > -2K/km
              zpm(jl,jk) <= zplimb(jl)) THEN   ! zpm not too low

          ! dtdz is valid - calculatye p_WMO by interpolation between 
          ! actual and layer above linear in p^kappa (Dieters new method)

          zag = (zdtdz(jl,jk)-zdtdz(jl,jk+1))/(zpmk(jl,jk)-zpmk(jl,jk+1))  
          zbg = zdtdz(jl,jk+1)-zag*zpmk(jl,jk+1)                           
          IF ((zag/=0.0_wp).AND.((zgwmo-zbg)*zag > 0.0_wp)) THEN
            zptph = ((zgwmo-zbg)/zag)**zzkap
          ELSE
            zptph = 0.0_wp
          END IF
          IF (zdtdz(jl,jk+1) >= zgwmo) zptph = zpm(jl,jk)
          zp2km   = zptph+zdeltaz*zfaktor*zpm(jl,jk)/ztm(jl,jk)   ! p at ptph + 2km
          zasum   = 0.0_wp                                        ! zdtdz above
          kcount  = 0.0_wp                                        ! number of levels above
          ! 2nd test: dt/dz above 2 km must not be lower than -2 K/km
          IF (zptph >= zplimt(jl)) THEN
            vertical_sub_loop: DO jj = jk, iplimt-1, -1
              
              IF (zpm(jl,jj) <= zptph .AND. zpm(jl,jj) >= zp2km )THEN
                zasum = zasum+zdtdz(jl,jj)
                kcount= kcount+1
                zamean = zasum/REAL(kcount,wp)            ! dt/dz mean
                IF (zamean <= zgwmo) CYCLE vertical_loop  ! dt/dz above < 2K/km
              ELSE IF (zpm(jl,jj) < zp2km) THEN           ! ptropo valid
                ztropo(jl) = zptph
                EXIT vertical_loop
              ENDIF

            END DO vertical_sub_loop
            
          ENDIF
        ENDIF
      END DO vertical_loop
    END DO nproma_loop

#endif



    ! if tropopause not find use previous value in time


#ifdef VECTOR
    WHERE (ztropov(1:kproma) > 0.0_wp)
      ptropo(1:kproma) = ztropov(1:kproma)
   END WHERE
#else
    WHERE (ztropo(1:kproma) > 0.0_wp)
      ptropo(1:kproma) = ztropo(1:kproma)
   END WHERE
#endif

!LK    DO jl = 1, kproma
!LK       IF (ztropo(jl) /= ztropov(jl)) &
!LK            write (0,*) 'Inconsistent tropopause value: ', &
!LK                        jl, ztropo(jl), ztropov(jl)
!LK    ENDDO

    ! calculation of tropopause on model levels

    ktrpwmop1(:) = iplimt

    DO jk = iplimt-1, iplimb+1
      DO jl = 1, kproma

        zdp = 0.5_wp*(papm1(jl,jk+1)-papm1(jl,jk))

        lo1 =           papm1(jl,jk)       < ptropo(jl)
        lo2 = lo1 .AND. papm1(jl,jk+1)     > ptropo(jl)
        lo3 = lo2 .AND. papm1(jl,jk+1)-zdp > ptropo(jl)

        ktrpwmop1(jl)  = MERGE(jk+1, ktrpwmop1(jl), lo2)
        ktrpwmo(jl)    = MERGE(jk,   ktrpwmop1(jl), lo3)

      END DO
    END DO
    
  END SUBROUTINE WMO_tropopause
  
END MODULE mo_tropopause
