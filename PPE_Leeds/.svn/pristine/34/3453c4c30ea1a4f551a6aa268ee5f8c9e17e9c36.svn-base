!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_sedimentation.f90
!!
!! \brief
!! Calculate sedimentation rates for aerosol tracers in the HAM model
!!
!! \author Michael Schulz (LSCE)
!! \author P. Stier (MPI-Met)
!!
!! \responsible_coder
!! [ John Doe, john.doe@blabla.com -Compulsory- ]
!!
!! \revision_history
!!   -# Michael Schulz (LSCE) - original code (2001-02-05)
!!   -# P. Stier (MPI-Met) - adapted to ECHAM5-HAM (2002)
!!   -# Martin grav. Schultz (FZ Juelich) - adapted to ECHAM6-HAMMOZ (integrated in modular structure) (2009)
!!
!! \limitations
!! None
!!
!! \details
!! This module is derived from the former xt_sedimentation routine (ECHAM5-HAM).
!!
!! \bibliographic_references
!!    - Seinfeld &Pandis (1998)
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

MODULE mo_ham_sedimentation

  USE mo_kind,             ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ham_prep_sedi, ham_sedimentation


  CONTAINS

!! The subroutine ham_prep_sedi calculates tracer-independent physical parameters 
!! needed in the ham_sedimentation routine

  SUBROUTINE ham_prep_sedi(kproma, kbdim, klev,       &
                           pt,     pq,    pap,  paph, & 
                           ptempc,                    &
                           pvis,   plair, prho,       &
                           pdpg,   pdz                )

  USE mo_physical_constants, ONLY: rd, vtmpc1, grav, tmelt

  !--- parameters
  INTEGER, INTENT(in)    :: kproma, kbdim, klev           ! indices
  REAL(dp), INTENT(in)   :: pt(kbdim, klev),       &      ! temperature 
                            pq(kbdim, klev),       &      ! specific humidity 
                            pap(kbdim, klev),      &      ! full level pressure
                            paph(kbdim, klev+1)           ! half level pressure
  REAL(dp), INTENT(out)  :: ptempc(kbdim, klev),     &    ! temp. above melting
                            pvis(kbdim, klev),       &    ! air viscosity
                            plair(kbdim, klev),      &    ! mean free path length
                            prho(kbdim, klev),       &    ! air density
                            pdpg(kbdim, klev),       &    ! layer thickness (pressure)
                            pdz(kbdim, klev)              ! layer thickness (length)

  !--- local variables
  INTEGER  :: jk
  REAL(dp) :: ztmp1(kbdim,klev) !SF #458   

  !--- code

  !--- Temperature at t=t+1:
  ptempc(1:kproma,:) = pt(1:kproma,:) - tmelt

  !---  Dynamic viscosity of air after Prup.Klett: pvis [Pa s]
  !>>SF #458 (replacing where statements)
  ztmp1(1:kproma,:) = 1.718_dp + 0.0049_dp*ptempc(1:kproma,:)

  pvis(1:kproma,:)  = MERGE( &
                           ztmp1(1:kproma,:)*1.E-5_dp, &
                           (ztmp1(1:kproma,:) - 1.2E-05_dp*(ptempc(1:kproma,:)**2))*1.E-5_dp, &
                           (ptempc(1:kproma,:) >= 0._dp))
  !<<SF #458 (replacing where statements)

  !--- Mean free path of air after Prupp. Klett: plair [10^-6 m]
  plair(1:kproma,:) = 0.066_dp *(1.01325E+5_dp/pap(1:kproma,:))       &
                      * (pt(1:kproma,:)/293.15_dp)*1.E-06_dp

  !--- Density of air:
  prho(1:kproma,:)=pap(1:kproma,:)/(pt(1:kproma,:)*rd*(1._dp+vtmpc1*pq(1:kproma,:)))

  !--- Air mass auxiliary variable:
  pdpg(1:kproma,1)=2._dp*(paph(1:kproma,2)-pap(1:kproma,1))/grav
  DO jk=2, klev
     pdpg(1:kproma,jk)=(paph(1:kproma,jk+1)-paph(1:kproma,jk))/grav
  END DO

  !--- Layer thickness:
  pdz(1:kproma,:)=pdpg(1:kproma,:)/prho(1:kproma,:)

  END SUBROUTINE ham_prep_sedi


!! The subroutine ham_sedimentation calculates the sedimentation rate for one
!! tracer. The tracer loop is contained in sedi_interface.

  SUBROUTINE ham_sedimentation(kproma, kbdim, klev, krow, & 
                               kt, pvis, plair, prho,     & 
                               pdpg, pdz,                 & 
                               pxtp1, pxtte,              &
                               pvsedi, psediflux)

  USE mo_ham_streams,   ONLY: rwet, densaer
  USE mo_time_control,  ONLY: time_step_len
  USE mo_tracdef,       ONLY: trlist, AEROSOLNUMBER, AEROSOLMASS
  USE mo_physical_constants, ONLY: grav
  USE mo_ham_m7ctl,     ONLY: sigma, sigmaln, cmedr2mmedr
  USE mo_exception,     ONLY: message_text, message, em_error
  USE mo_ham,           ONLY: nham_subm,HAM_M7,HAM_SALSA,HAM_BULK

  !--- parameters
  INTEGER, INTENT(in)    :: kproma, kbdim, klev, krow, kt ! indices
  REAL(dp), INTENT(in)   :: pvis(kbdim, klev),    &       ! air viscosity
                            plair(kbdim, klev),   &       ! mean free path
                            prho(kbdim, klev),    &       ! air density
                            pdpg(kbdim, klev),    &       ! layer thickness (pressure)
                            pdz(kbdim, klev),     &       ! layer thickness (length)
                            pxtp1(kbdim, klev)            ! updated tracer(kt) concentration
  REAL(dp), INTENT(inout):: pxtte(kbdim, klev)            ! tracer(kt) tendency
  REAL(dp), INTENT(out)  :: pvsedi(kbdim, klev),  &       ! sedimentation velocity
                            psediflux(kbdim, klev)        ! sedimentation flux


  !--- local variables
  INTEGER            :: imod

  LOGICAL :: ll1(kbdim,klev) !SF #458 dummy logical

  REAL(dp)           :: slinnfac
  REAL(dp)           :: zmd(kbdim, klev),   & !
                        zsedtend(kbdim, klev) ! tracer tendency
  REAL(dp), POINTER  :: rwet_p(:,:,:), densaer_p(:,:,:)

  REAL(dp) :: ztmp1(kbdim,klev), ztmp2(kbdim,klev) !SF #458 dummy temporary variables

  !--- code
  !--- initialize output
  pvsedi(:,:)    = 0._dp
  psediflux(:,:) = 0._dp

  !--- lookup mode of the tracer and get wet radius and density
  imod=trlist%ti(kt)%mode
  rwet_p      => rwet(imod)%ptr
  densaer_p   => densaer(imod)%ptr

  !--- Select diameter and limit it to maximal 50 um:
  !       If tracer is:
  !       aerosol number mixing ratio: use number median radius
  !       aerosol mass   mixing ratio: use mass   median radius

  !TB
  ! With SALSA bins are considered monodisperse so no need for cmedr2mmedr

  SELECT CASE(nham_subm)
      CASE(HAM_M7)
         IF (trlist%ti(kt)%nphase==AEROSOLNUMBER) THEN
            zmd(1:kproma,:)=MIN( rwet_p(1:kproma,:,krow)*2._dp ,  50.E-6_dp)
         ELSE IF (trlist%ti(kt)%nphase==AEROSOLMASS) THEN
            zmd(1:kproma,:)=MIN( rwet_p(1:kproma,:,krow)*cmedr2mmedr(imod)*2._dp , 50.E-6_dp)
         ELSE
            WRITE(message_text,'(a)') 'unexpected tracer phase'
            CALL message('ham_sedimentation', message_text, level=em_error)
         END IF
    
         !--- Slinn correction for sedimentation velocity of a
         !    size distribution with a given sigma:
         slinnfac=sigma(imod)**(2._dp*sigmaln(imod))     
    
      CASE(HAM_SALSA)
         zmd(1:kproma,:)=MIN( rwet_p(1:kproma,:,krow)*2._dp ,  50.E-6_dp)
         
         !TB
         !slinnfac not needed
         slinnfac=1.0_dp
  END SELECT

!>>SF #458 (replacing where statements)
  ll1(1:kproma,:) = (zmd(1:kproma,:) > 0._dp)

  ztmp1(1:kproma,:) = MERGE(zmd(1:kproma,:), 1._dp, ll1(1:kproma,:)) !SF protection against division by 0 further below
                                                                     !   1._dp is a dummy value
  !--- Stokes-velocity (S&P, Equation 8.42):
  ztmp2(1:kproma,:) = 2._dp/9._dp*(densaer_p(1:kproma,:,krow)-prho(1:kproma,:))  &
                       * grav/pvis(1:kproma,:)*(ztmp1(1:kproma,:)/2._dp)**2._dp

  !--- With Cunnigham slip- flow correction (S&P, Equation 8.34):
  SELECT CASE(nham_subm)
      CASE(HAM_M7)
          ztmp2(1:kproma,:) = ztmp2(1:kproma,:)                                    &
!>>DN #328
                             * (slinnfac +                                         &
                                1.246_dp*2._dp*plair(1:kproma,:)                   &
                                /ztmp1(1:kproma,:)*exp((0.5_dp*sigmaln(imod)**2._dp)))
!<<DN #328
      CASE(HAM_SALSA) !SF #328 (keep previous formulation for SALSA)
          ztmp2(1:kproma,:) = ztmp2(1:kproma,:)                                           &
                            * (1._dp+ 1.257_dp*plair(1:kproma,:)/ztmp1(1:kproma,:)*2._dp  &
                              + 0.4_dp*plair(1:kproma,:)/ztmp1(1:kproma,:)*2._dp          &
                                *EXP(-1.1_dp/(plair(1:kproma,:)/ztmp1(1:kproma,:)*2._dp)) )
  END SELECT

  !--- Calculate sedimentation in terms of mixing ratio tendencies:
  !--- Limit loss to the content of box:
  !    Multilayer crossing is not realised here as sedimentation
  !    velocity is in effect limited to dz/dt (grid velocity).
  !    (no multilayer crossing)
  ztmp2(1:kproma,:)  = MIN( ztmp2(1:kproma,:) , pdz(1:kproma,:)/time_step_len )

  pvsedi(1:kproma,:) = MERGE(ztmp2(1:kproma,:), 0._dp, ll1(1:kproma,:))

  !--- Loss in terms of mixing ratio tendency:
  ztmp1(1:kproma,:) = MAX(0._dp, pxtp1(1:kproma,:)*pvsedi(1:kproma,:)/pdz(1:kproma,:) )
  !--- Limit tendency to pxtp1/dt to avoid negative tracer concentrations
  ztmp1(1:kproma,:) = MIN(ztmp1(1:kproma,:), pxtp1(1:kproma,:)/time_step_len)

  zsedtend(1:kproma,:) = MERGE(ztmp1(1:kproma,:), 0._dp, ll1(1:kproma,:))
  
  !--- Apply loss throughout the column of mixing ratio tendencies:
  !    note: pxtte in this routine is already indexed for one tracer
  !    (normally this would be pxtte(1:kproma, :, kt) )

  pxtte(1:kproma,:)= pxtte(1:kproma,:)-zsedtend(1:kproma,:)

  !--- Transfer loss from mixing ratio tendency [kg/kg]
  !    to sedimentation flux [kg m-2 s-1]:

  psediflux(1:kproma,:) = zsedtend(1:kproma,:) * pdpg(1:kproma,:)
!<<SF #458 (replacing where statements)

  !--- Re-convert sedimentatation flux and add it to the mixing ratio tendency
  !    of the box below (conversion with zdpg of the box below):
  pxtte(1:kproma,2:klev)=pxtte(1:kproma,2:klev)                           &
                         + (psediflux(1:kproma,1:(klev-1))/pdpg(1:kproma,2:klev))

  END SUBROUTINE ham_sedimentation

END MODULE mo_ham_sedimentation

