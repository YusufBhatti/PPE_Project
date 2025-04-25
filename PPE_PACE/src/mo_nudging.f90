!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!BOI
!
! !TITLE:        Programmers Guide for Nudging
! !AUTHORS:      Ingo Kirchner
! !AFFILIATION:  Max Planch Institute for Meteorology, Hamburg
! !DATE:         August 2002
! !INTRODUCTION: External interface and structure of the implementation

! The nudging packages provides procedures for the simple data
! assimilation called "nudging". The nudging is activated with
! $LNUDGE=true$ in the namellist {\it RUNCTL} and all features of the
! nudging are controled by the namelist {\it NDGCTL}.
! 
! As minimum the user have to define the names of the input data
! files. The files contain the reference fields separated for each
! variable and the corresponding sea surface fields. Such data files can
! be prepared e.g. from the ERA15 data set using the interpolation
! program "intera". 

!\begin{verbatim}
! *** Nudging Interface to ECHAM4 ***
!
! *ECHAM*    *NUDGING*
!
!
! (1) initialization procedure
!
! INICTL ---> NudgingInit(NDG_INI_TIME)
!             NudgingInit(NDG_INI_STREAM)
! CONTROL --> NudgingInit(NDG_INI_MEM)
!             NudgingInit(NDG_CLEAN_MEM)
!
!
! (2) nudging procedure
!
! STEPON ---> NudgingReadSST                read new sst field
!         |
!       ....
!
!              GPC ------> NudgingSSTnew    reset sst at latitude circle
!
!       ....
!         |
!         +-> Nudging                       perform nudging
!         |   +--->GetNudgeData             read nudging data sets
!         |         +-->ReadOneBlock        read one data time step
!         |              +-->OpenOneBlock   open new data block
!         |              +-->CloseBlock     close data block
!         |
!         +-> NdgCorrBuffer                 correct accumulated nudging diagnostics
!
!
! --------------------------------------------------------------------
! April-2001 reorganization of modules
!
! layer1:
!          mo_nudgig
!          mo_nudging_init
!
! layer2: 
!          mo_nudging_sst
!          mo_nudging_io
!          mo_nudging_pattern
!
! layer3:
!          mo_nudging_utils
!
! layer4: 
!          mo_nudging_buffer
!          mo_nudging_constants
!
! --------------------------------------------------------------------
!\end{verbatim}

!EOI

!======================================================================
MODULE mo_nudging
!BOP

  ! !MODULE: mo_nudging (layer 1)

  ! !DESCRIPTION: 
  !         implementation of the nudging kernel procedure

  ! !REVISION HISTORY: 

  ! J. Feichter    Uni Hamburg    Jan 91 and Apr 93
  ! W. May         DMI-Copenhagen Mar 97
  ! M. Stendel     MPI-Hamburg    Sep 96 and May 97
  ! H.-S. Bauer    MPI-Hamburg    Jul 98
  ! I. Kirchner    MPI-Hamburg    Nov 98 and Oct/Dec 99, March 2000
  ! I. Kirchner    MPI-Hamburg    May 2000
  ! I. Kirchner    MPI-Hamburg    Aug 2000, upgrade echam5 release 5
  ! I. Kirchner/H.-S. Bauer MPI   Sep 2000, nudging sst setup
  ! I. Kirchner    MPI-Hamburg    Nov 2000, date/time control
  ! I. Kirchner    MPI-Hamburg    Apr 2001, revision
  ! R. Johanni,    IPP Garching,  May 2002, parallel version
  ! I. Kirchner    MPI-Hamburg    Aug 2002, revision

!EOP

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: Nudging

  CHARACTER(len=256) :: mess

CONTAINS

  !======================================================================

!BOP

  ! !IROUTINE:  Nudging

  ! !INTERFACE:

  SUBROUTINE Nudging

    ! !DESCRIPTION: 

    ! The 'nudging' routine for ECHAM4/5 to adjust the following 
    ! meteorological variables to observations by means of
    ! Newtonian relaxation:
    !
    !\begin{itemize}
    !\item divergence
    !\item vorticity
    !\item temperature
    !\item log surface pressure
    !\end{itemize}
    !
    ! The features were documented in {\it mo\_nudging\_constants.f90}
    ! 
    ! The method following Krishnamurti et al. (1991), Tellus 43AB, 53-81.
    ! This is "pure" Krishnamurti, which means that there are no other
    ! relaxation terms on rhs! see Eq. 7.3:
    !
    ! $X_{t+dt} = (X1_{t+dt} +2dt*N * X2_{t+dt}) / (1+2dt*N)$
    !
    ! where  $X1_{t+dt}$ is a predicted value of $X_{t+dt}$ prior to nudging.
    ! $X2_{t+dt}$ is a future value which the nudging is aimed at
    ! (in the twin case that is the observed value at the later time).
    ! $N$ is the nudging coefficient.
    ! The implementation of the method was modified (I.Kirchner) using the 
    ! following linear combination with modified coeffients A and B.
    !
    ! $ X = A * X1 + B * X2$
    !
    !\begin{description}
    !\item[implicit]
    !     $A = \frac{1}{1-2*dt*N}$
    !     $B = \frac{2*dt*N}{1-2*dt*N}$
    !
    !\item[explicit]
    !     $A = 1 - 2*dt*N$
    !     $B = 2*dt*N$
    !\end{description}

    ! !USES:

    !* basic functions
    USE mo_kind,              ONLY: dp
    USE mo_exception,         ONLY: message

    !* nudging features
    USE mo_nudging_constants, ONLY: &
         lnudgdbx, lnudgpat, ltintlin, lnudgfrd, lsite, lnudgwobs, &
         linp_div, linp_vor, linp_tem, linp_lnp, &
         nudgdsize, nudgdamp, ldamplin, lnudg_run, &
         nudgd, nudgv, nudgt, nudgp

    !* special features
    USE mo_nudging_pattern,   ONLY: Nudg_Update_Alpha
    USE mo_nmi,               ONLY: NMI_Make, NMI_MAKE_FMMO

    !* memory management
    USE mo_decomposition,     ONLY: ldc => local_decomposition
    USE mo_control,           ONLY: lnmi, nlev, nlevp1
    USE mo_memory_sp,         ONLY: sd, svo, stp
    USE mo_nudging_buffer,    ONLY: &
         sdten, svten, stten, ssten, sdtac, svtac, sttac, sstac, nssfc, &
         isano, itano, idano, ivano, asano, atano, adano, avano, &
         sdobs1, svobs1, stobs1, sdobs2, svobs2, stobs2, &
         sdobs,  svobs,  stobs, sdsite, svsite, stsite, &
         sd_o_n0, sv_o_n0, st_o_n0, sd_o_n1, sv_o_n1, st_o_n1, &
         sd_m_n0, sv_m_n0, st_m_n0, sd_m_n1, sv_m_n1, st_m_n1, &
         b_pat_div, b_pat_vor, b_pat_tem, b_pat_lnp, &
         nudgdo, nudgvo, nudgto, nudgpo, flagn, &
         nudgda, nudgva, nudgta, nudgsa, &
         nudgdb, nudgvb, nudgtb, nudgsb, &
         buf_ref0, buf_ref1, buf_ref2, buf_ref3, buf_ref, &
         isiteaccu, lsite_n0, lsite_n1, indg_accu

    !* i/o functions and synchronisation
    USE mo_time_control,      ONLY: get_time_step, out_convert_date,&
         next_date, delta_time, skip_nudging_read
    USE mo_nudging_io,        ONLY:  GetNudgeData, &
         lndgstep3, ndgstep0, ndgstep1, ndgstep2, ndgstep3

!EOP
!BOC
!BOX
    REAL(kind=dp)    :: wobs, wmod, rtwodt, tfact0, tfact1, tfact2, tfact3
    REAL(kind=dp)    :: tu0, ts0, tu1, ts1
    REAL(kind=dp)    :: dt0, dt1, dt2, dt3, dt4, dtx1, dtx2, dtx3
    INTEGER :: istep, my_day, jk, sheadd, sheadt
    LOGICAL :: lbetaprn

    INTEGER, EXTERNAL :: util_cputime

    istep = get_time_step()

    IF (.NOT. lnudg_run) THEN
       IF (lnudgdbx) THEN
          WRITE(mess,*) 'skip nudging for Step = ',istep
          CALL message('Nudging',mess)
       END IF
       RETURN
    END IF

    CALL out_convert_date(next_date,sheadd,sheadt)

    IF (lnudgdbx) THEN
       IF (util_cputime(tu0,ts0) == -1) THEN
          CALL message('Nudging','WARNING: Cannot determine used CPU time')
       END IF
    END IF

    ! *** initalize instantaneous tendencies
    sdten(:,:,:) = 0.0_dp
    svten(:,:,:) = 0.0_dp
    stten(:,:,:) = 0.0_dp
    ssten(:,:,:) = 0.0_dp
   
    rtwodt = 0.5_dp/delta_time

!EOX
    ! *** read nudging data arrays, minimum of four time steps *****************
!BOX
    IF (.NOT. (lndgstep3 .AND. skip_nudging_read())) CALL GetNudgeData


    ! --------------------------------------------------------------------------
    !
!EOX
    ! *** time interpolation and damping of nudging strength *******************
!BOX
    dt0 = REAL(ndgstep2 - ndgstep0,dp)
    dt1 = REAL(ndgstep2 - ndgstep1,dp)
    dt2 = REAL(ndgstep3 - ndgstep1,dp)

    dt3 = dt0/dt1
    dt4 = dt2/dt1

    dtx1 = REAL(istep+1 - ndgstep1,dp)/dt1

    ! define weights for linear interpolation
    tfact2 = dtx1
    tfact1 = 1._dp-tfact2

    IF (lnudgpat) THEN

       ! pattern correlation ***************************************************
       ! calculates b_pat_xxx = a_pat_xxx + NORM(model,obs)

       CALL Nudg_Update_Alpha
 
       ! set lbetaprn, print at first time in the month at 12:00 UTC
       my_day   = MOD(sheadd,100)
       lbetaprn = (sheadt == 120000)
 
       IF (linp_div) THEN
          DO jk=1,nlev
             sdobs(jk,:,:) = sdobs1(jk,:,:) * tfact1 * b_pat_div(jk,1) &
                           + sdobs2(jk,:,:) * tfact2 * b_pat_div(jk,2)
          END DO
          IF (lbetaprn) THEN
             WRITE(mess,'(a,i3,40f12.4)') 'BETA(div)1= ',my_day,b_pat_div(:,1)
             CALL message('Nudging',mess)
             WRITE(mess,'(a,i3,40f12.4)') 'BETA(div)2= ',my_day,b_pat_div(:,2)
             CALL message('Nudging',mess)
          END IF
       END IF

       IF (linp_vor) THEN
          DO jk=1,nlev
             svobs(jk,:,:) = svobs1(jk,:,:) * tfact1 * b_pat_vor(jk,1) &
                           + svobs2(jk,:,:) * tfact2 * b_pat_vor(jk,2)
          END DO
          IF (lbetaprn) THEN
             WRITE(mess,'(a,i3,40f12.4)') 'BETA(vor)1= ',my_day,b_pat_vor(:,1)
             CALL message('Nudging',mess)
             WRITE(mess,'(a,i3,40f12.4)') 'BETA(vor)2= ',my_day,b_pat_vor(:,2)
             CALL message('Nudging',mess)
          END IF
       END IF
 
       IF (linp_tem) THEN
          DO jk=1,nlev
             stobs(jk,:,:) = stobs1(jk,:,:) * tfact1 * b_pat_tem(jk,1) &
                           + stobs2(jk,:,:) * tfact2 * b_pat_tem(jk,2)
          END DO
          IF (lbetaprn) THEN
             WRITE(mess,'(a,i3,40f12.4)') 'BETA(tem)1= ',my_day,b_pat_tem(:,1)
             CALL message('Nudging',mess)
             WRITE(mess,'(a,i3,40f12.4)') 'BETA(tem)2= ',my_day,b_pat_tem(:,2)
             CALL message('Nudging',mess)
          END IF
       END IF
       
       IF (linp_lnp) THEN
          stobs(nlevp1,:,:) = stobs1(nlevp1,:,:) * tfact1 * b_pat_lnp(1) &
                            + stobs2(nlevp1,:,:) * tfact2 * b_pat_lnp(2)
          IF (lbetaprn) THEN
             WRITE(mess,'(a,i3,40f12.4)') 'BETA(lnp)1= ',my_day,b_pat_lnp(1)
             CALL message('Nudging',mess)
             WRITE(mess,'(a,i3,40f12.4)') 'BETA(lnp)2= ',my_day,b_pat_lnp(2)
              CALL message('Nudging',mess)
           END IF
        END IF

        IF (lnudgdbx) THEN
           WRITE(mess,'(a,2(f8.4,1x),a,i10)') &
               'time interpolation done, factors ... ',tfact1,tfact2,&
               ' at NSTEP+1= ',istep+1
           CALL message('Nudging',mess)
        END IF



    ELSE IF (ltintlin) THEN
       ! linear time interpolation *********************************************

       buf_ref(:,:,:,:) = buf_ref1(:,:,:,:)*tfact1 &
                        + buf_ref2(:,:,:,:)*tfact2

       IF (lnudgdbx) THEN
          WRITE(mess,'(a,2(f8.4,1x),a,i10)') &
               'time interpolation done, factors ... ',tfact1,tfact2,&
               ' at NSTEP+1= ',istep+1
          CALL message('Nudging',mess)
       END IF



    ELSE
       ! non-linear time interpolation, use polynom 3rd order solution *********

       dtx2   = dtx1 * dtx1
       dtx3   = dtx2        * dtx1
       tfact0 = ( -dtx2 +2._dp*dtx1 -1._dp )                              * dtx1 / dt3
       tfact1 = ( (2._dp*dt4-1._dp)*dtx3 +(1._dp-3._dp*dt4)*dtx2 +dt4 )          / dt4
       tfact2 = ( (1._dp-2._dp*dt3)*dtx2 +(3._dp*dt3-2._dp)*dtx1 +1._dp ) * dtx1 / dt3
       tfact3 = ( dtx1 -1._dp )                                           * dtx2 / dt4

       buf_ref(:,:,:,:) = buf_ref0(:,:,:,:)*tfact0 &
                        + buf_ref1(:,:,:,:)*tfact1 &
                        + buf_ref2(:,:,:,:)*tfact2 &
                        + buf_ref3(:,:,:,:)*tfact3

       IF (lnudgdbx) THEN
          WRITE(mess,'(a,4(f8.4,1x),a,i10)') &
               'time interpolation done, factors ... ',tfact0,tfact1,tfact2,tfact3,&
               ' at NSTEP+1= ',istep+1
          CALL message('Nudging',mess)
       END IF

    END IF



    ! --------------------------------------------------------------------------
    !
!EOX
    ! calculates weights inside the time window between two reference dates ****
!BOX
    IF ((nudgdsize < dtx1) .AND. (dtx1 < (1._dp-nudgdsize)) ) THEN
       wobs = nudgdamp
    ELSE
       IF (dtx1 > 0.5_dp ) dtx1 = 1._dp-dtx1    ! left side mirrored
       IF (nudgdsize <= 0.0_dp) THEN
          wobs = 0.0_dp
       ELSE IF (ldamplin) THEN
          ! linear time damping weight
          wobs = 1._dp - dtx1 * (1._dp-nudgdamp)/nudgdsize
       ELSE
          ! non -linear damping
          dtx1 = dtx1/nudgdsize             ! rescale x-axis
          wobs = dtx1*dtx1*(nudgdamp-1._dp)*(3._dp-2._dp*dtx1)+1._dp
       END IF
    END IF

    wobs = MAX(MIN(wobs,1.0_dp),0.0_dp)
    wmod = 1._dp-wobs
    IF (lnudgdbx) THEN
       WRITE(mess,'(a,f7.4,a,f7.4)') &
            'recombination weights ... Reference = ',wobs,' Model = ',wmod
       CALL message('Nudging',mess)
    END IF


    ! --------------------------------------------------------------------------
    !
!EOX
    ! *** Normal Mode Filter with time interpolated data ***********************
!BOX
    IF (lnmi .AND. (.NOT.lnudgfrd)) CALL NMI_Make(NMI_MAKE_FMMO)


    ! --------------------------------------------------------------------------
    !
    ! *** setup coefficient matrix *********************************************
    !
    ! inside nudging region
    !
    !   A = (nudg? + wmod*nudg?o)      B = wobs*nudg?o
    !
    ! outside nudging region
    !
    !   A = 1.0_dp                     B = 0.0_dp

    DO jk=1,nlev    ! 3-d fields
       nudgta(jk,:,:) = 1 + flagn(jk,:,:)*(nudgt(jk) + wmod*nudgto(jk) - 1)
       nudgtb(jk,:,:) =     flagn(jk,:,:)*             wobs*nudgto(jk)

       nudgda(jk,:,:) = 1 + flagn(jk,:,:)*(nudgd(jk) + wmod*nudgdo(jk) - 1)
       nudgdb(jk,:,:) =     flagn(jk,:,:)*             wobs*nudgdo(jk)

       nudgva(jk,:,:) = 1 + flagn(jk,:,:)*(nudgv(jk) + wmod*nudgvo(jk) - 1)
       nudgvb(jk,:,:) =     flagn(jk,:,:)*             wobs*nudgvo(jk)
    END DO
    DO jk=1,nssfc   ! surface fields
       nudgsa(jk,:,:) = 1 + flagn(nlevp1,:,:)*(nudgp + wmod*nudgpo - 1)
       nudgsb(jk,:,:) =     flagn(nlevp1,:,:)*         wobs*nudgpo
    END DO


    ! --------------------------------------------------------------------------
    !   Set the complex part of the global means to 0.

    IF (ldc%sm(1) == 0 .AND. ldc%snn0(1) == 0) buf_ref(:,2,1,:) = 0.0_dp

    !---------------------------------------------------------------------------
!EOX
    !  Now perform Newtonian relaxation.
!BOX
    !
    !  Apply full nudging for all scales. For T106 run, it may be
    !  necessary to decrease weight for smaller scales.
    !  Originally, there was full nudging up to T42 (NSP<=946), 50%
    !  at T63 (NSP=2080) and no nudging for NSP>2628 (T71).
    !
    ! --------------------------------------------------------------------------
    !   Tendencies are defined as the difference between the 
    !   nudged and the original fields
    !
    !        T_tendency := T_nudging  - T_original,    ergo 
    !   >>>  T_nudging   = T_original + T_tendency  <<<
    ! --------------------------------------------------------------------------
    
    ! store old model values
    sdten( :, :,:) = sd (      :,      :,:)
    svten( :, :,:) = svo(      :,      :,:)
    stten( :, :,:) = stp(      :nlev,  :,:)
    ssten(1:1,:,:) = stp(nlevp1:nlevp1,:,:)

    ! calculate new model value as linear combination of old model value
    ! and observations
    ! calculate in nudging window (nudglmin:nudglmax,nudgmin:nudgmax)
    !
!EOX
    ! NEW = A * OLD + B * OBS
!BOX
    stp(nlevp1,:,:)   = nudgsa(1,:,:)*stp(nlevp1,:,:) + nudgsb(1,:,:)*stobs(nlevp1,:,:)
    stp(:nlev ,:,:)   = nudgta(:,:,:)*stp(:nlev, :,:) + nudgtb(:,:,:)*stobs(:nlev, :,:)
    sd (:,     :,:)   = nudgda(:,:,:)*sd (:,     :,:) + nudgdb(:,:,:)*sdobs(:,     :,:)
    svo(:,     :,:)   = nudgva(:,:,:)*svo(:,     :,:) + nudgvb(:,:,:)*svobs(:,     :,:)
    
    ! --------------------------------------------------------------------------
    !
    ! compose diagnostics with flag-field
    !
    ! INSIDE
    ! store nudging term as a tendency term, (NEW-OLD)/2DT in UNIT/sec
    !
    ! nudging region update
    sdten(:,:,:) = flagn(:nlev, :,:) *(sd (:,     :,:) - sdten(:,:,:))*rtwodt
    svten(:,:,:) = flagn(:nlev, :,:) *(svo(:,     :,:) - svten(:,:,:))*rtwodt
    stten(:,:,:) = flagn(:nlev, :,:) *(stp(:nlev, :,:) - stten(:,:,:))*rtwodt
    ssten(1,:,:) = flagn(nlevp1,:,:) *(stp(nlevp1,:,:) - ssten(1,:,:))*rtwodt

    sdtac(:,:,:) = sdtac(:,:,:) + sdten(:,:,:)*delta_time
    svtac(:,:,:) = svtac(:,:,:) + svten(:,:,:)*delta_time
    sttac(:,:,:) = sttac(:,:,:) + stten(:,:,:)*delta_time
    sstac(:,:,:) = sstac(:,:,:) + ssten(:,:,:)*delta_time

    ! OUTSIDE
    ! calculate difference outside the nudging range, absolut value
    ! MODEL - OBSERVED in UNIT
    !
    ! additional diagnostics
    IF (lnudgwobs) THEN
       idano(:,:,:) = (1-flagn(:nlev, :,:))*(sd (:,     :,:) - sdobs(:,     :,:))
       ivano(:,:,:) = (1-flagn(:nlev, :,:))*(svo(:,     :,:) - svobs(:,     :,:))
       itano(:,:,:) = (1-flagn(:nlev, :,:))*(stp(:nlev, :,:) - stobs(:nlev, :,:))
       isano(1,:,:) = (1-flagn(nlevp1,:,:))*(stp(nlevp1,:,:) - stobs(nlevp1,:,:))

       adano(:,:,:) = adano(:,:,:) + idano(:,:,:)*delta_time
       avano(:,:,:) = avano(:,:,:) + ivano(:,:,:)*delta_time
       atano(:,:,:) = atano(:,:,:) + itano(:,:,:)*delta_time
       asano(:,:,:) = asano(:,:,:) + isano(:,:,:)*delta_time

    END IF

    indg_accu = indg_accu + 1   !ik temporary
!EOX    
    IF (lsite) THEN

      ! accumulate SITEs *****************************************************
!BOX
       IF (lsite_n0) THEN

          !   fraction is  (X[n+1]-X[n-1])/2dt - Residue - (O[n+1]-O[n-1])/2dt

          sdsite(:,     :,:) = sdsite(:,     :,:) - flagn(:nlev, :,:)*sdten(:,:,:) &
               + ( sd   (:,:,:) - sd_m_n0(:,:,:) )*rtwodt &
               - ( sdobs(:,:,:) - sd_o_n0(:,:,:) )*rtwodt

          svsite(:,     :,:) = svsite(:,     :,:) - flagn(:nlev, :,:)*svten(:,:,:) &
               + ( svo  (:,:,:) - sv_m_n0(:,:,:) )*rtwodt &
               - ( svobs(:,:,:) - sv_o_n0(:,:,:) )*rtwodt

          stsite(:nlev, :,:) = stsite(:nlev, :,:) - flagn(:nlev, :,:)*stten(:,:,:) &
               + ( stp  (:nlev,:,:) - st_m_n0(:nlev,:,:) )*rtwodt &
               - ( stobs(:nlev,:,:) - st_o_n0(:nlev,:,:) )*rtwodt

          stsite(nlevp1,:,:) = stsite(nlevp1,:,:) - flagn(nlevp1,:,:)*ssten(1,:,:) &
               + ( stp  (nlevp1,:,:) - st_m_n0(nlevp1,:,:) )*rtwodt &
               - ( stobs(nlevp1,:,:) - st_o_n0(nlevp1,:,:) )*rtwodt

          isiteaccu = isiteaccu + 1

       END IF

       ! store in SITE detection mode observed/modelled values from the last step
       IF (lsite_n1) lsite_n0 = .TRUE.
       sd_o_n0 = sd_o_n1;  sv_o_n0 = sv_o_n1;  st_o_n0 = st_o_n1
       sd_m_n0 = sd_m_n1;  sv_m_n0 = sv_m_n1;  st_m_n0 = st_m_n1

       ! store present time level
       lsite_n1 = .TRUE.
       sd_o_n1 = sdobs; sv_o_n1 = svobs; st_o_n1 = stobs
       sd_m_n1 = sd;    sv_m_n1 = svo;   st_m_n1 = stp

    END IF

    IF (lnudgdbx) THEN
       IF (util_cputime(tu1, ts1) == -1) THEN
          CALL message('Nudging','WARNING: Cannot determine used CPU time')
       ELSE
          WRITE (mess,'(a,i10,f10.3,a)') 'Performance at NSTEP ',&
               istep,(tu1+ts1)-(tu0+ts0),'s (user+system)'
          CALL message('Nudging',mess)
       END IF
    END IF

  END SUBROUTINE Nudging

  !======================================================================
!EOX
!EOC
END MODULE mo_nudging
