!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_chemistry
!!
!! \brief
!!  This module connects the submodel interface of ECHAM to the MOZART chemistry driver
!!
!! \author S. Rast (MPI Hamburg)
!! \author Martin G. Schultz (FZ-Juelich)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - first original version for ECHAM5-MOZ (s.rast, 2004)
!!  - cleaned merge with HAMMONIA version from H. Schmidt (m.schultz, Aug 2008)
!!  - draft ECHAM6-HAMMOZ version (m.schultz, 2010-05-03)
!!  - cleanup of diagnostics (m.schultz, 2014-03-19)
!!
!! \limitations
!!  none
!!
!! \details
!!  In moz_chemistry we prepare everything for the actual chemistry driver (chemdr) routine
!!  of MOZART and we handle the diagnostics resulting from the chemical transformations.
!!  MOZ_CHEMISTRY is the handover from the ECHAM world (kproma, etc.) to the MOZ world (plonl, etc.).
!!  MOZ-internal diagnostics such as rate coefficients, production and loss rates, OH concentration
!!  are calculated in chemdr (this is the only reason why kproma is passed down into chemdr).
!!  MOZ_CHEMISTRY calls the more sophisticated diagnostics which require, for example, global averaging.
!!
!! \bibliographic_references
!!  none
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_moz_chemistry


USE mo_kind,       ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: moz_chemistry

  !-- private module variables


CONTAINS


SUBROUTINE moz_chemistry (kproma, kbdim, klev, klevp1, ktrac, krow,          &
                          paph,   pap,   pt,   ptte,                         &
                          pq,     pqte,  pxl,  pxlte,  pxi,   pxite,         &
                          pxt,    pxtte,                                     &
                          pgeom1, pgeohm1,     ptsurf, paclc,                &
                          pfrl,   pfrw,  pfri, pfrg                         )
!
! Description:
!
! calls the MOZART chemistry driver CHEMDR after some necessary preparations
! and calls chemistry diagnostics which go beyond archiving of variabes from within chemdr.
!
! *moz_chemistry* is called from *physc_subm_3 in mo_submodel_interface.f90
!


  USE mo_kind,            ONLY: dp
  USE mo_exception,       ONLY: finish, message
  USE mo_control,         ONLY: ltimer
  USE mo_timer,           ONLY: timer_start,          &
                                timer_stop
  USE mo_hammoz_timer,    ONLY: timer_moz_chem_prep,  &
                                timer_moz_chem_diag 
  USE mo_decomposition,   ONLY: ldc => local_decomposition
  USE mo_time_control,    ONLY: current_date, time_step_len
  USE mo_time_conversion, ONLY: TC_convert, TC_get, time_intern
  USE mo_memory_base,     ONLY: t_stream, get_stream, get_stream_element
  USE mo_tracdef,         ONLY: trlist
!++ost: these diagnostics are currently missing, introduce later
! USE mo_diag_sfc,        ONLY: set_sfcstream
! USE mo_diag_chemistry,  ONLY: set_chemstream
  USE mo_moz_diag,        ONLY: dpo3netchem
  USE mo_physical_constants, ONLY: rgrav, amd, cpd
  USE mo_memory_g3a,      ONLY: geospm
!++lp: add acdnc for fastj
  USE mo_memory_g3b,      ONLY: acdnc
!--lp
  USE mo_vphysc,          ONLY: vphysc
  USE mo_submodel,        ONLY: id_moz, submlist, lhammonia, lhmzhet, linterh2o
  USE mo_boundary_condition, ONLY: bc_apply
  USE mo_geoloc,          ONLY: ilat 
  USE mo_vphysc,          ONLY: vphysc
  USE mo_moz_uvalbedo,    ONLY: ibc_uvalbedo, moz_uvalbedo
  USE mo_moz_chemdr,      ONLY: chemdr
  USE mo_moz_util,        ONLY: Rgas, mwair, get_spc_ndx
  USE mo_moz,             ONLY: lstrathet, ibc_sad, lphotolysis, e2m
  USE mo_moz_mods,        ONLY: pcnstm1
!++ost: currently missing, introduce later, !baustelle!
! USE mo_hammoz,          ONLY: hammoz_het_chemistry
!--ost
  USE mo_moz_lbc_ubc,     ONLY: moz_lbc_ubc
  USE mo_debugs

  IMPLICIT NONE
  
  ! Scalar arguments
  INTEGER, INTENT(in)     :: kproma, kbdim, klev, klevp1, ktrac, krow
  
  ! Array arguments
  REAL(dp), INTENT(in)    :: paph(kbdim,klevp1),     &
                             pap(kbdim,klev),        &
                             pt(kbdim,klev),         &
                             pq(kbdim,klev),         &
                             pxl(kbdim,klev),        &
                             pxi(kbdim,klev),        &
                             pxt(kbdim,klev,ktrac),  &
                             pgeom1(kbdim,klev),     &
                             pgeohm1(kbdim,klevp1),  &
                             ptsurf(kbdim),          &
                             paclc(kbdim,klev),      &
                             pfrl(kbdim),            &   ! land fraction
                             pfrw(kbdim),            &   ! water fraction
                             pfri(kbdim),            &   ! (sea) ice fraction
                             pfrg(kbdim)                 ! glacier (land ice) fraction

  REAL(dp), INTENT(inout) :: ptte(kbdim,klev),       &   ! temperature tendency (HAMMONIA)
                             pqte(kbdim,klev),       &   ! humidity tendency (strat. chemistry)
                             pxlte(kbdim,klev),      &   ! cloud water tendency (not yet used!)
                             pxite(kbdim,klev),      &   ! cloud ice tendency (not yet used!)
                             pxtte(kbdim,klev,ktrac)     ! tracer tendency

  ! Local variables
!++mgs20140319: removed because diagnostics were moved
!!  REAL(dp), PARAMETER:: secperyear = 365._dp*24._dp*3600._dp  ! seconds per year
!--mgs

  REAL(dp):: zaph(kbdim,klevp1),            &
             zap(kbdim,klev),               &
             zt(kbdim,klev),                &
             zq(kbdim,klev),                &
             zq_save(kbdim,klev),           &
             zx(kbdim,klev),                &
             zxi(kbdim,klev),               &
             zmmr(kbdim,klev,pcnstm1),      &
             zmmr_save(kbdim,klev,pcnstm1), &
             ztsurf(kbdim),                 &
             zaclc(kbdim,klev),             &
             zpdel(kbdim,klev),             &
             zma(kbdim,klev),               &
             zmah(kbdim,klevp1),            &
             zamu(kbdim,klev),              &
             zcp(kbdim,klev),               &
             zsad_sage(kbdim,klev),         &  ! stratospheric aerosol density from SAGE climatology
             zalbedo(kbdim),                & 
             zhch(kbdim,klev),              &  ! chemical heating rate
             zhrad(kbdim,klev),             &  ! radiative heating rate
             zacdnc(kbdim,klev)                !++lp cloud droplet number concenration for fastj  

!++mgs20140319: removed because diagnostics were moved
!!  !--diagnostics
!!  REAL(dp):: zkch4(kbdim,klev),             &  ! reaction rate constant OH+CH4
!!             zkmcl(kbdim,klev),             &  ! reaction rate constant OH+CH3CCl3
!!             zrate(kbdim,klev)                 ! buffer for reaction turnover
!--mgs

  REAL(dp):: zdum,                          &
             zdum1(kbdim),                  &
             zdum2(kbdim,klev),             &
             zdum2p(kbdim,klevp1),          &
             zdumzamu, zdumdayl

  INTEGER :: jk, jkk, jlon, jt
  INTEGER :: imztrac, ndx_o3

  REAL(dp), POINTER       :: tmdiag_p(:,:,:)
  TYPE(t_stream), POINTER :: tmdiag

  ! Timer
  TYPE(time_intern) :: mydate
  INTEGER :: yymmdd, hhmmss, sec, day

  !-- Initialize 
  IF (ltimer) CALL timer_start( timer_moz_chem_prep )
  CALL TC_convert(current_date, mydate)
  CALL TC_get(mydate, yymmdd, hhmmss)
  CALL TC_get(current_date, day, sec)

  zhch(:,:) = 0._dp
  zhrad(:,:) = 0._dp

  ! Set absolute geopotential height
  DO jk=1,klev
     jkk=klev-jk+1
     zma(1:kproma,jk)=(pgeom1(1:kproma,jk)+geospm(1:kproma,krow))*rgrav
  ENDDO
  DO jk=1,klevp1
     jkk=klev-jk+1
     zmah(1:kproma,jk)=(pgeohm1(1:kproma,jk)+geospm(1:kproma,krow))*rgrav
  ENDDO

  !-- Get SAGE climatology 
  IF (lstrathet) THEN
    CALL bc_apply(ibc_sad, kproma, krow, zsad_sage)
  ELSE
    zsad_sage(:,:) = 0._dp
  ENDIF

  !-- obtain UV albedo
!++sschr #379: added wrapping lphotolysis IF
  IF (lphotolysis) THEN
    IF (ibc_uvalbedo > 0) THEN
      CALL moz_uvalbedo(kbdim, kproma, krow, pfrl, pfrw, pfri, pfrg)
      CALL bc_apply(ibc_uvalbedo, kproma, krow, zalbedo)
    ELSE
      CALL finish('moz_chemistry', 'No UV albedo boundary condition defined!')
    END IF
  END IF
!--sschr #379

  !-- Make local copy of physical fields (ap, aph, t, q, xl, xi)
  ! This is necessary, because MOZART always loops over plonl and is not "nproma-safe"
  ! Therefore, we make sure that all vectors are "complete" everywhere
  zaph(1:kproma,:) = paph(1:kproma,:)
  zap(1:kproma,:)  = pap(1:kproma,:)
  zt(1:kproma,:)   = pt(1:kproma,:)
!++mgs 20140404: also make sure that humidity is positive!
  zq(1:kproma,:)   = MAX(pq(1:kproma,:), 0._dp)
!--mgs
  zx(1:kproma,:)   = pxl(1:kproma,:) + pxi(1:kproma,:)
!++sschr 2014/04/02: xi from ECHAM has negative values!
  zxi(1:kproma,:)  = MAX(pxi(1:kproma,:),0._dp)
!--sschr 2014/04/02: xi from ECHAM has negative values!
  ztsurf(1:kproma) = ptsurf(1:kproma)
  zaclc(1:kproma,:)= paclc(1:kproma,:)
!++lp 
  zacdnc(1:kproma,:)=acdnc(1:kproma,:,krow)
!--lp
  !-- Assure positivity of tracer mixing ratios (copy only MOZ tracers!)
   DO jt = 1, ktrac
     imztrac = e2m(jt)
     IF (imztrac /= -1) THEN
       DO jk = 1, klev
         DO jlon = 1, kproma
           zmmr(jlon,jk,imztrac) = MAX(pxt(jlon,jk,jt),0.0_dp)
         END DO
       END DO
     END IF
   END DO
  !-- fill vectors for last block 
  IF (kproma < kbdim) THEN
    DO jlon=kproma+1,kbdim
      zaph(jlon,:)          = zaph(kproma,:)
      zap(jlon,:)           = zap(kproma,:)
      zt(jlon,:)            = zt(kproma,:)
      zq(jlon,:)            = zq(kproma,:)
      zx(jlon,:)            = zx(kproma,:)
      zxi(jlon,:)           = zxi(kproma,:)
      zmmr(jlon,:,:)        = zmmr(kproma,:,:)
      ztsurf(jlon)          = ztsurf(kproma)
      zaclc(jlon,:)         = zaclc(kproma,:)
      zma(jlon,:)           = zma(kproma,:)
      zmah(jlon,:)          = zmah(kproma,:)
      zsad_sage(jlon,:)     = zsad_sage(kproma,:)
      zalbedo(jlon)         = zalbedo(kproma)
      vphysc%rhoam1(jlon,:,krow) = vphysc%rhoam1(kproma,:,krow)      !++mgs
      zacdnc(jlon,:)        = zacdnc(kproma,:) !++lp
    ENDDO
  END IF

  ! from here on we can (and should) assign all arrays (1:kbdim) except for those 
  ! which are directly used from a module

  ! Make another copy of updated concentrations (for chemistry budget diagnostics)
!!mgs###!! IF (lbudgetdiag) THEN
  zmmr_save(:,:,:)=zmmr(:,:,:)
!!mgs###!! END IF
  zq_save(:,:) = zq(:,:)     ! save water vapour for tendency update (linterh2o)

  !-- Get air parameters (use updated tracer concs in HAMMONIA else constants)
  IF (lhammonia) THEN
!!mgs!!  CALL airparams( nproma, nbdim, nlev, ntrac, zmmr(:,:,:), zamu, zcp, n2_mmr)
    CALL finish('moz_chemistry', 'HAMMONIA needs subroutine airparams! Not implemented!')
  ELSE
    zamu(:,:) = amd
    zcp(:,:)  = cpd
!!mgs!! may need to set n2_mmr as well....
  END IF

  ! Set pressure increments
  DO jk=1,klev
     zpdel(:,jk)=zaph(:,jk+1)-zaph(:,jk)
  ENDDO
  
  !--- Call MOZART chemistry driver. This will initialize reaction rates and call the 
  !    chemical solver(s) needed to advance concentrations by one time step
  IF (ltimer) CALL timer_stop( timer_moz_chem_prep )    ! not quite true, but safer

  CALL moz_lbc_ubc(kproma,krow,zmmr)

       !  ECHAM names                                  ! MOZART names
  CALL chemdr( yymmdd,                  &              ! ncdate
       sec,                             &              ! ncsec
       kproma,                          &              ! kproma
       krow,                            &              ! lat
       zmmr(1:kbdim,:,:),               &              ! mmr     (inout)
       zq(1:kbdim,:),                   &              ! sh      (inout)
       zx(1:kbdim,:),                   &              ! cwat    (inout) 
       zxi(1:kbdim,:),                  &              ! cice    (inout)
       zhch(1:kbdim,:),                 &              ! cheat (temp. tendency due to chem. heating)
       zhrad(1:kbdim,:),                &              ! rheat (temp. tendency due to rad. heating)
       time_step_len,                   &              ! delt 
       zaph(1:kbdim,:),                 &              ! pfull
       zap(1:kbdim,:),                  &              ! pmid
       zpdel(1:kbdim,:),                &              ! pdel
       pgeom1(1:kbdim,:),               &              ! pgeom1
       ztsurf(1:kbdim),                 &              ! tsurf
       zma(1:kbdim,:),                  &              ! zma
       zmah(1:kbdim,:),                 &              ! zi
       zaclc(1:kbdim,:),                &              ! cldfr  
       zt(1:kbdim,:),                   &              ! tfld
       vphysc%rhoam1(1:kbdim,:,krow),   &              ! prho (air density)
       zsad_sage(1:kbdim,:),            &              ! stratospheric aerosol density (SAGE)
       zalbedo(1:kbdim),                &              ! albs   ! ++mgs: local albedo variable
       zamu(1:kbdim,:),                 &              ! zamu (molecular mass of air)
       zcp(1:kbdim,:),                  &              ! zcp (specific heat)
       zacdnc(1:kbdim,:) )                             ! cdnc (cloud droplet number conc. for fastj photolysis)
  ! Note on ozone field for photolysis: depending on lchemfeedback switch this will be either
  ! the climatological field from ECHAM or the chemistry field in mmr.

  IF (ltimer) CALL timer_start( timer_moz_chem_diag )

  !---------------------------------------------------------------------------------------------
  !     ... evaluate tendencies and interactions with other modules
  !---------------------------------------------------------------------------------------------

  ! calculate tracers' tendency
  DO jt = 1, ktrac
     imztrac = e2m(jt)
     IF (imztrac /= -1) THEN
       pxtte(1:kproma,:,jt) = pxtte(1:kproma,:,jt)                       &
                              + (zmmr(1:kproma,:,imztrac)-zmmr_save(1:kproma,:,imztrac))/time_step_len
     END IF
  END DO
  IF (associated(dpo3netchem)) THEN
    ndx_o3   = get_spc_ndx('O3')
    IF (ndx_o3 > 0) dpo3netchem(1:kproma, :, krow) = (zmmr(1:kproma,:,ndx_o3)-zmmr_save(1:kproma,:,ndx_o3))/time_step_len
  END IF

!++lp 13/11, call heterogeneous chemistry subroutines
  IF (lhmzhet) THEN
!    CALL hammoz_het_chemistry(kproma,     kbdim,      klev,     &
!                          krow,                             &
!		           ptm1,       pxtm1,      pxtte,    &
!       		   paphp1,     papp1,      rhoam1(:,:,krow))
  ENDIF
!--lp

!-- feedback of MOZART water to ECHAM
  IF (linterh2o) THEN 
    pqte(1:kproma,:) = pqte(1:kproma,:) + (zq(1:kproma,:)-zq_save(1:kproma,:))/time_step_len
    ! ### add diagnostics for water vapour change due to chemistry
    ! ### warning: if luse_p1_vars=false in mo_submodel_interface, then pqte=(zq-zq_save)/dt !!
    ! perhaps update tendencies for xl and xi as well in later versions??
    ! pxlte(...),   pxite(...)
  END IF

  !---------------------------------------------------------------------------------------------
  !    ... diagnostics
  !    Note that diagnostics which require knowledge of MOZART species index or reaction index
  !    are done in mo_moz_chemdr.
  !---------------------------------------------------------------------------------------------

!++mgs: make sure these diagnostics are "found" also in moz3 version !!
! surface ozone and NO2 diagnostic (if these species exist)
!++ost: diagnostics currently missing, introduce later
! CALL set_sfcstream ( kproma, kbdim, klev, ktrac, krow, pxtm1, pxtte )

! production and loss diagnostics
! CALL set_chemstream (zpdel, krow, kbdim, kproma, klev)
!--ost

! write mass diagnostics, in mass mixing ratio per second
  IF (ANY(trlist%ti(1:ktrac)%nbudget == 1)) THEN
!     zdum2(1:kproma,:)=(mwair*1.e-3_dp/Rgas)*papp1(1:kproma,:)/ptm1(1:kproma,:)
!     CALL get_stream(tmdiag,'tmdiag')
!     DO jt=1,ktrac
!          IF ( trlist%ti(jt)%nbudget == 1 ) THEN
!             CALL get_stream_element(tmdiag,TRIM(trlist%ti(jt)%basename)//'_chem',tmdiag_p)
!             tmdiag_p(1:kproma,1:klev,krow)=(zxt(1:kproma,1:klev,jt)-zxt_save(1:kproma,1:klev,jt))/time_step_len
!*zdum2(1:kproma,1:klev)
!          END IF
!       ENDDO
    END IF
  IF (ltimer) CALL timer_stop( timer_moz_chem_diag ) 

END SUBROUTINE moz_chemistry

END MODULE mo_moz_chemistry
