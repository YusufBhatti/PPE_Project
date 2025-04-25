#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE setphys

  !
  ! preset physical constants and process according namelist 
  !
  ! M. Jarraud, ECMWF, December 1982, original version
  ! A. Tompkins, MPI, April 2000, added LUTs for cover
  ! L. Kornblueh, MPI, April 2003, moved calculation of convection parameterization LUTs to mo_convect_tables
  ! U. Schlese, M.Esch June 2007, added volcanic aerosol
  !
  ! This subroutine preset ,modifies and derive
  ! constants in modules used in the parameterization
  ! of diabatic processes except the constants internal
  ! to the radiation "black box".
  !

  USE mo_kind,           ONLY: dp 
  USE mo_mpi,            ONLY: p_io, p_parallel, p_parallel_io, p_bcast  
  USE mo_control,        ONLY: nn, l_volc, ngl
  USE mo_param_switches, ONLY: lphys, lrad, lvdiff, lcond, lsurf,  &
                               lconv, lgwdrag, lice, lconvmassfix, &
                               lclmi_progn, nclmi_progn, ncd_activ, &
                               nic_cirrus, lsecprod, &
                               nauto, lorocirrus, ldyn_cdnc_min, & !SF #475
                               cdnc_min_fixed, & !SF #589
                               nactivpdf, & !ZK
                               lemuphase_nic_cirrus, eta_emu_nic_cirrus, & !UP emulator
                               lemuphase_riming, eta_emu_riming, & !UP emulator
                               lemuphase_wbf, eta_emu_wbf, & !UP emulator
                               lemuphase_icaggr, eta_emu_icaggr, & !UP emu. #757
                               lemuphase_icaccr, eta_emu_icaccr, & !UP emu. #775
                               lemuphase_cdaccr, eta_emu_cdaccr, & !UP emu. #805
                               lemuphase_cdautc, eta_emu_cdautc, & !UP emu. #805
                               lemuphase_sci, eta_emu_sci, & !UP emulator
                               lemuphase_dep, eta_emu_dep, & !UP emulator, #801
                               lemuphase_icnucl, eta_emu_icnucl, & !UP emulator, #801
                               lemuphase_sub_evp_falling, eta_emu_subfis, & !UP #802
                               eta_emu_evpr, & !UP #802
                               lemuphase_cdnuc, eta_emu_cdnuc, & !UP emulator, #805
                               lemuphase_sip, eta_emu_sip, & !UP emulator, #805
                               lemuphase_mlt, eta_emu_mlt, & !UP emulator, #802
                               ldetr_liquid, & !UP #782.1
                               ldetr_liquid, & !UP #782.1
                               convicfactor, ldetr_convicfactor, & !UP #782.2
                               lnewdiags, & !UP: new diags for CD activation correction
                               lccnclimdiags, & !UP: new diags for CCNclim
                               lpr_corr, & !UP: diagnostic process rates corrections
                               rcd_vol_max, & !UP: make this tunable in settings file
                               lsedfix, & !UP #797
                               lpr_corr, & !UP: diagnostic process rates corrections
                               lslf, & !UP #839
                               lcmpsimpl_prescr, & !UP #821
                               ncmpsimpl_prescr_rime, &
                               ncmpsimpl_prescr_icnucl, &
                               ncmpsimpl_prescr_subfis, &
                               ncmpsimpl_prescr_sci, &
                               ncmpsimpl_prescr_icaccr, &
                               ncmpsimpl_prescr_mlt, &
                               lemuphase_hetfrz, eta_emu_hetfrz, & !MA emulator #790
                               lemuphase_fvic, eta_emu_fvic, & !MA emulator #795
                               lctrl_WBF !MA switch WBF on/off
!davidn
  USE mo_param_switches, ONLY: tun47zinhomi,tun47ccraut,tun47ccsaut,tun47cdncmin,lepsfix
  USE mo_param_switches, ONLY: tun47entrpen,tun47cprcon,tun47entrscv,tun47cmfctop
  USE mo_param_switches, ONLY: tun47zinhoml1,tun47zinhoml2,tun47zinhoml3,tun47zinpar
!davidn
  USE mo_echam_conv_constants, ONLY: lmfpen, lmfscv, lmfmid, lmfdd, lmfdudv
  USE mo_echam_convect_tables, ONLY: init_convect_tables  
  USE mo_namelist,       ONLY: open_nml, position_nml, POSITIONED
  USE mo_exception,      ONLY: finish, message, message_text, em_error, em_warn, em_param, em_info
  USE mo_volc_data,      ONLY: jpd, aod, reff, extrat, ssa, asym, init_volc_tables, read_volc_data
  USE mo_surface_ice,    ONLY: init_albedo_ice

  USE mo_submodel,       ONLY: lccnclim, print_value !++mgs 

  IMPLICIT NONE

  INCLUDE 'physctl.inc'

  REAL(dp) :: zq, zx

  INTEGER :: it, iq
  INTEGER :: ierr, inml, iunit   ! error return value from position_nml

  REAL(dp), EXTERNAL :: betai

  !     ------------------------------------------------------------
  !
  !*        1.       preset constants.
  !                  ------ ----------
  !
  !*        1.2      preset constants in MO_PARAM_SWITCHES.
  !
  lphys   = .TRUE.
  lrad    = .TRUE.
  lvdiff  = .TRUE.
  lconv   = .TRUE.
  lcond   = .TRUE.
  lsurf   = .TRUE.
  lmfpen  = .TRUE.
  lmfscv  = .TRUE.
  lmfmid  = .TRUE.
  lmfdd   = .TRUE.
  lmfdudv = .TRUE.
  lgwdrag = .TRUE.
  lice          = .TRUE.
  lconvmassfix  = .TRUE.
  lclmi_progn   = .FALSE.
  ncd_activ     = 2        ! default scheme for cdnc activation (will be 0 if lclmi_progn=false)
  nactivpdf     = 0        ! default scheme for updraft PDF !ZK
  nic_cirrus    = 2        ! default scheme for cirrus scheme   (will be 0 if lclmi_progn=false)
  nauto         = 1        ! default scheme for autoconversion  (will be 0 if lclmi_progn=false)
  lsecprod      = .FALSE.  ! switch for computing secondary ice production !SF #251
  lorocirrus    = .FALSE.  ! switch for gravity waves updraft velocity for icnc (orographic cirrus clouds)
  nclmi_progn   = 3        ! parameter to fine-tune the 2-moment calculations in mo_cloud_micro_2m 
  !UPcomment: this is the default, don't change here but in the settings file
  ldyn_cdnc_min = .FALSE.  ! switch to turn on the dynamical setting of the min CDNC !SF #475
  cdnc_min_fixed = 40      ! fixed value for min CDNC, in cm-3 (used when ldyn_cdnc_min is FALSE)
                           ! Warning! So far only values of 40 or 10 are accepted.
                           ! SF #489
  rcd_vol_max = 19.e-6_dp  !maximum droplet volume radius [m]
                           !UP: moved this here to make it tunable via settings
                           !file 
!davidn
  tun47zinhomi = 0.7_dp
  tun47zinhoml1= 0.8_dp
  tun47zinhoml2= 0.4_dp
  tun47zinhoml3= 0.8_dp
  tun47zinpar  = 0.1_dp
  tun47ccraut  = 1._dp
  tun47ccsaut  = 900._dp
  !tun47cdncmin = 10.0E6_dp !UP: original line
  tun47cdncmin = 10.0_dp !UP: changed units to cm-3 here
  tun47entrpen = 2.0E-4_dp
  tun47cprcon  = 9.0E-4_dp
  tun47entrscv = 3.0E-3_dp
  tun47cmfctop = 0.2_dp
  lepsfix      = .FALSE. !#680
!davidn
  !
  !>>UP #782.1
  ldetr_liquid = .FALSE.
  !<<UP #782.1
  !>>UP #797
  lsedfix = .FALSE.
  !<<UP #797
  !>>UP #783.3; as this applies real changes to the model and results, I put
  !default to false
  lpr_corr = .FALSE.
  !<<UP #783.3
  !>>UP #782.2
  convicfactor = 5.0_dp
  ldetr_convicfactor = .FALSE.
  !<<UP #782.2
  !>>UP
  !                  Preset constants in mo_cmp_emulator
  lemuphase_nic_cirrus      = .FALSE. ! IMPORTANT: when set to true, nic_cirrus needs to be 2 so that everything necessary gets computed
  eta_emu_nic_cirrus        = 1.0_dp ! needs to be between 0 and 2, issue #766, #776
  lemuphase_riming          = .FALSE.
  eta_emu_riming            = 1.0_dp ! needs to be between 0 and 2, issue #766, #776
  lemuphase_wbf             = .FALSE.
  eta_emu_wbf               = 1.0_dp ! needs to be between 0 and 1, issue #766, #776
  lemuphase_icaggr          = .FALSE.
  eta_emu_icaggr            = 1.0_dp ! needs to be between 0 and 2, issue #766, #776
  lemuphase_icaccr          = .FALSE.
  eta_emu_icaccr            = 1.0_dp ! needs to be between 0 and 2, issue #766, #776
  lemuphase_cdaccr          = .FALSE.
  eta_emu_cdaccr            = 1.0_dp ! needs to be between 0 and 2, issue #805
  lemuphase_cdautc          = .FALSE.
  eta_emu_cdautc            = 1.0_dp ! needs to be between 0 and 2, issue #805
  lemuphase_sci             = .FALSE.
  eta_emu_sci               = 1.0_dp ! needs to be between 0 and 2, issue #766, #776
  lemuphase_dep             = .FALSE.
  eta_emu_dep               = 1.0_dp ! needs to be between 0 and 2, issue #801
  lemuphase_sub_evp_falling = .FALSE.
  eta_emu_subfis            = 1.0_dp ! needs to be between 0 and 2, issue #802
  eta_emu_evpr              = 1.0_dp ! needs to be between 0 and 2, issue #802
  lemuphase_icnucl          = .FALSE.
  eta_emu_icnucl            = 1.0_dp ! needs to be between 0 and 2, issue #801
  lemuphase_cdnuc           = .FALSE.
  eta_emu_cdnuc             = 1.0_dp ! needs to be between 0 and 2, issue #805
  lemuphase_sip             = .FALSE.
  eta_emu_sip               = 1.0_dp ! needs to be between 0 and 2, issue #805
  lemuphase_mlt             = .FALSE.
  eta_emu_mlt               = 1.0_dp ! needs to be between 0 and 2, issue #802
  lnewdiags                 = .FALSE.
  lccnclimdiags             = .FALSE.
  !<<UP
  !>>MA
  lemuphase_hetfrz          = .FALSE. ! issue #790
  eta_emu_hetfrz            = 1.0_dp ! needs to be between 0 and 2, issue #790
  lemuphase_fvic            = .FALSE. ! issue #795
  eta_emu_fvic              = 1.0_dp ! needs to be between 0 and 2, issue #795
  lctrl_WBF = .TRUE. ! WBF prozess
  !<<MA
  lslf                    = .FALSE. !UP #839

  !>>UP #821
  lcmpsimpl_prescr        = .FALSE.
  ncmpsimpl_prescr_rime   = 0
  ncmpsimpl_prescr_icnucl = 0
  ncmpsimpl_prescr_subfis = 0
  ncmpsimpl_prescr_sci    = 0
  ncmpsimpl_prescr_icaccr = 0
  ncmpsimpl_prescr_mlt = 0
  !<<UP
 
  !
  !*         1.3     Initialise lookup tables for CUADJTQ
  !
  CALL init_convect_tables
  !
  !*         1.5     Initialise lookup tables for aerosol radiation parameters
  !
  IF(l_volc) THEN
    ALLOCATE ( aod(ngl,0:jpd+1))
    ALLOCATE (reff(ngl,0:jpd+1))
    CALL message('','l_volc =.TRUE.  --> volcanic forcing on')
    IF (p_parallel_io) THEN    
      CALL init_volc_tables
    ENDIF
    IF (p_parallel) THEN
       CALL p_bcast (extrat, p_io)
       CALL p_bcast (ssa, p_io)
       CALL p_bcast (asym, p_io)
    ENDIF
    IF (p_parallel_io) THEN
       CALL read_volc_data
    ENDIF
    IF (p_parallel) THEN
       CALL p_bcast (aod, p_io)
       CALL p_bcast (reff, p_io)
    ENDIF
  ELSE
    CALL message('','l_volc =.FALSE.  --> volcanic forcing off')
  ENDIF
  !
  !     ------------------------------------------------------------
  !
  !*        2.       READ NAMELIST.
  !                  ---- ---------
  !
  IF (p_parallel_io) THEN
    inml = open_nml('namelist.echam')
     iunit = position_nml ('PHYSCTL', inml, status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
       READ (iunit, physctl)
     END SELECT
  ENDIF
  IF (p_parallel) THEN
     CALL p_bcast (lphys, p_io)
     CALL p_bcast (lrad, p_io)
     CALL p_bcast (lvdiff, p_io)
     CALL p_bcast (lcond, p_io)
     CALL p_bcast (lsurf, p_io)
     CALL p_bcast (lconv, p_io)
     CALL p_bcast (lmfpen, p_io)
     CALL p_bcast (lgwdrag, p_io)
     CALL p_bcast (lice, p_io)
     CALL p_bcast (lconvmassfix, p_io)
     CALL p_bcast (lclmi_progn, p_io)
     CALL p_bcast (nclmi_progn, p_io)
     CALL p_bcast (ncd_activ, p_io)
     CALL p_bcast (nactivpdf, p_io)
     CALL p_bcast (nic_cirrus, p_io)
     CALL p_bcast (lsecprod, p_io)
     CALL p_bcast (lorocirrus, p_io)
     CALL p_bcast (ldyn_cdnc_min, p_io)
     CALL p_bcast (cdnc_min_fixed, p_io)
     CALL p_bcast (nauto, p_io)
     !>>UP #782.1
     CALL p_bcast (ldetr_liquid, p_io)
     !<<UP #782.1
     !>>UP #797
     CALL p_bcast (lsedfix, p_io)
     !<<UP #797
     !>>UP #783.3
     CALL p_bcast (lpr_corr, p_io)
     !<<UP #783.3
     !>>UP #782.2
     CALL p_bcast (convicfactor, p_io)
     CALL p_bcast (ldetr_convicfactor, p_io)
     !<<UP #782.2
     CALL p_bcast (rcd_vol_max, p_io) !UP: make this tunable
     !>>UP emulator
     CALL p_bcast (lemuphase_nic_cirrus, p_io)
     CALL p_bcast (eta_emu_nic_cirrus, p_io)
     CALL p_bcast (lemuphase_riming, p_io)
     CALL p_bcast (eta_emu_riming, p_io)
     CALL p_bcast (lemuphase_wbf, p_io)
     CALL p_bcast (eta_emu_wbf, p_io)
     CALL p_bcast (lemuphase_icaggr, p_io)
     CALL p_bcast (eta_emu_icaggr, p_io)
     CALL p_bcast (lemuphase_icaccr, p_io)
     CALL p_bcast (eta_emu_icaccr, p_io)
     CALL p_bcast (lemuphase_cdaccr, p_io)
     CALL p_bcast (eta_emu_cdaccr, p_io)
     CALL p_bcast (lemuphase_cdautc, p_io)
     CALL p_bcast (eta_emu_cdautc, p_io)
     CALL p_bcast (lemuphase_sci, p_io)
     CALL p_bcast (eta_emu_sci, p_io)
     CALL p_bcast (lemuphase_dep, p_io)
     CALL p_bcast (eta_emu_dep, p_io)
     CALL p_bcast (lemuphase_icnucl, p_io)
     CALL p_bcast (eta_emu_icnucl, p_io)
     CALL p_bcast (lemuphase_sub_evp_falling, p_io)
     CALL p_bcast (eta_emu_subfis, p_io)
     CALL p_bcast (eta_emu_evpr, p_io)
     CALL p_bcast (lemuphase_cdnuc, p_io)
     CALL p_bcast (eta_emu_cdnuc, p_io)
     CALL p_bcast (lemuphase_sip, p_io)
     CALL p_bcast (eta_emu_sip, p_io)
     CALL p_bcast (lemuphase_mlt, p_io)
     CALL p_bcast (eta_emu_mlt, p_io)
     !<<UP emulator
     !>>UP #821
     CALL p_bcast (lcmpsimpl_prescr, p_io)
     CALL p_bcast (ncmpsimpl_prescr_rime, p_io)
     CALL p_bcast (ncmpsimpl_prescr_icnucl, p_io)
     CALL p_bcast (ncmpsimpl_prescr_subfis, p_io)
     CALL p_bcast (ncmpsimpl_prescr_sci, p_io)
     CALL p_bcast (ncmpsimpl_prescr_icaccr, p_io)
     CALL p_bcast (ncmpsimpl_prescr_mlt, p_io)
     !<<UP
     !>>UP: new diags for activation corrections
     CALL p_bcast (lnewdiags, p_io)
     ! Additional CCNclim diags for investigation
     CALL p_bcast (lccnclimdiags, p_io)
     !<<UP
     !>>MA #790 and #795
     CALL p_bcast (lemuphase_hetfrz, p_io)
     CALL p_bcast (eta_emu_hetfrz, p_io) 
     CALL p_bcast (lemuphase_fvic, p_io)
     CALL p_bcast (eta_emu_fvic, p_io)
     !<<MA   
     !>>MA switch WBF
     CALL p_bcast(lctrl_WBF, p_io)
     !<<MA
     !davidn
     CALL p_bcast (tun47zinhomi, p_io)   
     CALL p_bcast (tun47zinhoml1, p_io)   
     CALL p_bcast (tun47zinhoml2, p_io)   
     CALL p_bcast (tun47zinhoml3, p_io)   
     CALL p_bcast (tun47zinpar, p_io)   
     CALL p_bcast (tun47ccraut, p_io)   
     CALL p_bcast (tun47ccsaut, p_io)   
     CALL p_bcast (tun47cdncmin, p_io)   
     CALL p_bcast (tun47entrpen, p_io)   
     CALL p_bcast (tun47entrscv, p_io)   
     CALL p_bcast (tun47cmfctop, p_io)   
     CALL p_bcast (tun47cprcon, p_io)   
     CALL p_bcast (lepsfix, p_io)   
    !davidn
    !>>UP FOR-ICE
     CALL p_bcast (lslf, p_io)
    !<<UP FOR-ICE
  ENDIF
!
!     ------------------------------------------------------------
!
!*        3.       MODIFY CONSTANTS.
!                  ------ ----------
!
!*        3.1      MODIFY CONSTANTS IN *mo_param_switches*.
!
  IF(.NOT.lphys) THEN
     lrad=.FALSE.
     lvdiff=.FALSE.
     lgwdrag=.FALSE.
     lconv=.FALSE.
     lcond=.FALSE.
     lsurf=.FALSE.
     lice=.FALSE.
  ELSE
     CALL message('',' Cloud cover scheme: diagnostic (Sundqvist)')
  END IF
!
  IF(.NOT.lconv) THEN
     lmfpen=.FALSE.
     lmfscv=.FALSE.
     lmfmid=.FALSE.
     lmfdd=.FALSE.
     lmfdudv=.FALSE.
  ELSE
     CALL message('',' Convection: Nordeng (default)')
  ENDIF
!
  IF(.NOT.lmfpen) THEN
     lconv=.FALSE.
     lmfscv=.FALSE.
     lmfmid=.FALSE.
     lmfdd=.FALSE.
     lmfdudv=.FALSE.
  END IF

  CALL message('','---')
  WRITE(message_text,'(a,a,1x,l1)') 'Prognostic eq. for cloud droplet- and ice crystal-number', &
                                    ' concentration (lclmi_progn)', lclmi_progn
  CALL message('',message_text, level=em_param)

  IF (.NOT. lclmi_progn) THEN
    ncd_activ  = 0
    nic_cirrus = 0
    nauto      = 0       !SF added for security only
    lsecprod    = .FALSE. !SF added for security only (#251)
    lorocirrus = .FALSE. !GF added for security only
    ldyn_cdnc_min = .FALSE. !SF security
    nclmi_progn = 0       !SF added for security only
    CALL message('setphys', 'ncd_activ, nic_cirrus and nauto set to 0 because lclmi_progn=.FALSE.',  &
                 level=em_warn)
    CALL message('setphys', 'lsecprod set to false because lclmi_progn=.FALSE.', level=em_warn)
    IF (lccnclim) THEN
       CALL message('setphys', 'lclmi_progn must be .TRUE. when lccnclim is active!', level=em_error)
    ENDIF
  ELSE !SF lclmi_progn=true
    IF (ncd_activ <= 0) THEN
      CALL message('setphys', 'ncd_activ must be >0 for lclmi_progn=.TRUE.!', level=em_error)
    ELSE !SF ncd_activ > 0
      IF (nauto <= 0) THEN
        CALL message('setphys', 'nauto must be >0 for ncd_activ>0!', level=em_error)
      ENDIF
      IF (nic_cirrus <= 0) THEN
        CALL message('setphys', 'nic_cirrus must be >0 for ncd_activ>0!', level=em_error)
      END IF
    END IF
  END IF

!gf-security check
  IF (lorocirrus .AND. nic_cirrus /= 2) THEN
     CALL message('setphys', 'lorocirrus=.TRUE. only possible with nic_cirrus=2', level=em_error)
  END IF
!gf-security check

!>>UP #782.1
  CALL message('','---')
  CALL print_value('Only liquid detrainment at mixed-phase temperatures (ldetr_liquid) = ', ldetr_liquid)
!<<UP #782.1
!>>UP #783.3
  CALL message('','---')
  CALL print_value('Correct process rate diagnostic for realistic correction terms (lpr_corr) = ', lpr_corr)
!<<UP #783.3
!>>UP #782.2
  CALL message('','---')
  CALL print_value('Cirrus tuning of immediate aggregation of detrained ice with updraft (ldetr_convicfactor) = ',&
         ldetr_convicfactor)
  CALL print_value('Cirrus tuning of immediate aggregation of detrained ice with updraft (convicfactor) = ',&
         convicfactor)
!<<UP #782.2
!>>UP #797
  CALL message('','---')
  CALL print_value('Restrict ice crystal sedimentation amount (lsedfix) = ', lsedfix)
!<<UP #797
!>>UP emulator
  CALL message('','---')
  CALL print_value('Emulator phase in/out (lemuphase_nic_cirrus) [outdated, needs checking] = ', lemuphase_nic_cirrus)
  ! WRITE(message_text,'(a,a,1x,l1)') 'Emulator phase in/out (lemuphase_nic_cirrus)', lemuphase_nic_cirrus
  !CALL message('',message_text, level=em_param)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_nic_cirrus) = ', eta_emu_nic_cirrus)
  !WRITE(message_text,'(a,a,1x,l1)') 'Emulator phase in/out (eta_emu_nic_cirrus)', eta_emu_nic_cirrus
  !CALL message('',message_text, level=em_param)
  CALL message('','---')
  CALL print_value('nic_cirrus = ', nic_cirrus)
  IF (lemuphase_nic_cirrus .AND. nic_cirrus /= 2) THEN
     CALL message('setphys', 'lemuphase_nic_cirrus=.TRUE. only possible with nic_cirrus=2', level=em_error)
  END IF
  CALL message('','---')
  CALL print_value('Emulator phase in/out (lemuphase_riming) = ', lemuphase_riming)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_riming) = ', eta_emu_riming)
  !>>UP #765
  CALL message('','---')
  CALL print_value('Emulator phase in/out (lemuphase_wbf) = ', lemuphase_wbf)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_wbf) = ', eta_emu_wbf)
  CALL message('','---')
  IF (lemuphase_wbf .AND. eta_emu_wbf > 1) THEN
     CALL message('setphys', 'eta_emu_wbf must be 0<eta_emu_wbf<=1', level=em_error)
  ENDIF
  !<<UP #765
  CALL print_value('Emulator phase in/out aggregation of ice crystals (lemuphase_icaggr) = ', lemuphase_icaggr)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_icaggr) = ', eta_emu_icaggr)
  CALL message('','---')
  CALL print_value('Emulator phase in/out accretion of ice crystals to snow (lemuphase_icaccr) = ', lemuphase_icaccr)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_icaccr) = ', eta_emu_icaccr)
  CALL message('','---')
  CALL print_value('Emulator phase in/out accretion of cloud droplets (lemuphase_cdaccr) = ', lemuphase_cdaccr)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_cdaccr) = ', eta_emu_cdaccr)
  CALL message('','---')
  CALL print_value('Emulator phase in/out autoconversion of cloud droplets (lemuphase_cdautc) = ', lemuphase_cdautc)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_cdautc) = ', eta_emu_cdautc)
  CALL message('','---')
  CALL print_value('Emulator phase in/out self-collection of ice crystals to ice (lemuphase_sci) = ', lemuphase_sci)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_sci) = ', eta_emu_sci)
  CALL message('','---')
  CALL print_value('Emulator phase in/out deposition (lemuphase_dep) = ', lemuphase_dep)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_dep) = ', eta_emu_dep)
  CALL message('','---')
  CALL print_value('Emulator phase in/out nucleation in cirrus scheme (lemuphase_icnucl) = ', lemuphase_icnucl)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_icnucl) = ', eta_emu_icnucl)
  CALL message('','---')
  CALL print_value('Emulator phase in/out subl/evp of falling i/s/r (lemuphase_sub_evp_falling) = ', lemuphase_sub_evp_falling)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_subfis) = ', eta_emu_subfis)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_evpr) = ', eta_emu_evpr)
  CALL message('','---')
  CALL print_value('Emulator phase in/out CD nucleation (lemuphase_cdnuc) = ', lemuphase_cdnuc)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_cdnuc) = ', eta_emu_cdnuc)
  CALL message('','---')
  CALL print_value('Emulator phase in/out secondary ice production (lemuphase_sip) = ', lemuphase_sip)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_sip) = ', eta_emu_sip)
  CALL message('','---')
  CALL print_value('Emulator phase in/out melting of sed. s/i (lemuphase_mlt) = ', lemuphase_mlt)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_mlt) = ', eta_emu_mlt)
  CALL message('','---')
  CALL print_value('New diagnostics on (lnewdiags) = ', lnewdiags)
  CALL message('','---')
  CALL print_value('New CCNclim diagnostics on (lccnclimdiags) = ', lccnclimdiags)
!<<UP
!>>UP #821
  CALL message('','---')
  CALL print_value('Prescribe tend. for simplified CMPs (lcmpsimpl_prescr) = ', lcmpsimpl_prescr)
  CALL message('','---')
  CALL message('','Prescribe tend. for specific CMP processes: 0 - no, 1 - constant, 2 - latlev profile')
  CALL message('','---')
  CALL print_value('Prescribe tend. for simplified CMPs (ncmpsimpl_prescr_rime) = ', ncmpsimpl_prescr_rime)
  CALL message('','---')
  CALL print_value('Prescribe tend. for simplified CMPs (ncmpsimpl_prescr_icnucl) = ', ncmpsimpl_prescr_icnucl)
  CALL message('','---')
  CALL print_value('Prescribe tend. for simplified CMPs (ncmpsimpl_prescr_subfis) = ', ncmpsimpl_prescr_subfis)
  CALL message('','---')
  CALL print_value('Prescribe tend. for simplified CMPs (ncmpsimpl_prescr_sci) = ', ncmpsimpl_prescr_sci)
  CALL message('','---')
  CALL print_value('Prescribe tend. for simplified CMPs (ncmpsimpl_prescr_icaccr) = ', ncmpsimpl_prescr_icaccr)
  CALL message('','---')
  CALL print_value('Prescribe tend. for simplified CMPs (ncmpsimpl_prescr_mlt) = ', ncmpsimpl_prescr_mlt)
!<<UP
!>>MA
  CALL message('','---')
  CALL print_value('Emulator phase in/out herterogeneous freezing (lemuphase_hetfrz) = ', lemuphase_hetfrz)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_hetfrz) = ', eta_emu_hetfrz)
  CALL message('','---')
  CALL print_value('Emulator phase in/out fall velocity of ice crystal mass (lemuphase_fvic) = ', lemuphase_fvic)
  CALL message('','---')
  CALL print_value('Emulator phase in/out (eta_emu_fvic) = ', eta_emu_fvic)
  CALL message('','---')
  CALL print_value('Switch WBF process on/off (lctrl_WBF) = ', lctrl_WBF)
!<<MA

    CALL message('','---')
    CALL print_value('Activation scheme (ncd_activ)                    = ', ncd_activ)
    SELECT CASE(ncd_activ)
      CASE(0)
       CALL message('','                            -> no activation',level=em_param)
      CASE(1)
       CALL message('','                            -> Lohmann et al. (1999) + Lin and Leaitch (1997)',level=em_param)
      CASE(2)
       CALL message('','                            -> Lohmann et al. (1999) + Abdul-Razzak and Ghan (2000)',level=em_param)
      CASE DEFAULT
        WRITE (message_text,*) '   ncd_activ = ',ncd_activ,' not supported.'
        CALL message('setphys',message_text, level=em_error)
    END SELECT

!>>UP FOR-ICE
  CALL print_value('SLF diagnostics on (lslf) = ', lslf)
!<<UP FOR-ICE

!>>ZK
    CALL message('','---')
    CALL print_value('Sub-grid updraft PDF (nactivpdf)                 = ', nactivpdf)
    SELECT CASE(nactivpdf)
      CASE(0)
       CALL message('','                            -> Mean updraft from TKE scheme without PDF',level=em_param)
      CASE(1)
       CALL message('','                            -> Coupling of updraft PDF with 20 bins to TKE scheme (West et al., 2013)', &
                    level=em_param)
      CASE(-1)
       CALL message('','                            -> Coupling of updraft PDF with 20 bins to ' &
                       // 'TKE scheme (West et al., 2013) with per-bin diagnostics',             &
                    level=em_param)
      CASE DEFAULT
       IF (nactivpdf > 1) THEN
        WRITE (message_text,*) '                            -> Coupling of updraft PDF with ', &
                               nactivpdf,' bins to TKE scheme (West et al., 2013)'
       ELSE
        WRITE (message_text,*) '                            -> Coupling of updraft PDF with ', &
                               -nactivpdf,' bins to TKE scheme (West et al., 2013) with per-bin diagnostics'
       END IF
       CALL message('', message_text, level=em_param)
    END SELECT

    IF (nactivpdf /= 0 .AND. ncd_activ /= 2) THEN
       CALL message('setphys', 'nactivpdf/=0 only possible with ncd_activ=2', level=em_error)
    END IF
!<<ZK

    CALL message('','---')
    CALL print_value('Cirrus scheme (nic_cirrus)                    = ', nic_cirrus)
    SELECT CASE(nic_cirrus)
      CASE(0)
       CALL message('','                            -> no cirrus scheme',level=em_param)
      CASE(1)
       CALL message('','                            -> Lohmann, JAS 2002',level=em_param)
      CASE(2)
       CALL message('','                            -> Kaercher & Lohmann, JGR 2002',level=em_param)
      CASE DEFAULT
        WRITE (message_text,*) '   nic_cirrus = ',nic_cirrus,' not supported.'
        CALL message('setphys',message_text, level=em_error)
    END SELECT

    CALL message('','---')
    CALL print_value('Autoconversion (nauto)                    = ', nauto)
    SELECT CASE(nauto)
      CASE (0)
        CALL message('setphys','nauto=0. Are you sure?',level=em_warn)
      CASE (1)
        CALL message('','                            -> Beheng (1994) - ECHAM5 Standard',level=em_param)
      CASE (2)
        CALL message('','                            -> Khairoutdinov and Kogan (2000)',level=em_param)
      CASE DEFAULT
        WRITE (message_text,*) '   nauto = ',nauto,' not supported.'
        CALL message('setphys',message_text, level=em_error)
    END SELECT

!>>gf
    CALL message('','---')
    CALL print_value('Orographic cirrus cloud (lorocirrus)                    = ', lorocirrus)
!<<gf

!>>SF #251
    CALL message('','---')
    CALL print_value('Secondary ice production (lsecprod)                    = ', lsecprod)
!<<SF

!>>SF #475
    CALL message('','---')
    CALL print_value('Dynamical minimum CDNC (ldyn_cdnc_min)                 = ', ldyn_cdnc_min)
!<<SF

!>>SF #489
    CALL message('','---')
    SELECT CASE(cdnc_min_fixed)
       CASE(10, 40)
         CALL print_value('Fixed minimum CDNC (cdnc_min_fixed) = ', cdnc_min_fixed)
       CASE DEFAULT
         WRITE (message_text,*) '   cdnc_min_fixed = ',cdnc_min_fixed,' not supported.'
         CALL message('setphys',message_text, level=em_error)
    END SELECT
!<<SF
!>>UP: make rcd_vol_max tunable via settings file
    CALL message('','---')
    CALL print_value('maximum droplet volume radius [m] (rcd_vol_max) = ', rcd_vol_max)
!<<UP

!>>SF
    CALL message('','---')
    CALL print_value('Detailed 2-moment scheme (nclmi_progn)                  = ', nclmi_progn)
!<<SF

    CALL message('','---')
!
!*        3.2      SET UP CONSTANTS IN *mo_physc2*.
!
  CALL iniphy
!
!*        3.3      Define albedo parameters depending on resolution
!
  CALL init_albedo_ice
!
END SUBROUTINE setphys

!-----------------------------------------------------------------------

