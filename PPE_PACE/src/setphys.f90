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
                               lcdnc_progn, ncd_activ, nic_cirrus, lsecprod, &
                               nauto, lorocirrus, ldyn_cdnc_min, & !SF #475
                               cdnc_min_fixed, & !SF #589
                               nactivpdf, & !ZK
                               t_ccsaut, t_ccraut, t_csecfrl, t_entrmid, t_entrpen, &
                               t_cprcon, t_zinhomi, &
                               ac_scale_activation, ac_scale_autoconversion, ac_scale_accretion,& !DN
                               ac_scale_KK_LWP_exponent !DN
  USE mo_p3_fields,      ONLY: iprog, ccclpwr, l2moment, ldisableupdraftcond, lconstnc, lsnow_sed, lfalling_ice, &
                               lrain_sed, lact, lctrl_frz_below_238K, lctrl_het_mxphase_frz, levap_rain, &
                               lmelting, lriming, lself_collection, ladjustment, inumberadjustment, &
                               lcover_adjustment, lmass_transport, lsubsat_cnd_dep, lcirrus, lwbf, lconvice, &
                               isublimation, nconstnc, csubw, xismall, xlsmall, ccnislf, ccqccol, iintscheme, &
                               ldiabheat, ihomfrz, ccsupci, lallow_si, rcmax, adjsupsat, &
                               ccftau, ccrsnow, nmicro, nmicro_max, nsedi, microprct, lprescribeaerosols, &
                               lprescribeaerosols_every, lprescribe_numbers, lprescribe_numbers_every, &
                               lfull_diag, lspeedrun
!>>DN ambient AOD per species
  USE mo_param_switches, ONLY: lnoh2oaod
!<<DN
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
  lcdnc_progn   = .FALSE.
  ncd_activ     = 2        ! default scheme for cdnc activation (will be 0 if lcdnc_progn=false)
  nactivpdf     = 0        ! default scheme for updraft PDF !ZK
  nic_cirrus    = 2        ! default scheme for cirrus scheme   (will be 0 if lcdnc_progn=false)
  nauto         = 1        ! default scheme for autoconversion  (will be 0 if lcdnc_progn=false)
  lsecprod      = .FALSE.  ! switch for computing secondary ice production !SF #251
  lorocirrus    = .FALSE.  ! switch for gravity waves updraft velocity for icnc (orographic cirrus clouds)
  ldyn_cdnc_min = .FALSE.  ! switch to turn on the dynamical setting of the min CDNC !SF #475

  !>> RD
  t_ccsaut     = 1200._dp
  t_ccraut     = 1.3_dp
  t_csecfrl    = 5e-7_dp
  t_entrmid    = 2.0E-4_dp ! Average entrainment rate for midlevel convection
  t_entrpen    = 1.E-4_dp
  t_cprcon     = 9.E-04_dp
  t_zinhomi    = 0.80_dp
  !>> RD
  !>>DN
  ac_scale_activation     = 1.0_dp
  ac_scale_autoconversion = -1.79_dp
  ac_scale_KK_LWP_exponent = 2.47_dp
  ac_scale_accretion = 1.0_dp
  !<<DN
!>>DN ambient AOD per species
  lnoh2oaod=.FALSE.
!<<DN

  cdnc_min_fixed = 40      ! fixed value for min CDNC, in cm-3 (used when ldyn_cdnc_min is FALSE)
                           ! Warning! So far only values of 40 or 10 are accepted.
                           ! SF #489
  !
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
     CALL p_bcast (lcdnc_progn, p_io)
     CALL p_bcast (ncd_activ, p_io)
     CALL p_bcast (nactivpdf, p_io)
     CALL p_bcast (nic_cirrus, p_io)
     CALL p_bcast (lsecprod, p_io)
     CALL p_bcast (lorocirrus, p_io)
     CALL p_bcast (ldyn_cdnc_min, p_io)
     CALL p_bcast (cdnc_min_fixed, p_io)
     CALL p_bcast (nauto, p_io)

!>>RD
     CALL p_bcast (t_ccsaut, p_io)
     CALL p_bcast (t_ccraut, p_io)
     CALL p_bcast (t_csecfrl, p_io)
     CALL p_bcast (t_entrmid, p_io)
     CALL p_bcast (t_entrpen, p_io)
     CALL p_bcast (t_cprcon, p_io)
     CALL p_bcast (t_zinhomi, p_io)
     CALL p_bcast (iprog, p_io)
     CALL p_bcast (ccclpwr, p_io)
     CALL p_bcast (ccsupci, p_io)
     CALL p_bcast (lallow_si, p_io)
     CALL p_bcast (l2moment, p_io)
     CALL p_bcast (lfull_diag, p_io)
     CALL p_bcast (lspeedrun, p_io)
     CALL p_bcast (ldisableupdraftcond, p_io)
     CALL p_bcast (lconstnc, p_io)
     CALL p_bcast (lsnow_sed, p_io)
     CALL p_bcast (lfalling_ice, p_io)
     CALL p_bcast (lrain_sed, p_io)
     CALL p_bcast (lact, p_io)
     CALL p_bcast (lctrl_frz_below_238K, p_io)
     CALL p_bcast (lctrl_het_mxphase_frz, p_io)
     CALL p_bcast (levap_rain, p_io)
     CALL p_bcast (lmelting, p_io)
     CALL p_bcast (lriming, p_io)
     CALL p_bcast (lself_collection, p_io)
     CALL p_bcast (ladjustment, p_io)
     CALL p_bcast (inumberadjustment, p_io)
     CALL p_bcast (lcover_adjustment, p_io)
     CALL p_bcast (lmass_transport, p_io)
     CALL p_bcast (lsubsat_cnd_dep, p_io)
     CALL p_bcast (lcirrus, p_io)
     CALL p_bcast (lwbf, p_io)
     CALL p_bcast (lconvice, p_io)
     CALL p_bcast (isublimation, p_io)
     CALL p_bcast (ihomfrz, p_io)
     CALL p_bcast (iintscheme, p_io)
     CALL p_bcast (ldiabheat, p_io)
     CALL p_bcast (nconstnc, p_io)
     CALL p_bcast (rcmax, p_io)
     CALL p_bcast (adjsupsat, p_io)
     CALL p_bcast (csubw, p_io)
     CALL p_bcast (xismall, p_io)
     CALL p_bcast (xlsmall, p_io)
     CALL p_bcast (ccnislf, p_io)
     CALL p_bcast (ccqccol, p_io)
     CALL p_bcast (ccftau, p_io)
     CALL p_bcast (ccrsnow, p_io)
     CALL p_bcast (nmicro, p_io)
     CALL p_bcast (nsedi, p_io)
     CALL p_bcast (nmicro_max, p_io)
     CALL p_bcast (microprct, p_io)
     CALL p_bcast (lprescribeaerosols, p_io)
     CALL p_bcast (lprescribeaerosols_every, p_io)
     CALL p_bcast (lprescribe_numbers, p_io)
     CALL p_bcast (lprescribe_numbers_every, p_io)
!>>RD
     !>>DN
     CALL p_bcast (ac_scale_activation, p_io)
     CALL p_bcast (ac_scale_autoconversion, p_io)
     CALL p_bcast (ac_scale_KK_LWP_exponent, p_io)
     CALL p_bcast (ac_scale_accretion, p_io)
     !<<DN
     !>>DN ambient AOD per species
     CALL p_bcast (lnoh2oaod, p_io)   
     !<<DN
  ENDIF

! ADJUST FOR CONSISTENCY
  IF(iprog==2) THEN
     l2moment = .TRUE.
  endif !iprog == 2

! PRINT OUT THE P3 NAMELIST
     CALL message ('','---------------- GLOBAL PARAMS ----------------')
     CALL print_value ('t_entrmid             ', t_entrmid)
     CALL print_value ('t_entrpen             ', t_entrpen)
     CALL print_value ('t_cprcon              ', t_cprcon)
     CALL print_value ('t_zinhomi             ', t_zinhomi)
     CALL message ('','----------------- P3 NAMELIST -----------------')
     CALL print_value ('t_ccsaut              ', t_ccsaut)
     CALL print_value ('t_ccraut              ', t_ccraut)
     CALL print_value ('t_csecfrl             ', t_csecfrl)
     CALL print_value ('iprog                 ', iprog)
     CALL print_value ('ccclpwr               ', ccclpwr)
     CALL print_value ('ccsupci               ', ccsupci)
     CALL print_value ('lallow_si             ', lallow_si)
     CALL print_value ('l2moment              ', l2moment)
     CALL print_value ('lfull_diag            ', lfull_diag)
     CALL print_value ('lspeedrun             ', lspeedrun)
     CALL print_value ('ldisableupdraftcond   ', ldisableupdraftcond)
     CALL print_value ('lconstnc              ', lconstnc)
     CALL print_value ('lsnow_sed             ', lsnow_sed)
     CALL print_value ('lfalling_ice          ', lfalling_ice)
     CALL print_value ('lrain_sed             ', lrain_sed)
     CALL print_value ('lact                  ', lact)
     CALL print_value ('lctrl_frz_below_238K  ', lctrl_frz_below_238K)
     CALL print_value ('lctrl_het_mxphase_frz ', lctrl_het_mxphase_frz)
     CALL print_value ('levap_rain            ', levap_rain)
     CALL print_value ('lmelting              ', lmelting)
     CALL print_value ('lriming               ', lriming)
     CALL print_value ('lself_collection      ', lself_collection)
     CALL print_value ('ladjustment           ', ladjustment)
     CALL print_value ('inumberadjustment     ', inumberadjustment)
     CALL print_value ('lcover_adjustment     ', lcover_adjustment)
     CALL print_value ('lmass_transport       ', lmass_transport)
     CALL print_value ('lsubsat_cnd_dep       ', lsubsat_cnd_dep)
     CALL print_value ('lcirrus               ', lcirrus)
     CALL print_value ('lwbf                  ', lwbf)
     CALL print_value ('lconvice              ', lconvice)
     CALL print_value ('isublimation          ', isublimation)
     CALL print_value ('ihomfrz               ', ihomfrz)
     CALL print_value ('iintscheme            ', iintscheme)
     CALL print_value ('ldiabheat             ', ldiabheat)
     CALL print_value ('nconstnc              ', nconstnc)
     CALL print_value ('rcmax                 ', rcmax)
     CALL print_value ('adjsupsat             ', adjsupsat)
     CALL print_value ('csubw                 ', csubw)
     CALL print_value ('xismall               ', xismall)
     CALL print_value ('xlsmall               ', xlsmall)
     CALL print_value ('ccnislf               ', ccnislf)
     CALL print_value ('ccqccol               ', ccqccol)
     CALL print_value ('ccftau                ', ccftau)
     CALL print_value ('ccrsnow               ', ccrsnow)
     CALL print_value ('nmicro                ', nmicro)
     CALL print_value ('nmicro_max            ', nmicro_max)
     CALL print_value ('nsedi                 ', nsedi)
     CALL print_value ('microprct             ', microprct)
     CALL print_value ('lprescribeaerosols    ', lprescribeaerosols)
     CALL print_value ('lprescribeaerosols_every', lprescribeaerosols_every)
     CALL print_value ('lprescribe_numbers    ', lprescribe_numbers)
     CALL print_value ('lprescribe_numbers_every', lprescribe_numbers_every)

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
                                    ' concentration (lcdnc_progn)', lcdnc_progn
  CALL message('',message_text, level=em_param)

  IF (.NOT. lcdnc_progn) THEN
    ncd_activ  = 0
    nic_cirrus = 0
    nauto      = 0       !SF added for security only
    lsecprod    = .FALSE. !SF added for security only (#251)
    lorocirrus = .FALSE. !GF added for security only
    ldyn_cdnc_min = .FALSE. !SF security
    CALL message('setphys', 'ncd_activ, nic_cirrus and nauto set to 0 because lcdnc_progn=.FALSE.',  &
                 level=em_warn)
    CALL message('setphys', 'lsecprod set to false because lcdnc_progn=.FALSE.', level=em_warn)
    IF (lccnclim) THEN
       CALL message('setphys', 'lcdnc_progn must be .TRUE. when lccnclim is active!', level=em_error)
    ENDIF
  ELSE !SF lcdnc_progn=true
    IF (ncd_activ <= 0) THEN
      CALL message('setphys', 'ncd_activ must be >0 for lcdnc_progn=.TRUE.!', level=em_error)
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
      CASE(3)
       CALL message('','                            -> Steffen',level=em_param)
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

