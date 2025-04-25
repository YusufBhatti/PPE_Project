!*    NAMELIST *PHYSCTL* TO CONTROL OPTIONS HAVING A BEARING ON THE     
!                        PARAMETERISATIONS OF DIABATIC PROCESSES.       
!                                                                       
!     R.G AND M.J        E.C.M.W.F.     23/11/82.                       
!                                                                       
!      NAME     TYPE      PURPOSE                                       
!      ----     ----      -------                                       
!                                                                       
!     *LPHYS*   LOGICAL   *TRUE FOR PARAMETERISATION OF DIABATIC        
!                          PROCESSES.                            
!     *LRAD*    LOGICAL   *TRUE FOR RADIATION SCHEME        
!     *LVDIFF*  LOGICAL   *TRUE FOR VERTICAL DIFFUSION.                 
!     *LCONV*   LOGICAL   *TRUE* TO ALLOW CONVECTION                   
!     *LMFPEN*  LOGICAL   *TRUE IF PENETRATIVE CONVECTION IS SWITCHED ON 
!     *LCOND*   LOGICAL   *TRUE FOR LARGE SCALE CONDENSATION SCHEME.    
!     *LGWDRAG* LOGICAL   *TRUE FOR GRAVITY WAVE DRAG SCHEME  
!     *LSURF*   LOGICAL   *TRUE FOR SURFACE EXCHANGES.                  
!     *LICE*    LOGICAL   *TRUE* FOR SEA-ICE TEMPERATURE CALCULATION    
!                                                                       
!  !!mgs: new switches for cloud activation
!     lclmi_progn logical  true for prognostic cloud activation scheme (former ncdnc and nicnc)
!     nclmi_progn integer allows to fine-tune the 2-moment cloud microphysics calculations in
!                         mo_cloud_micro_2m.f90
!                         WARNING: this is meant for special sensitivity studies!!
!                                  Do not modify if you don't know about it
!                         0: both CDNC and ICNC are diagnostic variables
!                            WARNING!! this is *not* equivalent to lclmi_progn = .false. because
!                            this still uses cloud_cdnc_icnc.f90 (with a different autoconversion
!                            scheme as compared to using cloud.f90).
!                         1: CDNC is a prognostic var, ICNC is a diagnostic var
!                         2: CDNC is a diagnotic var,  ICNC is a prognostic var
!                         3: CDNC and ICNC are prognostic vars [DEFAULT, normal case]         
!     ncd_activ   integer  select cloud droplet activation scheme:
!                         0: off (lclmi_progn=false)
!                         1: Lohmann et al. (1999) + Lin and Leaitch (1997)
!                         2: Lohmann et al. (1999) + Abdul-Razzak and Ghan (2000)
!>>ZK
!     nactivpdf   integer  Sub-grid scale pdf of updraft velocities in activation scheme:
!                         0: Mean updraft from TKE scheme without pdf
!                         1: Coupling of updraft pdf with 20 bins to TKE scheme (West et al., 2013)
!                         >1: as 1, but use the specified number of bins
!                         <0: as positive value, but output per-bin supersaturation diagnostics etc.
!<<ZK
!>>SF changed the former lice_supersat logical into an integer (nic_cirrus)
!     nic_cirrus  integer select ice crystal cirrus scheme:
!                         0: off (lclmi_progn=false)
!                         1: Lohmann, JAS 2002
!                         2: Kaercher & Lohmann, JGR 2002
!<<SF
!     nauto       integer   select autoconversion scheme for clouds:
!                         1  Beheng (1994) - ECHAM5 Standard
!                         2  Khairoutdinov and Kogan (2000)
!
!     lsecprod    logical turn on secondary ice production !SF #251
!     lorocirrus  logical turn on orographic cirrus clouds
!     ldyn_cdnc_min logical turn on dynamical setting of the min cloud droplet number concentration
!     cdnc_min_fixed  fixed value for min CDNC in cm-3 (used when ldyn_cdnc_min is FALSE)
!                     Warning! So far only values of 40 or 10 are accepted.

NAMELIST /PHYSCTL/ LPHYS,  LRAD,  LVDIFF, LCOND,   LSURF,               &
                   LCONV, LMFPEN, LGWDRAG, LICE, LCONVMASSFIX,          &
                   lclmi_progn, nclmi_progn, ncd_activ, nic_cirrus, nauto, &
                   lsecprod, lorocirrus, ldyn_cdnc_min, nactivpdf, cdnc_min_fixed, &!davidn
		   tun47zinhomi,tun47zinhoml1,tun47zinhoml2,tun47zinhoml3,&
		   tun47zinpar,tun47ccraut,tun47ccsaut,tun47cdncmin,&
		   tun47entrpen,tun47cprcon,tun47entrscv,tun47cmfctop,&!davidn
	 	   lemuphase_nic_cirrus, eta_emu_nic_cirrus,& !UP
		   lemuphase_riming, eta_emu_riming,& !UP
		   lemuphase_wbf, eta_emu_wbf,& !UP #765
		   lemuphase_icaggr, eta_emu_icaggr,& !UP #757
		   lemuphase_icaccr, eta_emu_icaccr,& !UP #775
		   lemuphase_cdaccr, eta_emu_cdaccr,& !UP #805
		   lemuphase_cdautc, eta_emu_cdautc,& !UP #805
		   lemuphase_sci, eta_emu_sci,& !UP
		   lemuphase_dep, eta_emu_dep,& !UP #801
		   lemuphase_icnucl, eta_emu_icnucl,& !UP #801
		   lemuphase_sub_evp_falling, eta_emu_subfis,& !UP #802
		   eta_emu_evpr, & !UP #802
		   lemuphase_cdnuc, eta_emu_cdnuc,& !UP #805
		   lemuphase_sip, eta_emu_sip,& !UP #805
		   lemuphase_mlt, eta_emu_mlt,& !UP #802
		   lcmpsimpl_prescr, & !UP #821
	           ncmpsimpl_prescr_rime, & !UP #821
	           ncmpsimpl_prescr_icnucl, & !UP #821
	           ncmpsimpl_prescr_subfis, & !UP #821
	           ncmpsimpl_prescr_sci, & !UP #821
	           ncmpsimpl_prescr_icaccr, & !UP #821
	           ncmpsimpl_prescr_mlt, & !UP #821
		   ldetr_liquid, & !UP #782.1
                   ldetr_convicfactor, convicfactor,& !UP #782.2
                   lnewdiags,lepsfix, &
		   lccnclimdiags, &
                   lpr_corr, & !UP #783.3
                   rcd_vol_max, & !UP
		   lsedfix, & !UP #797
		   lslf, & !UP FOR-ICE #839
		   lemuphase_hetfrz, eta_emu_hetfrz, & !MA #790
		   lemuphase_fvic, eta_emu_fvic, & !MA #795
		   lctrl_WBF !MA WBF switch
