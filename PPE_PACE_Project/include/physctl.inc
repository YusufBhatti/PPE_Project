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
!     lcdnc_progn logical  true for prognostic cloud activation scheme (former ncdnc and nicnc)
!     ncd_activ   integer  select cloud droplet activation scheme:
!                         0: off (lcdnc_progn=false)
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
!                         0: off (lcdnc_progn=false)
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
                   lcdnc_progn, ncd_activ, nic_cirrus, nauto, lsecprod, &
                   lorocirrus, ldyn_cdnc_min, nactivpdf, cdnc_min_fixed,&
		   ac_scale_activation, ac_scale_autoconversion, lnoh2oaod,&
		   ac_scale_accretion, ac_scale_KK_LWP_exponent
