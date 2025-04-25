!*    NAMELIST *DYNCTL* TO CONTROL OPTIONS HAVING A BEARING ON THE DY-  
!                       NAMICS AND THE GENERAL CONTROL OF THE FORECAST. 
!                                                                       
!     R.G AND M.J        E.C.M.W.F.     23/11/82.                       
!                                                                       
!      NAME     TYPE      PURPOSE                                       
!      ----     ----      -------                                       
!                                                                       
!     *DIAGDYN*           interval for diagnostics of dynamical part
!     *DIAGVERT*          interval for diagnostics for all levels
!     *NTRN*    INTEGER   *CRITICAL WAVENUMBERS FOR STRONG              
!                         STRATOSPHERIC DAMPING.                        
!     *NLVSTD1* INTEGER   *LAST (UPPERMOST) LAYER AT WHICH              
!                         STRATOSPHERIC HORIZONTAL DIFFUSION IS         
!                         ENHANCED.                                     
!     *NLVSTD2* INTEGER   *FIRST (LOWEST) LAYER AT WHICH                
!                         STRATOSPHERIC HORIZONTAL DIFFUSION IS         
!                         ENHANCED.                                     
!     *VCRIT*   REAL      *CRITICAL VELOCITY ABOVE WHICH                
!                         HORIZONTAL DIFFUSION IS ENHANCED FOR          
!                         T63 WITH DT=20MIN.                            
!     *VCHECK*  REAL      THRESHOLD VALUE FOR CHECK OF HIGH WINDSPEED   
!     *HDAMP*   REAL      *DAMPING FACTOR FOR STRONG                    
!                         STRATOSPHERIC DAMPING.                        
!     *ENSTDIF* REAL      *FACTOR BY WHICH STRATOSPHERIC                
!                         HORIZONTAL DIFFUSION IS INCREASED FROM ONE    
!                         LEVEL TO NEXT LEVEL ABOVE.                    
!     *LUMAX*   LOGICAL   *TRUE TO COMPUTE AND PRINT INFORMATION ON     
!                         MAXIMUM WIND.                                 
!     *EPS*     REAL      *COEFFICIENT FOR THE TIME FILTERING.          
!     *DAMPTH*  REAL      *DAMPING TIME (N=NK) IN HOURS FOR THE         
!                         HORIZONTAL                                    
!                         DIFFUSION OF VORTICITY (LINEAR SQUARE         
!                         LAPLACIAN), DIVERGENCE, TEMPERATURE.          
!     *APSURF*  REAL      FIXED GLOBAL MEAN OF SURFACE PRESSURE         
!         
!     *damhih*  REAL                                                    
!     *spdrag*  REAL      upper sponge layer coefficient (sec)-1
!     *enspodi* REAL      factor by which   upper sponge layer 
!                         coefficient is increased from one        
!                         level to next level above. 
!     *nlvspd1* INTEGER   last (uppermost) layer of upper sponge  
!     *nlvspd2* INTEGER   first (lowest) layer of upper sponge  

NAMELIST /DYNCTL/                                                      &
  NTRN, NLVSTD1, NLVSTD2, LUMAX, LZONDIA,                              &
  LDIAHDF, VCRIT, HDAMP, ENSTDIF, APSURF, VCHECK, EPS, DAMPTH,         &
  damhih, spdrag, enspodi, nlvspd1, nlvspd2, diagdyn, diagvert
!                                                                       
!     ------------------------------------------------------------------
