!*    NAMELIST *DYNVEG_CTL* To control JSBACH dynamic vegetation.
!                                                                 
!      
!      NAME                TYPE       PURPOSE                                       
!      ----                ----       -------                                       
!                                                                       
!     *READ_FPC*           LOGICAL    read initial cover fraction from an 
!                                     external file
!     *FPC_FILE_NAME*      CHARACTER  name of the external vegetation file
!     *DYNVEG_ALL*         LOGICAL    activate competition between woody types
!                                     and grass
!     *DYNVEG_FEEDBACK*    LOGICAL    activate feedback of dynveg on jsbach
!     *ACCELERATE_DYNVEG*  REAL       factor by which vegetation dynamics are
!                                     accelerated, should be 1 for production
!                                     runs

NAMELIST /DYNVEG_CTL/ &
      READ_FPC, FPC_FILE_NAME, DYNVEG_ALL, DYNVEG_FEEDBACK, ACCELERATE_DYNVEG
!                                                                       
!     ------------------------------------------------------------------
