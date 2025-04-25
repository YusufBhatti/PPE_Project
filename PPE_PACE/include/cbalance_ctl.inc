!*    NAMELIST *CBALANCE_CTL* To control the carbon pools.
!                                                                 
!      
!   NAME            TYPE      PURPOSE                                       
!   ----            ----      -------                                       
!
!  *READ_CPOOLS*     LOGICAL    Initialize carbon pools with data from an 
!                               external file 
!  *CPOOL_FILE_NAME* CHARACTER  Name of the file containing initial data for
!                               carbon pools (only used if read_cpools=true)
!  *READ_NPOOLS*     LOGICAL    Initialize nitrogen pools with data from an 
!                               external file 
!  *NPOOL_FILE_NAME* CHARACTER  Name of the file containing initial data for
!                               nitrogen pools (only used if read_npools=true)
!  *READ_NDEPO*      LOGICAL    Read in nitrogen deposition fields from file.
!  *NDEPO_FILE_NAME* CHARACTER  Name of the file containing the nitrogen deposition data.

NAMELIST /CBALANCE_CTL/ &
      READ_CPOOLS, CPOOL_FILE_NAME, &
      READ_NPOOLS, NPOOL_FILE_NAME, &
      READ_NDEPO,NDEPO_FILE_NAME 
!                                                                       
!     ------------------------------------------------------------------
