!*    NAMELIST *CLIMBUF_CTL* To control JSBACH climate buffer.
!                                                                 
!      
!      NAME                TYPE       PURPOSE                                       
!      ----                ----       -------                                       
!                                                                       
!     *INIT_RUNNING_MEANS* LOGICAL    initialize running means of climate in
!                                     the climate buffer
!     *READ_CLIMBUF*       LOGICAL    read initial climate variables from an
!                                     external file
!     *CLIMBUF_FILE_NAME*  CHARACTER  name of the external climate file

NAMELIST /CLIMBUF_CTL/ &
      INIT_RUNNING_MEANS, READ_CLIMBUF, CLIMBUF_FILE_NAME
!                                                                       
!     ------------------------------------------------------------------
