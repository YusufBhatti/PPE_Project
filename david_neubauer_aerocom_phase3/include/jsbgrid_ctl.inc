!*    NAMELIST *JSBGRID_CTL* To control the grid in jsbach stand alone runs
!                                                                 
!      
!      NAME            TYPE        PURPOSE                                       
!      ----            ----        -------                                       
!
!     *LON_EAST*      REAL        Eastern boundary of region
!     *LON_WEST*      REAL        Western boundary of region
!     *LAT_NORTH*     REAL        Northern boundary of region
!     *LAT_SOUTH*     REAL        Southern boundary of region
!     *SET_LON_ZERO*  LOGICAL     Set longitude to zero
!     *TIME_OFFSET*    REAL        shift between local time according to 
!                                 longitude and time implied by forcing

NAMELIST /JSBGRID_CTL/ &
      LON_EAST, LON_WEST, LAT_NORTH, LAT_SOUTH, SET_LON_ZERO, TIME_OFFSET
!                                                                       
!     ------------------------------------------------------------------
