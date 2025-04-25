!  ---------------
!   HYDROLOGY_CTL controls the hydrlology module
!  ---------------
!      
!    NAME           TYPE       PURPOSE                                       
!    ----           ----       -------                                       
!                                                                       
!    ldebughd      LOGICAL     debug HD model  
!    diag_water_budget LOGICAL additional prints to diagnose the global water budget 
!    lbase         LOGICAL     switch for baseflow calculations
!    locean        LOGICAL     closure of water budget for ocean coupling
!    nhd_diag      INTEGER     region number for outflow diagnostic (formerly isolog)
!                                 0   none
!                                 1   Bothnian Bay/Sea
!                                 2   Torneaelven
!                                 4   St.Lawrence
!                                 5   Paraguay
!                                 6   Oder
!                                 7   Elbe
!                                 8   Oranje
!                                 9   Amudarya
!                                10   Lena
!                                99   user defined (FBLOG1, FLLOG1, FBLOG2, FLLOG2)
!    lhd_highres   LOGICAL     swich for outflow diagnostic on HD model grid
!    puthd         TYPE        coupling interval - put data to   HD-model
!    gethd         TYPE        coupling interval - get data from HD-model
!    fllog1        REAL(dp)    user defined grid cells for diagnostics (with nhd_diag=99)
!    fblog1        REAL(dp)      fllog1, fblog1: longitude, latitude of grid cell 1
!    fllog2        REAL(dp)      fllog2, fblog2: longitude, latitude of grid cell 2
!    fblog2        REAL(dp)
!    lhd_rout      LOGICAL     switch for routing scheme
!                                 true: via direction arrays;  false: via index arrays

NAMELIST /HYDROLOGY_CTL/ &
     LDEBUGHD, DIAG_WATER_BUDGET, LBASE, LOCEAN, NHD_DIAG, LHD_HIGHRES, &
     GETHD, PUTHD, FLLOG1, FBLOG1, FLLOG2, FBLOG2, LHD_ROUT

!     ------------------------------------------------------------------
