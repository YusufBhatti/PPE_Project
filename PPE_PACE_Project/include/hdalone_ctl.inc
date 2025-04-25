!  -----------
!  HDALONE_CTL controls stand alone hydrlology moldel runs
!  -----------
!
!  In echam runs, some of these parameters are set in RUNCTL.
!                                                                 
!      
!    NAME           TYPE       PURPOSE                                       
!    ----           ----       -------                                       
!                                                                       
!    iswrit        INTEGER     restart time step (0: no restart file)
!    out_expname   CHARACTER   experiment name
!    out_datapath  CHARACTER   path to where the output data shall be written
!    year1         INTEGER     initial year of the run
!    nstep         INTEGER     number of time steps within the run
!    delta_time    REAL(dp)    model time step lenght in seconds
!    ufakru        REAL(dp)    unit factor for runoff and drainage input data
!    runoff_file   CHARACTER   file with input runoff data
!    drainage_file CHARACTER   file with input drainage data
!    forcing_freq  INTEGER     frequency of the forcing data (0: stepwise, 1: daily)
!    iout          INTEGER     averaging period of the HD output
!                                 1   30-day averages
!                                 2   decadal averages
!                                 3   weekly averages
!                                 4   monthly averages without leap years
!                                 5   monthly averages with leap years
!                                 6   daily averages
!                                 7   stepwise output
!    do_remapping  INTEGER     type of interpolation from input (atmospheric) grid to HD grid
!                                 0   no interpolation: input data grid equals HD grid
!                                 1   remapping in routine hd_remap (default)

NAMELIST /HDALONE_CTL/ &
     ISWRIT, OUT_EXPNAME, OUT_DATAPATH, YEAR1, NSTEP, DELTA_TIME, &
     UFAKRU, RUNOFF_FILE, DRAINAGE_FILE, FORCING_FREQ, IOUT, DO_REMAPPING

!     ------------------------------------------------------------------
