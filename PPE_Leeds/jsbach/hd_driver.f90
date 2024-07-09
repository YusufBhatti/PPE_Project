!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
PROGRAM hd_driver

  !------------------------------------------------------------------------------------
  !
  !   ******* Globale/Regionale Abfluss-Simulation: Hauptprogramm, dessen
  !           Hauptaufgabe darin besteht, die HDMODEL-Subroutine 
  !           HDMODEL.f aufzurufen und die Inputfelder pro Zeitschritt zu 
  !           uebergeben.
  !
  !
  !   ******** Version 1.0 - Oktober 1999 
  !            Programmierung und Entwicklung: Stefan Hagemann 
  !            Programmcode basiert auf der Offline-Version der HDModels
  !            die auch regional anwendbar ist (regsim.f)
  !            Service-Outputfile ist nun simple binaer statt direct access
  !
  !   ******** Version 1.1 - March 2000
  !            Steuerungsroutine: Inputdaten ueber Commonblock pinp.for
  !            gesteuert 
  !
  !   ******** Version 1.1 -- January 2001
  !            Ursprung Ocean-Longitude korrigiert.  0. statt -1.40625
  !
  !            Anmerkung: Input-Daten von Runoff und Drainage sollten
  !                       Einheit m/s haben.
  !
  !   ******** Version 1.2 -- October 2005
  !            Implementation of Unit Factor UFAKRU, which is applied to Input arrays
  !            of Surface ruonoff and drainage if UFAK.ne.1. Init factor is necessary
  !            if runoff/drainage input unit is not m/s.
  !
  !   ******** Version 1.2.1 -- August 2008
  !            Initialisierung von UFAKRU mit 0
  !            Einbau Zusaetzlicher Logpunkte fuer ISOLOG
  !            ISOLOG = 100 -> Koordinaten von Fluss 1 und 2 werden aus Input 
  !            file hdini.inp gelesen.
  !
  !   ******** May 2011, Veronika Gayler
  !            fortran 90 version to run with ECHAM6 mo_hydrology
  !
  !   ******** Feb 2015 Stefan Hagemann
  !            Correction of bug in HD offline version
  !            Remapping alternative implemented (no remapping)
  !
  !------------------------------------------------------------------------------------

  USE mo_kind,             ONLY: dp 
  USE mo_jsbach_grid,      ONLY: grid_type, domain_type
  USE mo_hydrology,        ONLY: init_hydrology, hydrology_model, &
                                 hydrology_restart, hydrology_offline, &
                                 cleanup_hydrology, &
                                 locean, lbase, hd_steps_per_day, &
                                 nl, nb, friv, water_to_ocean, lhd_highres, &
                                 do_remapping, runoff_s, runoff_dr
  USE mo_hd_highres_io,    ONLY: hd_highres_open
  USE mo_mpi,              ONLY: p_start, p_stop, p_init_communicators
  USE mo_machine,          ONLY: machine_setup
  USE mo_time_control,     ONLY: l_trigfiles, dt_start, init_manager, delta_time, &
                                 init_times, ec_manager_init, time_reset, &
                                 no_steps, time_set, time_reset, current_date, &
                                 write_date, lresume
  USE mo_exception,        ONLY: finish, message, message_text
  USE mo_jsbach_interface, ONLY: get_dates
  USE mo_control,          ONLY: nlon, ngl
  USE mo_filename,         ONLY: out_expname, out_datapath, find_next_free_unit

  IMPLICIT NONE

  ! local variables
  INTEGER  :: IO_timestep
  INTEGER  :: luoce          ! file unit
  INTEGER  :: idum           ! dummy integer
  INTEGER  :: step           ! hd model time step counter
  INTEGER  :: istep          ! initial time step

  ! local variables defined in namelist HDALONE_CTL
  INTEGER  :: nstep
  REAL(dp) :: ufakru
  INTEGER  :: iswrit
  INTEGER  :: year1
  INTEGER  :: forcing_freq
  INTEGER  :: iout

  ! local arrays that have to do with the forcing, defined on the atmosphere grid
  REAL(dp), ALLOCATABLE :: slm(:,:)
  REAL(dp), ALLOCATABLE :: runoff(:,:), drain(:,:)
  REAL(dp), ALLOCATABLE :: disch(:,:)
  REAL(dp), ALLOCATABLE :: awfre(:,:), apmecal(:,:)

  REAL,     ALLOCATABLE :: dummy_sp(:,:)

  ! grid and domain (HD offline only runs with one process)
  TYPE(grid_type)   :: grid
  TYPE(domain_type) :: domain

  ! local arrays needed with service format (.srv)
  INTEGER :: IHEAD(8)

  ! logical I/O units
  INTEGER lurun, lubas, luout

  ! file names
  CHARACTER(80) :: runoff_file    ! input file with runoff data
  CHARACTER(80) :: drainage_file  ! input file with drainage data
  CHARACTER(80) :: dnout          ! averaged output data (nml parameter iout)

  ! local parameters
  INTEGER, PARAMETER :: DAILY = 1    ! used for frequency of forcing data
  INTEGER, PARAMETER :: STEPWISE = 0

  !----------------
  ! initialization
  !----------------

  hydrology_offline = .TRUE.

  ! MPI initialization

  CALL p_start
  CALL p_init_communicators(1, 1, 0)    ! single process: nproca=1, nprocb=1, nprocio=0
  CALL machine_setup

  ! read the namelist

  CALL config_hd

  ! initialization of the echam time manager

  dt_start = (/year1,1,1,0,0,0/)                ! start date: Jan 1st of year1
  no_steps = nstep                              ! number of time steps within the run
  lresume = .FALSE.

  CALL get_dates(IO_timestep, istep)            ! get timestep
  CALL ec_manager_init(IO_timestep, istep)      ! time manager initialization
  CALL init_manager
  CALL init_times                               ! initialize all dates
  l_trigfiles = .FALSE.

  ! hd model initializations
 
  CALL hd_init_dims                             ! init HD model grid dimensions
  CALL hd_init_forcing                          ! allocate memory of the forcing arrays
  CALL init_hydrology(grid, domain, jsbach_offline=.FALSE.)
  CALL hd_init_io                               ! open in- and output files

  !----------------
  ! Time step loop
  !----------------

  DO step = istep+1, istep+nstep

     ! set time for echam time manager
     CALL time_set 

     IF ((forcing_freq == STEPWISE) .OR. &
          forcing_freq == DAILY .AND. MOD(step-1,hd_steps_per_day) == 0) THEN
        ! read forcing data
        CALL write_date(current_date, 'Read forcing data of :')

        CALL hd_update_forcing         ! read runoff and drainage from file

     END IF

     ! discharge calculations
     CALL hydrology_model((slm == 1._dp), jsbach_offline=.FALSE., &
          aros_offline=runoff, adrain_offline=drain, disch_offline=disch, &
          awfre_offline=awfre, apmecal_offline=apmecal)

     ! write the output
     CALL hd_write_output

     ! write restart file
     IF (step == istep+nstep) CALL hydrology_restart 

     ! update model time step
     CALL time_reset

  ENDDO

  !   *** Schreiben der letzten Mittelwert-Daten fuer Inflow per Gridbox 
  IDUM=-1
  CALL GMITWRI(luout, IHEAD, FRIV + water_to_ocean, IDUM)

  !-------------
  ! cleaning up
  !-------------

  CALL cleanup_hydrology

  ! close files

  CLOSE(lurun)
  CLOSE(luout)

  IF (locean) CLOSE(luoce)
  IF (lbase) CLOSE(lubas)

  ! MPI finalization

   call p_stop

CONTAINS
!
!****************************************************************************
  SUBROUTINE GMITWRI(LUOUT, IHEAD, FOUT, IOUT)
!****************************************************************************
!
!     ******** Routine zur Mittelung ueber NT Zeitschritte und Ausgabe
!              in einer Serviceformat-Binaerdatei.c 
!
!     ******** Version 1.0 - Oktober 1995
!              Programmierung und Entwicklung: Stefan Hagemann 
!
!     ******** Version 1.1 - Oktober 1996
!              Nun auch Monatsmittelwerte mit Schaltjahren moeglich (IOUT=5)
!
!     ******** Version 1.2 - Oktober 1999
!              Serviceformat: simple binaer statt direct access.
!
!     ******** Version 1.3 - May 2002
!              Daily Output possible --> IOUT = 6
!
!     ********* Variablenliste
!     ***
!     ***  LUOUT = Logical Unit der binaeren DIRECT Access Datei
!     ***  IHEAD = Header-Array fuer Serviceformat-Dateien
!     ***   FOUT = Globales Output-Array des Timesteps IREC
!     ***     NL = Anzahl der Laengengrade
!     ***     NB = Anzahl der Breitengrade
!     ***   LDEBUGHD = Kommentarvariable ( 0 = Kein Kommentar )
!     ***   IOUT = Mittelungsartvariable, d.h. ueber wieviel Zeitschritte
!     ***          1   30-Day Averages   --> NT = 30 * hd_steps_per_day
!     ***          2   Decadal Averages  --> NT = 10 * hd_steps_per_day
!     ***          3   Weekly Averages   --> NT = 7  * hd_steps_per_day
!     ***          4   Monthly Averages ohne Schaltjahre
!     ***          5   Monthly Averages inklusive Schaltjahre
!     ***          6   Daily Output
!     ***          -1  Mittelung ueber restliche Zeitschritte und Return
!     ***
!     ***   IREC = Timestep-Zaehler des gemittelten Feldes FSUM
!     ***   IMON = Monatsnummer bzw. New Record Ja/Nein-Variable 
!     ***          -1 ==> Beginn eines New Records
!     ***   it = Lokaler Timestep-Zaehler des Records IREC+1
!     ***
!     ******** Include of Parameter NL, NB 
!     ***       NL, NB = Globale/Regionale Feldgrenzen

    USE mo_time_control,     ONLY: get_date_components, get_month_len, current_date 

    INTEGER,  INTENT(in)  :: luout
    INTEGER,  INTENT(in)  :: iout
    REAL(dp), INTENT(in)  :: fout(nl,nb)
    INTEGER,  INTENT(out) :: ihead(8)

    REAL(dp), SAVE :: fsum(nl,nb)
    INTEGER,  SAVE :: nt                          ! number of time steps
    INTEGER        :: yr, mon, day, hr, min, sec  ! current date components
    INTEGER        :: date, time
    INTEGER,  SAVE :: irec = 0, nday(12)
    INTEGER,  SAVE :: it
    LOGICAL,  SAVE :: start_accumulation = .TRUE.


    ! initializations
    nday = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    ! define the header (service format)

    CALL get_date_components(current_date, yr, mon, day, hr, min, sec)
    date = yr*10000 + mon*100 + day
    time = hr*10000 + min*100 + sec

    ihead(1) = 0    ! no code number for river discharge available
    ihead(2) = 1    ! level
    ihead(3) = date
    ihead(4) = time
    ihead(5) = nl
    ihead(6) = nb

    IF (iout == -1) THEN        ! only at the very end
       IF (start_accumulation) RETURN   ! nothing to do, averaging had just been done

       ! averaging the remaining time steps
       fsum(:,:) = fsum(:,:) / nstep

       irec = irec+1
       ihead(7) = irec

       WRITE(luout) ihead
       WRITE(luout) fsum

       start_accumulation = .TRUE.
       RETURN
    ENDIF


    IF (start_accumulation) THEN

       ! initializations
       it = 0
       fsum(:,:) =0.
       start_accumulation = .FALSE.

       ! find number of time_steps for the averaging

       IF (iout == 1) THEN

          ! 30-day averages
          nt = 30 * hd_steps_per_day

       ELSE IF (iout == 2) THEN

          ! 10-day averages
          nt = 10 * hd_steps_per_day

       ELSE IF (iout == 3) THEN

          ! weekly
          nt = 7 * hd_steps_per_day

       ELSE IF (iout == 4) THEN

          ! monthly averages, without leap year
          nt = nday(mon) * hd_steps_per_day

       ELSE IF (iout == 5) THEN

          ! monthly averages, with leap year
          nt = get_month_len(yr,mon) * hd_steps_per_day

       ELSE IF (iout == 6) THEN

          ! daily averages
          nt = hd_steps_per_day

       ELSE IF (iout == 7) THEN

          ! stepwise output (each time routine is called)
          nt = 1

       ELSE
          WRITE (message_text,*) 'namelist parameter IOUT out of range: iout = ', iout
          CALL finish ('hd_gmitwri', message_text)
       ENDIF
    ENDIF

    ! accumulate data for the selected time period

    it = it + 1
    fsum(:,:) = fsum(:,:) + fout(:,:)
    IF (it == nt) THEN

       ! do the averaging

       fsum(:,:) = fsum(:,:) / nt
       irec = irec+1

       ihead(7) = irec
       WRITE(luout) ihead
       WRITE(luout) fsum

       start_accumulation = .TRUE.
    ENDIF

  END SUBROUTINE GMITWRI

  !------------------------------------------------------------------------------------
  ! Basic configuration of HD offline simulations
  !------------------------------------------------------------------------------------
  SUBROUTINE config_hd

    !----------------------------------------------------------------------------------
    ! parameters of namelist HDALONE_CTL
    !
    !        iswrit    restart time step (0: no restart file)
    !   out_expname    experiment name
    !  out_datapath    path to where the output data shall be written
    !         year1    initial year of the run
    !         nstep    number of time steps within the run
    !    delta_time    model time step lenght in seconds
    !   runoff_file    file with input runoff data
    ! drainage_file    file with input drainage data
    !        ufakru    unit factor for runoff and drainage input data
    !  forcing_freq    data frequency of the forcing data (STEPWISE/DAILY)
    !          iout    averaging period of some HD output
    !  do_remapping    type of interpolation from input (atmospheric) grid to HD grid
    !                         0   no interpolation: input data grid equals HD grid
    !                         1   remapping in routine hd_remap (default)
    !----------------------------------------------------------------------------------

    USE mo_namelist,         ONLY: open_nml, position_nml, POSITIONED
    USE mo_io_units,         ONLY: nout
 
    ! local variables

    INTEGER                :: read_status, inml, iunit

    INCLUDE 'hdalone_ctl.inc'

    ! set default values of the namelist parmeters

    iswrit = 0             ! time step to write restart file (0: no restart file)
    out_expname = 'hd'
    out_datapath = './'
    year1 = 1900
    nstep = 365
    delta_time = 86400._dp        ! time step in seconds (one day) 
    ufakru = 1._dp
    runoff_file = "runoff.nc"
    drainage_file = "drainage.nc"
    forcing_freq = STEPWISE
    iout = 5               ! averaging of the HD output: 1=30d, 2=10d, 3=7d, 5=monthly
    do_remapping = 1

    ! read namelist hdalone_ctl

    inml = open_nml ('namelist.hd')
    iunit = position_nml ('HDALONE_CTL', inml, status=read_status)
    SELECT CASE (read_status)
    CASE (POSITIONED)
       READ (iunit, hdalone_ctl)
       CALL message('config_hd', 'Namelist HDALONE_CTL: ')
       WRITE(nout, hdalone_ctl)
    END SELECT

  END SUBROUTINE config_hd

  !------------------------------------------------------------------------------------
  ! Read forcing coordinates and masks 
  !
  ! Array allocations for HD model simulations. In coupled echam/HD runs
  ! these variables are provided by echam.
  !------------------------------------------------------------------------------------
  SUBROUTINE hd_init_forcing

    USE mo_gaussgrid,     ONLY: gauaw, gridarea, philat, philon
    USE mo_physical_constants, ONLY: earth_radius
    USE mo_math_constants,     ONLY: pi
    USE mo_netcdf,        ONLY: file_info, nf_max_name, io_inq_dimid, &
                                io_inq_dimlen, io_inq_varid, io_get_var_double
    USE mo_io,            ONLY: io_open, io_close, io_read

    TYPE (file_info)  :: gridfile

    INTEGER dimid, varid, fileid

    CHARACTER(nf_max_name) :: filename
    LOGICAL           :: lex
    REAL(dp), ALLOCATABLE :: zgw(:), zgmu(:)
    INTEGER           :: jgl

    ! read land sea masks

    filename = 'jsbach.nc'

    INQUIRE (file=filename, exist=lex)
    IF (.NOT. lex) THEN
      WRITE (message_text,*) 'Could not open file <',TRIM(filename),'>'
      CALL message('read_hd_forcing', message_text)
      CALL finish ('read_hd_forcing', 'run terminated.')
    ENDIF

    gridfile%opened = .FALSE.
    CALL IO_open (filename, gridfile, IO_READ)
    WRITE (message_text,*) 'Reading grid dimensions from file ', TRIM(filename)
    CALL message('read_hd_forcing', message_text)

    fileID = gridfile%file_id

    CALL IO_inq_dimid (fileID, 'lon', dimid)
    CALL IO_inq_dimlen (fileID, dimid, nlon)
    ALLOCATE (philon(nlon))
    CALL IO_inq_varid (fileID, 'lon', varid)
    CALL IO_get_var_double (fileID, varid, philon)

    CALL IO_inq_dimid (fileID, 'lat', dimid)
    CALL IO_inq_dimlen (fileID, dimid, ngl)
    ALLOCATE (philat(ngl))
    CALL IO_inq_varid (fileID, 'lat', varid)
    CALL IO_get_var_double (fileID, varid, philat)

    CALL IO_inq_varid (fileID, 'slm', varid)
    ALLOCATE(slm(nlon,ngl))
    CALL IO_get_var_double (fileID, varid, slm)

    CALL IO_close(gridfile)

    !-- definition of grid and domain types

    grid%nlon  = nlon
    grid%nlat  = ngl
    grid%nland = NINT(SUM(slm))

    domain%ndim  = nlon
    domain%nblocks = ngl
    domain%nland = NINT(SUM(slm))

    !-- Compute Gaussian weights

    ALLOCATE (zgw(ngl))
    ALLOCATE (zgmu(ngl))
    ALLOCATE (gridarea(ngl))

    CALL gauaw(zgmu, zgw, ngl)
    DO jgl = 1, ngl
      gridarea(jgl)  = 0.5_dp * zgw(jgl)/nlon * 4 * pi * earth_radius**2
    END DO

    DEALLOCATE (zgw, zgmu)

    !-- Memory allocation

    ALLOCATE(runoff(nlon,ngl))
    ALLOCATE(drain(nlon,ngl))
    ALLOCATE(disch(nlon,ngl))

    ALLOCATE(awfre(nlon,ngl))
    awfre(:,:) = 0._dp
    ALLOCATE(apmecal(nlon,ngl))
    apmecal(:,:) = 0._dp

    IF (do_remapping == 1) THEN
       ALLOCATE(dummy_sp(nlon,ngl))
    ELSE
       ALLOCATE(dummy_sp(nl,nb))
    ENDIF

  END SUBROUTINE hd_init_forcing

  !------------------------------------------------------------------------------------
  ! init netcdf dimensions. This is needed with IO_open (hd_init_forcing). The
  ! dimension lenghts given here (HD grid) do not matter. 
  !------------------------------------------------------------------------------------
  SUBROUTINE hd_init_dims
    USE mo_netcdf, ONLY: add_dim

    CALL add_dim ("lon", nl, "longitude", "degrees_east" )
    CALL add_dim ("lat", nb, "latitude", "degrees_north" )

  END SUBROUTINE hd_init_dims

  !------------------------------------------------------------------------------------
  ! open hydrology model in- and output files
  !------------------------------------------------------------------------------------
  SUBROUTINE hd_init_io

    INTEGER  :: ios          ! I/O status

    ! open file for outflow diagnostics on the hd model grid
    IF (lhd_highres) CALL hd_highres_open

    ! open file for outflow data on the forcing data grid

    IF (locean) THEN
       luoce = find_next_free_unit (51,100)
       OPEN(luoce, FILE="discharge.srv", FORM='unformatted', IOSTAT=ios)

       IF (ios /= 0) THEN
          WRITE (message_text,*) 'Error opening file discharge.srv'
          CALL finish ('hd_init_io', message_text)
       ENDIF
    END IF

    ! open input file with overlandflow (= runoff) data

    lurun = find_next_free_unit (51,100)
    OPEN(lurun, FILE=runoff_file, FORM='unformatted', STATUS='old', IOSTAT=ios)
    IF (ios /= 0) THEN
       WRITE (message_text,*) 'Error opening file ', runoff_file
       CALL finish ('hd_init_io', message_text)
    ENDIF

    ! open input file with baseflow (= drainage) data

    IF (lbase) THEN
       lubas = find_next_free_unit (51,100)
       OPEN(lubas, FILE=drainage_file, FORM='unformatted', STATUS='old', IOSTAT=ios)
       IF (ios /= 0) THEN
          WRITE (message_text,*) 'Error opening file ', drainage_file
          CALL finish ('hd_init_io', message_text)
       ENDIF
    ENDIF

    ! open output file for time averaged outflow data

    dnout = "meanflowbin.srv"
    luout = find_next_free_unit (51,100)
    OPEN(luout, FILE=dnout, form='unformatted', IOSTAT=IOS)

    IF (IOS.NE.0) THEN
       WRITE (message_text,*) 'Error opening file ', dnout
       CALL finish ('hd_init_io', message_text)
    ENDIF

  END SUBROUTINE hd_init_io

  !------------------------------------------------------------------------------------
  ! read forcing data for the current time step
  !------------------------------------------------------------------------------------
  SUBROUTINE hd_update_forcing
    
    READ(lurun) IHEAD
    READ(lurun) dummy_sp
    IF (do_remapping == 1) THEN
       runoff(:,:) = REAL(dummy_sp,dp)
    ELSE
       runoff(:,:) = 0._dp
       runoff_s(:,:) = REAL(dummy_sp,dp)
    END IF
    IF (ABS(UFAKRU-1.) .GT. 1.E-6) RUNOFF=RUNOFF*UFAKRU
    IF (lbase) THEN
       READ(lubas) IHEAD
       READ(lubas) dummy_sp
       IF (do_remapping == 1) THEN
          drain(:,:) = REAL(dummy_sp,dp)
       ELSE
          drain(:,:) = 0._dp
          runoff_dr(:,:) = REAL(dummy_sp,dp)
       END IF
       IF (ABS(UFAKRU-1.) .GT. 1.E-6) DRAIN=DRAIN*UFAKRU
    ENDIF

  END SUBROUTINE hd_update_forcing

  !------------------------------------------------------------------------------------
  ! write hydrology output of the current time step
  !------------------------------------------------------------------------------------
  SUBROUTINE hd_write_output

    ! write inflow data for each gridbox. 
    ! Note: The time steps are shifted backwards (time step 2 -> 1). The current inflow
    ! is flowing at the end of the current time step, it thus corresponds to the
    ! beginning of the following time step.

    CALL GMITWRI(luout, ihead, friv + water_to_ocean, iout)

    ! write the discharge aray (on the grid of the forcing data)

    IF (locean) THEN
       ihead(1) = 219
       ihead(5) = nlon
       ihead(6) = ngl
       WRITE(luoce) ihead
       WRITE(luoce) disch
    ENDIF

  END SUBROUTINE hd_write_output

END PROGRAM hd_driver
