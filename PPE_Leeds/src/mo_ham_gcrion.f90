!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham_gcrion.f90
!!
!! \brief
!! mo_ham_gcrion provides routines and quantities for the calculation of the
!! galactic cosmic ray ionization rate.
!!
!! \author Jan Kazil (MPI-Met)
!!
!! \responsible_coder
!! Jan Kazil, jan.kazil@zmaw.de  
!!
!! \revision_history
!!   -# J. Kazil (MPI-Met) - original code (2008-02-21)
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ham_gcrion
  
  USE mo_io_units,   ONLY: nerr        !### replace with message from mo_exception
  USE mo_exception,  ONLY: finish
  USE mo_filename,   ONLY: find_next_free_unit
  USE mo_kind,       ONLY: dp
  USE mo_mpi,        ONLY: p_parallel,p_parallel_io,p_io,p_bcast
  USE mo_geopack,    ONLY: geo2mag,vertical_cutoff_rigidity
  
  IMPLICIT NONE
  
  !
  ! Public elements:
  !
  
  PUBLIC :: gcr_initialize, gcr_ionization, gcr_cleanup
  
  !
  ! Private elements:
  !
  
  PRIVATE :: read_obrien_gcr_ipr
  PRIVATE :: gcr_ionization_profile
  
  INTEGER, PRIVATE :: vcr_n ! Number of vertical cutoff rigidities
  INTEGER, PRIVATE :: mcd_n ! Number of mass column densities
  
  ! Vertical cutoff rigidity lookup table:
  REAL(dp), PRIVATE, ALLOCATABLE :: vcr_table(:) ! GV
  
  ! Mass column density lookup table:
  REAL(dp), PRIVATE, ALLOCATABLE :: mcd_table(:) ! g cm-2
  
  ! Solar minimum and maximum ion pair production rate lookup tables:
  REAL(dp), PRIVATE, ALLOCATABLE :: gcr_ipr_solmin_table(:,:) ! Ion pairs cm-3 s-1 @ 273.15 K, 1013.25 hPa
  REAL(dp), PRIVATE, ALLOCATABLE :: gcr_ipr_solmax_table(:,:) ! Ion pairs cm-3 s-1 @ 273.15 K, 1013.25 hPa
  
CONTAINS
  
  SUBROUTINE gcr_initialize()
    
    ! *gcr_initialize* reads Keran O'Brien's solar minimum and maximum galactic
    ! cosmic ray ionization rates from file and broadcasts them to all processors.
    
    IMPLICIT NONE
    
    vcr_n = 15
    mcd_n = 110
    
    ALLOCATE(mcd_table(mcd_n))
    ALLOCATE(vcr_table(vcr_n))
    
    ALLOCATE(gcr_ipr_solmin_table(vcr_n,mcd_n))
    ALLOCATE(gcr_ipr_solmax_table(vcr_n,mcd_n))
    
    IF (p_parallel) THEN ! In parallel mode
      
      IF (p_parallel_io) THEN ! We are on the I/O processor
        
        CALL read_obrien_gcr_ipr()
        
      ENDIF
      
      CALL p_bcast(vcr_n,p_io)
      CALL p_bcast(mcd_n,p_io)
      
      CALL p_bcast(vcr_table,p_io)
      CALL p_bcast(mcd_table,p_io)
      
      CALL p_bcast(gcr_ipr_solmin_table,p_io)
      CALL p_bcast(gcr_ipr_solmax_table,p_io)
      
    ELSE ! In serial mode
      
      CALL read_obrien_gcr_ipr()
      
    ENDIF
    
  END SUBROUTINE gcr_initialize
  
  !=============================================================================
  
  SUBROUTINE read_obrien_gcr_ipr()
    
    ! *read_obrien_gcr_ipr* reads Keran O'Brien's solar minimum and maxium ion
    ! pair production rates from files.
    
    IMPLICIT NONE
    
    !
    ! Local variables:
    !
    
    CHARACTER(32) :: sol_min_file,sol_max_file ! File names
    CHARACTER(32) :: string
    
    INTEGER :: iunit ! File number
    LOGICAL :: lex ! Flag to test if file exists
    
    INTEGER :: ji,jj ! Loop variables
    
    sol_min_file = 'gcr_ipr_solmin.txt'
    sol_max_file = 'gcr_ipr_solmax.txt'
    
    ! Test if files exist:
    
    INQUIRE(file=sol_min_file,exist=lex)
    
    IF (.not. lex) THEN
      CALL finish('read_obrien_gcr_ipr','File not found : '//sol_min_file)
    ENDIF
    
    INQUIRE(file=sol_max_file,exist=lex)
    
    IF (.not. lex) THEN
      CALL finish('read_obrien_gcr_ipr','File not found : '//sol_max_file)
    ENDIF
    
    ! Read the solar minimum data:
    
    iunit = find_next_free_unit(20,99)
    
    OPEN (iunit,file=sol_min_file,form='FORMATTED')
    
    READ(iunit,*)
    
    DO ji = 1, vcr_n
      
      READ(iunit,*)
      READ(iunit,'(F11.6,A)') vcr_table(ji),string
      READ(iunit,*)
      READ(iunit,*)
      
      DO jj=1, mcd_n
        
        READ(iunit,*) mcd_table(jj),gcr_ipr_solmin_table(ji,jj)
        
      ENDDO
      
    ENDDO
    
    close(iunit)
    
    ! Read the solar maximum data:
    
    iunit = find_next_free_unit(20,99)
    
    OPEN (iunit,file=sol_max_file,form='FORMATTED')
    
    READ(iunit,*)
    
    DO ji=1, vcr_n
      
      READ(iunit,*)
      READ(iunit,'(F11.6,A)') vcr_table(ji),string
      READ(iunit,*)
      READ(iunit,*)
      
      DO jj=1, mcd_n
        
        READ(iunit,*) mcd_table(jj),gcr_ipr_solmax_table(ji,jj)
        
      ENDDO
      
    ENDDO
    
    close(iunit)
    
  END SUBROUTINE read_obrien_gcr_ipr
  
  !=============================================================================
  
  SUBROUTINE gcr_cleanup()
    
    ! *gcr_cleanup* deallocates the memory of the module mo_gcr.
    
    IMPLICIT NONE
    
    DEALLOCATE(mcd_table)
    DEALLOCATE(vcr_table)
    
    DEALLOCATE(gcr_ipr_solmin_table)
    DEALLOCATE(gcr_ipr_solmax_table)
    
  END SUBROUTINE gcr_cleanup
  
  !=============================================================================
  
  SUBROUTINE gcr_ionization(krow,kproma,kbdim,klev,psolact,ptemp,ppress,pgcripr)
    
    ! *gcr_ionization* calculates the galactic cosmic ray ionization rate.
    
    USE mo_geoloc,       ONLY: philat_2d,philon_2d
    USE mo_time_control, ONLY: current_date,get_date_components,get_year_day
    USE mo_geopack,      ONLY: recalc,geo2mag
    
    IMPLICIT NONE
    
    !
    ! Input variables:
    !
    
    INTEGER,  INTENT(in) :: krow   ! ECHAM5 geographic block number
    INTEGER,  INTENT(in) :: kproma ! Number of locations in the ECHAM5 geographic block
    INTEGER,  INTENT(in) :: kbdim  ! Maximum number of locations in the ECHAM5 geographic block
    INTEGER,  INTENT(in) :: klev   ! Number of levels
    
    REAL(dp), INTENT(in) :: psolact ! Solar activity parameter [-1,1], where -1 represents solar minimum and 1 solar maximum
    
    REAL(dp), INTENT(in) :: ptemp(kbdim,klev)  ! Temperature (K)
    REAL(dp), INTENT(in) :: ppress(kbdim,klev) ! Pressure (Pa)
    
    !
    ! Output variables:
    !
    
    REAL(dp), INTENT(out) :: pgcripr(kbdim,klev) ! ion pairs cm-3 s-1
    
    !
    ! Local variables:
    !
    
    INTEGER :: jj
    
    INTEGER :: iyr,imo,idy,ihr,imn,ise,idoy ! Year, month, day of month, hour, minute, second ,day of year
    
    REAL(dp) :: zmaglon(kbdim),zmaglat(kbdim) ! Geomagnetic longitude and latitude (deg)
    
    REAL(dp) :: zpres(klev) ! Pressure profile [hPa]
    REAL(dp) :: ztemp(klev) ! Temperature profile [K]
    REAL(dp) :: zgcripr(klev) ! Ionization profile [cm-3 s-1]
    
    ! Get date information for the current time step:
    
    CALL get_date_components(current_date,iyr,imo,idy,ihr,imn,ise)
    
    idoy = aint(get_year_day(current_date))
    
    ! Make sure the year is in [1965,2010] for the purpose of calculating the
    ! transformation between geomagnetic and geographic coordinates (The
    ! orientation of the Earth's magnetic dipole is not known in the present
    ! implementation for years outside of this range). The year used here has
    ! *no* effect regarding the calculation of the GCR ionization rate as
    ! function of solar activity.
    
    iyr = max(iyr,1965)
    iyr = min(iyr,2010)
    
    ! Prepare the elements of rotation matrices for transformations between
    ! geographic and magnetic coordinates:
    CALL recalc(iyr,idoy,12,0,0)

    ! Calculate geomagnetic coordinates:
    DO jj=1,kproma ! Geographic locations loop
      CALL geo2mag(philat_2d(jj,krow),philon_2d(jj,krow),zmaglat(jj),zmaglon(jj))
    ENDDO
    
    ! NOTE: Some execution time could be saved by calling recalc and geo2mag
    ! only once a month or once a year, as the Earth's geomagnetic axis moves
    ! very slowly with respect to its rotational axis. This would, however,
    ! require that a way is found to determine when recalc and geo2mag need to
    ! be called, and to store the geomagnetic coordinates until the next call.
    
    ! Calculate the ionization rate:
    
    DO jj=1,kproma ! Geographic locations loop
      
      zpres = 0.01_dp*ppress(jj,1:klev) ! hPa
      ztemp = ptemp(jj,1:klev) ! K
      
      CALL gcr_ionization_profile(psolact,zmaglat(jj),klev,zpres,ztemp,zgcripr)
      
      pgcripr(jj,1:klev) = zgcripr(1:klev)
      
    ENDDO
    
  END SUBROUTINE gcr_ionization
  
  !============================================================================

  SUBROUTINE gcr_ionization_profile(psolact,pmaglat,klev,ppress,ptemp,pgcripr)
    
    ! *gcr_ionization_profile* returns a galactic cosmic ray ionization rate
    ! profile based on Keran O'Brien's solar minimum and maximum model
    ! ionization rates. The parameter psolact determines the solar cycle phase.
    
    USE mo_physical_constants, ONLY: grav
    
    IMPLICIT NONE
    
    !
    ! Input variables:
    !
    
    ! Solar activity parameter: Value in [-1,1], where -1 represents solar
    ! minimum, 1 solar maximum:
    REAL(dp), INTENT(in) :: psolact
    
    ! Geomagnetic coordinates:
    REAL(dp), INTENT(in) :: pmaglat ! degrees
    
    ! Number of levels:
    INTEGER,  INTENT(in) :: klev
    
    ! The pressures which the ion pair production rates are calculated for:
    REAL(dp), INTENT(in) :: ppress(klev) ! hPa
    
    ! The temperatures which the ion pair production rates are calculated for:
    REAL(dp), INTENT(in) :: ptemp(klev) ! K
    
    !
    ! Output variables:
    !
    
    ! Ion pair production rate profile:
    REAL(dp), INTENT(out) :: pgcripr(klev) ! ion pairs cm-3 s-1
    
    !
    ! Local variables:
    !
    
    ! Mass column density profile corresponding to the input pressures:
    REAL(dp) :: zmcd(klev) ! (g/cm2)
    
    ! Vertical cutoff rigidity for the input location:
    REAL(dp) :: zvcr ! (GV)
    
    ! Ion pair production rates at the input location and the input pressures,
    ! for solar minimum and maximum:
    REAL(dp) :: zgcriprsolmin(klev)
    REAL(dp) :: zgcriprsolmax(klev)
        
    REAL(dp) :: zv,zw
    
    ! Loop indices:
    
    INTEGER :: ji,jj
    
    INTEGER :: ivcr0,ivcr1,imcd0,imcd1
    
    ! Compute the vertical cutoff rigidity for the current location:
    zvcr = vertical_cutoff_rigidity(pmaglat)
    
    ! Check if the vertical cutoff rigidity for the current location is within
    ! the vertical cutoff rigidities of the lookup table:
    
    IF ((zvcr < vcr_table(1)).or.(zvcr > vcr_table(vcr_n))) THEN
      
      ! Inform the user:
      
!      WRITE(nerr,*) ''
!      WRITE(nerr,*) 'gcr_ionization_profile: Cutoff rigidity out of range.'
!      WRITE(nerr,*) ''
!      WRITE(nerr,*) 'Current cutoff rigidity                             (GV) :', zvcr
!      WRITE(nerr,*) ''
!      WRITE(nerr,*) 'Lookup table top level cutoff rigidity              (GV) :', vcr_table(1)
!      WRITE(nerr,*) 'Lookup table bottom level cutoff rigidity           (GV) :', vcr_table(vcr_n)
      
      IF (zvcr.lt.vcr_table(1)) THEN
        zvcr = vcr_table(1)
        ivcr0 = 1
        ivcr1 = 2
      ENDIF
      
      IF (zvcr.gt.vcr_table(vcr_n)) THEN
        zvcr = vcr_table(vcr_n)
        ivcr0 = vcr_n - 1
        ivcr1 = vcr_n
      ENDIF
      
      ! Inform the user:
      
!      WRITE(nerr,*) ''
!      WRITE(nerr,*) 'Resetting cutoff rigidity to                        (GV) :', zvcr
!      WRITE(nerr,*) ''
        
    ENDIF
    
    ! Determine the indices in the vertical cutoff rigidity array containing the
    ! current vertical cutoff rigidity:
    
    ivcr0 = 1
    ivcr1 = vcr_n
    
    DO WHILE (ivcr1-ivcr0.gt.1)
      
      ji = (ivcr1+ivcr0)/2
      
      IF (vcr_table(ji).gt.zvcr) THEN
        ivcr1 = ji
      ELSE
        ivcr0 = ji
      ENDIF
      
    ENDDO
    
    ! Compute the mass column densities for the given pressure levels, assuming
    ! that gravity is independent of altitude:
    
    zmcd(1:klev) = 10.0_dp*ppress(1:klev)/grav ! g/cm2
    
    ! Interpolate the solar minimum and maximum data for the given cutoff
    ! rigidity and column mass densities:
    
    DO ji = 1, klev
      
      ! Determine the indices in the mass column density array containing the
      ! current mass column density:
      
      imcd0 = 0
      imcd1 = 0
      
      DO jj = 2, mcd_n
        IF ((zmcd(ji).ge.mcd_table(jj-1)).and.(zmcd(ji).le.mcd_table(jj))) THEN
          imcd0 = jj - 1
          imcd1 = jj
        EXIT
        ENDIF
      ENDDO
      
      ! If the mass column density for the current location is outside the
      ! mass column densities of the lookup table then we replace it with the
      ! mass column density of the closer end of the lookup table:
      
      IF (imcd0*imcd1 == 0) THEN
        
        ! Inform the user:
        
!        WRITE(nerr,*) ''
!        WRITE(nerr,*) 'gcr_ionization_profile: Mass column density out of range.'
!        WRITE(nerr,*) ''
!        WRITE(nerr,*) 'Current mass column density                   (g cm-2) :', zmcd(ji)
!        WRITE(nerr,*) ''
!        WRITE(nerr,*) 'Lookup table top level mass column density    (g cm-2) :', mcd_table(1)
!        WRITE(nerr,*) 'Lookup table bottom level mass column density (g cm-2) :', mcd_table(mcd_n)
        
        IF (zmcd(ji).lt.mcd_table(1)) THEN
          zmcd(ji) = mcd_table(1)
          imcd0 = 1
          imcd1 = 2
        ENDIF
        
        IF (zmcd(ji).gt.mcd_table(mcd_n)) THEN
          zmcd(ji) = mcd_table(mcd_n)
          imcd0 = mcd_n - 1
          imcd1 = mcd_n
        ENDIF
        
        ! Inform the user:
        
!        WRITE(nerr,*) ''
!        WRITE(nerr,*) 'Resetting mass column density to              (g cm-2) :', zmcd(ji)
!        WRITE(nerr,*) ''
        
      ENDIF
      
      ! 2D interpolation:
      
      zv = (zvcr-vcr_table(ivcr0))/(vcr_table(ivcr1)-vcr_table(ivcr0))
      zw = (zmcd(ji)-mcd_table(imcd0))/(mcd_table(imcd1)-mcd_table(imcd0))
      
      zgcriprsolmin(ji) = (zv-1.0_dp)*(zw-1.0_dp)*gcr_ipr_solmin_table(ivcr0,imcd0) &
       + (zw-zv*zw)*gcr_ipr_solmin_table(ivcr0,imcd1) &
       + (zv-zv*zw)*gcr_ipr_solmin_table(ivcr1,imcd0) &
       + zv*zw*gcr_ipr_solmin_table(ivcr1,imcd1)
      
      zgcriprsolmax(ji) = (zv-1.0_dp)*(zw-1.0_dp)*gcr_ipr_solmax_table(ivcr0,imcd0) &
       + (zw-zv*zw)*gcr_ipr_solmax_table(ivcr0,imcd1) &
       + (zv-zv*zw)*gcr_ipr_solmax_table(ivcr1,imcd0) &
       + zv*zw*gcr_ipr_solmax_table(ivcr1,imcd1)
      
    ENDDO
    
    ! Interpolate the solar minimum and maximum profiles
    ! for the given solar activity phase:
    
    pgcripr(1:klev) = 0.5_dp* &
     ((1.0_dp-psolact)*zgcriprsolmin(1:klev) &
     +(1.0_dp+psolact)*zgcriprsolmax(1:klev))
    
    ! Transform the ionization rate from normal to local conditions:
    
    pgcripr(1:klev) = pgcripr(1:klev)*ppress(1:klev)/1013.25_dp*273.15_dp/ptemp(1:klev)
    
    END SUBROUTINE gcr_ionization_profile
    
    !=============================================================================
  
    REAL FUNCTION solar_activity()
    
    ! *solar_activity* parameterizes the solar activity with a cosine peaking 1
    ! January 1991 (solar maximum), and dipping 2 July 1996 (solar minimum),
    ! corresponding to a solar cycle length of 11 years. Hence a value of 1 is
    ! returned at solar maximum, a value of -1 at solar minimum, and values
    ! ]-1,1[ in between. The parameterized solar activity will follow the
    ! observed solar activity only approximately, as the  actual solar cycle
    ! length and the activity at solar minimum and maximum change from cycle to
    ! cycle.
    
    USE mo_math_constants, ONLY: pi
    USE mo_time_control, ONLY: current_date, &
                               get_date_components, &
                               get_year_day, &
                               get_year_len
    
    IMPLICIT NONE
    
    INTEGER :: iyr,imo,idy,ihr,imn,ise ! Year, month, day of month, hour, minute, second
    REAL(dp):: zdoy ! Day of year
    REAL(dp):: ztime
    
    CALL get_date_components(current_date,iyr,imo,idy,ihr,imn,ise)
    
    zdoy = get_year_day(current_date) ! Day of the year (1 January = 1)
    
    ztime = iyr + (zdoy-1.0_dp)/get_year_len(iyr) ! Time in years + fraction of the year
    
    solar_activity = cos(2.0_dp*pi*(1991.0_dp-ztime)/11.0_dp)
    
    END FUNCTION solar_activity
    
END MODULE mo_ham_gcrion
