!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_kappa.f90
!!
!! \brief
!! Defines data for ham_kappa routine in mo_m7
!!
!! \author Declan O'Donnell (MPI-M)
!!
!! \responsible_coder
!! Declan O'Donnell, declan.odonnell@fmi.fi
!!
!! \revision_history
!!   -# Declan O'Donnell (MPI-Met) - original code (2007)
!!   -# Declan O'Donnell (FMI) - lookup table version 2 (2014)
!!
!! \limitations
!! None
!!
!! \details
!! ham_kappa implements the parameterisation of the hygroscopic
!! growth of aerosols, as measured by the growth factor (gf), as a
!! function of  the hygroscopicity parameter kappa and 
!! the ambient temperature and relative humidity. The relation is:
!!            [    A   ]      (gf**3 -1)
!!      RH*exp[- ----- ] = ----------------
!!            [  Rd*gf ]   gf**3 - (1-kappa)
!! 
!! (eqn 11 of M.D. Petters and S.M. Kreidenweis, ACP 7, 2007, see 
!! bibliographic references below)
!! 
!! The above-quoted equation from the cited Petters & Kreidenweis paper expresses
!! the hygroscopic growth factor (gf) as a function of aerosol dry radius (Rd),
!! temperature (T), relative humidity (rh) and a substance property denoted kappa
!! that encapsulates hygroscopic properties of that substance. This transcendental
!! equation is solved offline for various Rd, T, rh and kappa of atmospheric
!! relevance and the results stored in a lookup table. 
!! The m7 subroutine ham_kappa uses the ambient T and rh, uses the mode count median 
!! dry radius and the volume-weighted average kappa to find the entry point into
!! the GF lookup table. Using the growth factor the wet median radius and water
!! uptake of the mode are calculated.
!!   
!! This module contains relevant constants and a routine to read the lookup table
!!
!! NEW FOR VERSION 2 LOOKUP TABLE
!! ==============================
!! In version 1, the relative humidity axis consists of evenly-spaced increments starting 
!! at 0.15 and ending at 0.95. This has several problems. Firstly, the maximum value is
!! too low which causes a significant underestimate of hygroscopic growth at high RH. 
!! Secondly, the growth factor (GF) is highly sensitive to RH at high RH
!! values, but not at low RH. This means that the coarse, even spacing of the RH axis values
!! results in inaccurate GF estimates at high RH. Thirdly, GF is very low at low RH and 
!! in practice we can neglect it below rh of 0.3 (or perhaps 0.4). 
!! Version 2 of the lookup table implements a RH axis that has uneven value spacing, 
!! monotonically decreasing from the minimum RH (now 0.3) to the maximum (now 0.995).
!!
!! The nth RH axis value is now given by
!! RH(n) = RH(n-1) + deltaRH*c**(n-1),      n=1,2,3...
!! where c is a constant (c < 1) and deltaRH is the initial spacing factor (i.e. RH(1)-RH(0)).
!! c is chosen to be 0.96, RH(0) = 0.30 and deltaRH = 0.03 for practical reasons. 
!! c and deltaRH are obtained from the lookup table by:
!!
!! deltaRH = RH(1) - RH(0)             (actually RH(2)-RH(1) since fortran does not use zero indexing)
!! c = (RH(2) - RH(1))/deltaRH         (again add 1 to index for fortran implementation)
!!
!! Other changes:
!! The water uptake is almost totally insensitive to the temperature value (in the kappa-koehler
!! equation considered as an independent variable from RH). Accordingly, the temperature axis has been
!! reduced in range and more coarsely spaced. Taking advantage of the resulting size (=memory)
!! saving, the dry radius axis is more finely spaced.
!! The kappa axis is slightly reduced by removing values above the range in the model (max. now 1.12).
!! Lastly, version number is included in the file as a string "x.y". x is intended to denote a major
!! revision (which I define as any revision that requires an update to the HAMMOZ fortran code)
!! and y is a minor revision (anything else). Note this means that minor revisions in the file can 
!! result in major changes in model results and vice-versa.
!!
!! NB This implementation is compatible with both versions 1 and 2 of the lookup table. The same
!!    file name is expected (lut_kappa.nc). However, version 1 is DEPRECATED and should be phased out.
!! 
!! \bibliographic_references
!! M.D. Petters and S.M. Kreidenweis, A single parameter representation of hygroscopic 
!! growth and cloud condensation nucleus activity, ACP 7, 1961-1971, 2007
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

MODULE mo_ham_kappa

  !--- inherited types, functions and data
  USE mo_kind,           ONLY: dp
  
  IMPLICIT NONE

  !--- public member functions
  PUBLIC :: start_kappa, term_kappa
  
  !--- module parameters 
  !>>dod redmine #260 ... also removed unused parameter sigma_sa
  INTEGER,  PUBLIC :: lut_kappa_version                 ! Lookup table version
  !<<dod

  !--- public module data 
  !    array dimensions
  REAL(dp), ALLOCATABLE, PUBLIC :: Rd(:)
  REAL(dp), ALLOCATABLE, PUBLIC :: T(:)
  REAL(dp), ALLOCATABLE, PUBLIC :: rh(:)
  REAL(dp), ALLOCATABLE, PUBLIC :: kappa(:)

  !   aerosol hygroscopic growth factor lookup table
  REAL(dp), ALLOCATABLE, PUBLIC :: gf(:,:,:,:)

  !    max and min values in each coordinate
  REAL(dp), PUBLIC :: Rd_min, Rd_max, ln_Rd_min, ln_Rd_max
  REAL(dp), PUBLIC :: T_min, T_max
  REAL(dp), PUBLIC :: rh_min, rh_max
  REAL(dp), PUBLIC :: kappa_min, kappa_max

  !    number of increments on each coordinate axis
  INTEGER, PUBLIC :: N_Rd
  INTEGER, PUBLIC :: N_T
  INTEGER, PUBLIC :: N_rh
  INTEGER, PUBLIC :: N_kappa

  !>>dod redmine #260
  !   For non-linear spacing of the RH axis in version 2: initial rh, initial step size
  !   (plus its reciprocal) and step compression factor (plus its log)
  REAL(dp), PUBLIC :: rh_init_step
  REAL(dp), PUBLIC :: inv_rh_init_step
  REAL(dp), PUBLIC :: rh_step_compress
  REAL(dp), PUBLIC :: log_rh_step_compress
  !<<dod

  !--- local data
  ! 
  !--- parameters
  !
  !    file name for lookup table
  CHARACTER(LEN=32), PARAMETER, PRIVATE :: lut_fn = "lut_kappa.nc"

  !    dimension names in lookup table
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: cdimn_Rd = "DryRadius"
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: cdimn_T  = "Temperature"
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: cdimn_rh = "RelativeHumidity"
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: cdimn_k = "kappa"

  !    variable names in lookup table
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: cvarn_Rd = "DryRadius"
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: cvarn_T = "Temperature"
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: cvarn_rh = "RelativeHumidity"
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: cvarn_k = "kappa"
  CHARACTER(LEN=16), PARAMETER, PRIVATE :: cvarn_gf= "GF"

  !    attribute names in lookup table
  CHARACTER(LEN=7),  PARAMETER, PRIVATE :: cattr_version = "Version"

  !--- member functions
CONTAINS
  
  SUBROUTINE start_kappa

    ! start_kappa reads the growth factor lookup table
    ! start_kappa is called from init_subm in mo_submodel_interface

    USE mo_mpi,           ONLY: p_parallel, p_parallel_io, p_bcast, p_io
    USE mo_exception,     ONLY: finish, message, em_info, em_param, em_warn
    USE mo_submodel,      ONLY: print_value
    USE mo_read_netcdf77, ONLY: read_var_nf77_1d, read_var_nf77_4d
    USE mo_netcdf,        ONLY: nf_check, nf__open, nf_close, nf_nowrite, nf_global, chunksize, &
                                IO_GET_ATT_TEXT, IO_INQ_DIMID, IO_INQ_DIMLEN
                                
    IMPLICIT NONE
    
    !---local data
    INTEGER :: zncid                 ! netcdf file identifier
    INTEGER :: zdimid_Rd, zdimid_T, &! netcdf dimension identifiers
               zdimid_rh, zdimid_k
    INTEGER :: idum

    LOGICAL :: lfilex                ! file existence flag

    CHARACTER(LEN=32) :: attrname    ! needed for reasons of mo_netcdf implementation of IO_GET_ATT_TEXT
    CHARACTER(LEN=3) :: versval      ! version number, as string (x.y)
    CHARACTER(LEN=1) :: versmaj      ! major version number (x), as string

    !--- executable procedure

    !--- check that the lookup table exists
    INQUIRE(file=TRIM(lut_fn), exist=lfilex)

    IF (.NOT. lfilex) CALL finish('mo_ham_kappa.start_kappa', 'Missing lookup table file '//TRIM(lut_fn))
    
    IF (p_parallel_io) THEN

       !---open netcdf file
       CALL nf_check(nf__open(TRIM(lut_fn), nf_nowrite, chunksize, zncid), fname=TRIM(lut_fn))

       !>>dod redmine #260
       !--- read file version (not included in version 1, so default to that version if
       !    the version attribute is missing)
       !    ECHAM's implementation is poor: instead of writing "not available" in the attribute 
       !    value string when the requested attribute is not found, it overwrites the 
       !    attribute name
       attrname = cattr_version
       CALL IO_GET_ATT_TEXT(zncid, nf_global, attrname, versval)
       
       IF (TRIM(attrname) == "not available") THEN
          lut_kappa_version = 1
          versval="1.0"
          CALL message('mo_ham_kappa.start_kappa', &
                       'Deprecated lookup table version! Consider updating to the latest version', level=em_warn)
       ELSE
          versmaj = versval(1:1)
          IF (versmaj == "2") THEN
             lut_kappa_version = 2
          ELSE
             CALL finish('mo_ham_kappa.start_kappa', &
                         'Lookup table '//TRIM(lut_fn)//' version '//versval//' not supported')
          END IF
       END IF
       !<<dod
       
       !---get dimensions
       !   dry radius
       CALL IO_INQ_DIMID(zncid, TRIM(cdimn_Rd), zdimid_Rd)
       CALL IO_INQ_DIMLEN(zncid, zdimid_Rd, N_Rd)

       !    temperature
       CALL IO_INQ_DIMID(zncid, TRIM(cdimn_T), zdimid_T)
       CALL IO_INQ_DIMLEN(zncid, zdimid_T, N_T)

       !    relative humidity
       CALL IO_INQ_DIMID(zncid, TRIM(cdimn_rh), zdimid_rh)
       CALL IO_INQ_DIMLEN(zncid, zdimid_rh, N_rh)

       !    kappa
       CALL IO_INQ_DIMID(zncid, TRIM(cdimn_k), zdimid_k)
       CALL IO_INQ_DIMLEN(zncid, zdimid_k, N_kappa)

    END IF

    !---broadcast version and dimension information over processors
    IF (p_parallel) THEN
       !>>dod redmine #260
       CALL p_bcast(lut_kappa_version, p_io)
       !<<dod
       CALL p_bcast(N_Rd, p_io)
       CALL p_bcast(N_T,  p_io)
       CALL p_bcast(N_rh, p_io)
       CALL p_bcast(N_kappa,  p_io)
    END IF

    !---allocate arrays
    !   1.coordinate variables
    !   In runtime, we use only the minimum and maximum value of each
    !   coordinate variable. So we read these from the lookup file 
    !   on the I/O processor, from which we find Tmin, Tmax, etc. and 
    !   then we broadcast the min/max values over the processors.
    
    IF (p_parallel_io) THEN
       ALLOCATE(Rd(N_Rd))
       ALLOCATE(T(N_T))
       ALLOCATE(rh(N_rh))
       ALLOCATE(kappa(N_kappa))
    END IF

    !---GF lookup. This is needed on every processor
    ALLOCATE(gf(N_Rd, N_T, N_rh, N_kappa))

    !---read variables: start with the coordinate variables, then read the GF 
    !   lookup table

    IF (p_parallel_io) THEN

       !---coordinate variables
       CALL read_var_nf77_1d(lut_fn, TRIM(cdimn_Rd), TRIM(cvarn_Rd), Rd, idum)
       CALL read_var_nf77_1d(lut_fn, TRIM(cdimn_T), TRIM(cvarn_T), T, idum)
       CALL read_var_nf77_1d(lut_fn, TRIM(cdimn_rh), TRIM(cvarn_rh), rh, idum)
       CALL read_var_nf77_1d(lut_fn, TRIM(cdimn_k), TRIM(cvarn_k), kappa, idum)

       !---GF lookup table
       CALL read_var_nf77_4d(lut_fn, TRIM(cdimn_Rd), TRIM(cdimn_T), TRIM(cdimn_rh), &
                                     TRIM(cdimn_k), TRIM(cvarn_gf), gf, idum)   
          
       !---close file
       CALL nf_check(nf_close(zncid), fname=TRIM(lut_fn))

       !---get min and max values of the coordinate variables
       Rd_min = Rd(1)
       Rd_max = Rd(N_Rd)
       ln_Rd_min = LOG(Rd_min)
       ln_Rd_max = LOG(Rd_max)
       T_min = T(1)
       T_max = T(N_T)
       rh_min = rh(1)
       rh_max = rh(N_rh)
       kappa_min = kappa(1)
       kappa_max = kappa(N_kappa)
       !<<dod redmine #260
       rh_init_step = rh(2)-rh(1)
       inv_rh_init_step = 1._dp/rh_init_step
       rh_step_compress = (rh(3)-rh(2))/rh_init_step  ! compression ratio 
       log_rh_step_compress = LOG(rh_step_compress)
       !<<dod
    END IF

    !---Broadcast over processors: min/max values, number of coordinate points
    !   and GF lookup table
    IF (p_parallel) THEN
       CALL p_bcast(Rd_min,     p_io)
       CALL p_bcast(Rd_max,     p_io)
       CALL p_bcast(ln_Rd_min,  p_io)
       CALL p_bcast(ln_Rd_max,  p_io)
       CALL p_bcast(N_Rd,       p_io)
       CALL p_bcast(T_min,      p_io)
       CALL p_bcast(T_max,      p_io)
       CALL p_bcast(N_T,        p_io)
       CALL p_bcast(rh_min,     p_io)
       CALL p_bcast(rh_max,     p_io)
       CALL p_bcast(N_rh,       p_io)
       !>>dod redmine #260
       CALL p_bcast(rh_init_step, p_io)
       CALL p_bcast(inv_rh_init_step, p_io)
       CALL p_bcast(rh_step_compress, p_io)
       CALL p_bcast(log_rh_step_compress, p_io)
       !<<dod
       CALL p_bcast(kappa_min,  p_io)
       CALL p_bcast(kappa_max,  p_io)
       CALL p_bcast(N_kappa,    p_io)

       CALL p_bcast(gf,      p_io) 
    END IF

    !--- report lookup table data
    CALL message('','',level=em_param)
    CALL message('mo_ham_kappa.start_kappa','Aerosol hygroscopic growth lookup table data:',level=em_info)
    CALL message('','',level=em_param)
    CALL message(' Table version: ', versval,level=em_param)
    CALL print_value(' Size of radius dimension : ', N_Rd)
    CALL print_value(' Min. radius: ', Rd_min)
    CALL print_value(' Max. radius: ', Rd_max)
    CALL print_value(' Size of temperature dimension : ', N_T) 
    CALL print_value(' Min. temperature: ', T_min)
    CALL print_value(' Max. temperature: ', T_max)
    CALL print_value(' Size of relative humidity dimension : ', N_rh) 
    CALL print_value(' Min. RH : ', rh_min)
    CALL print_value(' Max. RH : ', rh_max)
    CALL print_value(' Size of hygroscopicity dimension : ', N_kappa )
    CALL print_value(' Min. kappa: ', kappa_min)
    CALL print_value(' Max. kappa: ', kappa_max)
    CALL message('','',level=em_param)

    !---Release memory holding the coordinate variables
    IF (p_parallel_io) THEN
       DEALLOCATE(Rd)
       DEALLOCATE(T)
       DEALLOCATE(rh)
       DEALLOCATE(kappa)
    END IF

  END SUBROUTINE start_kappa


  SUBROUTINE term_kappa

    ! term_kappa frees memory allocated for the growth fator lookup table
    ! term_kappa is called from free_subm_memory in mo_submodel_interface

    IMPLICIT NONE

    ! --- executable procedure

    DEALLOCATE(gf)

  END SUBROUTINE term_kappa

END MODULE mo_ham_kappa
