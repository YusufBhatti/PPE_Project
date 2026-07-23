!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_srtm_netcdf

  USE mo_kind, ONLY : wp
  USE mo_mpi,       ONLY: p_parallel_io, p_io, p_bcast
  USE mo_netcdf,    ONLY: io_inq_varid, io_get_vara_double
  USE mo_io,        ONLY: io_open, io_close, io_read, file_info  
  USE mo_exception, ONLY: finish

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: srtm_read

  INTEGER, PARAMETER :: &
       keylower      = 9,  &
       keyupper      = 5,  &
       Tdiff         = 5,  &
       ps            = 59, &
       plower        = 13, &
       pupper        = 47, &
       Tself         = 10, &
       Tforeignlower = 3,  &
       Tforeignupper = 2,  &
       pforeign      = 4,  &
       T             = 19, &
       band          = 14, &
       GPoint        = 16, &
       GPointSet     = 2

  INTEGER, PARAMETER ::  &
       maxAbsorberNameLength   = 5,  &
       Absorber                = 12, &
       maxKeySpeciesNameLength = 3,  &
       maxKeySpeciesNames      = 2

  CHARACTER(len = maxAbsorberNameLength), DIMENSION(Absorber), PARAMETER :: &
       AbsorberNames = (/        &
       'N2   ',  &
       'CCL4 ',  &
       'CFC11',  &
       'CFC12',  &
       'CFC22',  &
       'H2O  ',  &
       'CO2  ',  &
       'O3   ',  &
       'N2O  ',  & 
       'CO   ',  &
       'CH4  ',  &
       'O2   '  /)

  CHARACTER(len = maxKeySpeciesNameLength),   &
       DIMENSION(band,maxKeySpeciesNames), PARAMETER :: &
       KeySpeciesNamesLower = RESHAPE( SOURCE = (/      &
       'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', 'H2O', &
       'H2O', 'H2O', 'H2O', '   ', 'O3 ', 'O3 ', 'H2O', &
       'CH4', 'CO2', 'CH4', 'CO2', '   ', 'CO2', 'O2 ', &
       '   ', 'O2 ', '   ', '   ', '   ', 'O2 ', '   '  /), &
       SHAPE = (/ band, maxKeySpeciesNames /) )

  CHARACTER(len = maxKeySpeciesNameLength),  &
       DIMENSION(band,maxKeySpeciesNames), PARAMETER :: &
       KeySpeciesNamesUpper = RESHAPE( SOURCE = (/      &
       'CH4', 'H2O', 'CH4', 'CO2', 'H2O', 'H2O', 'O2 ', &
       '   ', 'O2 ', '   ', '   ', 'O3 ', 'O3 ', 'CO2', &
       '   ', 'CO2', '   ', '   ', '   ', 'CO2', '   ', &
       '   ', '   ', '   ', '   ', '   ', 'O2 ', '   '  /), &
       SHAPE = (/ band, maxKeySpeciesNames /) )

  INTEGER, PARAMETER :: gPointSetNumber = 1
  INTEGER :: varid

  TYPE(file_info) :: rrtmg_sw

CONTAINS 

  SUBROUTINE getAbsorberIndex(AbsorberName,AbsorberIndex)

    CHARACTER(len = *), INTENT(in) :: AbsorberName
    INTEGER, INTENT(out)           :: AbsorberIndex

    INTEGER :: m

    AbsorberIndex = -1
    DO m = 1, Absorber
      IF (TRIM(AbsorberNames(m)) == TRIM(AbsorberName)) THEN
        AbsorberIndex = m
      END IF
    END DO

    IF (AbsorberIndex == -1) THEN
      CALL finish('Absorber name index lookup failed.')
    END IF
  END SUBROUTINE getAbsorberIndex

  !=============================================================================

  SUBROUTINE srtm_read

    IF (p_parallel_io) THEN
      CALL io_open ('rrtmg_sw.nc', rrtmg_sw, io_read)
    END IF

      CALL sw_kgb16  ! molecular absorption coefficients
      CALL sw_kgb17
      CALL sw_kgb18
      CALL sw_kgb19
      CALL sw_kgb20
      CALL sw_kgb21
      CALL sw_kgb22
      CALL sw_kgb23
      CALL sw_kgb24
      CALL sw_kgb25
      CALL sw_kgb26
      CALL sw_kgb27
      CALL sw_kgb28
      CALL sw_kgb29

    IF (p_parallel_io) THEN
      CALL io_close(rrtmg_sw)
    ENDIF

  END SUBROUTINE srtm_read

  SUBROUTINE sw_kgb16

    USE rrsw_kg16, ONLY: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no16

    INTEGER, PARAMETER :: bandNumber = 1, numGPoints = no16

    REAL(kind=wp) :: ncrayl(1)

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/1,1,1,1/), ncrayl)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OSelfAbsorptionCoefficients',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tself,numGPoints,1,1/), selfrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tforeignlower,numGPoints,1,1/), forrefo)

    rayl = ncrayl(1)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(selfrefo, p_io) 
    CALL p_bcast(forrefo, p_io) 
    CALL p_bcast(rayl, p_io) 

  END SUBROUTINE sw_kgb16
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb17

    USE rrsw_kg17, ONLY: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no17

    INTEGER, PARAMETER :: bandNumber = 2
    INTEGER, PARAMETER :: numGPoints = no17

    REAL(kind=wp) :: ncrayl(1)

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,keyupper,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/1,1,1,1/), ncrayl)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keyupper,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OSelfAbsorptionCoefficients',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tself,numGPoints,1,1/), selfrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tforeignlower,numGPoints,1,1/), forrefo(1:3,:))

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/2,1,bandNumber,gPointSetNumber/), &
         (/1,numGPoints,1,1/), forrefo(4,:))

    rayl = ncrayl(1)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(selfrefo, p_io) 
    CALL p_bcast(forrefo, p_io) 
    CALL p_bcast(rayl, p_io) 

  END SUBROUTINE sw_kgb17
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb18

    USE rrsw_kg18, ONLY: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no18

    INTEGER, PARAMETER :: bandNumber = 3
    INTEGER, PARAMETER :: numGPoints = no18

    REAL(kind=wp) :: ncrayl(1)

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,keylower,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/1,1,1,1/), ncrayl)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OSelfAbsorptionCoefficients',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tself,numGPoints,1,1/), selfrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tforeignlower,numGPoints,1,1/), forrefo)

    rayl = ncrayl(1)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(selfrefo, p_io) 
    CALL p_bcast(forrefo, p_io) 
    CALL p_bcast(rayl, p_io) 

  END SUBROUTINE sw_kgb18
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb19

    USE rrsw_kg19, ONLY: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no19

    INTEGER, PARAMETER :: bandNumber = 4
    INTEGER, PARAMETER :: numGPoints = no19

    REAL(kind=wp) :: ncrayl(1)

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,keylower,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/1,1,1,1/), ncrayl)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OSelfAbsorptionCoefficients',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tself,numGPoints,1,1/), selfrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tforeignlower,numGPoints,1,1/), forrefo)

    rayl = ncrayl(1)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(selfrefo, p_io) 
    CALL p_bcast(forrefo, p_io) 
    CALL p_bcast(rayl, p_io) 

  END SUBROUTINE sw_kgb19
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb20

    USE rrsw_kg20, ONLY: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, absch4o, no20

    INTEGER :: ab
    INTEGER, PARAMETER :: bandNumber = 5
    INTEGER, PARAMETER :: numGPoints = no20

    REAL(kind=wp) :: ncrayl(1)

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/1,1,1,1/), ncrayl)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OSelfAbsorptionCoefficients',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tself,numGPoints,1,1/), selfrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tforeignlower,numGPoints,1,1/), forrefo(1:3,:))

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/2,1,bandNumber,gPointSetNumber/), &
         (/1,numGPoints,1,1/), forrefo(4,:))

    CALL getAbsorberIndex('CH4',ab)
    CALL io_inq_varid(rrtmg_sw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,ab,bandNumber,gPointSetNumber/), &
         (/1,1,numGPoints,1,1,1/), absch4o)

    rayl = ncrayl(1)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(selfrefo, p_io) 
    CALL p_bcast(forrefo, p_io) 
    CALL p_bcast(rayl, p_io) 
    CALL p_bcast(absch4o, p_io) 

  END SUBROUTINE sw_kgb20
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb21

    USE rrsw_kg21, ONLY: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no21

    INTEGER, PARAMETER :: bandNumber = 6
    INTEGER, PARAMETER :: numGPoints = no21

    REAL(kind=wp) :: ncrayl(1)

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,keylower,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/1,1,1,1/), ncrayl)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keyupper,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OSelfAbsorptionCoefficients',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tself,numGPoints,1,1/), selfrefo )

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tforeignlower,numGPoints,1,1/),  forrefo(1:3,:))

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/2,1,bandNumber,gPointSetNumber/), &
         (/1,numGPoints,1,1/), forrefo(4,:))

    rayl = ncrayl(1)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(selfrefo, p_io) 
    CALL p_bcast(forrefo, p_io) 
    CALL p_bcast(rayl, p_io) 

  END SUBROUTINE sw_kgb21
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb22	

    USE rrsw_kg22, ONLY: sfluxrefo, kao, kbo, selfrefo, forrefo, rayl, no22

    INTEGER, PARAMETER :: bandNumber = 7
    INTEGER, PARAMETER :: numGPoints = no22

    REAL(kind=wp) :: ncrayl(1)

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,keylower,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/1,1,1,1/), ncrayl)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OSelfAbsorptionCoefficients',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tself,numGPoints,1,1/), selfrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tforeignlower,numGPoints,1,1/), forrefo)

    rayl = ncrayl(1)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(selfrefo, p_io) 
    CALL p_bcast(forrefo, p_io) 
    CALL p_bcast(rayl, p_io) 

  END SUBROUTINE sw_kgb22
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb23		

    USE rrsw_kg23, ONLY: sfluxrefo, kao, selfrefo, forrefo, raylo, no23

    INTEGER, PARAMETER :: bandNumber = 8
    INTEGER, PARAMETER :: numGPoints = no23

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), raylo)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OSelfAbsorptionCoefficients',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tself,numGPoints,1,1/), selfrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tforeignlower,numGPoints,1,1/), forrefo)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(selfrefo, p_io) 
    CALL p_bcast(forrefo, p_io) 
    CALL p_bcast(raylo, p_io) 

  END SUBROUTINE sw_kgb23
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb24	

    USE rrsw_kg24, ONLY: sfluxrefo, kao, kbo, selfrefo, forrefo, &
         raylao, raylbo, abso3ao, abso3bo, no24

    INTEGER :: ab
    INTEGER, PARAMETER :: bandNumber = 9
    INTEGER, PARAMETER :: numGPoints = no24

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,keylower,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,keylower,1,1/), raylao)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), raylbo)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OSelfAbsorptionCoefficients',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tself,numGPoints,1,1/), selfrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tforeignlower,numGPoints,1,1/), forrefo)

    CALL getAbsorberIndex('O3',ab)
    CALL io_inq_varid(rrtmg_sw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,ab,bandNumber,gPointSetNumber/), &
         (/1,1,numGPoints,1,1,1/), abso3ao)
    CALL io_inq_varid(rrtmg_sw%file_id,'AbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,ab,bandNumber,gPointSetNumber/), &
         (/1,1,numGPoints,1,1,1/), abso3bo)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(selfrefo, p_io) 
    CALL p_bcast(forrefo, p_io) 
    CALL p_bcast(raylao, p_io) 
    CALL p_bcast(raylbo, p_io) 
    CALL p_bcast(abso3ao, p_io) 
    CALL p_bcast(abso3bo, p_io) 

  END SUBROUTINE sw_kgb24
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb25		

    USE rrsw_kg25, ONLY: sfluxrefo, kao, raylo, abso3ao, abso3bo, no25

    INTEGER :: ab	
    INTEGER, PARAMETER :: bandNumber = 10
    INTEGER, PARAMETER :: numGPoints = no25

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), raylo)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,plower,numGPoints,1,1/), kao)

    CALL getAbsorberIndex('O3',ab)
    CALL io_inq_varid(rrtmg_sw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,ab,bandNumber,gPointSetNumber/), &
         (/1,1,numGPoints,1,1,1/), abso3ao)

    CALL io_inq_varid(rrtmg_sw%file_id,'AbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,ab,bandNumber,gPointSetNumber/), &
         (/1,1,numGPoints,1,1,1/), abso3bo)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(raylo, p_io) 
    CALL p_bcast(abso3ao, p_io) 
    CALL p_bcast(abso3bo, p_io) 

  END SUBROUTINE sw_kgb25
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb26

    USE rrsw_kg26, ONLY: sfluxrefo, raylo, no26

    INTEGER, PARAMETER :: bandNumber = 11
    INTEGER, PARAMETER :: numGPoints = no26

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), raylo)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(raylo, p_io) 

  END SUBROUTINE sw_kgb26
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb27

    USE rrsw_kg27, ONLY: sfluxrefo, kao, kbo, raylo, no27

    INTEGER, PARAMETER :: bandNumber = 12
    INTEGER, PARAMETER :: numGPoints = no27

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), raylo)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(raylo, p_io) 

  END SUBROUTINE sw_kgb27
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb28		

    USE rrsw_kg28, ONLY: sfluxrefo, kao, kbo, rayl, no28

    INTEGER, PARAMETER :: bandNumber = 13
    INTEGER, PARAMETER :: numGPoints = no28  

    REAL(kind=wp) :: ncrayl(1)

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,keyupper,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/1,1,1,1/), ncrayl)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keylower,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/keyupper,Tdiff,pupper,numGPoints,1,1/), kbo)

    rayl = ncrayl(1)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(rayl, p_io) 

  END SUBROUTINE sw_kgb28
  !*******************************************************************************

  !*******************************************************************************
  SUBROUTINE sw_kgb29	

    USE rrsw_kg29, ONLY: sfluxrefo, kao, kbo, selfrefo, forrefo, &
         absh2oo, absco2o, rayl, no29

    INTEGER :: ab
    INTEGER, PARAMETER :: bandNumber = 14
    INTEGER, PARAMETER :: numGPoints = no29

    REAL(kind=wp) :: ncrayl(1)

    IF (p_parallel_io) THEN
    CALL io_inq_varid(rrtmg_sw%file_id,'SolarSourceFunctionLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/numGPoints,1,1,1/), sfluxrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'RayleighExtinctionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/1,1,1,1/), ncrayl)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,plower,numGPoints,1,1/), kao)

    CALL io_inq_varid(rrtmg_sw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,1,bandNumber,gPointSetNumber/), &
         (/1,Tdiff,pupper,numGPoints,1,1/), kbo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OSelfAbsorptionCoefficients',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tself,numGPoints,1,1/), selfrefo)

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,bandNumber,gPointSetNumber/), &
         (/Tforeignlower,numGPoints,1,1/), forrefo(1:3,:))

    CALL io_inq_varid(rrtmg_sw%file_id,'H2OForeignAbsorptionCoefficientsUpperAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/2,1,bandNumber,gPointSetNumber/), &
         (/1,numGPoints,1,1/), forrefo(4,:))	

    CALL getAbsorberIndex('H2O',ab)
    CALL io_inq_varid(rrtmg_sw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,ab,bandNumber,gPointSetNumber/), &
         (/1,1,numGPoints,1,1,1/), absh2oo)

    CALL getAbsorberIndex('CO2',ab)
    CALL io_inq_varid(rrtmg_sw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
    CALL io_get_vara_double(rrtmg_sw%file_id, varid, &
         (/1,1,1,ab,bandNumber,gPointSetNumber/), &
         (/1,1,numGPoints,1,1,1/), absco2o)

    rayl = ncrayl(1)

    ENDIF
    CALL p_bcast(sfluxrefo, p_io) 
    CALL p_bcast(kao, p_io) 
    CALL p_bcast(kbo, p_io) 
    CALL p_bcast(selfrefo, p_io) 
    CALL p_bcast(forrefo, p_io) 
    CALL p_bcast(rayl, p_io) 
    CALL p_bcast(absh2oo, p_io) 
    CALL p_bcast(absco2o, p_io) 

  END SUBROUTINE sw_kgb29
  !*******************************************************************************

END MODULE mo_srtm_netcdf

